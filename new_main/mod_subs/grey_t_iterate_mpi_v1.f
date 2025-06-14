!
! Subroutine designed to compute a new estimate of the Temperature structure
! through iterating on the grey temperature structure. Written for use
! with the hydrostatic structure iteration. Based on code in set_new_model_estmates.f
!
	SUBROUTINE GREY_T_ITERATE_MPI_V1(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
	1            LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
	USE SET_KIND_MODULE
	USE MPI
	USE ANG_QW_MOD
	USE MOD_CMFGEN
	USE OPAC_MOD
	USE CONTROL_VARIABLE_MOD
	USE LINE_VEC_MOD
	USE LINE_MOD
	IMPLICIT NONE
!
	INTEGER NC
	INTEGER ND
	INTEGER NP
	INTEGER NT
!
	INTEGER NCF				!Number of continuum frequencies
	INTEGER N_LINE_FREQ			!Number of lines
	INTEGER LUER				!Unit for error messages
	INTEGER LUIN				!Unit for input
	INTEGER MAX_SIM				!Maximum number of lines that can be treated simultaneously.
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) Z_POP(NT)			!Vector containing Z of atom/ion (not core)
!
	REAL(KIND=LDP) FQW(NCF)
	REAL(KIND=LDP) NU_EVAL_CONT(NCF)
	REAL(KIND=LDP) NU(NCF)
!
! These are set in CMFGEN.
!
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) OPLIN,EMLIN
!
! Arrays for improving on the initial T structure --- partition functions.
! Need one for each atomic species.
!
        REAL(KIND=LDP), ALLOCATABLE :: U_PAR_FN(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: PHI_PAR_FN(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: HIGH_POP(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: Z_PAR_FN(:)
!
	REAL(KIND=LDP) TGREY(ND)			!Grey temperature structure
	REAL(KIND=LDP) T_SAVE(ND)
	REAL(KIND=LDP) ROSSMEAN(ND)			!Rosseland mean opacity
	REAL(KIND=LDP) PLANCKMEAN(ND)                   !Planck mean opacity
!
! These are all work vectors.
!
	REAL(KIND=LDP) RJ(ND)
	REAL(KIND=LDP) DTAU(ND)
	REAL(KIND=LDP) Z(ND)
	REAL(KIND=LDP) dCHIdR(ND)
!
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) TC(ND)
	REAL(KIND=LDP) QH(ND)
	REAL(KIND=LDP) Q(ND)
	REAL(KIND=LDP) GAM(ND)
	REAL(KIND=LDP) GAMH(ND)
	REAL(KIND=LDP) H(ND)
	REAL(KIND=LDP) SOB(ND)
	REAL(KIND=LDP) XM(ND)
	REAL(KIND=LDP) FEDD(ND)
!
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) HBC_J
	REAL(KIND=LDP) NU_DOP
	REAL(KIND=LDP) FL		!Current frequency
	REAL(KIND=LDP) CONT_FREQ	!Frequency at which current ETA/CHI was evaluated
!
	INTEGER FREQ_INDX 	!Index of current frequency in NU
	INTEGER ML		!Same as FREQ_INDX
	INTEGER LAST_LINE	!Next line to be accessed
	INTEGER GREY_IOS	!Used to return error if GREY_SCL_FAC_IN can't be read.
!
	INTEGER I,J,K,L
	INTEGER ISPEC
	INTEGER ICNT
	INTEGER ID
	INTEGER ID_SAV
	INTEGER NL,NUP
	INTEGER MNL_F,MNUP_F
	INTEGER MNL,MNUP
	INTEGER MAIN_COUNTER
!
	LOGICAL LST_DEPTH_ONLY
	LOGICAL FIRST
	LOGICAL COMPUTED
	LOGICAL TMP_LOG
!
	CHARACTER*80 TMP_STRING
	CHARACTER*20 SECTION
	CHARACTER(LEN=12) SAVED_TWO_PHOTON_METHOD

!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
	LST_DEPTH_ONLY=.FALSE.
	SECTION='CONTINUUM'
	GREY_IOS=0
	SAVED_TWO_PHOTON_METHOD=TWO_PHOTON_METHOD
	TWO_PHOTON_METHOD='LTE'
!
! temperature distribution at depth corresponds to the GREY solution.
! We use the Rosseland mean opacities to evaluate the GREY temperature
! distribution. We then compute non-LTE partition functions which
!
! This page computes the Rosseland mean opacity from the temperature
! distribution and the population levels. TA is a working vector. The
! Rosseland opacity is given in ROSSMEAN.
!
	IF(MYPE .EQ. 0)WRITE(6,*)'Executing barrier statement before T iterate)'
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	CALL TUNE(1,'T_ITERATE')
	MAIN_COUNTER=1
	DO WHILE (MAIN_COUNTER .LE. MAX_NO_GREY_ITERATIONS)
!
!	    WRITE(6,*)'Begining T iteration -- NEED TO FIX CODE for GRID option'; FLUSH(UNIT=6)
	    IF(.NOT. ALLOCATED(U_PAR_FN))THEN
	      ALLOCATE (U_PAR_FN(DST:DEND,NUM_IONS),STAT=IOS)
	      IF(IOS .EQ. 0)ALLOCATE (PHI_PAR_FN(DST:DEND,NUM_IONS),STAT=IOS)
	      IF(IOS .EQ. 0)ALLOCATE (HIGH_POP(DST:DEND,NUM_SPECIES),STAT=IOS)
	      IF(IOS .EQ. 0)ALLOCATE (Z_PAR_FN(NUM_IONS),STAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Unable to allocate PHI_PAR_FN in SET_NEW_MODEL_ESTIMATES'
	        STOP
	      END IF
	    END IF
! 
!
! Set 2-photon data with current atomic models and populations.
!
	    DO ID=1,NUM_IONS-1
	       ID_SAV=ID
	       CALL SET_TWO_PHOT_ATM_MPI_V1(ION_ID(ID), ID_SAV,
	1          ATM(ID)%XzVLTE,          ATM(ID)%NXzV,
	1          ATM(ID)%XzVLTE_F_ON_S,   ATM(ID)%XzVLEVNAME_F,
	1          ATM(ID)%EDGEXzV_F,       ATM(ID)%GXzV_F,
	1          ATM(ID)%F_TO_S_XzV,      ATM(ID)%NXzV_F, DST, DEND,
	1          ATM(ID)%ZXzV,            ATM(ID)%EQXzV,  ATM(ID)%XzV_PRES)
	    END DO
!
! 
!
! We ensure that LAST_LINE points to the first LINE that is going to
! be handled in the BLANKETING portion of the code.
!
	    LAST_LINE=0	    		!Updated as each line is done
	    DO WHILE(LAST_LINE .LT. N_LINE_FREQ .AND.
	1             VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	            LAST_LINE=LAST_LINE+1
	    END DO
!
! ROSSMEAN is initially used to accumulate the integral of 1/chi (weighted
! by dB/DT). After the frequency loop it is corrected so that it contains
! Rosseland mean opacity.
!
	    CALL DP_ZERO(ROSSMEAN,ND)
	    CALL DP_ZERO(PLANCKMEAN,ND)
	    TSTAR=T(ND)			!Required for IC in OPACITIES
	    CONT_FREQ=0.0_LDP
	    DO ML=1,NCF
	      FREQ_INDX=ML
	      FL=NU(ML)
!
	      IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	        COMPUTE_NEW_CROSS=.TRUE.
	        CONT_FREQ=NU_EVAL_CONT(ML)
	      ELSE
	        COMPUTE_NEW_CROSS=.FALSE.
	      END IF
	      CALL SET_PHOT_CROSS_SECTIONS_V1(CONT_FREQ,NU,NU_EVAL_CONT,FREQ_INDX,NCF)
	      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
	      CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
!
! Compute the line opacity. This initializes the storage locations on
! entry when FREQ_INDX=ML=1.
!
              CALL SET_LINE_OPAC(POPS,NU,FREQ_INDX,LAST_LINE,N_LINE_FREQ,
	1            LST_DEPTH_ONLY,LUER,ND,NT,NCF,MAX_SIM)

!
! Now add in line opacity to continuum opacity.
!
	      DO SIM_INDX=1,MAX_SIM
	        IF(RESONANCE_ZONE(SIM_INDX))THEN
	          DO I=DST,DEND
	            CHI(I)=CHI(I) + CHIL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	            ETA(I)=ETA(I) + ETAL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	          END DO
	        END IF
	      END DO
!
! CHECK for negative line opacities.
!	        CHI_NOSCAT(I)=MAX(0.0_LDP,CHI(I)-ESEC(I))
!	        IF(CHI(I) .LT. 0.1_LDP*ESEC(I))CHI(I)=0.1_LDP*ESEC(I)
!
	      DO I=DST,DEND
	        CHI_NOSCAT(I)=MAX(0.0_LDP,CHI(I)-CHI_SCAT(I))
	        IF(CHI(I) .LT. 0.1_LDP*CHI_SCAT(I))CHI(I)=0.1_LDP*CHI_SCAT(I)
	      END DO
!
! Note division by T**2 is included with Stefan-Boltzman constant.
!
	      T1=-HDKT*NU(ML)
	      T2=FQW(ML)*TWOHCSQ*(NU(ML)**3)
	      T3=-T1*FQW(ML)*TWOHCSQ*(NU(ML)**3)
	      DO I=DST,DEND
	        PLANCKMEAN(I)=PLANCKMEAN(I) + T2*CHI_NOSCAT(I)*EMHNUKT(I)/(1.0_LDP-EMHNUKT(I))
	        ROSSMEAN(I)=ROSSMEAN(I) + T3*EMHNUKT(I)/CHI(I)/(1.0_LDP-EMHNUKT(I))**2
	      END DO
	    END DO
!
! Compute CHI, and then optical depth scale.
! Stefan-Boltzman constant *1D-15*1D+16/PI (T**4/PI). NB --- T1 is a factor
! of 10^15 larger than in MAINGEN as FQW has already been multiplied by
! 10^15 for dv integrations.
!
! If clumping is important, we need to correct the Rosseland mean opacity
! for clumping. Since it is a simple scale factor at each depth, we can do
! it here, rather than adjust CHI for each frequency.
!
	    T1=1.8047E+11_LDP
	    DO I=DST,DEND
	      ROSSMEAN(I)=4.0_LDP*CLUMP_FAC(I)*T1*(T(I)**5)/ROSSMEAN(I)
	      PLANCKMEAN(I)=CLUMP_FAC(I)*PLANCKMEAN(I)/T1/(T(I)**4)
	    END DO
!
	    CALL MPI_ALLREDUCE(ROSSMEAN,TA,ND,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,IERR)
	    ROSSMEAN(1:ND)=TA(1:ND)
            CALL MPI_ALLREDUCE(PLANCKMEAN,TA,ND,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,IERR)
	    PLANCKMEAN(1:ND)=TA(1:ND)
!
	    IF(MYPE .EQ. 0)THEN
	      CALL WRITV(T,ND,'Current temperature',88)
	      CALL WRITV(ROSSMEAN,ND,'Rosseland Mean Opacity',88)
	      CALL WRITV(PLANCKMEAN,ND,'Planck Mean Opacity',88)
	      TA(1:ND)=1.0E-10_LDP*ROSSMEAN(1:ND)/DENSITY(1:ND)
	      TB(1:ND)=1.0E-10_LDP*PLANCKMEAN(1:ND)/DENSITY(1:ND)
	      CALL WRITV(TA,ND,'Rosseland mean mass absorption coefficient',88)
	      CALL WRITV(TB,ND,'Planck mean mass absorption coefficient',88)
	      FLUSH(UNIT=88)
	    END IF
!
	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
! 
!
! Check that inner boundary is deep enough so that LTE can be fully recovered. SOURCE and
! TC are used as temporary vectors.
!
	    IF(MYPE .EQ. 0 .AND. MAIN_COUNTER .EQ. 1)THEN
	      CALL TORSCL(TA,ROSSMEAN,R,TB,TC,ND,METHOD,' ')
	      CALL ESOPAC(ESEC,ED,ND)
	      CALL TORSCL(TB,ESEC,R,SOURCE,TC,ND,METHOD,' ')
	      WRITE(LUER,*)' '
	      WRITE(LUER,'(A,ES10.3)')' Thompson scattering optical depth at inner boundary is:',TB(ND)
	      WRITE(LUER,'(A,ES10.3)')' Rosseland optical depth at inner boundary is:          ',TA(ND)
	      WRITE(LUER,'(A,ES10.3)')' Rosseland optical depth at outer boundary is:          ',TA(1)
	      WRITE(LUER,*)' '
	      IF(TA(ND) .LT. 10.0_LDP)THEN
	        WRITE(LUER,*)('*',I=1,70)
	        WRITE(LUER,*)('*',I=1,70)
	        WRITE(LUER,*)' '
	        WRITE(LUER,*)'Warning --- your core optical depth is probably too low'
	        WRITE(LUER,*)'You should use a value in excess of 10'
	        WRITE(LUER,*)' '
	        WRITE(LUER,*)('*',I=1,70)
	        WRITE(LUER,*)('*',I=1,70)
	      END IF
	    END IF
!
! Compute the grey temperature structure and the Rosseland optical depth scale.
! ROSSMEAN already includes the effect of clumping.
!
	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	    CHI(1:ND)=ROSSMEAN(1:ND)
	    CALL COMP_GREY_V4(POPS,TGREY,TA,ROSSMEAN,PLANCKMEAN,COMPUTED,LUER,NC,ND,NP,NT)
	    IF(.NOT. COMPUTED)THEN
	      WRITE(LUER,*)'Unable to compute grey temperature structure'
	      WRITE(LUER,*)'As this is needed to provide initial T estimate, stopping code'
	      STOP
	    END IF
!
! SCALE_GREY modifies the computed grey temperature distribution according
! to that computed in a previous model.
!
! i.e. TGREY = TGREY . (T/TGREY)_old
!
	   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	   IF(GREY_IOS .EQ. 0)THEN
	      CALL SCALE_GREY(TGREY,TA,GREY_IOS,LUIN,ND)
	      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	   END IF
!
! Now correct T distribution towards grey value. As we don't require the
! old T, we can overwrite it straight away. If we multiplied T1 by a number
! less than  unity, this would be equivalent to only a partial correction
! of T towards TGREY:
!          GREY_PAR=0 set T=TGREY
!          GREY_PAR=INFINITY leaves T=T.
!
! T3 and T2 are used to determine the current largest correction.
!
	    TC(1:ND)=T(1:ND)
	    IF( MAIN_COUNTER .EQ. 1)THEN
	      DO I=1,ND
	        T_SAVE(I)=T(I)		!Save original T for use when
	      END DO                    !correcting T towards TGREY.
	    END IF
 	    T2=0.0_LDP
	    DO I=DST,DEND
	      IF(GREY_PAR .LE. 0)then
	        T1=1.0_LDP
	      ELSE
	        T1=1.0_LDP-EXP(-TA(I)/GREY_PAR)
	      END IF
	      IF(TA(I) .LT. 0.1_LDP*GREY_PAR)T1=0.0_LDP     !Changed T1 to TA(I) [14-Jan-2009]
	      T3=ABS( T1*(TGREY(I)-T(I)) )
	      T(I)=T1*TGREY(I)+(1.0_LDP-T1)*T_SAVE(I)
	      T(I)=MAX(T(I),0.95_LDP*T_MIN)
	      T2=MAX(T3/T(I),T2)
	    END DO
	    CALL MPI_ALLREDUCE(MPI_IN_PLACE,T2,IONE,MY_MPI_DP,MPI_MAX,MPI_COMM_WORLD,IERR)
	    CALL ALL_GATHERV_VEC_MPI_V1(T,ND)
	    IF(MYPE .EQ. 0)THEN
	      WRITE(LUER,'('' Largest correction to T in GREY initialization loop is '',ES9.2,'' %'')')100.0*T2
	    END IF
!
! Now compute non-LTE partition functions. These assume that the
! departure coefficients are independent of Temperature. This
! is a good assumption at depth where b is approximately unity.
!
! GAM_SPECIES is used as a storage location for the population of the
! highest ionization stage. Must be done in forward direction.
!
! NB: After calling PAR_FUN_V2, ROOT(ID)%XzV_F will contain DCs -
!       NOT populations.
!
! TMP_STRING is used to indicate whether we interpolate in depature
! coefficients (DC) or excitation temperatures (TX).
!
	    TMP_STRING='DC'
	    IF(DC_INTERP_METHOD .EQ. 'RTX')TMP_STRING='TX'
	    DO ID=1,NUM_IONS
	      J=ID
	      ISPEC=SPECIES_LNK(ID)
	      CALL PAR_FUN_MPI_V1(U_PAR_FN, PHI_PAR_FN, Z_PAR_FN, HIGH_POP,
	1          ATM(ID)%XzV_F,     ATM(ID)%LOG_XzVLTE_F,  ATM(ID)%W_XzV_F,
	1          ATM(ID)%DXzV_F,    ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,
	1          ATM(ID)%GIONXzV_F, ATM(ID)%ZXzV,T, TC, ED,
	1          ATM(ID)%NXzV_F, DST, DEND, ND, 
	1          ISPEC, NUM_SPECIES, J, NUM_IONS,
	1          ATM(ID)%XzV_PRES,ION_ID(ID),TMP_STRING)
	   END DO
	   J=DEND-DST+1
!
! The non-LTE partition functions are density independent, provided
! we assume the departure coefficients remain fixed.
!
! We now evaluate the contribution to the electron density by each
! species, using the non-LTE partition functions.
!
! We use H for ED(est)
! We use QH for dED(est)/dT.
!
	    T1=1.0_LDP
	    ICNT=0
	    DO WHILE (T1 .GT. 1.0E-04_LDP)
	      FIRST=.TRUE.
!
! Recall GAM_SPECIES is set to be the population of the highest ionization
! stage.
!
	      DO ISPEC=1,NUM_SPECIES
	        ID=SPECIES_BEG_ID(ISPEC)
	        J=SPECIES_END_ID(ISPEC)
	        IF(SPECIES_PRES(ISPEC))THEN
	          CALL EVAL_ED_MPI_V1(H,QH,U_PAR_FN,PHI_PAR_FN,HIGH_POP,
	1                  Z_PAR_FN,ED,POP_SPECIES(1,ISPEC),
	1                  XM,TB,TC,ISPEC,NUM_SPECIES,
	1                  ID,J,NUM_IONS,DST,DEND,ND,FIRST)
	        END IF
	      END DO
!	    WRITE(MYPE+230,*)'Aft EVAL_ED',MYPE,DST,DEND
!	    CALL WRITE_VEC(ED(DST),DEND-DST+1,'ED',230+MYPE); FLUSH(UNIT=230+MYPE)
!	    CALL WRITE_VEC(H(DST),DEND-DST+1,'H',230+MYPE); FLUSH(UNIT=230+MYPE)
!
	      T1=0.0
	      DO I=DST,DEND
	        TA(I)=-(H(I)-ED(I))/(QH(I)-1.0_LDP)/ED(I)
	        T1=MAX(T1,ABS(TA(I)))
	        IF(TA(I) .LT. -0.9_LDP)TA(I)=-0.9_LDP
	        IF(TA(I) .GT. 9.0_LDP)TA(I)=9.0_LDP
	        ED(I)=ED(I)*(1.0_LDP+TA(I))
	      END DO
	      ICNT=ICNT+1
	      IF(ICNT .GT. 20)THEN
	        WRITE(LUER,*)'Error --- Computation of ED in EVAL_ED section'//
	1                 ' has taken more than 20 iterations'
	        WRITE(LUER,*)'Current error T1 is',T1
	        STOP
	      END IF
	    END DO
	    CALL ALL_GATHERV_VEC_MPI_V1(ED,ND)
!	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
! 
!
! Now need to compute LTE populations, and populations.
! Since T and Ne have altered, we revise the vectors for evaluating the
! level dissolution. These constants are the same for all species. These are
! stored in a common block, and are required by SUP_TO_FULL and LTE_POP_WLD.
!
! NB: POPION will also alter but in W-R and LBV's all species will be ionized,
! and hence POPION will not change from iteration to iteration. In any event,
! it has a smaller effect than changes in Ne.
!
	    CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! We do low ionization species second, as first need DION.
!
	    DO ISPEC=1,NUM_SPECIES
	      FIRST=.TRUE.
	      DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	        IF(ATM(ID)%XzV_PRES)THEN
	          CALL LTEPOP_WLD_V2(ATM(ID)%XzVLTE_F, ATM(ID)%LOG_XzVLTE_F,  ATM(ID)%W_XzV_F,
	1               ATM(ID)%EDGEXzV_F,  ATM(ID)%GXzV_F,  ATM(ID)%ZXzV,
	1               ATM(ID)%GIONXzV_F,  ATM(ID)%NXzV_F,  ATM(ID)%DXzV_F,
	1               ED,T, DST, DEND, ND)
	          CALL CNVT_FR_DC_V2(ATM(ID)%XzV_F, ATM(ID)%LOG_XzVLTE_F,
	1               ATM(ID)%DXzV_F,   ATM(ID)%NXzV_F,
	1               TB,               TA, DST, DEND, ND, FIRST,      ATM(ID+1)%XzV_PRES)
	          IF(ID .NE. SPECIES_BEG_ID(ISPEC))ATM(ID-1)%DXzV_F(DST:DEND)=TB(DST:DEND)
	        END IF
	      END DO
	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
! We need to scale the populations to ensure that the change in temperature
!   has not causes some population to blow up. We always do this --- the
! DO_POP_SCALE option has no effect.
!
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        CALL SCALE_POPS_MPI_V1(ATM(ID)%XzV_F, ATM(ID)%DXzV_F,
	1           POP_SPECIES(1,SPECIES_LNK(ID)),TA, ATM(ID)%NXzV_F, DST, DEND, ND)
	      END DO
	    END DO
!
	    MAIN_COUNTER=MAIN_COUNTER+1
!
! We now need to compute the populations for the model atom with Super-levels.
! We do this in reverse order (i.e. highest ionization stage first) in order
! that we the ion density for the lower ionization stage is available for
! the next call.
!
! For 1st call to FULL_TO_SUP, Last line contains FeX etc as FeXI not installed.
!
	    DO ID=NUM_IONS-1,1,-1
	      CALL FULL_TO_SUP_MPI_V1(
	1      ATM(ID)%XzV,   ATM(ID)%NXzV,       ATM(ID)%DXzV,      ATM(ID)%XzV_PRES,
	1      ATM(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F,    ATM(ID)%DXzV_F,
	1      ATM(ID+1)%XzV, ATM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES, DST, DEND)
	    END DO
	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
!	  IF(MYPE .EQ. 0)THEN
!	    WRITE(6,*)MYPE,'B-Calling FULL_TO_SUP'
!	    DO ID=NUM_IONS-1,1,-1
!	       CALL FULL_TO_SUP(
!	1          ROOT(ID)%XzV,   ATM(ID)%NXzV,      ROOT(ID)%DXzV,   ATM(ID)%XzV_PRES,
!	1          ROOT(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ROOT(ID)%DXzV_F,
!	1          ROOT(ID+1)%XzV, ATM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES, DST, DEND)
!	    END DO
!	  END IF
!
! Store all quantities in POPS array. This is done here (rather than
! after final iteration) as it enable POPION to be readily computed.
!
	    POPS=0.0_LDP
	    DO ID=1,NUM_IONS-1
	      CALL IONTOPOP(POPS, ATM(ID)%XzV, ATM(ID)%DXzV, ED,T,
	1         ATM(ID)%EQXzV, ATM(ID)%NXzV, NT, DST, DEND,  ND,
	1         ATM(ID)%XzV_PRES)
	    END DO
	    J=ND*NT
	    CALL MPI_ALLREDUCE(MPI_IN_PLACE,POPS,J,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,IERR)
	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
	    IF(MYPE .EQ. 0)THEN
	      DO ID=1,NUM_IONS-1
	        CALL POPTOION(POPS, ROOT(ID)%XzV, ROOT(ID)%DXzV,ED,T,
	1              ATM(ID)%EQXzV, ATM(ID)%NXzV, NT, IONE, ND, ND, ATM(ID)%XzV_PRES)
	      END DO
	      DO ID=NUM_IONS-1,1,-1
	         CALL FULL_TO_SUP_MPI_V1(
	1          ROOT(ID)%XzV,   ATM(ID)%NXzV,       ROOT(ID)%DXzV,      ATM(ID)%XzV_PRES,
	1          ROOT(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F,     ROOT(ID)%DXzV_F,
	1          ATM(ID+1)%XzV, A TM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES, DST, DEND)
	       END DO
	    END IF
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	    DO J=1,ND
	      POPION(J)=0.0_LDP
	      DO I=1,NT
	        IF(Z_POP(I) .GT. 0.01_LDP)POPION(J)=POPION(J)+POPS(I,J)
	      END DO
	    END DO
!
! While the following may seem superfolous, it ensures absolute consistency. Its possible
! that, for H, He etc that the upper levels with interpolating sequences may not be fully
! consistent (due to rounding errors, from another model, etc).
!
	    CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!	    CALL MPI_FINALIZE (ierr)
!	    STOP
!
! Revise ALL LTE populations.
!
	    CALL EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
	    IF(MYPE .EQ. 0)THEN
	      CALL EVAL_ROOT_LTE_MPI_V1(DO_LEV_DISSOLUTION,ND)
	    END IF
!
	END DO	
	CALL TUNE(2,'T_ITERATE')
	IF(MYPE .EQ. 0)THEN
	   WRITE(6,*)'Finished T iteration in GREY_T_ITERATE_MPI_V1'
	   FLUSH(UNIT=6)
	END IF
!
	IF(ALLOCATED(U_PAR_FN))THEN
	  DEALLOCATE (U_PAR_FN,STAT=IOS)
	  DEALLOCATE (PHI_PAR_FN,STAT=IOS)
	  DEALLOCATE (HIGH_POP,STAT=IOS)
	  DEALLOCATE (Z_PAR_FN,STAT=IOS)
	  CALL GATHER_ATM_MPI_V1(ND)
	END IF
!
!	TMP_LOG=.TRUE.
!	CALL WR2D_MPI_V1(ATM(1)%XzVLTE,ATM(1)%NXzV,DST,DEND,ND,'Hyd SL LTEPOP',' ',TMP_LOG,410)
!	IF(MYPE .EQ. 0)WRITE(6,*)'Called WR2D_MPI',MYPE
!
! Restore two photon method option.
!
	TWO_PHOTON_METHOD=SAVED_TWO_PHOTON_METHOD
	IF(MYPE .EQ. 0)WRITE(6,*)'Leaving GREY_T_ITERATE_MPI_V1'; FLUSH(UNIT=6)
!
	RETURN
	END
