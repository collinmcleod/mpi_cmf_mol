! 
!
! Associate charge exchange reactions with levels in the model atoms.
!
	CALL SET_CHG_LEV_ID_V4(ND,LUMOD)
	CALL VERIFY_CHG_EXCH_V3()
!
! Determine the number of important variables for each species, and
! set the links.
!
	CALL DETERMINE_NSE(NION,XRAYS)
        CALL CREATE_IV_LINKS_V2(NT,NION)
	CALL WR_LEVEL_LINKS
!
! Read in data for treating non-thermal ionization.
!
	IF(TREAT_NON_THERMAL_ELECTRONS)THEN
	  CALL READ_ARNAUD_ION_DATA(ND)
	  CALL READ_NT_OMEGA_DATA()
	END IF
!
! Allocate memory for STEQ and BA arrays.
!
        CALL SET_BA_STORAGE_MPI_V1(NT,NUM_BNDS,ND,NION)
	DO ID=1,NION
	  IF(SE(ID)%XzV_PRES)THEN
	     SE(ID)%STEQ=0.0_LDP
	  END IF
	END DO
!	CALL WR_ASCI_STEQ_MPI_V1(NION,DST,DEND,ND,'testin STEQ',330)
!
! Read in BA and STEQ arrays. We only attempt this if we have an existing
! model.
!
	CHK=.FALSE.
	IF(.NOT. NEWMOD)THEN
          CALL READ_BA_DATA_MPI_V1(LU_BA,NION,NUM_BNDS,CHK,FIXED_T,SUCCESS,'BAMAT')
	END IF
	IF(.NOT. SUCCESS .OR. LAMBDA_ITERATION)THEN
	  TMP_LOGICAL=SUCCESS
	  CALL MPI_BCAST(TMP_LOGICAL,IONE,MPI_LOGICAL,IZERO,MPI_COMM_WORLD,IERR)
	  IF(TMP_LOGICAL .NEQV. SUCCESS)THEN
	    WRITE(6,*)'Error -- inconsistent SUCCESS option for reading BAMAT access'
	    WRITE(6,*)'MYPE =',MYPE
	    STOP
	  END IF
	  NLBEGIN=0
          COMPUTE_BA=.TRUE.
	  WRBAMAT=.FALSE.
	  IF(N_ITS_TO_FIX_BA .GT. 0)WRBAMAT=.TRUE.
	ELSE IF(NLBEGIN .EQ. -999)THEN		!Indicate completed iteration
	  I=NLBEGIN
	  CALL MPI_BCAST(I,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD,IERR)
	  IF(I .NE. NLBEGIN)THEN
            WRITE(6,*)'Error -- inconsistent NLBEGIN for reading BAMAT access'
            WRITE(6,*)'MYPE =',MYPE
            STOP
          END IF
	  NLBEGIN=0				!hence BA matrix available.
	  COMPUTE_BA=COMPUTE_BARDIN
	  WRBAMAT=.FALSE.
	ELSE
	  WRBAMAT=WRBAMAT_RDIN
 	END IF
	IF(FLUX_CAL_ONLY)THEN
	   COMPUTE_BA=.FALSE.
	   WRBAMAT=.FALSE.
	   LAMBDA_ITERATION=.FALSE.
	   MAXCH=0.0_LDP
!
! For coherent electon scattering, only need 1 iteration.
!
           IF(RD_COHERENT_ES)NUM_ITS_TO_DO=1
	END IF
!
! Removed TMIN consistency check since I now use LOG(LTE pops).
!
	IF(MYPE .EQ. 0)THEN
	  CALL CHECK_IONS_PRESENT(ND,NUM_IONS)
	END IF
!
! Temporary check
!
!	CALL WRITE_SEQ_TIME_FILE_V1(SN_AGE_DAYS,ND,LUIN)
!	CALL TST_RD_EQ_FILE(POPS,ND,NT,LUIN)
!
!************************************************************************************************
!  Running the gamma-ray code for SNe. This calculation only needs to be run one, since it
!  (generally) only depends on the total electron density.
!
!************************************************************************************************
!
	INQUIRE(FILE='GAMRAY_ENERGY_DEP',EXIST=CHK)
	IF(GAMRAY_TRANS .EQ. 'RAD_TRANS' .AND. .NOT. CHK)THEN
	  WRITE(6,'(/,A)')' Running the gamma-ray routine GAMRAY_SUB_V3'
	  CALL TUNE(IONE,'FULL_GAMMA')
	  CALL GAMRAY_SUB_V3(ND,NC,NP,P,R,V,SIGMA,VDOP_VEC,CLUMP_FAC,
	1         MU_AT_RMAX,HQW_AT_RMAX,DELV_FRAC_FG,REXT_FAC,METHOD,
	1         INSTANTANEOUS_ENERGY_DEPOSITION,SN_AGE_DAYS)
	  CALL TUNE(ITWO,'FULL_GAMMA')
	  IF(MYPE .EQ. 0)CALL TUNE(3,' ')
	ELSE IF(GAMRAY_TRANS .EQ. 'RAD_TRANS')THEN
	  WRITE(6,'(/,A)')' Using previosuly computed GAMRAY_ENERGY_DEP file'
	END IF
        WRITE(STRING,*)ML; STRING='1ML='//ADJUSTL(TRIM(STRING))
        CALL WRITV(T,ND,TRIM(STRING),852+MYPE); FLUSH(852+MYPE)

! 
!
!**************************************************************************
!**************************************************************************
!
!                    MAIN ITERATION LOOP
!
!**************************************************************************
!**************************************************************************
!
! MAIN_COUNTER is an integer variable which keeps track of the TOTAL number
! of iterations performed. NB - A NG acceleration is counted as a single
! acceleration.
!
! NUM_ITS_TO_DO indicates the number of iterations left to do. For the
! last iteration this will be zero in the "2000" LOOP.
!
! LST_ITERATION is a logical variable which indicate that the current
! iteration is the last one, and hence DEBUGING and INTERPRETATION data
! should be written out. Its equivalent to NUM_ITS_TO_DO=0 in the loop.
!
	MAIN_COUNTER=NITSF			!Initialize main loop counter
	LAST_LAMBDA=NITSF
	LAST_AV=NITSF
	NEXT_AV=0
!
	IF(DST .EQ. 31)THEN
          WRITE(125,'(10ES14.4)')(ATM(12)%XzV(I,DST),I=4,11)
          WRITE(125,'(10ES14.4)')(ATM(12)%XzV(I,DST+1),I=4,11)
	  FLUSH(UNIT=125)
        END IF
!
20000	CONTINUE
	CALL TUNE(IONE,'GIT')
	NUM_ITS_TO_DO=NUM_ITS_TO_DO-1
	IF(NUM_ITS_TO_DO .EQ. 0)LST_ITERATION=.TRUE.
	MAIN_COUNTER=MAIN_COUNTER+1
!
	IF(MYPE .EQ. 0)THEN
	  WRITE(LUER,*)' Start of GIT loop -- MYPE =',MYPE
	  WRITE(LUER,'(A)')' Variable summary in each threadfollows:'
	  WRITE(LUER,'(A)')' '
	  WRITE(LUER,'(A)')' '
	  WRITE(LUER,'(8A10)')'MYPE','IT_COUNT','RD_LAM','LAMBDA','FIXED_T',
	1                       'COMP._BA','COH._ES','SN_MODEL'
	  FLUSH(LUER)
	END IF
!
	DO I=0,NTHREAD-1
	  IF(MYPE .EQ. I)THEN
	    WRITE(LUER,'(2I10,6(9X,L1))')MYPE,MAIN_COUNTER,RD_LAMBDA,LAMBDA_ITERATION,FIXED_T,
	1                   COMPUTE_BA,COHERENT_ES,SN_MODEL
	    FLUSH(LUER)
	  END IF
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	END DO
!
! Used as a initializing switch for COMP_OBS.
!
	FIRST_OBS_COMP=.TRUE.
!
	IF(LST_ITERATION .AND. WRITE_RATES .AND. MYPE .EQ. 0)THEN
	  CALL GEN_ASCI_OPEN(LU_NET,'NETRATE','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL GEN_ASCI_OPEN(LU_DR,'TOTRATE','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL GEN_ASCI_OPEN(LU_EW,'EWDATA','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL GEN_ASCI_OPEN(LU_HT,'LINEHEAT','UNKNOWN',' ',' ',IZERO,IOS)
	ELSE IF(LST_ITERATION)THEN
	  CALL GEN_ASCI_OPEN(LU_NEG,'NEG_OPAC','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL SET_LINE_BUFFERING(LU_NEG)
	END IF
	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	  OPEN(UNIT=199,FILE='CHECK_ON_BA_UPDATE',STATUS='UNKNOWN',ACTION='WRITE')
	END IF
!
	IF(IMPURITY_CODE)THEN
	  I=WORD_SIZE*(4*ND+1)/UNIT_SIZE
	  OPEN(UNIT=LU_EDD,FILE='IMPURITYJ',FORM='UNFORMATTED',
	1      ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening IMPURITYJ in CMFGEN. MYPE=',MYPE
	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	    CALL MPI_FINALIZE(IERR)
	    STOP
	  END IF
	ELSE
	
	  IF(ACCURATE .OR. EDD_CONT .OR. EDD_LINECONT)THEN
!
! NB: If not ACCURATE, NDEXT was set to ND. The +1 arises since we write
! NU on the same line as RJ. J is used to get the REC_LENGTH, while string
! will contain the date.
!
	    CALL OPEN_RW_EDDFACTOR(R,V,LANG_COORD,ND,
	1     REXT,VEXT,LANG_COORDEXT,NDEXT,
	1     ACCESS_F,NEWMOD,COMPUTE_EDDFAC,USE_FIXED_J,'EDDFACTOR',LU_EDD)
!
	  END IF
	END IF
!
! Now open file containing the electron scatterin J (i.e. the convolution of
! J with the e.s. redistribution function.)
!
! If we don't have EDDFACTOR file it is assumed that we don't have
! J_CONV also.
!
	IF(COMPUTE_EDDFAC)COHERENT_ES=.TRUE.
	IF(.NOT. COHERENT_ES .AND. MYPE .EQ. 0)THEN
	   I=WORD_SIZE*(ND+1)/UNIT_SIZE
	   OPEN(UNIT=LU_ES,FILE='ES_J_CONV',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	     IF(IOS .NE. 0)THEN
	       IF(MYPE .EQ. 0)WRITE(LUER,*)'Error opening ES_J_CONV - will compute new J'
	       COHERENT_ES=.TRUE.
	     END IF
	END IF
!
! 
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
! 
!
! This routine not only evaluates the LTE populations of both model atoms, but
! it also evaluates the dln(LTE Super level Pop)/dT.
!
	CALL EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
	IF(MYPE .EQ. 0)CALL EVAL_ROOT_LTE_MPI_V1(DO_LEV_DISSOLUTION,ND)
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
! 
!
! Set 2-photon data with current atomic models and populations.
!
	DO ID=1,NUM_IONS-1
	  ID_SAV=ID
	  CALL SET_TWO_PHOT_ATM_MPI_V1(ION_ID(ID), ID_SAV,
	1       ATM(ID)%XzVLTE,          ATM(ID)%NXzV,
	1       ATM(ID)%XzVLTE_F_ON_S,   ATM(ID)%XzVLEVNAME_F,
	1       ATM(ID)%EDGEXzV_F,       ATM(ID)%GXzV_F,
	1       ATM(ID)%F_TO_S_XzV,      ATM(ID)%NXzV_F,  DST,   DEND,
	1       ATM(ID)%ZXzV,            ATM(ID)%EQXzV,   ATM(ID)%XzV_PRES)
	END DO
!
	CALL WR2D_MPI_V1(ATM(1)%XzVLTE,ATM(1)%NXzV,DST,DEND,ND,'Hyd SL LTEPOP',' ',.TRUE.,411)
	IF(MYPE .EQ. 0)CLOSE(UNIT=411)
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
! 
!
! Section to compute DT/DR and D(DT/DR)/D? for use in the
! diffusion approximation. DIFFW is used to store the second
! derivatives.
!
! This was based on the subroutine DTSUB. No longer a subroutine as
! too many variables need to be included.
!
! Zero variation of DTDR vector
!
	DO I=1,NT
	  DIFFW(I)=0.0_LDP
	END DO
!
! We ensure that LAST_LINE points to the first LINE that is going to
! be handled in the portion of the code that computes dTdR.
!
	LAST_LINE=0			!Updated as each line is done
	DO WHILE(LAST_LINE .LT. N_LINE_FREQ .AND.
	1             VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	        LAST_LINE=LAST_LINE+1
	END DO
	DO SIM_INDX=1,MAX_SIM
	  LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	END DO
        WRITE(STRING,*)ML; STRING='2ML='//ADJUSTL(TRIM(STRING))
        CALL WRITV(T,ND,TRIM(STRING),852+MYPE); FLUSH(852+MYPE)
!
	CALL TUNE(IONE,'DTDR')
	DTDR=0.0_LDP
	SECTION='DTDR'
	IF(IMPURITY_CODE .OR. USE_FIXED_J .OR. FLUX_CAL_ONLY .OR. (RD_LAMBDA .AND. NEWMOD .AND. .NOT. SN_MODEL))THEN
	  DTDR=(T(ND)-T(ND-1))/(R(ND-1)-R(ND))
	  DIFFW(1:NT)=0.0_LDP
	ELSE 
	  IF(DEND .EQ. ND)THEN
!
! We only need to compute the opacity at the innermost depth, but to save
! programing we will compute it at all depths. As this is only done once
! per iteration, not much time will be wasted.
!
! Setting LST_DEPTH_ONLY to true limits the computation of CHI, ETA, and
! dCHI and dETA to the inner boundary only (in some cases).
!
	    LST_DEPTH_ONLY=.TRUE.
!
! RJ is used in VARCONT to compute the varaition of ETA. In this section
! we only want the variation of CHI, so we initialize its value to zero.
! This prevents a floating point exception.
!
	    RJ(1:ND)=0.0_LDP
	    CONT_FREQ=0.0_LDP
	    FL=NU(1)
	    DO ML=1,NCF
	      FREQ_INDX=ML
!
	      FL_OLD=FL
	      FL=NU(ML)
	      IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	        COMPUTE_NEW_CROSS=.TRUE.
	        CONT_FREQ=NU_EVAL_CONT(ML)
	      ELSE
	        COMPUTE_NEW_CROSS=.FALSE.
	      END IF
!
	      CALL TUNE(IONE,'DTDR_OPAC')
	      CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
	      CALL TUNE(ITWO,'DTDR_OPAC')
!
! 
!
! Compute variation of opacity/emissivity. Store in VCHI and VETA.
!
	      IF(.NOT. LAMBDA_ITERATION .AND. COMPUTE_BA)THEN
	        CALL TUNE(IONE,'DTDR_VOPAC')
	         CALL COMP_VAR_OPAC_MPI_V1(POPS,RJ,FL,CONT_FREQ,FREQ_INDX,
	1                  SECTION,NUM_BNDS,ND,NT,LST_DEPTH_ONLY)
	        CALL TUNE(ITWO,'DTDR_VOPAC')
	      END IF
! 
!
! Compute contribution to CHI and VCHI by lines.
!
! Section to include lines automatically with the continuum.
! Only computes line opacity at final depth point. This is used in the
! computation of dTdR.
!
! NB: Care must taken to ensure that this section remains consistent
!      with that in continuum calculation section.
!
	      CALL TUNE(IONE,'SET_LINE_OPAC')
	        CALL SET_LINE_OPAC(POPS,NU,FREQ_INDX,LAST_LINE,N_LINE_FREQ,
	1            LST_DEPTH_ONLY,LUER,ND,NT,NCF,MAX_SIM)
	      CALL TUNE(ITWO,'SET_LINE_OPAC')
!
! Add in line opacity.
!
	      DO SIM_INDX=1,MAX_SIM
	        IF(RESONANCE_ZONE(SIM_INDX))THEN
	          CHI(ND)=CHI(ND)+CHIL_MAT(ND,SIM_INDX)*LINE_PROF_SIM(ND,SIM_INDX)
	        END IF
	      END DO
!
! Now do the line variation. This presently ignores the effect of a
! temperature variation.
!
	      IF(.NOT. LAMBDA_ITERATION .AND. COMPUTE_BA)THEN
	        DO SIM_INDX=1,MAX_SIM
	          IF(RESONANCE_ZONE(SIM_INDX))THEN
	            NL=SIM_NL(SIM_INDX)
	            NUP=SIM_NUP(SIM_INDX)
	            VCHI(NL,ND)=VCHI(NL,ND)+LINE_PROF_SIM(ND,SIM_INDX)*
	1                LINE_OPAC_CON(SIM_INDX)*L_STAR_RATIO(ND,SIM_INDX)
	            VCHI(NUP,ND)=VCHI(NUP,ND)-LINE_PROF_SIM(ND,SIM_INDX)*
	1               LINE_OPAC_CON(SIM_INDX)*U_STAR_RATIO(ND,SIM_INDX)*GLDGU(SIM_INDX)
	          END IF
	        END DO
	      END IF
!
! 
!
! Update DTDR. Ordering changed to fix issues when BA not computed.
!
	      T1=HDKT*NU(ML)/T(ND)
	      T3=FQW(ML)*TWOHCSQ*( NU(ML)**3 )*T1*EMHNUKT(ND)/
	1         CHI(ND)/T(ND)/(1.0_LDP-EMHNUKT(ND))**2
	      DTDR=DTDR+T3
!
! Set TA = to the variation vector at the inner boundary.
!
	      CALL TUNE(IONE,'DTDR_VEC')
	      IF(.NOT. LAMBDA_ITERATION .AND. COMPUTE_BA)THEN
!
! Increment Parameters
!
	        T1=HDKT*NU(ML)/T(ND)
	        T3=FQW(ML)*TWOHCSQ*( NU(ML)**3 )*T1*EMHNUKT(ND)/CHI(ND)/T(ND)/(1.0_LDP-EMHNUKT(ND))**2
	        DO I=1,NT-1
	          DIFFW(I)=DIFFW(I)+T3*VCHI(I,ND)/CHI(ND)
	        END DO
	        DIFFW(NT)=DIFFW(NT)+T3*(VCHI(NT,ND)/CHI(ND)-(T1*(1.0_LDP+EMHNUKT(ND))
	1           /(1.0_LDP-EMHNUKT(ND))-2.0_LDP)/T(ND))
	      END IF
	      CALL TUNE(ITWO,'DTDR_VEC')
	      WRITE(STRING,*)ML; STRING='ML='//ADJUSTL(TRIM(STRING))
	      WRITE(840,'(A,ES20.12)')TRIM(STRING),CHI(ND)
!
	    END DO
!
! The luminosity of the Sun is 3.826D+33 ergs/sec. For convenience
! DTDR will have the units  (D+04K)/(D+10cm) .
!
	    T1=LUM*7.2685E+11_LDP/R(ND)/R(ND)
	    DTDR=T1/DTDR
	    IF(LAMBDA_ITERATION .OR. .NOT. COMPUTE_BA)THEN
	      DIFFW(1:NT)=0.0_LDP
	    ELSE
	      T1=( DTDR**2 )/T1
	      DO I=1,NT
	        DIFFW(I)=DIFFW(I)*T1
	      END DO
	    END IF
	  END IF
	  I=NTHREAD-1
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
! NB: As the last thread compute DTDR we have to use NYTHREAD-1 as the broadcaster.
!
	  IF(MYPE .EQ. 0)WRITE(6,*)'About to receive DTDR'
	  CALL MPI_BCAST(DTDR,IONE,MPI_DOUBLE_PRECISION,I,MPI_COMM_WORLD,IERR)
	END IF
	CALL TUNE(ITWO,'DTDR')
!
	LST_DEPTH_ONLY=.FALSE.
	IF(MYPE .EQ. 0)THEN
	  WRITE(LUER,'(/,1X,A,ES16.8)')'The value of DTDR is:',DTDR
	  WRITE(LUER,*)'We will now zero the BA matrices'
	  WRITE(LUER,'(A)')' '
	END IF

	IF(DST .EQ. 31)THEN
          WRITE(125,'(A,10ES14.4)')'B ',(ATM(12)%XzV(I,DST),I=4,11)
          WRITE(125,'(A,10ES14.4)')'B ',(ATM(12)%XzV(I,DST+1),I=4,11)
	  FLUSH(UNIT=125)
	END IF
!
! Zero STEQ and BA arrays.
!
	CALL TUNE(IONE,'ZBA')
	FORALL (ID=1:NION)
	  SE(ID)%STEQ   =0.0_LDP
	  SE(ID)%BA     =0.0_LDP
	  SE(ID)%BA_PAR =0.0_LDP
	  SE(ID)%T_EHB =0.0_LDP
	END FORALL
	CALL TUNE(ITWO,'ZBA')
!
	STEQ_ED=0.0_LDP
	STEQ_T=0.0_LDP
	STEQ_T_EHB=0.0_LDP
	STEQ_T_NO_SCL=0.0_LDP
        BA_ED   = 0.0_LDP
        BA_T    = 0.0_LDP
        BA_T_PAR=0.0_LDP
        BA_T_EHB    = 0.0_LDP
        BA_T_PAR_EHB=0.0_LDP
!
	DO ID=1,NUM_IONS
	  ATM(ID)%DIERECOM=0.0_LDP
	  ATM(ID)%ADDRECOM=0.0_LDP
	  ATM(ID)%DIECOOL=0.0_LDP
	  ATM(ID)%X_RECOM=0.0_LDP
	  ATM(ID)%X_COOL=0.0_LDP
	END DO
!
	DIELUM(:)=0.0_LDP
	IF(MYPE .EQ. 0)THEN
	  CALL TUNE(3,' ')
	  WRITE(6,*)'Zeroed BA matrices --- MYPE=',MYPE
	  FLUSH(UNIT=6)
	END IF
!
	  IF(DST .EQ. 31)THEN
            WRITE(125,'(A,10ES14.4)')'F ',(ATM(12)%XzV(I,DST),I=4,11)
            WRITE(125,'(A,10ES14.4)')'F ',(ATM(12)%XzV(I,DST+1),I=4,11)
	    FLUSH(UNIT=125)
	  END IF
! 
!
! Compute the value of the S.E. equations and compute the variation
! matrix for terms that are independent of Jv .
!
! DST and DEND can be adjusted to that we can avoid reading in the entire
! diagonal of the BA array for each call to STEQ_MULTI.
!
! Assume all BA mtarix is in memory/
!
	  IF(MYPE .EQ. 0)WRITE(6,*)'Call STEQ routines',MYPE,DST,DEND
	  CALL TUNE(IONE,'STEQ')
          DO ID=1,NUM_IONS-1
            LOC_ID=ID
	    IF(ATM(ID)%XzV_PRES)THEN
	      TMP_STRING=TRIM(ION_ID(ID))//'_COL_DATA'
              CALL STEQ_MULTI_MPI_V1(ED,T,
	1         ATM(ID)%XzV,            ATM(ID)%XzVLTE,         ATM(ID)%dlnXzVLTE_dlnT,
	1         AVE_ENERGY(ATM(ID)%EQXzV),
	1         ATM(ID)%NXzV,           ATM(ID)%DXzV,           ATM(ID)%XzV_F,
	1         ATM(ID)%XzVLTE_F_ON_S,  ATM(ID)%W_XzV_F,        ATM(ID)%AXzV_F,
	1         ATM(ID)%EDGEXzV_F,      ATM(ID)%GXzV_F,         ATM(ID)%XzVLEVNAME_F,
	1         ATM(ID)%NXzV_F,         ATM(ID)%F_TO_S_XzV,
	1         POP_SPECIES(1,SPECIES_LNK(ID)), ATM(ID+1)%XzV_PRES, ATM(ID)%ZXzV,
	1         LOC_ID,TMP_STRING,OMEGA_GEN_V3,
	1         ATM(ID)%EQXzV,NUM_BNDS,DST,DEND,ND,NION,NT,
	1         COMPUTE_BA,FIXED_T,LST_ITERATION)
!
	  IF(DST .EQ. 31)THEN
	    WRITE(125,*)ID
            WRITE(125,'(A,10ES14.4)')'G ',(ATM(12)%XzV(I,DST),I=4,11)
            WRITE(125,'(A,10ES14.4)')'G ',(ATM(12)%XzV(I,DST+1),I=4,11)
	    FLUSH(UNIT=125)
	  END IF
!
! Handle states which can partially autoionize.
!
	      TMP_STRING=TRIM(ION_ID(ID))//'_AUTO_DATA'
              CALL STEQ_AUTO_MPI_V1(ED,T,
	1         ATM(ID)%XzV,        ATM(ID)%NXzV,         ATM(ID)%DXzV,
	1         ATM(ID)%XzV_F,      ATM(ID)%XzVLTE_F,     ATM(ID)%EDGEXzV_F,
	1         ATM(ID)%GXzV_F,     ATM(ID)%XzVLEVNAME_F, ATM(ID)%NXzV_F,
	1         ATM(ID)%F_TO_S_XzV, LOC_ID,
	1         ATM(ID)%DIERECOM,   ATM(ID)%DIECOOL,
	1         TMP_STRING,NUM_BNDS,DST,DEND,ND,COMPUTE_BA)
	    END IF
	  IF(DST .EQ. 31)THEN
	    WRITE(125,*)ID
            WRITE(125,'(A,10ES14.4)')'H ',(ATM(12)%XzV(I,DST),I=4,11)
            WRITE(125,'(A,10ES14.4)')'H ',(ATM(12)%XzV(I,DST+1),I=4,11)
	    FLUSH(UNIT=125)
	  END IF
	  END DO
	  CALL TUNE(ITWO,'STEQ')
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
! Update charge equation. No longer done in STEQHEII
!
          CALL STEQNE_V4(ED,NT,DIAG_INDX,ND,COMPUTE_BA,DST,DEND)
!
	  IF(DST .EQ. 31)THEN
            WRITE(125,'(A,10ES14.4)')'BC ',(ATM(12)%XzV(I,DST),I=4,11)
            WRITE(125,'(A,10ES14.4)')'BC ',(ATM(12)%XzV(I,DST+1),I=4,11)
	    FLUSH(UNIT=125)
	  END IF
!
	IF(VERBOSE_OUTPUT .AND. LST_ITERATION)THEN
	  CALL GET_LU(LU_T_EHB,'Electron energy balance')
	  OPEN(UNIT=LU_T_EHB,FILE='CHECK_EHB_BALANCE',STATUS='UNKNOWN',ACTION='WRITE')
	  CALL WR2D_MPI_V1(TA,IONE,DST,DEND,ND,'After collison terms','&',L_TRUE,LU_T_EHB)
	END IF
!
	IF(TREAT_NON_THERMAL_ELECTRONS)THEN
	  IF(MYPE .EQ. 0)WRITE(6,*)'Beginning calculation of non-thermal electron spectrum: ED next'
	  CALL TUNE(IONE,'NON_THERM')
	    CALL ELECTRON_NON_THERM_SPEC_MPI_V1(ND)
	  CALL TUNE(ITWO,'NON_THERM')
	  CALL TUNE(IONE,'SE_NON_THERM')
	    CALL SE_BA_NON_THERM_MPI_V1(dE_RAD_DECAY,dE_SHOCK_POWER,COMPUTE_BA,NT,ND,DEC_NRG_SCL_FAC)
	  CALL TUNE(ITWO,'SE_NON_THERM')
	ELSE
	  dE_RAD_DECAY=0.0_LDP; dE_SHOCK_POWER=0.0_LDP
	END IF
!
!	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	  CALL WR_ASCI_STEQ_MPI_V1(NION,DST,DEND,ND,'STEQ ARRAY- Collisional Terms',LU_DR)
!	END IF
! 
!
! Compute the collisional cooling terms for digestion.
!
	  IF(MYPE .EQ. 0)WRITE(6,*)'Calling COLCOOL routines',MYPE,DST,DEND,NUM_IONS,NION
	  CALL TUNE(IONE,'COL_DIGEST')
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      TMP_STRING=TRIM(ION_ID(ID))//'_COL_DATA'
	      CALL COLCOOL_SL_MPI_V1(
	1        ATM(ID)%CPRXzV, ATM(ID)%CRRXzV,  ATM(ID)%COOLXzV,
	1        ATM(ID)%XzV,    ATM(ID)%XzVLTE,  ATM(ID)%dlnXzVLTE_dlnT,
	1        ATM(ID)%NXzV,   ATM(ID)%XzV_F,   ATM(ID)%XzVLTE_F_ON_S,
	1        ATM(ID)%AXzV_F, ATM(ID)%W_XzV_F, ATM(ID)%EDGEXzV_F,
	1        ATM(ID)%GXzV_F, ATM(ID)%XzVLEVNAME_F,
	1        ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ATM(ID)%ZXzV,
	1        ID,TMP_STRING,OMEGA_GEN_V3,ED,T,DST,DEND,ND)
	    END IF
	  END DO
	  CALL TUNE(ITWO,'COL_DIGEST')
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
! 
!
	EDDINGTON=EDD_LINECONT
	IF(ACCURATE .OR. EDDINGTON)THEN
	  IF(COMPUTE_EDDFAC)THEN
	    IF(MYPE .EQ. 0)WRITE(LU_EDD,REC=DIE_CONT_REC)ACCESS_F
	  ELSE
	    READ(LU_EDD,REC=DIE_CONT_REC)ACCESS_F
	  END IF
	END IF
!
	DO ML=1,NDIETOT
	  SECTION='DIELECTRONIC'
	  DO ID=1,NUM_IONS-1
	    IF(SPECDIE(ML) .EQ. ION_ID(ID))THEN
	      MNL_F=LEVDIE(ML)			!Level in full atom
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)		!Level in small_atom atom
	      NL=MNL+ATM(ID)%EQXzV-1			!Level in pops
	      GLOW=ATM(ID)%GXzV_F(LEVDIE(ML))
	      GION=ATM(ID)%GIONXzV_F
	      EQBAL=ATM(ID)%EQXzVBAL
	      EQION=ATM(ID+1)%EQXzV
	      EQSPEC=EQ_SPECIES(SPECIES_LNK(ID))
	      FL=ATM(ID)%EDGEXzV_F(LEVDIE(ML))-EDGEDIE(ML)
	      DO K=1,ND
		LOW_OCC_PROB(K)=ATM(ID)%W_XzV_F(MNL_F,K)		!Occupation prob.
		L_STAR_RATIO(K,1)=ATM(ID)%XzVLTE_F(MNL_F,K)/ATM(ID)%XzVLTE(MNL,K)
		dL_RAT_dT(K,1)=L_STAR_RATIO(K,1)*
	1         (-1.5_LDP-HDKT*ATM(ID)%EDGEXzV_F(MNL_F)/T(K)
	1               -ATM(ID)%dlnXzVLTE_dlnT(MNL,K))/T(K)
	      END DO
	      ID_SAV=ID
	      EXIT
	    END IF
	  END DO
!
	  COMPUTE_NEW_CROSS=.TRUE.
	  CONT_FREQ=FL
!
! Determine which method will be used to compute continuum intensity.
! Present form is temporary measure for consistency with SAO.
!
	  IF(ALL_FREQ)THEN
	    THIS_FREQ_EXT=.TRUE.
	  ELSE
	    THIS_FREQ_EXT=.FALSE.
	  END IF
!
	  GLDGU(1)=GLOW/GUPDIE(ML)
	  EINA(1)=EINADIE(ML)
	  OSCIL(1)=EINA(1)*EMLIN/( GLDGU(1)*OPLIN*TWOHCSQ*(FL**2) )
!
	  DO I=1,ND
	    DION(I)=POPS(EQION,I)
	  END DO
!
! Compute the LTE population for upper autoionizing state. We multiply
! NUST by Occupation probability of the lower state to correct for the
! fact that some transitions effectively keep the atom ionized.
!
	  CALL LTEPOP(NUST,ED,DION,GUPDIE(ML),EDGEDIE(ML),T,GION,IONE,ND)
	  DO I=1,ND
	    NUST(I)=NUST(I)*LOW_OCC_PROB(I)
	  END DO
!
! Compute line opacity and emissivity.
!
	  T1=OSCIL(1)*OPLIN
	  T2=FL*EINA(1)*EMLIN
	  DO I=1,ND
	    CHIL(I)=T1*( POPS(NL,I)*L_STAR_RATIO(I,1)-GLDGU(1)*NUST(I) )
	    ETAL(I)=T2*NUST(I)
	  END DO
! 
!
	  IF(IMPURITY_CODE)THEN
!
! Obtain previously compute continuum opacities, and mean intensities.
!
	    INCLUDE 'GET_J_CHI.INC'
	  ELSE
!
! Compute continuum opacity and emissivity at the line frequency.
!
	     CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
!
! Solve for the continuous radiation field.
!
!	    INCLUDE 'COMP_JCONT_V4.INC'
            CALL COMP_J_BLANK_MPI_V1(SECTION,EDDINGTON,FL,FREQ_INDX,FIRST_FREQ,LST_ITERATION,
	1                              MAXCH,LUER,LU_ES,LU_JCOMP,LU_EDD,ACCESS_F,
	1                              ND,NC,NP,NCF,NDEXT,NCEXT,NPEXT)
	  END IF
!
! SOURCE is used by SOBJBAR and in VARCONT. Note that SOURCE is corrupted
! (i.e. set to line source function) in CMFJBAR. Also note that SOURCEEXT
! has previously been computed (not for impurity species).
!
	  DO I=1,ND
	    SOURCE(I)=ZETA(I)+RJ(I)*THETA(I)
	  END DO
! 
!
! We define VC(I)=N* d(N* Z)/dN* and VB(I)=d(N* Z)/dNL.
!
	  T3=1.0_LDP/( TWOHCSQ*(FL**3) )
	  T2=T3/GLDGU(1)
	  DO I=DST,DEND
	    ZNET(I)=1.0_LDP-RJ(I)*CHIL(I)/ETAL(I)
	    VC(I)=(1.0_LDP+RJ(I)*T3)*NUST(I)
	    VB(I)=-RJ(I)*T2
	  END DO
!
! Evaluate contribution to statistical equilibrium equation, and
! and increment variation matrices.
! To convert cooling to ergs/cm**3/s, we use EDGEDIE(ML) and not FL for the
! electron cooling component since this is the energy spent in exciting the
! autoionizing state. Can show from statistical equilibrium equations.
! Note that EDGEDIE(ML) is negative.
!
	  T2=-1.256637E-09_LDP*EINA(1)*EMLIN*EDGEDIE(ML)
	  T3=EINA(1)*FL*EMLIN
!
	  DO K=DST,DEND
	    SE(ID)%STEQ(NL,K)=SE(ID)%STEQ(NL,K)+EINA(1)*NUST(K)*ZNET(K)
	    STEQ_T(K)=STEQ_T(K)-T3*NUST(K)*ZNET(K)
	    DIELUM(K)=DIELUM(K)+ETAL(K)*ZNET(K)
	    ATM(ID)%DIERECOM(K)=ATM(ID)%DIERECOM(K) + EINA(1)*NUST(K)*ZNET(K)
	    ATM(ID)%DIECOOL(K)=ATM(ID)%DIECOOL(K)   + T2*NUST(K)*ZNET(K)
	  END DO
!
! Update ionization balance equations if required.
!
	  MNUP=ATM(ID)%NXzV+1
	  DO K=DST,DEND
	    SE(ID)%STEQ(MNUP,K)=SE(ID)%STEQ(MNUP,K)-EINA(1)*NUST(K)*ZNET(K)
	  END DO
! 
!
! Update BA matrix - this section must be done even if we are
! performing a LAMBDA iteration. We define a LAMBDA iteration by assuming
! that the variation of J is zero.
!
	  IF(COMPUTE_BA)THEN
!
	    MNUP=ATM(ID)%NXzV+1
	    NUP=ATM(ID)%EQXzV+ATM(ID)%NXzV
	    NIV=SE(ID)%N_IV
	    DO K=DST,DEND
	      L=GET_DIAG(K)
	      SE(ID)%BA(MNL,MNL,L,K) =SE(ID)%BA(MNL,MNL,L,K)  +EINA(1)*VB(K)
	      SE(ID)%BA(MNL,MNUP,L,K)=SE(ID)%BA(MNL,MNUP,L,K) +EINA(1)*VC(K)/DION(K)
	      SE(ID)%BA(MNL,NIV-1,L,K)=SE(ID)%BA(MNL,NIV-1,L,K) +EINA(1)*VC(K)/ED(K)
	      SE(ID)%BA(MNL,NIV,L,K)  =SE(ID)%BA(MNL,NIV,L,K)   -
	1          EINA(1)*VC(K)*(1.5_LDP+HDKT*EDGEDIE(ML)/T(K))/T(K) +
	1          EINA(1)*VB(K)*POPS(NL,K)*dL_RAT_dT(K,1)/L_STAR_RATIO(K,1)
!
	      BA_T(NL,L,K) =BA_T(NL,L,K)-T3*VB(K)
	      BA_T(NUP,L,K) =BA_T(NUP,L,K)-T3*VC(K)/DION(K)
	      BA_T(NT-1,L,K)=BA_T(NT-1,L,K)-T3*VC(K)/ED(K)
	      BA_T(NT,L,K)  =BA_T(NT,L,K)+
	1                   T3*VC(K)*(1.5_LDP+HDKT*EDGEDIE(ML)/T(K))/T(K)
	    END DO
!
! Update ionization equation.
!
	    DO K=DST,DEND
	      L=GET_DIAG(K)
	      SE(ID)%BA(MNUP,MNL,L,K)  =SE(ID)%BA(MNUP,MNL,L,K)   - EINA(1)*VB(K)
	      SE(ID)%BA(MNUP,MNUP,L,K) =SE(ID)%BA(MNUP,MNUP,L,K)  - EINA(1)*VC(K)/DION(K)
	      SE(ID)%BA(MNUP,NIV-1,L,K)=SE(ID)%BA(MNUP,NIV-1,L,K) - EINA(1)*VC(K)/ED(K)
	      SE(ID)%BA(MNUP,NIV,L,K    )=SE(ID)%BA(MNUP,NIV,L,K) +
	1         EINA(1)*VC(K)*(1.5_LDP+HDKT*EDGEDIE(ML)/T(K))/T(K) -
	1         EINA(1)*VB(K)*POPS(NL,K)*dL_RAT_dT(K,1)/L_STAR_RATIO(K,1)
	    END DO
	  END IF
! 
!
! Allow for the variation of the continuous radiation field.
!
	  IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION .AND.
	1                      .NOT. IMPURITY_CODE)THEN
!
	    DO I=DST,DEND
	      BETAC(I)=CHIL(I)/ETAL(I)
	    END DO
!
!	    INCLUDE 'VARCONT.INC'
            CALL DO_VAR_CONT_MPI_V1(POPS,SECTION,EDDINGTON,
	1                  FL,CONT_FREQ,FREQ_INDX,FIRST_FREQ,TX_OFFSET,
	1                  ND,NC,NP,NUM_BNDS,DIAG_INDX,NT,NM,
	1                  NDEXT,NCEXT,NPEXT,MAX_SIM,NM_KI)
!
! Increment the large simultaneous perturbation matrix due to a variation
! in the continuum. This is only incremented if ZNET is not within 1% of
! 1.0, which indicates that the continuum term is important.
!
	    CALL TUNE(IONE,'DIECONTBA')
	    DO L=DST,DEND   	  		  	!S.E. equation depth
	      T1=EINA(1)*NUST(L)*BETAC(L)
	      T2=ETAL(L)*BETAC(L)
	      DO K=BNDST(L),BNDEND(L)	  		!Variable depth.
	        LS=BND_TO_FULL(K,L)
   	        DO J=1,SE(ID)%N_IV	 	   		!Variable
	          JJ=SE(ID)%LNK_TO_F(J)
	          SE(ID)%BA(MNL,J,K,L)=SE(ID)%BA(MNL,J,K,L) - T1*VJ(JJ,K,L)
	          BA_T(JJ,K,L)=BA_T(JJ,K,L) + T2*VJ(JJ,K,L)
	        END DO
	      END DO
	    END DO
!
	    DO L=DST,DEND	  	  		  	!S.E. equation depth
	      T1=EINA(1)*NUST(L)*BETAC(L)
	      DO K=BNDST(L),BNDEND(L)	  		!Variable depth.
	        LS=BND_TO_FULL(K,L)
   	        DO J=1,NT	 	   		!Variable
	          JJ=SE(ID)%LNK_TO_F(J)
	          SE(ID)%BA(MNUP,J,K,L)=SE(ID)%BA(MNUP,J,K,L) + T1*VJ(JJ,K,L)
	        END DO
	      END DO
	    END DO
	    CALL TUNE(ITWO,'DIECONTBA')
	  END IF			!BA Matrix computed (compute_ba).
!
	  IF(LST_ITERATION .AND. WRITE_RATES)THEN
!
! Estimate the line EW using a Modified Sobolev approximation.
!
! We use TA as a temporary vector which indicates the origin
! of the line emission. Not required in this code as used only
! for display purposes. Variable after THK_CONT is true as we
! want to assume the line opacity is zero --- since dielectronic
! transition.
!
	    CALL SOBEW(SOURCE,CHI,CHI_SCAT,CHIL,ETAL,
	1              V,SIGMA,R,P,AQW,HQW,TA,EW,CONT_INT,
	1              FL,DIF,DBB,IC,THK_CONT,L_TRUE,NC,NP,ND,METHOD)
!
	    T1=LAMVACAIR(FL)			!Wavelength(Angstroms)
	    CALL EW_FORMAT(EW_STRING,DIENAME(ML),T1,CONT_INT,EW,L_TRUE)
	    L=ICHRLEN(EW_STRING)
	    WRITE(LU_NET,40002)EW_STRING(1:L)
	    WRITE(LU_DR,40002)EW_STRING(1:L)
	    WRITE(LU_EW,40005)EW_STRING(1:L)
	    WRITE(LU_HT,40002)EW_STRING(1:L)
	    WRITE(LU_NET,40009)GUPDIE(ML),FL,EINA(1)
	    WRITE(LU_DR,40009)GUPDIE(ML),FL,EINA(1)
	    WRITE(LU_HT,40009)GUPDIE(ML),FL,EINA(1)
	    WRITE(LU_NET,40003)( ZNET(I),I=1,ND )
	    WRITE(LU_DR,40003)(  ( ZNET(I)*NUST(I)*EINA(1) ),I=1,ND  )
	    WRITE(LU_HT,40003)(  ( ZNET(I)*ETAL(I) ),I=1,ND  )
40002	    FORMAT(//,A)		!From SOB section
40003	    FORMAT(3X,1P,5E16.5)
40005	    FORMAT(A)
40009	    FORMAT(1X,F5.0,2X,1P,2E12.4)
	    CLOSE(UNIT=LU_NET)
	    CALL GEN_ASCI_OPEN(LU_NET,'NETRATE','OLD','APPEND',' ',IZERO,IOS)
	  END IF
	END DO		!End of Dielectronic section [do ML]
	IF(MYPE .EQ. 0)WRITE(6,*)'Done dielectronic section',MYPE,DST,DEND
! 
!***************************************************************************
!***************************************************************************
!
!                         CONTINUUM LOOP
!
!***************************************************************************
!***************************************************************************
!
	EDDINGTON=EDD_CONT
	IF(ACCURATE .OR. EDDINGTON)THEN
	  IF(COMPUTE_EDDFAC .AND. MYPE .EQ. 0)THEN
	    WRITE(LU_EDD,REC=EDD_CONT_REC)ACCESS_F,NCF,NDEXT
	  ELSE 
	    READ(LU_EDD,REC=EDD_CONT_REC)ACCESS_F
	  END IF
	END IF
	DO I=0,NTHREAD-1
	  IF(MYPE .EQ. 0)THEN
	    FLUSH(UNIT=6); WRITE(6,'(A)')' '
	    WRITE(6,'(4A12)')'MYPE','EDDINGTON','ACCESS_F','COMP_F'
	    WRITE(6,'(9X,I3,11X,L1,10X,I2,11X,L1)')MYPE,EDDINGTON,ACCESS_F,COMPUTE_EDDFAC	
	    FLUSH(UNIT=6)
	  ELSE IF(MYPE .EQ. 1)THEN
	    WRITE(6,'(9X,I3,11X,L1,10X,I2,11X,L1)')MYPE,EDDINGTON,ACCESS_F,COMPUTE_EDDFAC	
	    FLUSH(UNIT=6)
	  END IF
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	END DO
!
! Decide whether to use an file with old J values to provide an initial
! estimate of J with incoherent electron scattering. The options
!
!             .NOT. RD_COHERENT_ES .AND. COHERENT_ES
!
! indicate that no ES_J_CONV file exists already.
! Note that TEXT and NDEXT contain T and ND when ACCURATE is FALSE.
!
	IF(USE_OLDJ_FOR_ES .AND. .NOT. RD_COHERENT_ES .AND. COHERENT_ES)THEN
	  COHERENT_ES=RD_COHERENT_ES
	  I=SIZE(VJ)
	  CALL COMP_J_CONV_V2(VJ,I,NU,TEXT,NDEXT,NCF,LUIN,'OLD_J_FILE',
	1           EDD_CONT_REC,L_FALSE,L_TRUE,LU_ES,'ES_J_CONV')
!
! Now open the file so it can be read in the CONTINUUM loop (read in COMP_JCONT).
!
	  OPEN(UNIT=LU_ES,FILE='ES_J_CONV',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening ES_J_CONV - will compute new J'
	      COHERENT_ES=.TRUE.
	    END IF
	END IF
!
! We ensure that LAST_LINE points to the first LINE that is going to
! be handled in the BLANKETING portion of the code.
!
	LAST_LINE=0			!Updated as each line is done
	DO WHILE(LAST_LINE .LT. N_LINE_FREQ .AND.
	1             VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	        LAST_LINE=LAST_LINE+1
	END DO
	DO SIM_INDX=1,MAX_SIM
	  LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	END DO
!
! Ensure none of the storage location for the variation of J with CHIL etc
! are being pointed at.
!
	DO SIM_INDX=1,MAX_SIM
	  LOW_POINTER(SIM_INDX)=0
	  UP_POINTER(SIM_INDX)=0
	END DO
!
	DO I=1,NM
	  VAR_IN_USE_CNT(I)=0
	  VAR_LEV_ID(I)=0
	  IMP_TRANS_VEC(I)=.FALSE.
	END DO
!
	NUM_OF_WEAK_LINES=0.0_LDP
	CONT_FREQ=0.0_LDP
!
	ETA(1:ND)=1.0E-10_LDP; CHI(1:ND)=1.0E-10_LDP; ESEC(1:ND)=1.0E-11_LDP
	ETA_CONT(1:ND)=1.0E-10_LDP; CHI_CONT(1:ND)=1.0E-10_LDP
	ZETA(1:ND)=1.0E-10_LDP; THETA(1:ND)=0.1_LDP; CHI_SCAT(1:ND)=1.0E-10_LDP
	SOURCE(1:ND)=1.0E-10_LDP; ETA_NOSCAT(1:ND)=1.0E-10_LDP; CHI_NOSCAT(1:ND)=1.0E-10_LDP
!
! Enter loop for each continuum frequency.
!
	SUM_BA=0.0_LDP
	FL=NU(1)
	FLUSH(LUER)
	CALL TUNE(IONE,'MLCF')
	CALL TUNE(IONE,'10000')
	IF(MYPE .EQ. 0)WRITE(6,*)'Starting 10000 loop', WRITE_JH; FLUSH(6)
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        WRITE(STRING,*)ML; STRING='ML='//ADJUSTL(TRIM(STRING))
        CALL WRITV(T,ND,TRIM(STRING),852+MYPE); FLUSH(852+MYPE)
	DO 10000 ML=1,NCF
	  FREQ_INDX=ML
	  FL=NU(ML)
	  IF(ML .EQ. 1)THEN
	    FIRST_FREQ=.TRUE.
	  ELSE
	    FIRST_FREQ=.FALSE.
	  END IF
	  SECTION='CONTINUUM'
!
	  IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	    COMPUTE_NEW_CROSS=.TRUE.
	    CONT_FREQ=NU_EVAL_CONT(ML)
	  ELSE
	    COMPUTE_NEW_CROSS=.FALSE.
	  END IF
	  FINAL_CONSTANT_CROSS=.TRUE.
	  IF(ML .EQ. NCF)THEN
	    FINAL_CONSTANT_CROSS=.TRUE.
	  ELSE
	    IF(NU_EVAL_CONT(ML+1) .EQ. CONT_FREQ)
	1                      FINAL_CONSTANT_CROSS=.FALSE.
	  END IF
!
! Compute quadrature weights for statistical equilibrium equations.
! TA is used as a work vector (dim ND)
!
	  CALL TUNE(IONE,'QUAD')
	  TMP_LOGICAL=FIXED_T .OR. LAMBDA_ITERATION
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES .AND. COMPUTE_NEW_CROSS)THEN
	       CALL QUAD_MULTI_MPI_V1(ATM(ID)%WSXzV, ATM(ID)%dWSXzVdT,
	1             ATM(ID)%WCRXzV,ATM(ID)%dWCRXzVdT,
	1             ATM(ID)%XzVLTE, ATM(ID)%dlnXzVLTE_dlnT, ATM(ID)%NXzV,
	1             ATM(ID)%XzVLTE_F_ON_S, ATM(ID)%EDGEXzV_F, ATM(ID)%NXzV_F,
	1             ATM(ID)%F_TO_S_XzV, CONT_FREQ,T,
	1             DST,DEND,ND,
	1             COMPUTE_BA,TMP_LOGICAL,LST_ITERATION,
	1             ION_ID(ID), ATM(ID)%ZXzV, ATM(ID)%N_XzV_PHOT, ID)
	    END IF
	  END DO
	  CALL TUNE(ITWO,'QUAD')
!
!	  ID=1; L=71
!	  DO K=1,ATM(ID)%NXzV
!	   WRITE(314,'(ES16.8,I6,4ES14.4)')FL,K,ATM(ID)%WSXzV(K,L,1),ATM(ID)%dWSXzVdT(K,L,1),
!	1             ATM(ID)%WCRXzV(K,L,1),ATM(ID)%dWCRXzVdT(K,L,1)
!	  END DO
!
! 
!
	IF(XRAYS .AND. COMPUTE_NEW_CROSS)THEN
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	      T1=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV		!Number of electrons
	      CALL QUAD_X_GEN_MPI_V1(AT_NO(SPECIES_LNK(ID)),T1,
	1           ATM(ID)%WSE_X_XzV,      ATM(ID)%WCR_X_XzV, CONT_FREQ,
	1           ATM(ID)%XzVLTE,         ATM(ID)%NXzV,
	1           ATM(ID)%XzVLTE_F_ON_S,  ATM(ID)%EDGEXzV_F,
	1           ATM(ID)%F_TO_S_XzV,     ATM(ID)%NXzV_F,
	1           ATM(ID+1)%EDGEXzV_F,    ATM(ID+1)%NXzV_F, DST, DEND, ND)
	    END IF
	  END DO
	END IF
!
! 
!
! Include lines
!
	CALL TUNE(IONE,'SET_LINE_OPAC')
        CALL SET_LINE_OPAC(POPS,NU,ML,LAST_LINE,N_LINE_FREQ,
	1         LST_DEPTH_ONLY,LUER,ND,NT,NCF,MAX_SIM)
	CALL TUNE(ITWO,'SET_LINE_OPAC')
!
	CALL INIT_LINE_OPAC_VAR_V2(LAST_LINE,LUER,ND,TX_OFFSET,MAX_SIM,NM)
!
! Determine which method will be used to compute continuum intensity.
!
	  IF(ACCURATE .AND. ALL_FREQ)THEN
	    THIS_FREQ_EXT=.TRUE.
	  ELSE IF( ACCURATE .AND. FL .GT. ACC_FREQ_END )THEN
	    THIS_FREQ_EXT=.TRUE.
	  ELSE
	    THIS_FREQ_EXT=.FALSE.
	  END IF
!
	  IF(IMPURITY_CODE)THEN
!
! Obtain previously compute continuum opacities, and mean intensities.
!
	    INCLUDE 'GET_J_CHI.INC'
	  ELSE
!
! Compute opacity and emissivity.
!
	    IF(USE_FIXED_J)THEN
	      IF(COMPUTE_NEW_CROSS)THEN
                T1=-HDKT*CONT_FREQ
                EMHNUKT_CONT(1:ND)=EXP(T1/T(1:ND))
                T1=-HDKT*FL
                EMHNUKT(1:ND)=EXP(T1/T(1:ND))
	        CHI=0.0_LDP; ETA=0.0_LDP
	        CHI(DST:DEND)=2.0*ED(DST:DEND)*6.65D-15
	        CHI_SCAT(DST:DEND)=ED(DST:DEND)*6.65D-15
	      END IF
	    ELSE
	      CALL TUNE(IONE,'C_OPAC')
	      CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONlY)
	      CALL TUNE(ITWO,'C_OPAC')
	    END IF
!
! Since resonance zones included, we must add the line opacity and
! emissivity to the raw continuum values. We first save the pure continuum
! opacity and emissivity. These are used in carrying the variation of J from
! one frequency to the next.
!
	    CALL TUNE(IONE,'RS_ZONE')
	    DO SIM_INDX=1,MAX_SIM
	      IF(RESONANCE_ZONE(SIM_INDX))THEN
	        DO I=DST,DEND
	          CHI(I)=CHI(I)+CHIL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	          ETA(I)=ETA(I)+ETAL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	        END DO
	      END IF
	    END DO
	    CALL TUNE(ITWO,'RS_ZONE')
!
! CHECK for negative line opacities. NEG_OPAC_FAC is the factor we
! multiply the line opacities by so that the total opacity is positive.
! We do not distinguish between lines. Two different options are possible.
! The first, 'ESEC_CHK' was in use for years, and is probably the preferred
! option. The second, 'SRCE_CHK', was introduced to overcome problems in
! O Star models. In particular, in some models a negative optical depth
! could occur on some iteartions at depths where ABS(TAUL) was still very
! large (primraily in far IT transitions [e.g. H(9-8)].
!
	    AT_LEAST_ONE_NEG_OPAC=.FALSE.
	    NEG_OPACITY(1:ND)=.FALSE.
	    NEG_OPAC_FAC(1:ND)=1.0_LDP
	    IF(NEG_OPAC_OPTION .EQ. 'SRCE_CHK')THEN
	      DO I=DST,DEND
	        IF(CHI(I) .LT. CHI_CONT(I) .AND.
	1            CHI(I) .LT. 0.1_LDP*ETA(I)*CHI_NOSCAT(I)/ETA_CONT(I) )THEN
	          CHI(I)=0.1_LDP*ETA(I)*CHI_NOSCAT(I)/ETA_CONT(I)
	          NEG_OPACITY(I)=.TRUE.
	          NEG_OPAC_FAC(I)=0.0_LDP
	          AT_LEAST_ONE_NEG_OPAC=.TRUE.
	        ELSE IF(CHI(I) .LT. 0.1_LDP*CHI_SCAT(I))THEN
	          CHI(I)=0.1_LDP*CHI_SCAT(I)
	          NEG_OPACITY(I)=.TRUE.
	          NEG_OPAC_FAC(I)=0.0_LDP
	          AT_LEAST_ONE_NEG_OPAC=.TRUE.
	        END IF
	      END DO
	    ELSE IF(NEG_OPAC_OPTION .EQ. 'ESEC_CHK')THEN
	      DO I=DST,DEND
	        IF(CHI(I) .LT. 0.1_LDP*CHI_SCAT(I))THEN
	          T1=CHI(I)
	          CHI(I)=0.1_LDP*CHI_SCAT(I)
	          NEG_OPACITY(I)=.TRUE.
!	          NEG_OPAC_FAC(I)=(CHI(I)-CHI_CONT(I))/(T1-CHI_CONT(I))
	          NEG_OPAC_FAC(I)=0.0_LDP
	          AT_LEAST_ONE_NEG_OPAC=.TRUE.
	        END IF
	      END DO
	    END IF
!
            CALL MPI_ALLREDUCE(CHI,TA,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR); CHI(1:ND)=TA(1:ND)
            CALL MPI_ALLREDUCE(ETA,TA,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR); ETA(1:ND)=TA(1:ND)
!
	    IF(LST_ITERATION .AND. 3 .EQ.4)THEN
	      CALL MPI_REDUCE(MPI_IN_PLACE,AT_LEAST_ONE_NEG_OPAC,IONE,MPI_LOGICAL,MPI_LAND,IZERO,MPI_COMM_WORLD,IERR)
	      IF(LST_ITERATION .AND. AT_LEAST_ONE_NEG_OPAC)THEN
	        CALL MPI_REDUCE(MPI_IN_PLACE,NEG_OPACITY,ND,MPI_LOGICAL,MPI_LAND,IZERO,MPI_COMM_WORLD,IERR)
	      END IF
	    END IF
!I
	    IF(LST_ITERATION .AND. AT_LEAST_ONE_NEG_OPAC .AND. MYPE .EQ. 0)THEN
	      WRITE(LU_NEG,'(A,1P,E14.6)')' Neg opacity for transition for frequency ',FL
	      DO SIM_INDX=1,MAX_SIM
	        IF(RESONANCE_ZONE(SIM_INDX))THEN
	          WRITE(LU_NEG,'(1X,A)')TRANS_NAME_SIM(SIM_INDX)
	        END IF
	      END DO
	      J=0
	      K=0
	      DO I=DST,DEND
	       IF(NEG_OPACITY(I) .AND. K .EQ. 0)K=I
	       IF(NEG_OPACITY(I))J=I
	      END DO
	      WRITE(LU_NEG,'(A,2X,I3,5X,A,2XI3)')' 1st depth',K,'Last depth',J
	    END IF
!
	    DO I=DST,DEND
	      ZETA(I)=ETA(I)/CHI(I)
	      THETA(I)=CHI_SCAT(I)/CHI(I)
	    END DO
!
	    IF(LST_ITERATION .AND. ML .NE. NCF .AND. MYPE .EQ. 0)THEN
	      DO I=1,N_TAU_EDGE
	        IF(NU(ML) .GE. TAU_EDGE(I) .AND.
	1                       NU(ML+1) .LT. TAU_EDGE(I))THEN
	          T1=LOG(CHI_CONT(5)/CHI_CONT(1))/LOG(R(1)/R(5))
	          IF(I .EQ. 1)WRITE(LUER,'(A)')' '
	          WRITE(LUER,'(A,1P,E11.4,A,E10.3)')' Tau(Nu=',NU(ML),
	1            ') at outer boundary is:',CHI_CONT(1)*R(1)/MAX(T1-1.0D0,1.0D0)
	          IF(I .EQ. N_TAU_EDGE)WRITE(LUER,'(A)')' '
	        END IF
	      END DO
	    END IF
!
! Compute continuum intensity.
!
	    CALL TUNE(IONE,'COMP_J')
!	    INCLUDE 'COMP_JCONT_V4.INC'	
	    CALL COMP_J_BLANK_MPI_V1(SECTION,EDDINGTON,FL,FREQ_INDX,FIRST_FREQ,LST_ITERATION,
	1                              MAXCH,LUER,LU_ES,LU_JCOMP,LU_EDD,ACCESS_F,
	1                              ND,NC,NP,NCF,NDEXT,NCEXT,NPEXT)
	    CALL TUNE(ITWO,'COMP_J')
	  END IF
! 
!
! Increment the RADIATIVE EQUILIBRIUM equation due to radiation field at this
! frequency. The correction for non-coherent electron scattering allows
! for the fact that the electron-scattering emissivity is ESEC*RJ_ES, not
! ESEC*RJ.
!
! ETA_CONT includes all emissivity sources, including X-ray emission produced
! by mechanical or magnetic energy deposition. This should not be included
! in the radiatively equilibrium equation, hence we subtract out the
! emissivity due to mechanical processes. NB: ETA_NOSCAT does not include
! mechanical term.
!
	  DO K=DST,DEND
	    STEQ_T(K)=STEQ_T(K)+ FQW(ML)*(CHI_NOSCAT(K)*RJ(K) - ETA_NOSCAT(K))
	    IF(.NOT. COHERENT_ES)STEQ_T(K)=STEQ_T(K)+FQW(ML)*ESEC(K)*(RJ(K)-RJ_ES(K))
	  END DO
!
	  CALL COMP_VAR_JREC_MPI_V1(JREC,dJRECdT,JPHOT,JREC_CR,dJREC_CRdT,JPHOT_CR,BPHOT_CR,
	1       RJ,EMHNUKT,T,NU(ML),FQW(ML),TWOHCSQ,HDKT,DST,DEND,ND,COMPUTE_NEW_CROSS)
!
! Increment the S.E. equations due to radiation field at this
! frequency.
!
! At the same time, we compute the quadrature weights associated with
! the intensity for each depth point and each equation. QFV must be zeroed
! before calling EVALSE_QWVJ. QFV is incremented - not set. Allows for
! bound-free processes to both the ground and excited states (necessary
! for CIII).
!
	IF(FINAL_CONSTANT_CROSS)THEN
	  DO ID=1,NION
	    IF(SE(ID)%XzV_PRES)THEN
	      SE(ID)%QFV_R(:,:)=0.0_LDP		!NT,ND
	      SE(ID)%QFV_P(:,:)=0.0_LDP
	      SE(ID)%QFV_R_EHB(:,:)=0.0_LDP		!NT,ND
	      SE(ID)%QFV_P_EHB(:,:)=0.0_LDP
	    END IF
	  END DO
	END IF
!
! Need to fix parallization issue with update of STEQ_T_ED inside  EVALSE_QWVJ_V8
!
	IF(FINAL_CONSTANT_CROSS)THEN
	  CALL TUNE(IONE,'EVALSE')
	  DO ID=1,NUM_IONS-1
	    ID_SAV=ID
	    IF(ATM(ID)%XzV_PRES)THEN
	      CALL EVALSE_QWVJ_MPI_V1(ID_SAV,
	1       ATM(ID)%WSXzV, ATM(ID)%WCRXzV,
	1       ATM(ID)%XzV, ATM(ID)%XzVLTE,
	1       ATM(ID)%NXzV, ATM(ID)%XzV_ION_LEV_ID,
	1       ATM(ID+1)%XzV, ATM(ID+1)%LOG_XzVLTE,
	1       ATM(ID+1)%NXzV,ATM(ID)%N_XzV_PHOT,
	1       JREC,JPHOT,JREC_CR,JPHOT_CR,NT,DST, DEND, ND)
	    END IF
	  END DO
	  CALL TUNE(ITWO,'EVALSE')
	END IF
!
	CALL TUNE(IONE,'LOWT')
	CALL EVALSE_LOWT_MPI_V1(RJ,NU(ML),FQW(ML),COMPUTE_BA,NT,ND,FIRST_FREQ)
	CALL TUNE(ITWO,'LOWT')
!
! 
! Note that ATM(ID+2)%EQXzV is the ion equation. Since Auger ionization,
! 2 electrons are ejected.
!
	IF(XRAYS .AND. FINAL_CONSTANT_CROSS)THEN
	  DO ID=1,NUM_IONS-1
	    ID_SAV=ID
	    IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	      CALL EVALSE_X_QWVJ_MPI_V1(ID_SAV,ATM(ID)%WSE_X_XzV,
	1          ATM(ID)%XzV,   ATM(ID)%XzVLTE,   ATM(ID)%NXzV,
	1          ATM(ID+1)%XzV, ATM(ID+1)%XzVLTE, ATM(ID+1)%NXzV, ATM(ID+2)%EQXzV,
	1          JREC,JPHOT,DST, DEND, ND,NION)
	    END IF
	  END DO
	END IF
!
! 
!
! Compute the recombination, photoionization and cooling rates.
!
	IF(LST_ITERATION .AND. FINAL_CONSTANT_CROSS)THEN ! .AND. .NOT.  USE_FIXED_J)THEN
	  CALL TUNE(IONE,'PRRRCOOL')
!
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      CALL PRRR_SL_MPI_V1(
	1          ATM(ID)%APRXzV,        ATM(ID)%ARRXzV,
	1          ATM(ID)%BFCRXzV,      ATM(ID)%FFXzV,
	1          ATM(ID)%WSXzV,        ATM(ID)%WCRXzV,
	1          ATM(ID)%XzV,          ATM(ID)%XzVLTE,
	1          ATM(ID)%NXzV,         ATM(ID)%ZXzV,
	1          ATM(ID+1)%XzV,        ATM(ID+1)%LOG_XzVLTE,
	1          ATM(ID+1)%NXzV,       ATM(ID)%XzV_ION_LEV_ID, ATM(ID)%N_XzV_PHOT,
	1          ED,T,JREC,JPHOT,JREC_CR,JPHOT_CR,BPHOT_CR,
	1          FL,CONT_FREQ,ZERO_REC_COOL_ARRAYS,DST,DEND,ND)
!	      IF(TRIM(ION_ID(ID)) .EQ. 'NI')THEN
!	        WRITE(184,'(3ES16.8,4ES14.4)')FL,CONT_FREQ,JREC(10),JPHOT(10),ATM(ID)%ARRXzV(1,10)
!	      END IF
	    END IF
	    IF(ID .EQ. 12 .AND. DST .EQ. 31)THEN
	      WRITE(120,'(10ES14.4)')(ATM(12)%APRXzV(I,DST),I=4,11),RJ(DST),JPHOT(DST)
	      WRITE(120,'(10ES14.4)')(ATM(12)%APRXzV(I,DST+1),I=4,11),RJ(DST+1),JPHOT(DST+1)
	      WRITE(120,'(10ES14.4)')(ATM(12)%WSXzV(I,DST,1),I=4,11)
	      WRITE(120,'(10ES14.4)')(ATM(12)%WSXzV(I,DST+1,1),I=4,11)
	      FLUSH(UNIT=120)
	    END IF
!
	    IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES .AND. XRAYS)THEN
	      CALL X_RRR_COOL_MPI_V1(ATM(ID)%X_RECOM,
	1            ATM(ID)%X_COOL, ATM(ID)%WSE_X_XzV, ATM(ID)%WCR_X_XzV,
	1            ATM(ID)%XzV,        ATM(ID)%LOG_XzVLTE,     ATM(ID)%NXzV,
	1            ATM(ID+1)%XzV_F,    ATM(ID+1)%LOG_XzVLTE_F, ATM(ID+1)%NXzV_F,
	1            JREC,JPHOT,JREC_CR,JPHOT_CR,
	1            ZERO_REC_COOL_ARRAYS,DST,DEND,ND,L_TRUE)
	      IF(ID .EQ. 12 .AND. DST .EQ. 31)THEN
	        WRITE(120,'(10ES14.4)')(ATM(12)%APRXzV(I,DST),I=4,11)
	        WRITE(120,'(10ES14.4)')(ATM(12)%APRXzV(I,DST+1),I=4,11)
	      END IF
	    END IF
	  END DO
	  ZERO_REC_COOL_ARRAYS=.FALSE.
!
	CALL TUNE(ITWO,'PRRRCOOL')
	END IF 			!Only evaluate if last iteration.
!
	IF(LST_ITERATION .AND. .NOT. USE_FIXED_J)THEN
	  CALL PRRR_LOWT_MPI_V1(RJ,NU(ML),FQW(ML),ND,FIRST_FREQ)
	END IF
!
! 
!
! Update line net rates, and the S.E. Eq. IFF we have finished a line
! transition.
!
	CALL TUNE(IONE,'JBAR_SIM')
	DO SIM_INDX=1,MAX_SIM
	  IF(RESONANCE_ZONE(SIM_INDX))THEN
	    DO I=DST,DEND
	      ZNET_SIM(I,SIM_INDX)=ZNET_SIM(I,SIM_INDX) + LINE_QW_SIM(I,SIM_INDX)*
	1          (1.0_LDP-RJ(I)*CHIL_MAT(I,SIM_INDX)/ETAL_MAT(I,SIM_INDX))
	      JBAR_SIM(I,SIM_INDX)=JBAR_SIM(I,SIM_INDX) + LINE_QW_SIM(I,SIM_INDX)*RJ(I)
	      LINE_QW_SUM(I,SIM_INDX)=LINE_QW_SUM(I,SIM_INDX) + LINE_QW_SIM(I,SIM_INDX)
	    END DO
	  END IF
	END DO
	CALL TUNE(ITWO,'JBAR_SIM')
!
! Update the S.E. Eq. IFF we have finished a line transition (i.e. are at
! the final point of the resonance zone.)
!
! NB: The line term in the RE equations is not needed since it is included
! directly with continuum integration.
!
	DO SIM_INDX=1,MAX_SIM
	  IF( END_RES_ZONE(SIM_INDX) )THEN
            T1=FL_SIM(SIM_INDX)*EMLIN
	    NL=SIM_NL(SIM_INDX)
	    NUP=SIM_NUP(SIM_INDX)
	    I=SIM_LINE_POINTER(SIM_INDX)
	    ID=VEC_ID(I)
	    MNL_F=VEC_MNL_F(I);     MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	    MNUP_F=VEC_MNUP_F(I);   MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	    SCL_FAC=1.0_LDP
	    IF(SCL_LINE_COOL_RATES)THEN
	      SCL_FAC=(AVE_ENERGY(NL)-AVE_ENERGY(NUP))/FL_SIM(SIM_INDX)
	      IF(ABS(SCL_FAC-1.0_LDP) .GT. SCL_LINE_HT_FAC)SCL_FAC=1.0_LDP
	    END IF
	    DO K=DST,DEND
	      T4=SCL_FAC
	      IF(POP_ATOM(K) .GE. SCL_LINE_DENSITY_LIMIT)T4=1.0_LDP
	      T2=EINA(SIM_INDX)*ATM(ID)%XzV_F(MNUP_F,K)*ZNET_SIM(K,SIM_INDX)
	      T3=T4*ETAL_MAT(K,SIM_INDX)*ZNET_SIM(K,SIM_INDX)
	      SE(ID)%STEQ(MNUP,K)=SE(ID)%STEQ(MNUP,K) - T2
	      SE(ID)%STEQ(MNL,K) =SE(ID)%STEQ(MNL,K) + T2
	      STEQ_T(K)=STEQ_T(K) - T3
	    END DO
	  END IF
	END DO
!
! 
!
! Allow for the variation of the continuous radiation field.
!
	  IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION .AND.
	1                     .NOT. IMPURITY_CODE)THEN
!
! Solve for the perturbations to J in terms of the perturbations
! to CHI and ETA.
!
	    CALL TUNE(IONE,'C_VARCONT')
!	      INCLUDE 'VARCONT.INC'
              CALL DO_VAR_CONT_MPI_V1(POPS,SECTION,EDDINGTON,
	1                    FL,CONT_FREQ,FREQ_INDX,FIRST_FREQ,TX_OFFSET,
	1                    ND,NC,NP,NUM_BNDS,DIAG_INDX,NT,NM,
	1                    NDEXT,NCEXT,NPEXT,MAX_SIM,NM_KI)
	    CALL TUNE(ITWO,'C_VARCONT')
!
! NB: VJ, VCHI, and VETA must not be modified until we have updated the
!     BA array.
!
	  END IF
!
! 
	  IF(USE_ELEC_HEAT_BAL .AND. COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION)THEN
	    CALL BA_EHB_BF_UPDATE_MPI_V1(VJ,ETA,CHI,POPS,RJ,
	1              FL,FQW(ML),COMPUTE_NEW_CROSS,FINAL_CONSTANT_CROSS,DO_SRCE_VAR_ONLY,
	1              NION,NT,NUM_BNDS,DST,DEND,ND)
	  END IF
!
! Modify the BA matrix for terms in the statistical equilibrium
! equations which are multiplied by RJ. NB. This is not for the
! variation of RJ - rather the multiplying factors. This section
! must be done for a LAMBDA iteration.
!
	CALL UPDATE_BA_FOR_LINE_MPI_V1(FL,FQW(ML),FREQ_INDX,
	1              POPS,JREC,dJRECdT,JPHOT,
	1              ND,NT,NUM_BNDS,NION,DIAG_INDX,
	1              TX_OFFSET,MAX_SIM,NM,NCF,NLF,
	1              LUER,FINAL_CONSTANT_CROSS,LST_ITERATION)
!
! 
!
! Free up LINE storage locations. Removal is done in 3 ways, but only 2 here.
! Recall that frequencies are ordered from highest to lowest, and that we
! integrate from blue to red.
!
! 1. Lines interact over at most 2Vinf from last point of resonance zone. Thus
!      when the current frequency is 2Vinf (converted to frequency units)
!      lower than the last frequency in the lines resonance zone it can safely
!      be removed.
!
! 2. Line is removed when the current frequency is lower by EXT_LINE_VAR*VINF
!      (converted to frequency units) than the last frequency in the resonance
!       zone. This is a control parameter, and may be used to speed up the code.
!       NB: EXT_LINE_VAR >= 0. Due to strong line overlap, it was found that
!       this method of line removal can cause issues when Vinf ~ 0 (e.g., in a
!       plane-parallel model). We thus put in a restriction of 300 km/s.
!
! 3. To make way for another line. This is only done when necessary, and is
!      done elsewhere. Only requirement is that the current frequency
!      is lower that the last frequency of the resonance zone.
!
	CALL TUNE(IONE,'CHK_L_FIN')
        T1=1.0_LDP-EXT_LINE_VAR*MAX(V(1),600.0_LDP)/2.998E+05_LDP
	DO SIM_INDX=1,MAX_SIM
	  IF(LINE_STORAGE_USED(SIM_INDX))THEN
!
! Check whether need storage location for net rate etc. We keep the storage
! until the line levels are removed the line variation setcion.
!
	    L=SIM_LINE_POINTER(SIM_INDX)
	    IF(NU(ML) .LT. NU(LINE_END_INDX_IN_NU(L))*T1)THEN
	      SIM_LINE_POINTER(SIM_INDX)=0
              LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	      LINE_LOC(L)=0
	    END IF
!
! Zero variation storage, if in use, and the storage location is not also
! being used by some other line.
!
	    IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION .AND.
	1        NU(ML) .LT. NU(LINE_END_INDX_IN_NU(L))*T1 .AND.
	1                                 .NOT. WEAK_LINE(SIM_INDX))THEN
	      NL=LOW_POINTER(SIM_INDX)
	      VAR_IN_USE_CNT(NL)=VAR_IN_USE_CNT(NL)-1
	      IF(VAR_IN_USE_CNT(NL) .EQ. 0)THEN
	        TX(:,:,NL)=0.0_LDP	!ND,ND,NM
	        TVX(:,:,NL)=0.0_LDP	!ND-1,ND,NM
	        dZ(NL,:,:,:)=0.0_LDP	!NM,NUM_BNDS,ND,MAX_SIM
	        VAR_LEV_ID(NL)=0
	      END IF
	      LOW_POINTER(SIM_INDX)=0
!
	      NUP=UP_POINTER(SIM_INDX)
	      VAR_IN_USE_CNT(NUP)=VAR_IN_USE_CNT(NUP)-1
	      IF(VAR_IN_USE_CNT(NUP) .EQ. 0)THEN
	        TX(:,:,NUP)=0.0_LDP	!ND,ND,NM
	        TVX(:,:,NUP)=0.0_LDP	!ND-1,ND,NM
	        dZ(NUP,:,:,:)=0.0_LDP	!NM,NUM_BNDS,ND,MAX_SIM
	        VAR_LEV_ID(NUP)=0
	      END IF
	      UP_POINTER(SIM_INDX)=0
	    END IF			!Outside region of influence by line?
	  END IF			!Line is in use.
	END DO				!Loop over line
	CALL TUNE(ITWO,'CHK_L_FIN')
!
	CALL TWO_PHOT_RATE_MPI_V1(T,RJ,FL,FQW(ML),DST,DEND,ND,NT)
!
! 
!
! Compute flux distribution and luminosity (in L(sun)) of star. NB: For
! NORDFLUX we always assume coherent scattering.
!
!
	CALL TUNE(IONE,'FLUX_DIST')
	IF(THIS_FREQ_EXT .AND. .NOT. CONT_VEL)THEN
!
! Since ETAEXT is not required any more, it will be used flux.
!
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  CALL MULTVEC(SOURCEEXT,ZETAEXT,THETAEXT,RJEXT,NDEXT)
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,
	1               SOURCEEXT,CHIEXT,dCHIdR,HQWEXT,ETAEXT,
	1               S1,THK_CONT,DIF,DBB,IC,NCEXT,NDEXT,NPEXT,METHOD)
	  CALL UNGRID(SOB,ND,ETAEXT,NDEXT,POS_IN_NEW_GRID)
	  SOB(2)=ETAEXT(2)				!Special case
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  N_OBS=NCF
	  OBS_FREQ(ML)=FL
	  OBS_FLUX(ML)=6.599341_LDP*SOB(1)*2.0_LDP		!2 DUE TO 0.5U
	ELSE IF(CONT_VEL)THEN
!
! TA is a work vector. TB initially used for extended SOB.
!
	   IF(ACCURATE)THEN
	     CALL REGRID_H(TB,REXT,RSQHNU,HFLUX_AT_OB,HFLUX_AT_IB,NDEXT,TA)
	     DO I=1,ND
	       SOB(I)=TB(POS_IN_NEW_GRID(I))
	     END DO
	   ELSE
	     CALL REGRID_H(SOB,R,RSQHNU,HFLUX_AT_OB,HFLUX_AT_IB,ND,TA)
	   END IF
	   IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NO_V)THEN
	     SOB(1)=HFLUX_AT_OB; SOB(ND)=HFLUX_AT_IB
	     SOB(1:ND)=SOB(1:ND)*R(ND)*R(ND)
	   END IF
	   IF(MYPE .EQ. 0)THEN
	     CALL COMP_OBS_V2(IPLUS,FL,
	1             IPLUS_STORE,NU_STORE,NST_CMF,
	1             MU_AT_RMAX,HQW_AT_RMAX,OBS_FREQ,OBS_FLUX,N_OBS,
	1             V_AT_RMAX,RMAX_OBS,'IPLUS','LIN_INT',DO_FULL_REL_OBS,
	1             FIRST_OBS_COMP,NP_OBS)
	   END IF
!
	ELSE
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  CALL MULTVEC(SOURCE,ZETA,THETA,RJ,ND)
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,R,Z,P,SOURCE,CHI,THETA,HQW,SOB,
	1               S1,THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  N_OBS=NCF
	  OBS_FREQ(ML)=FL
	  OBS_FLUX(ML)=6.599341_LDP*SOB(1)*2.0_LDP		!2 DUE TO 0.5U
	END IF
!
	IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NO_V)THEN
	  H_MOM(1:ND)=SOB(1:ND)/R(ND)/R(ND)
	ELSE
	  H_MOM(1:ND)=SOB(1:ND)/R(1:ND)/R(1:ND)
	END IF
!
! Evaluate (CMF) observd X-ray luminosities.
!
	  IF(ML .EQ. 1)THEN
	    OBS_XRAY_LUM_0P1=0.0_LDP
	    OBS_XRAY_LUM_1keV=0.0_LDP
	  END IF
	  T3=4.1274E-12_LDP
	  IF(NU(ML) .GT. 241.7988_LDP)OBS_XRAY_LUM_1keV=OBS_XRAY_LUM_1keV+T3*FQW(ML)*SOB(1)
	  IF(NU(ML) .GT. 24.17988_LDP)OBS_XRAY_LUM_0P1=OBS_XRAY_LUM_0P1+T3*FQW(ML)*SOB(1)
!
! Compute the luminosity, the FLUX mean opacity, and the ROSSELAND
! mean opacities.
!
	IF(ML .EQ. 1)THEN		!Need to move to main loop imit.
	  RLUMST(:)=0.0_LDP
	  J_INT(:)=0.0_LDP
	  H_INT(:)=0.0_LDP
	  DJDt_TERM(:)=0.0_LDP
	  dE_DJDt(:)=0.0_LDP
	  K_INT(:)=0.0_LDP
	  FLUX_MEAN(:)=0.0_LDP
	  ROSS_MEAN(:)=0.0_LDP
	  PLANCK_MEAN(:)=0.0_LDP
	  ABS_MEAN(:)=0.0_LDP
	  INT_dBdT(:)=0.0_LDP
	END IF
	T1=TWOHCSQ*HDKT*FQW(ML)*(NU(ML)**4)
	T3=TWOHCSQ*FQW(ML)*(NU(ML)**3)
	DO I=DST,DEND		              !(4*PI)**2*Dex(+20)/L(sun)
	  T2=SOB(I)*FQW(ML)*4.1274E-12_LDP
	  RLUMST(I)=RLUMST(I)+T2
	  J_INT(I)=J_INT(I)+RJ(I)*FQW(ML)*4.1274E-12_LDP
	  H_INT(I)=H_INT(I)+H_MOM(I)*FQW(ML)*4.1274E-12_LDP
	  K_INT(I)=K_INT(I)+K_MOM(I)*FQW(ML)*4.1274E-12_LDP
	  dE_DJDt(I)=dE_DJDt(I)+DJDt_TERM(I)*FQW(ML)*4.1274E-12_LDP
	  FLUX_MEAN(I)=FLUX_MEAN(I)+T2*CHI(I)
	  ABS_MEAN(I)=ABS_MEAN(I)+4.1274E-12_LDP*FQW(ML)*(CHI(I)-CHI_SCAT(I))*RJ(I)
	  T2=T1*EMHNUKT(I)/(  ( (1.0_LDP-EMHNUKT(I))*T(I) )**2  )
	  INT_dBdT(I)=INT_dBdT(I)+T2
	  ROSS_MEAN(I)=ROSS_MEAN(I)+T2/CHI(I)
!	  PLANCK_MEAN(I)=PLANCK_MEAN(I)+T3*CHI_NOSCAT(I)*EMHNUKT(I)/(1.0D0-EMHNUKT(I))
	  PLANCK_MEAN(I)=PLANCK_MEAN(I)+T3*(CHI(I)-CHI_SCAT(I))*EMHNUKT(I)/(1.0_LDP-EMHNUKT(I))
	END DO
	T1=SPEED_OF_LIGHT()*1.0E-05_LDP
	DO J=1,N_FLUX_MEAN_BANDS
	  IF(0.01_LDP*T1/FL .LT. LAM_FLUX_MEAN_BAND_END(J))THEN
	     BAND_FLUX_MEAN(DST:DEND,J)=FLUX_MEAN(DST:DEND)
	     BAND_FLUX(DST:DEND,J)=RLUMST(DST:DEND)
	     EXIT
	  END IF
	END DO
	CALL TUNE(ITWO,'FLUX_DIST')
!
! The current opacities and emissivities are stored for the variation of the
! radiation field at the next frequency.
!
	DO I=1,ND
	  CHI_PREV(I)=CHI_CONT(I)
	  CHI_NOSCAT_PREV(I)=CHI_NOSCAT(I)
	  CHI_SCAT_PREV(I)=CHI_SCAT(I)
	  ETA_PREV(I)=ETA_CONT(I)
	END DO
!
	IF(LST_ITERATION .AND. WRITE_RATES)THEN
	  DO SIM_INDX=1,MAX_SIM
	    IF(END_RES_ZONE(SIM_INDX))THEN
	      LS=SIM_LINE_POINTER(SIM_INDX)
	      IF(MYPE .EQ. 0)THEN
	        WRITE(LU_NET,'(/,1X,I6,2X,A,2X,F10.6,4(2X,I6))')
	1           LS,TRANS_NAME_SIM(SIM_INDX),VEC_FREQ(LS),
	1              VEC_NL(LS),VEC_NUP(LS),VEC_MNL_F(LS),VEC_MNUP_F(LS)
	        WRITE(LU_DR,'(/,1X,I6,2X,A,2X,F10.6,4(2X,I6))')
	1           LS,TRANS_NAME_SIM(SIM_INDX),VEC_FREQ(LS),
	1              VEC_NL(LS),VEC_NUP(LS),VEC_MNL_F(LS),VEC_MNUP_F(LS)
	        T3=(AVE_ENERGY(SIM_NL(SIM_INDX))-
	1              AVE_ENERGY(SIM_NUP(SIM_INDX)))/VEC_FREQ(LS)
	        WRITE(LU_HT,'(/,1X,I8,2X,A,2X,F10.6,2X,I6,2X,I6,ES14.5)')
	1           LS,TRANS_NAME_SIM(SIM_INDX),VEC_FREQ(LS),
	1              VEC_NL(LS),VEC_NUP(LS),T3
	    END IF
!
	      TA=0.0_LDP; TA(DST:DEND)=ZNET_SIM(DST:DEND,SIM_INDX)
	      CALL MPI_REDUCE(TA,TB,ND,MPI_DOUBLE_PRECISION,MPI_SUM,IZERO,MPI_COMM_WORLD,IERR)
	      IF(MYPE .EQ. 0)THEN
	        WRITE(LU_NET,'(1P,5E14.6)')(TB(I),I=1,ND)
	        FLUSH(UNIT=LU_NET)
	      END IF
!
	      TA=0.0_LDP
	      DO I=DST,DEND
	         TA(I)=ZNET_SIM(I,SIM_INDX)*POPS(SIM_NUP(SIM_INDX),I)*U_STAR_RATIO(I,SIM_INDX)*EINA(SIM_INDX)
	      END DO
	      CALL MPI_REDUCE(TA,TB,ND,MPI_DOUBLE_PRECISION,MPI_SUM,IZERO,MPI_COMM_WORLD,IERR)
	      IF(MYPE .EQ. 0)WRITE(LU_DR,40003)(TB(I),I=1,ND)
!
	      IF(SCL_LINE_COOL_RATES .OR. SCL_SL_LINE_OPAC)THEN
	        SCL_FAC=(AVE_ENERGY(SIM_NL(SIM_INDX))-
	1            AVE_ENERGY(SIM_NUP(SIM_INDX)))/VEC_FREQ(LS)
	        IF(ABS(SCL_FAC-1.0_LDP) .GT. SCL_LINE_HT_FAC)SCL_FAC=1.0_LDP
	      ELSE
	        SCL_FAC=1.0_LDP
	      END IF
	      T3=SCL_FAC; IF(SCL_SL_LINE_OPAC)T3=1.0_LDP
	      TA=0.0_LDP
	      DO I=DST,DEND
	        TA(I)=T3*ZNET_SIM(I,SIM_INDX)*ETAL_MAT(I,SIM_INDX)
	      END DO
	      CALL MPI_REDUCE(TA,TB,ND,MPI_DOUBLE_PRECISION,MPI_SUM,IZERO,MPI_COMM_WORLD,IERR)
	      IF(MYPE .EQ. 0)WRITE(LU_HT,'(1X,1P,5E12.4)')(TB(I),I=1,ND)
!
! As these are used only for diagnostic purposes, we scale them independent of
! the density.
!
	      T3=SCL_FAC
	      IF(SCL_LINE_COOL_RATES)THEN
	        DO K=DST,DEND
	          T2=ETAL_MAT(K,SIM_INDX)*ZNET_SIM(K,SIM_INDX)
	          STEQ_T_SCL(K)=STEQ_T_SCL(K) - T2*T3
	          STEQ_T_NO_SCL(K)=STEQ_T_NO_SCL(K) - T2
	        END DO
	      ELSE
	        DO K=DST,DEND
	          T2=ETAL_MAT(K,SIM_INDX)*ZNET_SIM(K,SIM_INDX)
	          STEQ_T_SCL(K)=STEQ_T_SCL(K) - T2
	          STEQ_T_NO_SCL(K)=STEQ_T_NO_SCL(K) - T2/T3
	        END DO
	      END IF
	      TA=0.0_LDP; TA(DST:DEND)=STEQ_T_SCL(DST:DEND)
	      CALL MPI_REDUCE(TA,TB,ND,MPI_DOUBLE_PRECISION,MPI_SUM,IZERO,MPI_COMM_WORLD,IERR)
	      WRITE(LU_HT,'(/,(1X,1P,5E12.4))')(TB(I), I=1,ND)
	      TA=0.0_LDP; TA(DST:DEND)=STEQ_T_NO_SCL(DST:DEND)
	      CALL MPI_REDUCE(TA,TB,ND,MPI_DOUBLE_PRECISION,MPI_SUM,IZERO,MPI_COMM_WORLD,IERR)
	      WRITE(LU_HT,'(/,(1X,1P,5E12.4))')(TB(I), I=1,ND)
	    END IF
	  END DO
	END IF
!
! Needs revising.
!
!
	IF( .NOT. USE_FIXED_J .AND. (USE_ELEC_HEAT_BAL .OR. COMP_STEQ_T_EHB) )THEN
	  CALL TUNE(IONE,'FF_EB_COR')
	  IF(ML .EQ. 1)WRITE(6,*)'Free-free correction may need revising'
	  CALL  COMP_FREE_FREE_MPI_V1(CHI,ETA,VCHI,VETA,CONT_FREQ,FL,
	1           FIRST_FREQ,USE_ELEC_HEAT_BAL,COMPUTE_BA,ND,NT)
	  T1=1.0E-10_LDP*16.0_LDP*ATAN(1.0_LDP)*FQW(ML)
	  DO I=DST,DEND
	    STEQ_T_EHB(I)=STEQ_T_EHB(I)+T1*(CHI(I)*RJ(I)-ETA(I))
	  END DO
	  CALL TUNE(ITWO,'FF_EB_COR')
	  IF(ML .EQ. NCF .AND. LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	    CALL WR2D_MPI_V1(STEQ_T_EHB,IONE,DST,DEND,ND,'After free-free','&',L_TRUE,LU_T_EHB)
	  END IF
	END IF
!
	IF(.NOT. USE_FIXED_J .AND. USE_ELEC_HEAT_BAL .AND. COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION)THEN
	  CALL BA_EHB_FF_UPDATE_MPI_V1(VJ,VCHI,VETA,
	1              ETA,CHI,T,POPS,RJ,
	1              FQW(ML),COMPUTE_NEW_CROSS,FINAL_CONSTANT_CROSS,DO_SRCE_VAR_ONLY,
	1              NION,NT,NUM_BNDS,ND,IONE,ND)
!
	  TB(1:ND)=BA_T_PAR_EHB(NT,1:ND)
	  DO ID=1,NUM_IONS
	    ID_SAV=ID
	    IF(ATM(ID)%XzV_PRES .AND. FINAL_CONSTANT_CROSS)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        CALL VEHB_BYJ_V1(ID_SAV,
	1             ATM(ID)%WSXzV(1,1,J), ATM(ID)%dWSXzVdT(1,1,J),
	1             ATM(ID)%WCRXzV(1,1,J), ATM(ID)%dWCRXzVdT(1,1,J),
	1             ATM(ID)%XzV, ATM(ID)%XzVLTE, ATM(ID)%dlnXzVLTE_dlnT,
	1             ATM(ID)%NXzV,ATM(ID)%EQXzV,
	1             ATM(ID+1)%XzV, ATM(ID+1)%LOG_XzVLTE,
	1             ATM(ID+1)%dlnXzVLTE_dlnT, ATM(ID+1)%NXzV,
	1             ATM(ID)%XzV_ION_LEV_ID(J),ED,T,
	1             JREC,dJRECdt,JPHOT,
	1             JREC_CR,dJREC_CRdt,JPHOT_CR,
	1             FIXED_T,ND,IONE,ND,NT)
	     END DO
	   END IF
	  END DO
!
	  IF( MOD(FREQ_INDX,N_PAR) .EQ. 0 .OR. FREQ_INDX .EQ. NCF )THEN
	    BA_T_EHB(:,DIAG_INDX,:)=BA_T_EHB(:,DIAG_INDX,:)+BA_T_PAR_EHB
	    BA_T_PAR_EHB=0.0_LDP
	  END IF
!	  WRITE(294,'(I5,6ES14.4)')ML,BA_T_EHB(1:3,DIAG_INDX,75),BA_T_EHB(NT-2:NT,DIAG_INDX,75)
!	  WRITE(295,'(I5,6ES14.4)')ML,BA_T_PAR_EHB(1:3,75),BA_T_PAR_EHB(NT-2:NT,75)
!	  FLUSH(UNIT=294); FLUSH(UNIT=295)
	END IF
!
	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	  T1=FQW(ML)*(CHI_NOSCAT(DPTH_INDX)*RJ(DPTH_INDX) - ETA_NOSCAT(DPTH_INDX))
	  IF(MOD(FREQ_INDX,N_PAR) .EQ. 0)THEN
	    WRITE(199,'(I10,12ES18.8)')ML,FL,STEQ_T(DPTH_INDX),T1,FQW(ML)*ETA_NOSCAT(DPTH_INDX),
	1                           STEQ_T_SCL(DPTH_INDX),STEQ_T_NO_SCL(DPTH_INDX),
	1                           BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),BA_T(NT,DIAG_INDX,DPTH_INDX)
	  ELSE
	    WRITE(199,'(I10,12ES18.8)')ML,FL,STEQ_T(DPTH_INDX),T1,FQW(ML)*ETA_NOSCAT(DPTH_INDX),
	1                           STEQ_T_SCL(DPTH_INDX),STEQ_T_NO_SCL(DPTH_INDX),
	1                           BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX)+BA_T_PAR(VAR_INDX,DPTH_INDX),
	1                           BA_T(NT,DIAG_INDX,DPTH_INDX)+BA_T_PAR(NT,DPTH_INDX)
	  END IF
	END IF
!
	IF(MYPE .EQ. 11)THEN
	  WRITE(524,*)ML,SE(1)%BA(1,1,DIAG_INDX,DST),SE(1)%BA(ATM(1)%NXzV:SE(1)%N_SE,1,DIAG_INDX,DST)
	  FLUSH(UNIT=524)
	END IF
!
	IF(MYPE .EQ. 0)THEN
          WRITE(STRING,*)ML; STRING='ML='//ADJUSTL(TRIM(STRING))
          CALL WRITV(CHI,ND,TRIM(STRING),841)
          CALL WRITV(ETA,ND,TRIM(STRING),842)
	END IF

!	IF(MOD(ML,5000) .EQ. 0)CALL TUNE(3,' ')
10000	CONTINUE
	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	  WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'End cont loop'
	  CALL WR2D_MPI_V1(STEQ_T_EHB,IONE,DST,DEND,ND,'After continuum loop','&',L_TRUE,LU_T_EHB)
	END IF
	CALL TUNE(ITWO,'10000')
!
	IF(MYPE .EQ. 0) WRITE(LUER,'(/,A,I10,/)')' Number of weak lines is:',NUM_OF_WEAK_LINES
!
! 
!
	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	  CALL WR2D_MPI_V1(STEQ_T,IONE,DST,DEND,ND,'Before two phot','&',L_TRUE,183)
	END IF
	CALL STEQ_BA_TWO_PHOT_RATE_MPI_V1(POPS,NT,DST,DEND,ND,
	1         DIAG_INDX,COMPUTE_BA,LUMOD,L_FALSE)                  !LST_ITERATION)
	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	  WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'two'
	  CALL WR2D_MPI_V1(STEQ_T,IONE,DST,DEND,ND,'After two phot','&',L_TRUE,183)
	END IF
!
! 
!
! Include influence of charge exchange reactions.
!
	DO ID=1,NUM_IONS-1
	  ID_SAV=ID
	  IF(ATM(ID)%XzV_PRES)THEN
	    CALL SET_CHG_EXCH_V4(ID_SAV, ATM(ID)%XzVLEVNAME_F,
	1       ATM(ID)%EDGEXzV_F,  ATM(ID)%GXzV_F,
	1       ATM(ID)%F_TO_S_XzV, ATM(ID)%GIONXzV_F,
	1       ATM(ID)%NXzV_F, ATM(ID)%NXzV, ND,
	1       ATM(ID)%EQXzV, EQ_SPECIES(SPECIES_LNK(ID)), T)
	  END IF
	END DO
!
	CALL STEQ_BA_CHG_EXCH_MPI_V1(POPS,T,NT,DST,DEND,ND,DIAG_INDX,COMPUTE_BA)
	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	  WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'chg_1'
	END IF
!
	DO ID=1,NUM_IONS-1
	  CALL EVAL_CHG_RATES_MPI_V1(ATM(ID)%CHG_PRXzV, ATM(ID)%CHG_RRXzV, ION_ID(ID),POPS,T,DST,DEND,ND,NT)
	END DO
	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	  WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'chg_2'
	  CALL WR2D_MPI_V1(STEQ_T,IONE,DST,DEND,ND,'After charge','&',L_TRUE,183)
	END IF
!
! Penning
!
!    HeI(1s_2s_3Se) + H(1s_2Se)  --->  HeI(1s2_1Se) +  H+
!
	IF(INCL_PENNING_ION)THEN
	  CALL DO_PENNING_ION_MPI_V1(HDKT,COMPUTE_BA,DIAG_INDX,ND)
	END IF
!
! 
!
! Output errors that have occurred in MOM_J_CMF
!
	CALL WRITE_J_CMF_ERR(MAIN_COUNTER)
!
! Allow for advection terms.
!
	IF(SN_MODEL .AND. DO_CO_MOV_DDT)THEN
          CALL STEQ_CO_MOV_DERIV_MPI_V1(POPS,ADVEC_RELAX_PARAM,LINEAR_ADV,
	1             DO_CO_MOV_DDT,LAMBDA_ITERATION,COMPUTE_BA,
	1             TIME_SEQ_NO,NUM_BNDS,ND,NT)
	  IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	    WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'SN_DDT'
	    CALL WR2D_MPI_V1(STEQ_T,IONE,DST,DEND,ND,'After comoving derivative.','&',L_TRUE,183)
	  END IF
	ELSE
	  DO ID=1,NION
	    SE(ID)%STEQ_ADV=0.0_LDP
	  END DO
	  IF(INCL_ADVECTION)THEN
	    IF(MYPE .EQ. 0)WRITE(6,*)'STE_ADVEC_V4 needs to be updated to USE POPS'; FLUSH(UNIT=6)
	    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	    CALL MPI_FINALIZE (ierr)
	    STOP  
	  END IF
	  CALL STEQ_ADVEC_MPI_V1(ADVEC_RELAX_PARAM,LINEAR_ADV,NUM_BNDS,ND,
	1            INCL_ADVECTION,LAMBDA_ITERATION,COMPUTE_BA)
	  IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	    WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'SN_DDT'
	    CALL WR2D_MPI_V1(STEQ_T,IONE,DST,DEND,ND,'After advection','&',L_TRUE,183)
	  END IF
	END IF
!
! Allow for adiabatic cooling, if requested.
!
	IF(SN_MODEL .AND. DO_CO_MOV_DDT)THEN
	  CALL EVAL_TEMP_DDT_V2(dE_WORK,AD_COOL_V,AD_COOL_DT,
	1                       POPS,AVE_ENERGY,HDKT,
	1                       COMPUTE_BA,INCL_ADIABATIC,
	1                       TIME_SEQ_NO,DIAG_INDX,NUM_BNDS,NT,ND)
	  WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'SN_ADD'
	ELSE IF(PLANE_PARALLEL_NO_V .OR.  PLANE_PARALLEL)THEN
	  AD_COOL_V(DST:DEND)=0.0_LDP; AD_COOL_DT(DST:DEND)=0.0_LDP
	ELSE
	  dE_WORK=0.0_LDP; TA(DST:DEND)=STEQ_T_EHB(DST:DEND)
	  CALL EVAL_ADIABATIC_MPI_V1(AD_COOL_V,AD_COOL_DT,
	1                       POPS,AVE_ENERGY,HDKT,
	1                       COMPUTE_BA,INCL_ADIABATIC,
	1                       DIAG_INDX,NUM_BNDS,NT,ND)
	  IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	    WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'Adiabatict'
	    CALL WR2D_MPI_V1(STEQ_T,IONE,DST,DEND,ND,'After adiabatic','&',L_TRUE,183)
	    CALL WR2D_MPI_V1(STEQ_T_EHB,IONE,DST,DEND,ND,'After adiabatic','&',L_TRUE,LU_T_EHB)
	    TA(DST:DEND)=STEQ_T_EHB(DST:DEND)-TA(DST:DEND)
	    CALL WR2D_MPI_V1(TA,IONE,DST,DEND,ND,'Adiabatic term','&',L_TRUE,LU_T_EHB)
	  END IF
	END IF
!
	IF(INC_SHOCK_POWER)THEN
	  CALL EVAL_SHOCK_POWER(dE_SHOCK_POWER,ND,SHOCK_POWER_FAC)
	  STEQ_T_EHB=STEQ_T_EHB+SHOCK_POWER_FAC*dE_SHOCK_POWER
	  IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	    CALL WR2D_MPI_V1(STEQ_T_EHB,IONE,DST,DEND,ND,'After shock power','&',L_TRUE,LU_T_EHB)
	  END IF
	ELSE
	  dE_SHOCK_POWER=0.0_LDP
	END IF
!
	IF(SN_MODEL .AND. INCL_RADIOACTIVE_DECAY)THEN
	  IF(TREAT_NON_THERMAL_ELECTRONS)THEN
	    STEQ_T_EHB=STEQ_T_EHB+dE_RAD_DECAY
	  ELSE IF(SN_MODEL .AND. INCL_RADIOACTIVE_DECAY)THEN
	    CALL EVAL_RAD_DECAY_V1(dE_RAD_DECAY,NT,ND)
	  END IF
	  IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	    WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'SN_Rad_Decay'
	    CALL WR2D_MPI_V1(STEQ_T_EHB,IONE,DST,DEND,ND,'After rad. decay','&',L_TRUE,LU_T_EHB)
	  END IF
	ELSE
	  dE_RAD_DECAY=0.0_LDP
	END IF
!
! Prevent T from becoming too small by adding a extra heating term.
!
	CALL PREVENT_LOW_T_MPI_V1(ARTIFICIAL_HEAT_TERM,T_MIN,COMPUTE_BA,LAMBDA_ITERATION,
	1                   T_MIN_BA_EXTRAP,DIAG_INDX,NUM_BNDS,ND,NT)
	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	  WRITE(199,'(I10,2ES18.8,3X,A)')ML,STEQ_T(DPTH_INDX),BA_T(VAR_INDX,DIAG_INDX,DPTH_INDX),'Artificial heat'
	  CALL WR2D_MPI_V1(STEQ_T_EHB,IONE,DST,DEND,ND,'After artificial heat','&',L_TRUE,LU_T_EHB)
	END IF
!
! Write pointer file and store BA, BA_ED and B_T matrices.
!
!	BA_T(:,:,ND)=0.0D0; STEQ_T(ND)=T(ND)-T(ND-1)
!	BA_T(NT,DIAG_INDX,ND)=1.0D0
!	IF(DIAG_INDX .NE. 1)BA_T(NT,DIAG_INDX-1,ND)=-1.0D0
	IF(COMPUTE_BA .AND. WRBAMAT .AND. .NOT. FLUX_CAL_ONLY .AND. .NOT. LAMBDA_ITERATION)THEN
	  CALL TUNE(IONE,'STORE_BA')
	    CALL STORE_BA_DATA_MPI_V1(LU_BA,NION,NUM_BNDS,COMPUTE_BA,FIXED_T,'BAMAT')
	  CALL TUNE(ITWO,'STORE_BA')
	END IF
!
! Store radiative equlibrium equation so we can check influence on radiation field.
!
	DEP_RAD_EQ(DST:DEND)=STEQ_T(DST:DEND)
! 
!
! Write out recombination, photoionization and cooling terms for digestion.
!
! Since X_RECOM (and X_COOL) were dimensioned (ND,0:NION), and since it was
! initialized, we can assume X_RECOM(1,0) is zero, and hence it can be used
! as a  dummy vector for the special case when ATM(ID-1)%INDX=0.
!
! NB: For the first ionization stage of a species, ATM(ID-1)%INDX will
!     refer to the highest level of the previous ion.  By convention in
!     CMFGEN this refers to a 1 level state with an INDEX of 0 (and
!     XzV_PRES for this species is FALSE).
!
	IF(LST_ITERATION)THEN   ! .AND. VERBOSE)THEN  ! .AND. .NOT. USE_FIXED_J)THEN
	  CALL WRITE_RECOM_MPI_V1(ND)
	  IF(MYPE .EQ. 0)THEN
	    WRITE(420,*)de_SHOCK_POWER; FLUSH(UNIT=420)
	  END IF
	  IF(MYPE .EQ. 1)THEN
	    WRITE(421,*)de_SHOCK_POWER; FLUSH(UNIT=421)
	  END IF
	  IF(MYPE .EQ. 0)THEN
	    WRITE(422,*)de_RAD_DECAY; FLUSH(UNIT=422)
	  END IF
	  IF(MYPE .EQ. 1)THEN
	    WRITE(423,*)de_RAD_DECAY; FLUSH(UNIT=423)
	  END IF
	  CALL WR_COOL_MPI_V1(AD_COOL_V,AD_COOL_DT,ARTIFICIAL_HEAT_TERM,
	1                dE_RAD_DECAY,dE_SHOCK_POWER,
	1                XRAY_LUM_TOT,INCL_ADIABATIC,ND)
	END IF		!Only output if last iteration.
! 
!
	TA=RLUMST;   CALL MPI_ALLREDUCE(TA,RLUMST,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	WRITE(STRING,'(I5)')MAIN_COUNTER; STRING=ADJUSTL(STRING)
	IF(MYPE .EQ. 0)THEN
	  STRING=' Luminosity of star (d=1,ND)(iteration '//TRIM(STRING)//') is:'
	  WRITE(LUER,'(A,2ES18.8,/)')TRIM(STRING),RLUMST(1),RLUMST(ND)
	  FLUSH(LUER)
	END IF
	IF(RLUMST(1) .LE. 0.0_LDP)RLUMST(1)=1.0E-20_LDP
	J=0
	DO I=DST,DEND
	  IF( ABS(RLUMST(I)) .LT. 1.0E-10)THEN
	     RLUMST= SIGN(1.0E-010_LDP,RLUMST); J=J+1
	  END IF
	END DO
	IF(J .NE. 0)WRITE(LUER,*)'|RLUMST| set to 1.0E-10 at some depths'
	RLUMST_BND=RLUMST(2)		!2 is used to avoid glitch at outer boundary in some models.
!
	IF(MYPE .EQ. 0)THEN
	  CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFLUX','UNKNOWN',' ',' ',IZERO,IOS)
	    WRITE(STRING,'(I10)')N_OBS
	    STRING=ADJUSTL(STRING)
            STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	    CALL WRITV_V2(OBS_FREQ,N_OBS,ISEV,TRIM(STRING),LU_FLUX)
	    CALL WRITV_V2(OBS_FLUX,N_OBS,IFOUR,'Observed intensity (Janskys)',LU_FLUX)
	    CALL WRITV(RLUMST,ND,'Luminosity',LU_FLUX)
	  CLOSE(UNIT=LU_FLUX)
!
! We don't check the last iteration, since the observed spectrum is computed wih a 
! different (and more accurate routine).
!
	  IF(.NOT. LAMBDA_ITERATION .AND. .NOT. LST_ITERATION .AND. .NOT. USE_FIXED_J)THEN
	    CALL CHECK_SPEC_CONV(OBS_FLUX,OBS_FREQ,T1,T2,CHK,N_OBS)
	  END IF
	END IF
	CALL TUNE(ITWO,'MLCF')
!
! Compute ROSSELAND and FLUX mean opacities. These MEAN opacities DO NOT
! include the effect of clumping. Compute the respective optical depth scales;
! TA is used for the FLUX mean optical depth scale, TB for the ROSSELAND mean
! optical depth scale, and DTAU for the electron scattering optical depth
! scale. The optical depth scale INCLDUES the effects of clumping. TCHI is used
! as a temporary work vector.
!
! T1=4 * [STEFAN BOLTZMAN CONS] * 1.0D+16 / pi
!
	T1=7.218771E+11_LDP
	DO I=DST,DEND
	  FLUX_MEAN(I)=FLUX_MEAN(I)/RLUMST(I)
	  ABS_MEAN(I)=ABS_MEAN(I)/J_INT(I)
	  INT_dBdT(I)=INT_dBdT(I)/ROSS_MEAN(I)		!Program rosseland opac.
	  ROSS_MEAN(I)=T1*( T(I)**3 )/ROSS_MEAN(I)
	  PLANCK_MEAN(I)=4.0_LDP*PLANCK_MEAN(I)/T1/(T(I)**4)
	END DO
	DO J=1,N_FLUX_MEAN_BANDS
	  BAND_FLUX_MEAN(DST:DEND,J)=BAND_FLUX_MEAN(DST:DEND,J)/RLUMST(DST:DEND)
	END DO
!
	TA=ROSS_MEAN;   CALL MPI_ALLREDUCE(TA,ROSS_MEAN,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	TA=ABS_MEAN;    CALL MPI_ALLREDUCE(TA,ABS_MEAN,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	TA=FLUX_MEAN;   CALL MPI_ALLREDUCE(TA,FLUX_MEAN,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	TA=PLANCK_MEAN; CALL MPI_ALLREDUCE(TA,PLANCK_MEAN,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	TA=INT_dBdT;    CALL MPI_ALLREDUCE(TA,INT_dBdT,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	K=ND*N_FLUX_MEAN_BANDS
	DO J=1,N_FLUX_MEAN_BANDS
	  CALL ROOT_GATHERV_VEC_MPI_V1(BAND_FLUX(:,J),ND)
	  CALL ROOT_GATHERV_VEC_MPI_V1(BAND_FLUX_MEAN(:,J),ND)
	END DO
!
	IF(MYPE .EQ. 0)THEN
	  TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(TA,TCHI,R,R,dCHIdR,ND)
	  TCHI(1:ND)=FLUX_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  IF(MINVAL(TCHI) .GT. 0)THEN
	    CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
	  ELSE
	    dCHIdR(1:ND)=0.0_LDP
	  END IF
          CALL NORDTAU(TB,TCHI,R,R,dCHIdR,ND)
	  TCHI(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
!
	  TA(ND)=0.0_LDP
	  TB(ND)=0.0_LDP
	  TC(ND)=0.0_LDP
	  DTAU(ND)=0.0_LDP
	  ROSS_PHOT_DPTH_INDX=ND
!
	  CALL GEN_ASCI_OPEN(LU_OPAC,'MEANOPAC','UNKNOWN',' ',' ',IZERO,IOS)
	    WRITE(LU_OPAC,
	1  '( ''         R          I   Tau(Ross)   /\Tau  Rat(Ross)'//
	1  ' Chi(Ross)  Chi(ross)  Chi(Flux)   Chi(es) '//
	1  '  Tau(Flux)  Tau(es)  Rat(Flux)  Rat(es)   Kappa(R)     V(km/s)'' )' )
	    IF(R(1) .GE. 1.0E+05_LDP)THEN
	      FMT='(ES17.10,I4,2ES10.3,ES10.2,4ES11.3,4ES10.2,2ES11.3)'
	    ELSE
	      FMT='( F17.10,I4,2ES10.3,ES10.2,4ES11.3,4ES10.2,2ES11.3)'
	    END IF
	    DO I=1,ND
	      IF(I .EQ. 1)THEN
	        T1=LOG(ROSS_MEAN(1)*CLUMP_FAC(1)/ROSS_MEAN(4)/CLUMP_FAC(4))/LOG(R(4)/R(1))
	        IF(T1 .LT. 2.0_LDP)T1=2.0_LDP
	        T1=ROSS_MEAN(1)*CLUMP_FAC(1)*R(1)/(T1-1.0_LDP)		!Rosseland optical depth scale
	        T2=LOG(ABS(FLUX_MEAN(1)*CLUMP_FAC(1)/FLUX_MEAN(4)/CLUMP_FAC(4)))/LOG(R(4)/R(1))
	        IF(T2 .LT. 2.0_LDP)T2=2.0_LDP
	        T2=FLUX_MEAN(1)*CLUMP_FAC(1)*R(1)/(T2-1.0_LDP)		!Flux optical depth scale
	        T3=LOG(ESEC(1)*CLUMP_FAC(1)/ESEC(4)/CLUMP_FAC(4))/LOG(R(4)/R(1))
	        IF(T3 .LT. 2.0_LDP)T3=2.0_LDP
	        T3=ESEC(1)*CLUMP_FAC(1)*R(1)/(T3-1.0_LDP)			!Electon scattering optical depth scale
	        TC(1:3)=0.0_LDP
	      ELSE
	        T1=T1+TA(I-1)                                             !Rosseland optical depth scale
	        IF(T1 .LT. 0.667_LDP)ROSS_PHOT_DPTH_INDX=I
	        T2=T2+TB(I-1)                                             !Flux optical depth scale
	        T3=T3+DTAU(I-1)                                           !Electon scattering optical depth scale
	        TC(1)=TA(I)/TA(I-1)
	        TC(2)=TB(I)/TB(I-1)
	        TC(3)=DTAU(I)/DTAU(I-1)
	      END IF
	      WRITE(LU_OPAC,FMT)R(I),I,T1,TA(I),TC(1),
	1        ROSS_MEAN(I),INT_dBdT(I),FLUX_MEAN(I),ESEC(I),
	1      T 2,T3,TC(2),TC(3),1.0D-10*ROSS_MEAN(I)/DENSITY(I),V(I)
	    END DO
	    IF(T1 .LT. 30.0_LDP .AND. .NOT. SN_MODEL .AND. .NOT. USE_FIXED_J)THEN
	      WRITE(6,*)'This model does not extend sufficiently deeply to guarantee the accuracy of'
	      WRITE(6,*)' the inner boundary condition. You need to adjust you model so that the'
	      WRITE(6,*)' minimum optical depth is 50, and preferably 100'
	      IF(STOP_IF_MAJOR_WARNING)STOP
	    END IF
	    WRITE(LU_OPAC,'(//,A,A)')
	1     'NB: Mean opacities do not include effect of clumping',
	1     'NB: Optical depth scale includes effect of clumping'
	  CLOSE(UNIT=LU_OPAC)
	END IF
!
	IF(MYPE .EQ. 0)WRITE(6,*)'Grey recompute has been compmented out'
	IF(MYPE .EQ. 50 .AND. LST_ITERATION .AND. .NOT. RD_FIX_T)THEN
!
! Compute the grey temperature distribution and the Rosseland optical
! depth scale (returned in TA). When CHK is TRUE, the grey temperature
! distribution has been successfully computed.
!
	  CHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  TCHI(1:ND)=PLANCK_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  IF(COMP_GREY_LST_IT)THEN
	    CALL COMP_GREY_V4(POPS,TGREY,TA,CHI,TCHI,CHK,LUER,NC,ND,NP,NT)
	    IF(CHK)THEN
	      WRITE(LUER,'(/,1X,A,/)')'Grey solution was successfully computed'
	    ELSE
	      WRITE(LUER,'(/,1X,A,/)')'Grey solution was NOT successfully computed'
	    END IF
	  END IF
!
	  IF(CHK .AND. COMP_GREY_LST_IT .AND. .NOT. USE_FIXED_J)THEN
	    OPEN(UNIT=LUIN,FILE='GREY_SCL_FACOUT',STATUS='UNKNOWN')
	      WRITE(LUIN,'(A)')'!'
	      WRITE(LUIN,'(A,8X,A,7X,A,7X,A,6X,A)')'!','Log(Tau)','T/T(grey)','T(10^4 K)','L'
	      WRITE(LUIN,'(A)')'!'
	      WRITE(LUIN,*)ND
	      DO I=1,ND
	        IF(TA(I) .GT. 0)THEN
	          WRITE(LUIN,'(2X,3ES16.6,4X,I3)')LOG10(TA(I)),T(I)/TGREY(I),T(I),I
	        ELSE
	          WRITE(LUER,'(A)')' Bad Roseeland optical depth scale for T/TGREY output'
	          WRITE(LUIN,'(A)')' Bad Roseeland optical depth scale for T/TGREY output'
	          EXIT
	        END IF
              END DO
	    END IF
	  CLOSE(LUIN)
	END IF
!
! Output hydrodynamical terms to allow check on radiation driving of the wind.
!
	IF(MYPE .EQ. 0 .AND. .NOT. SN_MODEL .AND. .NOT. USE_FIXED_J .AND. .NOT. PNT_SRCE_MOD)THEN
!	IF(.NOT. USE_FIXED_J)THEN
	  I=18
	  CALL WRITE_VEC(RLUMST,ND,'RLUMST',743)
	  CALL WRITE_VEC(FLUX_MEAN,ND,'FLUX_MEAN',743)
	  CALL WR2D(BAND_FLUX,ND,N_FLUX_MEAN_BANDS,'BAND_FLUX',743)
	  CALL WR2D(BAND_FLUX_MEAN,ND,N_FLUX_MEAN_BANDS,'BAND_FLUX_MEAN',743)
!
	  CALL HYDRO_TERMS_V5(POP_ATOM,R,V,T,SIGMA,ED,CLUMP_FAC,RLUMST,
	1                 LOGG,STARS_MASS,MEAN_ATOMIC_WEIGHT,
	1		  FLUX_MEAN,ROSS_MEAN,ESEC,
	1                 LAM_FLUX_MEAN_BAND_END,BAND_FLUX_MEAN,
	1		  PRESSURE_VTURB,PLANE_PARALLEL,PLANE_PARALLEL_NO_V,
	1                 LST_ITERATION,BAND_FLUX,N_FLUX_MEAN_BANDS,I,ND)
	END IF
!
	IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	     CALL WR_ASCI_STEQ_MPI_V1(NION,DST,DEND,ND,'STEQ ARRAY- Continuum Terms',LU_DR)
	END IF
!
! 
	CALL TUNE(IONE,'LINE LOOP')
!
! Zero LLUMST here, since total line luminosity output to OBSFLUX file
! even if the section of code is not executed.
!
	LLUMST(1:ND)=0.0_LDP
	IF(GLOBAL_LINE_SWITCH(1:5) .EQ. 'BLANK')THEN
!
! These lines have been treated with the continuum.
!
	ELSE IF(FLUX_CAL_ONLY .AND. .NOT. DO_SOBOLEV_LINES)THEN
!
! Don't waste time by computing rates and EW's for lines treated with the
! Sobolev approximation.
!
	ELSE
!
	END IF
!
! Write pointer file and then store BA and STEQ matrices.
!
! These are the final matrices.
!
	IF(COMPUTE_BA
	1          .AND. WRBAMAT
	1          .AND. .NOT. FLUX_CAL_ONLY
	1          .AND. .NOT. LAMBDA_ITERATION)THEN
	  CALL TUNE(IONE,'BAMAT_WR')
	  CALL STORE_BA_DATA_MPI_V1(LU_BA,NION,NUM_BNDS,COMPUTE_BA,FIXED_T,'BAMAT')
	  CALL TUNE(ITWO,'BAMAT_WR')
	END IF
!
! Ends check on GLOBAL_LINE.
!
	CALL TUNE(ITWO,'LINE LOOP')
!
! Compute the the total line luminosity, and the total dielectronic line
! emission emitted in each shell between i and i+1.
! Not all of this flux will be received by the observer due to continuum
! absorption. The factor of 0.5 arrises from the average of R*R*ZNET in the
! shell. [ 4.1274D-12=(4pi*1.0D+20)/Lsun * 4PI ]
!
	IF(.NOT. XRAYS)THEN
	  XRAY_LUM_TOT(1:ND)=0.0_LDP
	  XRAY_LUM_0P1(1:ND)=0.0_LDP
	  XRAY_LUM_1KEV(1:ND)=0.0_LDP
	END IF
	T1=4.1274E-12_LDP
	T2=1.0E+10_LDP*T1/4.0_LDP/ACOS(-1.0_LDP)
	IF(USE_FIXED_J)THEN
	  DIELUM=1.0E-20_LDP; DEP_RAD_EQ=1.0E-20_LDP; dE_WORK=1.0E-20_LDP
	  RAD_DECAY_LUM=1.0E-20_LDP
	  XRAY_LUM_TOT=1.0E-20_LDP; XRAY_LUM_0P1=0.0_LDP; XRAY_LUM_1KEV=0.0_LDP
	END IF
!
	DO I=DST,DEND
	  LLUMST(I)=LLUMST(I)*R(I)*R(I)*T1
	  DIELUM(I)=DIELUM(I)*R(I)*R(I)*T1
	  DEP_RAD_EQ_LUM(I)=DEP_RAD_EQ(I)*R(I)*R(I)*T1
	  XRAY_LUM_TOT(I)=dE_XRAY_TOT(I)*R(I)*R(I)*T1
	  XRAY_LUM_0P1(I)=dE_XRAY_0P1(I)*R(I)*R(I)*T1
	  XRAY_LUM_1KEV(I)=dE_XRAY_1KEV(I)*R(I)*R(I)*T1
	  dE_WORK_LUM(I)=dE_WORK(I)*R(I)*R(I)*T1*CLUMP_FAC(I)
	  RAD_DECAY_LUM(I)=dE_RAD_DECAY(I)*R(I)*R(I)*T2*CLUMP_FAC(I)		!As ergs/cm^3
	  SHOCK_POWER_LUM(I)=dE_SHOCK_POWER(I)*R(I)*R(I)*T2*CLUMP_FAC(I)
	  DJDT_LUM(I)=dE_DJDT(I)
	END DO
!
	IF(PLANE_PARALLEL_NO_V)THEN
	  MECH_LUM(DST:DEND)=0.0_LDP
	ELSE IF(USE_DJDT_RTE .AND. USE_Dr4JDt)THEN
	  MECH_LUM(DST:DEND)=0.0_LDP
	ELSE IF(PLANE_PARALLEL)THEN
	  T1=R(ND)*R(ND)*1.0E+05_LDP/SPEED_OF_LIGHT()	!As V in km/s, c in cgs units.
	  DO I=DST,DEND
	    MECH_LUM(I)=T1*V(I)*(1.0_LDP+SIGMA(I))*K_INT(I)/R(I)
	  END DO
	ELSE IF(USE_J_REL .AND. INCL_REL_TERMS)THEN
	  T1=1.0E+05_LDP/SPEED_OF_LIGHT()			!As V in km/s, c in cgs units.
	  DO I=DST,DEND
	    T2=1.0_LDP/(1.0_LDP-(T1*V(I))**2)		!Gamma^2
	    T3=SQRT(T2)					!Gamma
	    MECH_LUM(I)=T1*R(I)*V(I)*T3*( J_INT(I)-K_INT(I) +
	1          T2*(SIGMA(I)+1.0_LDP)*(K_INT(I)+T1*V(I)*RLUMST(I)/R(I)/R(I)) )
	    RLUMST(I)=T3*(RLUMST(I)+T1*V(I)*J_INT(I)*R(I)*R(I))
	  END DO
	ELSE
	  T1=1.0E+05_LDP/SPEED_OF_LIGHT()	!As V in km/s, c in cgs units.
	  DO I=DST,DEND
	    MECH_LUM(I)=T1*R(I)*V(I)*(J_INT(I)+SIGMA(I)*K_INT(I))
	  END DO
	END IF
!
	CALL WRITE_VEC(DEP_RAD_EQ_LUM,ND,'DEP_RAD_EQ_LUM',300+MYPE)
	CALL WRITE_VEC(DJDT_LUM,ND,'DJDT_LUM',300+MYPE)
	CALL WRITE_VEC(SHOCK_POWER_LUM,ND,'SHOCK_POWER_LUM',300+MYPE)
	CALL WRITE_VEC(RAD_DECAY_LUM,ND,'RAD_DECAY_LUM',300+MYPE)
	CALL WRITE_VEC(XRAY_LUM_TOT,ND,'XRAY_LUM_TOT',300+MYPE)
	FLUSH(UNIT=300+MYPE)
!
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
	CALL GATHER_SELF_VEC_MPI_V1(LLUMST,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(DIELUM,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(dE_WORK_LUM,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(DEP_RAD_EQ_LUM,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(SHOCK_POWER_LUM,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(XRAY_LUM_TOT,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(XRAY_LUM_0P1,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(XRAY_LUM_1keV,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(DJDT_LUM,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(MECH_LUM,DST,DEND,ND)
	CALL GATHER_SELF_VEC_MPI_V1(RAD_DECAY_LUM,DST,DEND,ND)
!
	IF(MYPE .EQ. 0)THEN
	  CALL LUM_FROM_ETA_V2(LLUMST,R,LUM_FROM_ETA_METHOD,ND)
	  CALL LUM_FROM_ETA_V2(DIELUM,R,LUM_FROM_ETA_METHOD,ND)
	  CALL LUM_FROM_ETA_V2(DEP_RAD_EQ_LUM,R,LUM_FROM_ETA_METHOD,ND)
	  CALL LUM_FROM_ETA_V2(dE_WORK_LUM,R,LUM_FROM_ETA_METHOD,ND)
	  CALL LUM_FROM_ETA_V2(RAD_DECAY_LUM,R,LUM_FROM_ETA_METHOD,ND)
	  CALL LUM_FROM_ETA_V2(SHOCK_POWER_LUM,R,LUM_FROM_ETA_METHOD,ND)
	  CALL LUM_FROM_ETA_V2(XRAY_LUM_TOT,R,LUM_FROM_ETA_METHOD,ND)
	  CALL LUM_FROM_ETA_V2(XRAY_LUM_0P1,R,LUM_FROM_ETA_METHOD,ND)
	  CALL LUM_FROM_ETA_V2(XRAY_LUM_1KeV,R,LUM_FROM_ETA_METHOD,ND)
	  CALL LUM_FROM_ETA(MECH_LUM,R,ND)
	  CALL LUM_FROM_ETA(DJDT_LUM,R,ND)
	  WRITE(6,*)'Done lum computation'; FLUSH(UNIT=6)
!
! Increment the continuum luminosity by the total line luminosity.
!
          T1=0.0_LDP
	  T2=0.0_LDP
	  DO I=1,ND-1
	    DO J=I,ND-1
	      RLUMST(I)=RLUMST(I)+LLUMST(J)+DIELUM(J)
	    END DO
            T1=T1+LLUMST(I)
	    T2=T2+DIELUM(I)
	  END DO
	  CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFLUX','OLD','APPEND',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening OBSFLUX to output Rec emission'
	    WRITE(LUER,*)'IOS=',IOS
	  END IF
	  IF(USE_J_REL)CALL WRITV(RLUMST,ND,'Luminosity [g.r^2.H + beta.g.r^2.J]',LU_FLUX)
	  IF(SUM(DIELUM) .NE. 0.0_LDP)CALL WRITV(DIELUM,ND,
	1    'Dielectronic and Implicit Recombination Line Emission',LU_FLUX)
	  IF(SUM(LLUMST) .NE. 0.0_LDP)CALL WRITV(LLUMST,ND,'Line Emission',LU_FLUX)
	  CALL WRITV(MECH_LUM,ND,'Mechanical Luminosity',LU_FLUX)
	  IF(SUM(de_WORK) .NE. 0.0_LDP)CALL WRITV(dE_WORK,ND,'Internal/adiabatic term',LU_FLUX)
	  IF(SUM(DJDt_LUM) .NE. 0.0_LDP)CALL WRITV(DJDt_LUM,ND,'Flux arrising from Dr^3J/Dt term',LU_FLUX)
          IF(SUM(SHOCK_POWER_LUM) .NE. 0.0_LDP)
	1      CALL WRITV(SHOCK_POWER_LUM,ND,'Energy deposited locally by the shock',LU_FLUX)
	  IF(SUM(RAD_DECAY_LUM) .NE. 0.0_LDP)
	1     CALL WRITV(RAD_DECAY_LUM,ND,'Energy deposited locally due to radioactive decay',LU_FLUX)
	  CALL WRITV(DEP_RAD_EQ_LUM,ND,'Departure from Rad Equilibrium Correction',LU_FLUX)
	  CALL WRITV(RLUMST,ND,'Total Radiative Luminosity',LU_FLUX)
	  IF(SUM(XRAY_LUM_TOT) .NE. 0.0_LDP)CALL WRITV(XRAY_LUM_TOT,ND,'Total Shock Luminosity (Lsun)',LU_FLUX)
!
! Include the machanical luminosity imparted to the wind by the radiation
! field in the total luminosity, and subtract out radioactive energy deposition..
!
! Altered: 28_Feb-2009: Changed J to  I in dE_RAD_DECAY.
! Altered: 06-Dec-2019: Normlaize in opposite directions for stars and SN.
!
	  T3=0.0_LDP
	  IF(SN_MODEL)THEN
	    DO I=1,ND-1
	      T3=T3+SHOCK_POWER_LUM(I)+RAD_DECAY_LUM(I)-MECH_LUM(I)-DJDT_LUM(I)-dE_WORK_LUM(I)
	      RLUMST(I+1)=RLUMST(I+1) + T3
	    END DO
	  ELSE
	    DO I=ND-1,1,-1
	      T3=T3+MECH_LUM(I)+DJDT_LUM(I)+dE_WORK_LUM(I)-RAD_DECAY_LUM(I)-SHOCK_POWER_LUM(I)
	      RLUMST(I)=RLUMST(I) + T3
	    END DO
	  END IF
	  CALL WRITV(RLUMST,ND,'Luminosity Check (not observed luminosity)',LU_FLUX)
!
	  IF(SN_MODEL)THEN
	    TA(1:ND)=RLUMST(1:ND)/RLUMST(2)
	    CALL WRITV(TA,ND,'Normalized luminosity check (normalized by L[2])',LU_FLUX)
	  ELSE
	    TA(1:ND)=RLUMST(1:ND)/RLUMST(ND)
	    CALL WRITV(TA,ND,'Normalized luminosity check (normalized by L[ND])',LU_FLUX)
	    IF(.NOT. USE_FIXED_J)THEN
	      DO I=ROSS_PHOT_DPTH_INDX,ND
	        IF(TA(I) .LT. 0.5_LDP .OR. TA(I) .GT. 2.0_LDP)THEN
	          WRITE(6,*)'You have poor flux accuracy in the inner region of the model'
	          WRITE(6,*)'You should restart you model changing the input options or input files'
	          WRITE(6,*)'Some normalized luminosities are out by a factor of 2'
	          WRITE(6,*)'You should check OBSFLUX'
	          IF(STOP_IF_MAJOR_WARNING)STOP
	        ELSE IF(TA(I) .LT. 0.8_LDP .OR. TA(I) .GT. 1.25_LDP)THEN
	          WRITE(6,*)'You have poor flux accuracy in the inner region of the model'
	          WRITE(6,*)'You may need to restart you model changing the input options or input files'
	          WRITE(6,*)'Some normalized luminosities are out by 25%'
	          WRITE(6,*)'You should check OBSFLUX'
	        END IF
	       END DO
	    END IF
      	  END IF             			!SN loop 
!
	  T3=0.0_LDP
	  IF(SN_MODEL)THEN
	    DO I=1,ND-1
	      T3=T3-DEP_RAD_EQ_LUM(I)
	      RLUMST(I+1)=RLUMST(I+1) + T3
	    END DO
	  ELSE
	    DO I=ND-1,1,-1
	      T3=T3+DEP_RAD_EQ_LUM(I)
	      RLUMST(I)=RLUMST(I) + T3
	    END DO
	  END IF
	  CALL WRITV(RLUMST,ND,'Consistency check (include dep. from rad. equil.)',LU_FLUX)
!
	  WRITE(LU_FLUX,'(A)')' '
	  WRITE(LU_FLUX,'(A,T60,1PE12.4)')'Total Line luminosity:',T1
	  WRITE(LU_FLUX,'(A,T60,1PE12.4)')'Total Dielectronic and Implicit Recombination Luminosity:',T2
	  WRITE(LU_FLUX,'(A,T60,1PE12.4)')'Total Mechanical Luminosity:',SUM(MECH_LUM)
	  WRITE(LU_FLUX,'(A,T60,1PE12.4)')'            Total DJDT_LUM:',SUM(DJDT_LUM)
	  WRITE(LU_FLUX,'(A,T60,1PE12.4)')'Total Rad. decay luminosity:',SUM(RAD_DECAY_LUM)
	  WRITE(LU_FLUX,'(A,T60,1PE12.4)')'              Total dE_WORK:',SUM(dE_WORK_LUM)
!
! The second XRAY flux printed is the OBSERVED XRAY luminosity. Its should be very similar
! to the earlier value for optically thin winds, and assuming that the star is not hot
! enough to act as its own source of X-rays.
!
	  T1=LUM
	  IF(SN_MODEL .AND. .NOT. USE_FIXED_J)T1=RLUMST_BND
	  WRITE(LU_FLUX,'(A,T60,ES12.4)')'Total Shock Luminosity (Lsun):',SUM(XRAY_LUM_TOT)
	  WRITE(LU_FLUX,'(A,T60,2ES12.4)')'Emitted & observed X-ray Luminosity (> 0.1 keV, Lsun) :',
	1                                     SUM(XRAY_LUM_0P1),OBS_XRAY_LUM_0P1
	  WRITE(LU_FLUX,'(A,T60,2ES12.4)')'Emitted & observed X-ray Luminosity (> 1 keV, Lsun):',
	1                                     SUM(XRAY_LUM_1KEV),OBS_XRAY_LUM_1KEV
	  WRITE(LU_FLUX,'(A,T60,2ES12.4,3X,A,ES11.4,A)')
	1                                 'Emitted & observed X-ray Luminosity (> 0.1 keV, Lstar) :',
	1                                     SUM(XRAY_LUM_0P1)/T1,OBS_XRAY_LUM_0P1/T1,
	1                                     '(Lstar[CMF]=',T1,')'
	  WRITE(LU_FLUX,'(A,T60,2ES12.4)')'Emitted & observed X-ray Luminosity (> 1 keV, Lstar):',
	1                                     SUM(XRAY_LUM_1KEV)/T1,OBS_XRAY_LUM_1KEV/T1
	  CLOSE(UNIT=LU_FLUX)
!
! Quick and dirty way of ensuring people don't take notic of OBSFLUX when USE_FIXED_J is TRUE.
!
	  IF(USE_FIXED_J)THEN
	    CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFLUX','OLD',' ',' ',IZERO,IOS)
	    WRITE(LU_FLUX,'(///,A)')'THIS FILE IS USELESS WHEN USE_FIXED_J is TRUE'
	    CLOSE(UNIT=LU_FLUX)
	  END IF
	END IF         			 ! MYPE=0 write loop
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
! Because of the use of SL's, it is possible that a error in Luminosity
! can occur at depth, particulary in SN models where line escape is
! enhanced because of the velocity field. TA provides an estimate of
! this error. This error is not relevant in the outer regions.
!
	IF(LST_ITERATION .AND. WRITE_RATES)THEN
	  DO I=DST,DEND
	    TA(I)=4.1274E-12_LDP*(STEQ_T_NO_SCL(I)-STEQ_T_SCL(I))*R(I)*R(I)
	  END DO
	  CALL MPI_REDUCE(TA,TB,ND,MPI_DOUBLE_PRECISION,MPI_SUM,IZERO,MPI_COMM_WORLD,IERR)
	  CALL LUM_FROM_ETA(TA,R,ND)
	  WRITE(LU_HT,'(//,A)')' Estimated error in L due to use of SLs'
	  WRITE(LU_HT,'(/,(1X,1P,5E12.4))')(TA(I), I=1,ND)
	END IF
!
! Output R, V and the Langrangian coordiante to the EDDACTOR file.
! If we adjust R, these values are consistent with the old grid --
! not the current grid.
!
	IF(MYPE .EQ. 0)THEN
          CALL OUT_RV_TO_EDDFACTOR(R,V,LANG_COORD,ND,
	1        REXT,VEXT,LANG_COORDEXT,NDEXT,
	1        ACCESS_F,'EDDFACTOR',LU_EDD)
	  IF(COMPUTE_EDDFAC)WRITE(LU_EDD,REC=FINISH_REC)T1
	  FLUSH(UNIT=LU_EDD)
	END IF
!
! Insure Eddington factor file is closed, and indicate all f's successfully
! computed. If COMPUTE_EDDFAC is true we can safely write to record 5
! as file must have new format.
!
	T1=1.0_LDP
	IF(.NOT. COHERENT_ES)CLOSE(UNIT=LU_ES)
	COMPUTE_JEW=.FALSE.
	COMPUTE_EDDFAC=.FALSE.
!
! IF we are doing a FLUX computation only we do not corrupt the SCTRMEP file
! or the BA matrix. We also do not output th populations. This can be done
! quickly by setting FLUX_CAL_ONLY=.FALSE. and putting N_ITS=0.
!
	IF(FLUX_CAL_ONLY .AND. RD_COHERENT_ES)THEN
	   WRITE(LUER,*)'Stopping CMFGEN as finished FLUX calculation.'
	   WRITE(LUER,*)'For a FLUX calculation we do 1 iteration only'
	   STOP
!
! Compute the convolution of J with the electron redistribution function.
! May need to compute a flux spectrum sveral times in order for e.s.
! redistributon to be correctly allowed for.
! RD_NU and ALLOW_UNEQUAL_FREQ are both set to FALSE. Note that TEXT and NDEXT
! contain T and ND when ACCURATE is FALSE.
!
	 ELSE IF(FLUX_CAL_ONLY .AND. .NOT. RD_COHERENT_ES)THEN
	   COHERENT_ES=RD_COHERENT_ES
	   I=SIZE(VJ)
	   CALL COMP_J_CONV_V2(VJ,I,NU,TEXT,NDEXT,NCF,LU_EDD,'EDDFACTOR',
	1             EDD_CONT_REC,L_FALSE,L_FALSE,LU_ES,'ES_J_CONV')
!
! Close units 2 and 16 to force writing of information.
!
	   CLOSE(UNIT=LUER)
	   CLOSE(UNIT=LU_SE)
	   CALL GEN_ASCI_OPEN(LUER,'OUTGEN','OLD','APPEND',' ',IZERO,IOS)
	   CALL GEN_ASCI_OPEN(LU_SE,'STEQ_VALS','OLD','APPEND',' ',IZERO,IOS)
!
! Now do another iteration. We continue to iterate so that the accuracy of
! RJ_ES is improved (i.e. to allow for multiple scattering). The number of
! iterations is set by NUM_ITS_TO_DO in the input file.
!
	   IF(.NOT. LST_ITERATION)GOTO 20000
	   STOP
	END IF
!
	CALL WR2D_MPI_V1(STEQ_T,IONE,DST,DEND,ND,'Radiative Equlibrium Equation','&',L_TRUE,LU_SE)
	IF(USE_ELEC_HEAT_BAL .OR. COMP_STEQ_T_EHB)THEN
	  DO ID=1,NION
	    STEQ_T_EHB=STEQ_T_EHB+SE(ID)%T_EHB
	  END DO
	  CALL WR2D_MPI_V1(STEQ_T_EHB,IONE,DST,DEND,ND,'Electron Energy Balance Equation','%',L_TRUE,LU_SE)
	  IF(LST_ITERATION .AND. VERBOSE_OUTPUT)THEN
	    CALL WR2D_MPI_V1(STEQ_T_EHB,IONE,DST,DEND,ND,'After BF update','%',L_TRUE,LU_T_EHB)
	  END IF
	ELSE
	  STEQ_T_EHB=0.0_LDP
	  CALL WR2D_MPI_V1(STEQ_T_EHB,IONE,DST,DEND,ND,'Electron Energy Balance Equation (not computed)','%',L_TRUE,LU_SE)
	END IF
!
! 
! Reread in BA array if we are not to compute it. There should be no problem
! with this read since it has previously been read in. This double reading
! is necessary to save a rewrite of the STEQ*** routines.
!
	IF(.NOT. COMPUTE_BA)THEN
	  CALL READ_BA_DATA_MPI_V1(LU_BA,NION,NUM_BNDS,COMPUTE_BA,FIXED_T,SUCCESS,'BAMAT')
	  IF(.NOT. SUCCESS)THEN
	    WRITE(LUER,*)'Major Error - cant read BA File'
	    WRITE(LUER,*)'Previously read successfully - before continuum loop'
	    STOP
	  END IF
	END IF
!
	CALL TUNE(IONE,'SOLVE_FOR_POPS')
	CALL SOLVE_FOR_POPS_MPI_V1(POPS,NT,NION,ND,NC,NP,NUM_BNDS,DIAG_INDX,
	1      MAXCH,MAIN_COUNTER,IREC,LU_SE,LUSCR,LST_ITERATION)
	IF(LAST_NG .EQ. MAIN_COUNTER)NUM_ITS_TO_DO=MAX(2,NUM_ITS_TO_DO)
	CALL TUNE(ITWO,'SOLVE_FOR_POPS')
!
! If we have changed the R grid, we need to recompute the angular quadrature weitghts,
! and put the atom density ect on the new radius grid.
!
	IF(REVISE_R_GRID .AND. R_GRID_REVISED)THEN
	  CALL SET_ANG_QW_V2(R,NC,ND,NP,REXT,NCEXT,NDEXT,NPEXT,
	1                 R_PNT_SRCE,NC_PNT_SRCE,TRAPFORJ,ACCURATE)
!
! Compute CLUM_FAC(1:ND) which allow for the possibility that the wind is
! clumped. At the sime time, we compute the vectors which give the density,
! the atom density, and the species density at each depth.
! The new call replaces the interpolation done in
!				  CALL ADJUST_DEN_VECS(R_OLD,ND)
!
	  CALL SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
!
! Need to update the electron non-thermal spectrum on the new grid.
!
	  IF(TREAT_NON_THERMAL_ELECTRONS)NT_ITERATION_COUNTER=0
	END IF
!
! 
!
! Store populations back into individual arrays. At present, program stops
! if an error occurs in the NG acceleration. Could be changed by doing
! the conversion below before the NG call. If the NG accelerate worked,
! we would do the conversion again. IF it failed, we would do the
! reverse conversion (as POPS might be corrupted).
!
	DO ID=1,NUM_IONS-1
	  CALL POPTOION(POPS, ATM(ID)%XzV, ATM(ID)%DXzV,ED,T,
	1    ATM(ID)%EQXzV, ATM(ID)%NXzV, NT, DST, DEND, ND, ATM(ID)%XzV_PRES)
	END DO
!
	IF(MYPE .EQ. 0)THEN
	  DO ID=1,NUM_IONS-1
	    CALL POPTOION(POPS, ROOT(ID)%XzV, ROOT(ID)%DXzV,ED,T,
	1      ATM(ID)%EQXzV, ATM(ID)%NXzV, NT, IONE, ND, ND, ATM(ID)%XzV_PRES)
	  END DO
	END IF
!
! We have now store the revised populations back in their individual
! storage locations. For some species we have 2 atomic models. For these
! species we need to take the super-level populations and compute:
!
! 1. The LTE population off all level ls in the FULL atom.
! 2. The population off all levels in the FULL atom.
! 3. The LTE population off all super-levels.
!
! This is done by the include file SUP_TO_FULL which calls the subroutine
! SUP_TO_FULL.FOR
!
! We first, however, need to compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	DO J=1,ND
	  POPION(J)=0.0_LDP
	  DO I=1,NT
	    IF(Z_POP(I) .GT. 0.01_LDP)POPION(J)=POPION(J)+POPS(I,J)
	  END DO
	END DO
!
! Revise vector constants for evaluating the level dissolution. These
! constants are the same for all species. These are stored in a common block,
! and are required by SUP_TO_FULL and LTE_POP_WLD.
!
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
	CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!	INCLUDE 'SUP_TO_FULL_V4.INC'
!
! Perform a revision to th hydrostatic structure of the atmosphere. If we have changed the
! R grid, we need to recomput the angular quadrature weitghts, and put the atom density etc
! on the new radius grid.
!
	DONE_HYDRO_REVISION=.FALSE.
	IF(DO_HYDRO)THEN
	  WRITE(6,*)'Need to fix DO_CMF_HYDRO_V2'
	  FLUSH(UNIT=6)
	  CALL MPI_FINALIZE (ierr)
	  STOP
	  K=MAIN_COUNTER
!	  CALL DO_CMF_HYDRO_V2(POPS,LUM,TEFF,LOGG,STARS_MASS,RP,RMAX,RMDOT,VINF,V_BETA1,
!	1          PRESSURE_VTURB,PLANE_PARALLEL,PLANE_PARALLEL_NO_V,
!	1          K,DONE_HYDRO_REVISION,NC,ND,NP,NT)
!	  IF(DON_HYDRO_REVISION)
!	    MAIN_COUNTER=MAIN_COUNTER+1
!	    CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
!	1                LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
!	  END IF
	END IF
	IF(.TRUE. .AND. .NOT. DONE_HYDRO_REVISION)THEN
	  R_OLD(1:ND)=R(1:ND)
	  CALL DO_WIND_VEL_V1(R,V,SIGMA,CLUMP_FAC,STARS_MASS,CHK,MAIN_COUNTER,ND)
	END IF
!
	IF(DONE_HYDRO_REVISION .OR. CHK)THEN
	    LAMBDA_ITERATION=.TRUE.
            FIXED_T=.TRUE.
	    FIX_IMPURITY=.FALSE.
	    COMPUTE_BA=.TRUE.
	    MAIN_COUNTER=MAIN_COUNTER+1
	    IF(LST_ITERATION .AND. WRITE_RATES)THEN
	      LST_ITERATION=.FALSE.
	      CLOSE(LU_NET); CLOSE(LU_DR); CLOSE(LU_EW)
	      CLOSE(LU_HT); CLOSE(LU_NEG)
	    END IF
	    IF(ACCURATE)THEN
	      I=ND-DEEP
	      CALL REXT_COEF_V2(REXT,COEF,INDX,NDEXT,R,POS_IN_NEW_GRID,
	1              ND,NPINS,L_TRUE,I,ST_INTERP_INDX,END_INTERP_INDX)
	      TA(1:ND)=1.0_LDP	!TEXT not required, T currently zero
	      CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,V,TA,SIGMA,ND)
              VDOP_VEC_EXT(1:NDEXT)=12.85_LDP*SQRT( TDOP/AMASS_DOP + (VTURB/12.85_LDP)**2 )
	    END IF
	    CALL SET_ANG_QW_V2(R,NC,ND,NP,REXT,NCEXT,NDEXT,NPEXT,
	1                      R_PNT_SRCE,NC_PNT_SRCE,TRAPFORJ,ACCURATE)
	    CALL SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
	    IF(DONE_HYDRO_REVISION)THEN
	       WRITE(6,*)'Need to edit GREY_T_ITERATE'
	       FLUSH(UNIT=6)
	       CALL MPI_FINALIZE (ierr)
	       STOP
!	      CALL GREY_T_ITERATE(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
!	1            LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
	    ELSE
	       WRITE(6,*)'Need to edit WIND_SCAL_POPS'
	       FLUSH(UNIT=6)
	       CALL MPI_FINALIZE (ierr)
	       STOP
!	      CALL WIND_SCALE_POPS_V1(POPS,R_OLD,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
	    END IF
	    CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1                LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	END IF
!
	IF(DO_CLUMP_MODEL .AND. MAXCH .LT. 50.0_LDP .AND. .NOT. FIXED_T .AND.
	1         LAST_LAMBDA .NE. MAIN_COUNTER)THEN
	  CALL AUTO_CLUMP_REV(POPS,CLUMP_LAW,CLUMP_PAR,N_CLUMP_PAR,CHK,ND,NT,LUIN)
	  IF(CHK)THEN
	    CALL SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
	    CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1                LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	    LAMBDA_ITERATION=.TRUE.
            FIXED_T=.TRUE.
	    RD_FIX_T=.TRUE.
	    FIX_IMPURITY=.FALSE.
	    COMPUTE_BA=.TRUE.
	    MAIN_COUNTER=MAIN_COUNTER+1
	    MAXCH=200
	    IF(LST_ITERATION .AND. WRITE_RATES)THEN
	      LST_ITERATION=.FALSE.
	      CLOSE(LU_NET); CLOSE(LU_DR); CLOSE(LU_EW)
	      CLOSE(LU_HT); CLOSE(LU_NEG)
	    END IF
	  END IF
	END IF
!
! Initialize pointer file for storage of BA matrix.
!
	I=-1000
	CALL INIT_BA_DATA_PNT_MPI_V1(LU_BA,NION,NUM_BNDS,COMPUTE_BA,FIXED_T,'BAMAT')
!
! 
!
! Compute the convolution of J with the electron redistribution function.
! Fist neeed to UPDATE TEXT because of the linearization. We need TEXT for
! convolving J with the electron scattering redistribution function.
! VEXT and SIGMAEXT have already been computed.
!
	TEXT(1:ND)=T(1:ND)
	IF(ACCURATE)THEN
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1        V,T,SIGMA,ND)
	END IF
	COHERENT_ES=RD_COHERENT_ES
	IF(.NOT. COHERENT_ES)THEN
	  I=SIZE(VJ)
	  CALL COMP_J_CONV_V2(VJ,I,NU,TEXT,NDEXT,NCF,LU_EDD,'EDDFACTOR',
	1             EDD_CONT_REC,L_FALSE,L_FALSE,LU_ES,'ES_J_CONV')
	END IF
!
! 
! Output brief summary of the model. This is to facilate the creation
! of compact model logs.
!
	IF(LST_ITERATION .AND. MYPE .EQ. 0)THEN
!
	  CALL GEN_ASCI_OPEN(LUMOD,'MOD_SUM','UNKNOWN',' ',' ',IZERO,IOS)
!
	  WRITE(LUMOD,'(/,''Model Started on:'',15X,(A))')TIME
	  CALL DATE_TIME(TIME)
	  WRITE(LUMOD,'(''Model Finalized on:'',13X,(A))')TIME
	  WRITE(LUMOD,'(''Main program last changed on:'',3X,(A))')PRODATE
	  WRITE(LUMOD,'()')
!
	  STRING=' '
	  NEXT_LOC=1
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'ND',ND)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NC',NC)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NP',NP)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NT',NT)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
	  STRING=' '
	  NEXT_LOC=1
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NUM_BNDS',NUM_BNDS)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NCF',NCF)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NLINES',N_LINE_FREQ)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
	  WRITE(LUMOD,'(A)')' '
!
! Output brief summary of atomic models.
!
	  DO ISPEC=1,NUM_SPECIES
	    STRING=' '
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	      IF(ATM(ID)%XzV_PRES)THEN
	       CALL WR_SL_INFO(STRING,ATM(ID)%NXzV,ATM(ID)%NXzV_F,
	1                        ATM(ID)%ZXzV,ION_ID(ID),LUMOD)
	      END IF
	    END DO
	    IF(STRING .NE. ' ')WRITE(LUMOD,'(A)')TRIM(STRING)
	  END DO
	  WRITE(LUMOD,'(A)')' '
!
! Output stellar parameters.
!
	  STRING=' '
	  NEXT_LOC=1
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'L*',LUM)
	  T1=RMDOT/3.02286E+23_LDP
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Mdot',T1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'R*  ',RP)
	  T1=RMAX/RP   ; CALL WR_VAL_INFO(STRING,NEXT_LOC,'RMAX/R*',T1)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
! Compute the Radius and Velocity at Tau=10, and at Tau=2/3.
!
	  TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
!          CALL WR2D_V2(CLUMP_FAC,IONE,ND,'YCHI','#',L_TRUE,LUER)
!          CALL WR2D_V2(ROSS_MEAN,IONE,ND,'YCHI','#',L_TRUE,LUER)
!          CALL WR2D_V2(TCHI,IONE,ND,'YCHI','#',L_TRUE,LUER)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	  TA(1:ND)=0.0_LDP ; DO I=2,ND ; TA(I) = TA(I-1)+DTAU(I-1) ; END DO
	  TB(1)=MIN(2.0_LDP/3.0_LDP,TA(ND))  ; TB(2)=MIN(10.0_LDP,TA(ND))
	  TB(3)=MIN(20.0_LDP,TA(ND))
	  CALL MON_INTERP(TC,ITHREE,IONE,TB,ITHREE,R,ND,TA,ND)
	  CALL MON_INTERP(AV,ITHREE,IONE,TB,ITHREE,V,ND,TA,ND)
!
	  LUM_FOR_TEFF=LUM
	  IF(SN_MODEL)LUM_FOR_TEFF=RLUMST_BND
!
! Output summary of Teff, R, and V at RSTAR, Tau=10, and TAU=2/3.
! For a plane-parallel atmosphere, R(ND) defines Teff.
!
	  NEXT_LOC=1  ;   STRING=' '
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TA(ND))
	  T1=1.0E+10_LDP*RP/RAD_SUN() ; CALL WR_VAL_INFO(STRING,NEXT_LOC,'R*/Rsun',T1)
	  T1=TEFF_SUN()*(ABS(LUM_FOR_TEFF)/T1**2)**0.25_LDP					!ABS for SN
	  IF(PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL)THEN
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff(K)',T1)
	  ELSE
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'T*(K)',T1)
	  END IF
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',V(ND))
	  IF(DO_HYDRO)THEN
	    T1=LOG10(1.0E-20_LDP*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()/RP/RP)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Log g',T1)
	  END IF
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
	  IF(TA(ND) .GT. 20.0_LDP .AND. .NOT. (PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL) )THEN
	    NEXT_LOC=1  ;   STRING=' '
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(3))		!20.0D0
	    T1=1.0E+10_LDP*TC(3)/RAD_SUN()
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	    T1=TEFF_SUN()*(ABS(LUM_FOR_TEFF)/T1**2)**0.25_LDP
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff(K)',T1)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(3))
	    IF(DO_HYDRO)THEN
	      T1=LOG10(1.0E-20_LDP*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()/TC(3)/TC(3))
	      CALL WR_VAL_INFO(STRING,NEXT_LOC,'Log g',T1)
	    END IF
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	  END IF
!
	  IF(TA(ND) .GT. 10.0_LDP .AND. .NOT. (PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL) )THEN
	    NEXT_LOC=1  ;   STRING=' '
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(2))		!10.0D0
	    T1=1.0E+10_LDP*TC(2)/RAD_SUN()
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	    T1=TEFF_SUN()*(ABS(LUM_FOR_TEFF)/T1**2)**0.25_LDP
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff(K)',T1)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(2))
	    IF(DO_HYDRO)THEN
	      T1=LOG10(1.0E-20_LDP*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()/TC(2)/TC(2))
	      CALL WR_VAL_INFO(STRING,NEXT_LOC,'Log g',T1)
	    END IF
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	  END IF
!
	  IF(TA(ND) .GT. 0.67_LDP .AND. .NOT. (PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL) )THEN
	    NEXT_LOC=1  ;   STRING=' '
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(1))		!0.67D0
	    T1=1.0E+10_LDP*TC(1)/RAD_SUN()
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	    T1=TEFF_SUN()*(ABS(LUM_FOR_TEFF)/T1**2)**0.25_LDP
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff(K)',T1)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(1))
	    IF(DO_HYDRO)THEN
	      T1=LOG10(1.0E-20_LDP*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()/TC(1)/TC(1))
	      CALL WR_VAL_INFO(STRING,NEXT_LOC,'Log g',T1)
	    END IF
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	  END IF
!
	  STRING=' '
	  NEXT_LOC=1
	  T1=4.9376E+07_LDP*(RMDOT/3.02286E+23_LDP)*V(1)/LUM
	  T2=8.235E+03_LDP*(RMDOT/3.02286E+23_LDP)*V(1)*V(1)/LUM
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Eta',T1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Ek/L(%)',T2)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
	  WRITE(LUMOD,'(A)')' '
!
! Velocity law information.
!
	  NEXT_LOC=1  ;   STRING=' '
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Vinf1',VINF1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Beta1',V_BETA1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'SCL_HT/RP',SCL_HT)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
	  NEXT_LOC=1  ;   STRING=' '
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'VCORE',VCORE)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'VPHOT',VPHOT)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
	  IF(VELTYPE .EQ. 6)THEN
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Vinf2',VINF2)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Beta2',V_BETA2)
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	    WRITE(LUMOD,'(A)')' '
	  END IF
!
! Output abundance information.
!
	  DO ISPEC=1,NUM_SPECIES
	    CALL WR_ABUND_INFO_V3(SPECIES(ISPEC),AT_MASS(ISPEC),
	1           AT_ABUND(ISPEC),ABUND_SUM,MEAN_ATOMIC_WEIGHT,
	1           SOL_MASS_FRAC(ISPEC),SOL_ABUND_REF_SET,LUMOD)
	  END DO
!
	  WRITE(LUMOD,'(A)')' '
	  IF(DO_CLUMP_MODEL)THEN
	    WRITE(LUMOD,'(A,A)')
	1            'Running clumped model: ',TRIM(CLUMP_LAW)
	    WRITE(LUMOD,'(A,1PE10.3)')
	1            'Filling factor at boundary is: ',CLUMP_FAC(1)
	    STRING=' '
	    NEXT_LOC=1
	    DO I=1,N_CLUMP_PAR
	      TEMP_CHAR='CL_P_'
	      WRITE(TEMP_CHAR(6:6),'(I1)')I
	      CALL WR_VAL_INFO(STRING,NEXT_LOC,TEMP_CHAR,CLUMP_PAR(I))
	      IF(NEXT_LOC .GT. 80)THEN
	        WRITE(LUMOD,'(A)')TRIM(STRING)
	        STRING=' '
	        NEXT_LOC=1
	      END IF
	    END DO
	    IF(STRING .NE. ' ')WRITE(LUMOD,'(A)')TRIM(STRING)
	    WRITE(LUMOD,'(A)')' '
	  END IF
!
	  WRITE(LUMOD,'(A,1PE10.3)')
	1            'Maximum correcion (%) on last iteration: ',MAXCH
	  CLOSE(LUMOD)
!
	END IF
!
! Check to see if corrections are reasonable.
!
	IF(MAXCH .GT. MAX_CHNG_LIM)THEN
	  WRITE(LUER,*)'Error - bad initial population guesses.'
	  WRITE(LUER,*)'Predicted changes are too large. '
	  WRITE(LUER,*)'New populations written to SCRTEMP file.'
	  WRITE(LUER,*)'Edit POINT1 file to recover older populations.'
	  STOP
	END IF
