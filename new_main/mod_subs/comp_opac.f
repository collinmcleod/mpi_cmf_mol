!
! Subroutine to compute the opacities and emissivities for CMFGEN.
!
	SUBROUTINE COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
	USE SET_KIND_MODULE
	USE MPI
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	USE MOD_J_TWO_PHOT_MPI_V1, ONLY : GET_J_FOR_TWO_PHOT
	USE OPAC_MOD
	IMPLICIT NONE
!
! Altered 26-Jun-2020 : Fixed bug related to Rayleigh scattering when H- is present.
! Altered 05-Apr-2011 : Now call GENOPAETA_V10 (6-Feb-2011)
! Altered 11-Jun-2006: Installed CHI_NOSCAT and ETA_NOSCAT. Scattering
!                        opacity now computed after bound-free, free-free,
!                        two photon, and X-ray opacity have been computed.
! Altered 11-Jun-2005: Bug fix: X-ray opacity/emissivity was not being
!                        added correctly.
! Altered 13-Sep-2004: Installed ETA_MECH to keep track of mechanical
!                        energy input.
! Altered 03-Mar-2004: Bug fix. EMHNUKT was not being computed when
!                        COMPUTE_NEW_CROSS=.FALSE. and NU=NU_CONT.
!                        Value at earlier frequncy was being used.
!                        Computation ox X-ray cross-sections now
!                        directly included (not INCLUDE file).
! Created 16-Feb-2004: Based on OPACITIES_V4
!
	INTEGER ND
	INTEGER NT
	INTEGER NCF
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) NU_EVAL_CONT(NCF)
	REAL(KIND=LDP) FQW(NCF)
	REAL(KIND=LDP) FL
	REAL(KIND=LDP) CONT_FREQ
!
	INTEGER FREQ_INDX
	CHARACTER*(*) SECTION
	LOGICAL LST_DEPTH_ONLY
!
! Constants for opacity etc.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Internally used variables
!
        REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
	REAL(KIND=LDP) XCROSS_V2
	REAL(KIND=LDP) GFF
	EXTERNAL XCROSS_V2,GFF
!
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) T1,T2,T3,T4
	INTEGER I,J
	INTEGER ID
	INTEGER MDST
	INTEGER PHOT_ID
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!	include 'mpif.h'
!
	CHI(1:ND)=0.0_LDP
	ETA(1:ND)=0.0_LDP
	MDST=DST
	IF(LST_DEPTH_ONLY)MDST=DEND
!
! Compute opacity and emissivity. This is a general include file
! provided program uses exactly the same variables. Can be achieved
! by copying declaration statements from CMFGEN. Always advisable
! to use ``IMPLICIT NONE''.
!
	IF(COMPUTE_NEW_CROSS)THEN
!
! Compute EXP(-hv/kT) and zero CHI, and ETA.
!
	  T1=-HDKT*CONT_FREQ
	  DO I=1,ND
	    EMHNUKT_CONT(I)=EXP(T1/T(I))
	    CHI_RAY(I)=0.0_LDP
	    CHI_SCAT(I)=0.0_LDP
	    ESEC(I)=0.0_LDP
	  END DO
	  CHI_C_EVAL=0.0_LDP; ETA_C_EVAL=0.0_LDP
!
! Compute continuum intensity incident from the core assuming a TSTAR
! blackbody.
!
          T1=EXP(-HDKT*CONT_FREQ/TSTAR)
	  IC=TWOHCSQ*T1*(CONT_FREQ**3)/(1.0_LDP-T1)
!
! Free-free and bound-free opacities.
!
	  CALL TUNE(1,'GEN_OPA_ETA')
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        PHOT_ID=J
	        CALL GENOPAETA_MPI_V1(ID,CHI,ETA,CONT_FREQ,
	1           ATM(ID)%XzV_F,      ATM(ID)%XzVLTE_F,     ATM(ID)%LOG_XzVLTE_F,  ATM(ID)%EDGEXzV_F,
	1           ATM(ID)%GIONXzV_F,  ATM(ID)%ZXzV,         ATM(ID)%NXzV_F,
	1           ATM(ID+1)%XzV,      ATM(ID+1)%LOG_XzVLTE, ATM(ID+1)%NXzV,
	1           PHOT_ID,            ATM(ID)%XzV_ION_LEV_ID(J),
	1           ED,T,EMHNUKT_CONT,L_TRUE,DST, DEND, ND,LST_DEPTH_ONLY)
	      END DO
	    END IF
	  END DO
	  CALL TUNE(2,'GEN_OPA_ETA')
!
	  IF(LST_DEPTH_ONLY .AND. DEND .NE. ND)RETURN
!
	  IF(ADD_ADDITIONAL_OPACITY)THEN
	     DO I=MDST,DEND
	       T1=ADD_OPAC_SCL_FAC*6.65E-15_LDP*POP_ATOM(I)
	       CHI(I)=CHI(I)+T1
	       ETA(I)=ETA(I)+T1*TWOHCSQ*(CONT_FREQ**3)*EMHNUKT_CONT(I)/(1.0_LDP-EMHNUKT_CONT(I))
	     END DO
	     IF(FREQ_INDX .EQ. 1)THEN
	       LUER=ERROR_LU()
	       WRITE(LUER,*)'Warning --- adding additional opacity to model'
	       WRITE(LUER,*)'Remember to remove additional opacity for final model'
	     END IF
	  END IF
!
! 
!
! Add in 2-photon emissivity and opacity.
!
	  IF(LST_DEPTH_ONLY)THEN
	    CALL TWO_PHOT_OPAC_MPI_V1(ETA,CHI,POPS,T,CONT_FREQ,'LTE',DST,DEND,ND,NT)
	  ELSE IF(COMPUTE_EDDFAC .AND. TWO_PHOTON_METHOD .EQ. 'USE_RAD')THEN
	    CALL TWO_PHOT_OPAC_MPI_V1(ETA,CHI,POPS,T,CONT_FREQ,'OLD_DEFAULT',DST,DEND,ND,NT)
	  ELSE
	    CALL TWO_PHOT_OPAC_MPI_V1(ETA,CHI,POPS,T,CONT_FREQ,TWO_PHOTON_METHOD,DST,DEND,ND,NT)
	  END IF
!
! Compute X-ray opacities and emissivities due to K (& L) shell ionization. In all cases
! it is assumed that 2 electrons are ejected. The K (& L) shell cross-sections are
! assumed to be independent of the level of the valence electron. In practice,
! ionizations will generally be determined by the population of the ground
! configuration.
!
! Since the cross-sections are level independent, we can sum over the levels before
! we add the contribution to the opacity/emissivity.
!
! This section was originaly XOPAC_V4.INC. See that file for erlier corrections.
!
	  CALL TUNE(1,'XRAY_FST_OPAC')
	  IF(XRAYS)THEN
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	        T2=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV
	        T1=XCROSS_V2(CONT_FREQ,AT_NO(SPECIES_LNK(ID)),T2,IZERO,IZERO,L_FALSE,L_FALSE)
	        IF(T1 .NE. 0.0_LDP)THEN
	          J=1
	          DO I=MDST,DEND
	            T2=0.0_LDP			!Temporary CHI
	            T3=0.0_LDP			!Temporary ETA
	            T4=(ATM(ID+1)%XzVLTE_F(1,I)*EMHNUKT_CONT(I))/ATM(ID+1)%XzV_F(1,I)
	            DO J=1,ATM(ID)%NXzV_F
		      T2=T2+ATM(ID)%XzV_F(J,I)
	              T3=T3+ATM(ID)%XzVLTE_F(J,I)
	            END DO
	            CHI(I)=CHI(I)+T1*(T2-T3*T4)
	            ETA(I)=ETA(I)+T1*T3*T4*TWOHCSQ*(CONT_FREQ**3)
	          END DO
	        END IF
	      END IF
	    END DO
	  END IF
	  CALL TUNE(2,'XRAY_FST_OPAC')
!
! Compute scattering opacity. ESEC is zeroed in ESOPAC.
!
	  CALL TUNE(1,'NEW_COMP_OAPC')
	  CALL ESOPAC(ESEC,ED,ND)		!Electron scattering emission factor.
!
! Add in Rayleigh scattering contribution.
!
	  CHI_RAY(MDST:DEND)=0.0_LDP
	  ID=1; IF(ION_ID(1) .EQ. 'HMI')ID=2
	  IF(SPECIES_PRES(1) .AND. INCL_RAY_SCAT)THEN
	    CALL RAYLEIGH_SCAT_MPI_V1(CHI_RAY,ATM(ID)%XzV_F,ATM(ID)%AXzV_F,ATM(ID)%EDGEXZV_F,
	1             ATM(1)%NXzV_F,CONT_FREQ,DST,DEND,ND)
	  END IF
!
	  CHI_NOSCAT(MDST:DEND)=CHI(MDST:DEND)
	  ETA_NOSCAT(MDST:DEND)=ETA(MDST:DEND)
	  CHI_SCAT(MDST:DEND)=ESEC(MDST:DEND)+CHI_RAY(MDST:DEND)
!
! Now compute total opacity --- scattering + non scattering.
!
	  CHI(MDST:DEND)=CHI(MDST:DEND)+CHI_SCAT(MDST:DEND)
!
	  IF(.NOT. LST_DEPTH_ONLY)THEN
	    CALL TUNE(1,'ALL_RED_COMP')
	    CHI_C_EVAL(DST:DEND)=CHI(DST:DEND)
	    CALL MPI_ALLREDUCE(MPI_IN_PLACE,CHI_C_EVAL,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	    CALL MPI_ERR_CHECK(IERR,'MPI_ALLREDUCE -- COMP_OPAC')
	    ETA_C_EVAL(DST:DEND)=ETA(DST:DEND)
	    CALL MPI_ALLREDUCE(MPI_IN_PLACE,ETA_C_EVAL,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	    CALL MPI_ERR_CHECK(IERR,'MPI_ALLREDUCE -- ETA_C_EVAL  -- COMP_OPAC')
	    CALL TUNE(2,'ALL_RED_COMP')
	    CALL TUNE(1,'IN_RED_COMP')
	    CALL MPI_ALLREDUCE(MPI_IN_PLACE,CHI_RAY,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	    CALL MPI_ERR_CHECK(IERR,'MPI_ALLREDUCE -- CHI_RAY -- COMP_OPAC')
	    CALL MPI_ALLREDUCE(MPI_IN_PLACE,CHI_SCAT,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
	    CALL MPI_ERR_CHECK(IERR,'MPI_ALLREDUCE -- SCAT -- COMP_OPAC')
	    CALL TUNE(2,'IN_RED_COMP')
	  ElSE
	    CHI_C_EVAL(ND)=CHI(ND)
	    ETA_C_EVAL(ND)=ETA(ND)
	  END IF
!	  
	  CHI_NOSCAT_EVAL(:)=CHI_NOSCAT(:)
	  ETA_NOSCAT_EVAL(:)=ETA_NOSCAT(:)
	  CALL TUNE(2,'NEW_COMP_OAPC')
!
	ELSE
	  IF(LST_DEPTH_ONLY .AND. DEND .NE. ND)RETURN
	END IF
!
! 
!
! Evaluate EXP(-hv/kT) for current frequency. This is needed by routines
! such as COMP_VAR_JREC etc.
!
	  DO J=MDST,ND
	    EMHNUKT(J)=EXP(-HDKT*FL/T(J))
	  END DO
!
! Section to revise continuum opacities etc so that they are computed at
! the correct frequency. We have stored the orginal continuum opacity and
! emissivity in CHI_C_EVAL and ETA_C_EVAL, which were computed at CONT_FREQ.
!
	IF(FL .NE. CONT_FREQ)THEN
!
! Compute continuum intensity incident from the core assuming a TSTAR
! blackbody.
!
          T1=EXP(-HDKT*FL/TSTAR)
	  IC=TWOHCSQ*T1*(FL**3)/(1.0_LDP-T1)
!
! We assume that the photoionization cross-section has not changed since the
! last iteration. Using the result that the stimulated emission occurs in
! LTE and is given by
!                     ETA/(2hv^3/c^2)
! we can adjust CHI and ETA so that the condition of constant photoionization
! cross-section is met. This adjustment automatically ensures that ETA/CHI
! gives the Planck function in LTE.
!
	  T1=(FL/CONT_FREQ)**3
	  T2=TWOHCSQ*(CONT_FREQ**3)
	  T3=TWOHCSQ*(FL**3)
	  DO J=MDST,DEND
	    T4=ETA_C_EVAL(J)*T1*EXP(-HDKT*(FL-CONT_FREQ)/T(J))
	    CHI(J)=CHI_C_EVAL(J)+(ETA_C_EVAL(J)/T2-T4/T3)
	    ETA(J)=T4
	    CHI_NOSCAT(J)=CHI_NOSCAT_EVAL(J)+(ETA_C_EVAL(J)/T2-T4/T3)
	    ETA_NOSCAT(J)=T4
	  END DO
!
	ELSE
!
! We reset CHI and ETA in case shock X-ray emission has been added to ETA,
! or CONT_FREQ was not the first frequency.
!
	  CHI(MDST:DEND)=CHI_C_EVAL(MDST:DEND)
	  ETA(MDST:DEND)=ETA_C_EVAL(MDST:DEND)
	  CHI_NOSCAT(MDST:DEND)=CHI_NOSCAT_EVAL(MDST:DEND)
	  ETA_NOSCAT(MDST:DEND)=ETA_NOSCAT_EVAL(MDST:DEND)
	END IF
!
! 
!
! The shock emission is added separately since it does not occur at the
! local electron temperature.
!
	CALL TUNE(1,'XRAY_COMP_OPAC')
	IF(XRAYS)THEN
!
	  ZETA=0.0_LDP
	  IF(FF_XRAYS)THEN
!
! Since T_SHOCK is depth indpendent, Z^2 * (the free-free Gaunt factors)
! are depth independent.
!
! We use T3 for the Electron density. We asume H, He, and C are fully ionized
! in the X-ray emitting plasma. All other species are assumed to have Z=6.0
!
	    IF(T_SHOCK_1 .NE. 0.0_LDP)THEN
	      T1=CHIFF*TWOHCSQ*(FILL_FAC_XRAYS_1**2)*
	1                    EXP(-HDKT*CONT_FREQ/T_SHOCK_1)/SQRT(T_SHOCK_1)
	      T2=1.0_LDP ; TA(1)=GFF(CONT_FREQ,T_SHOCK_1,T2)
	      T2=2.0_LDP ; TA(2)=4.0_LDP*GFF(CONT_FREQ,T_SHOCK_1,T2)
	      T2=6.0_LDP ; TA(3)=36.0_LDP*GFF(CONT_FREQ,T_SHOCK_1,T2)
	      DO I=MDST,DEND
	        T2=TA(1)*POP_SPECIES(I,1)+TA(2)*POP_SPECIES(I,2) +
	1         TA(3)*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        T3=POP_SPECIES(I,1)+POP_SPECIES(I,2) +
	1          6.0_LDP*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        ZETA(I)=T1*T2*T3*EXP(-V_SHOCK_1/V(I))	 !Zeta is temporary
              END DO
	    END IF
	    IF(T_SHOCK_2 .NE. 0.0_LDP)THEN
	      T1=CHIFF*TWOHCSQ*(FILL_FAC_XRAYS_2**2)*
	1                    EXP(-HDKT*CONT_FREQ/T_SHOCK_2)/SQRT(T_SHOCK_2)
	      T2=1.0_LDP ; TA(1)=GFF(CONT_FREQ,T_SHOCK_2,T2)
	      T2=2.0_LDP ; TA(2)=4.0_LDP*GFF(CONT_FREQ,T_SHOCK_2,T2)
	      T2=6.0_LDP ; TA(3)=36.0_LDP*GFF(CONT_FREQ,T_SHOCK_2,T2)
	      DO I=MDST,DEND
	        T2=TA(1)*POP_SPECIES(I,1)+TA(2)*POP_SPECIES(I,2) +
	1       TA(3)*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        T3=POP_SPECIES(I,1)+POP_SPECIES(I,2) +
	1          6.0_LDP*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        ZETA(I)=T1*T2*T3*EXP(-V_SHOCK_2/V(I))	 !Zeta is temporary
              END DO
	    END IF
	  ELSE
!
! Use X-ray emission as tubulated by a PLASMA code. Emission should  be
! tabulated per electron and per ion.
!
	    CALL GET_SCL_XRAY_FLUXES_V1(CONT_FREQ,
	1              XRAY_EMISS_1,XRAY_EMISS_2,
	1              NU_EVAL_CONT,NCF,FREQ_INDX,
	1              VSMOOTH_XRAYS,SECTION)
!
! We use T3 for the Electron density. We asume H, He, and C are fully ionized
! in the X-ray emitting plasma. All other species are assumed have Z=6.0
!
	    DO I=MDST,DEND
	      T1=EXP(-V_SHOCK_1/V(I))*(FILL_FAC_XRAYS_1)**2
	      T2=EXP(-V_SHOCK_2/V(I))*(FILL_FAC_XRAYS_2)**2
	      T3=POP_SPECIES(I,1)+2.0_LDP*POP_SPECIES(I,2)+
	1              6.0_LDP*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	      ZETA(I)=(T1*XRAY_EMISS_1+T2*XRAY_EMISS_2)*T3*POP_ATOM(I)
	    END DO
	  END IF
	  IF(XRAY_SMOOTH_WIND)ZETA(MDST:DEND)=ZETA(MDST:DEND)*CLUMP_FAC(MDST:DEND)		!Should be divided?
!
          ETA(MDST:DEND)=ETA(MDST:DEND)+ZETA(MDST:DEND)
	  ETA_MECH(MDST:DEND)=TA(MDST:DEND)
!
! Changed 06-Aug-2003: Clumping was not beeing allowed for when computing
! the shock luminosity.
!
	  IF(.NOT. LST_DEPTH_ONLY)THEN
	    CALL MPI_ALLREDUCE(MPI_IN_PLACE,ZETA,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	    CALL MPI_ERR_CHECK(IERR,'MPI_ALLREDUCE -- FST_ZETA -- COMP_OPAC')
	  END IF
	  IF(DST .EQ. 1)THEN
	    T1=0.241838_LDP		!eV to 10^15Hz
	    IF(SECTION .EQ. 'CONTINUUM')THEN
	      IF(FREQ_INDX .EQ. 1)THEN
	         dE_XRAY_TOT=0.0_LDP
	         dE_XRAY_0P1=0.0_LDP
	         dE_XRAY_1KEV=0.0_LDP
	      END IF
	      TA(DST:DEND)=ZETA(DST:DEND)*CLUMP_FAC(DST:DEND)*FQW(FREQ_INDX)
	      dE_XRAY_TOT(DST:DEND)=dE_XRAY_TOT(DST:DEND)+TA(DST:DEND)
	      IF(FL .GE. 100.0_LDP*T1)dE_XRAY_0P1(DST:DEND)=dE_XRAY_0P1(DST:DEND)+TA(DST:DEND)
	      IF(FL .GE. 1000.0_LDP*T1)dE_XRAY_1KEV(DST:DEND)=dE_XRAY_1KEV(DST:DEND)+TA(DST:DEND)
	    END IF
	  END IF
	ELSE
	  ETA_MECH=0.0_LDP
	END IF
	CALL TUNE(2,'XRAY_COMP_OPAC')
!
! Set a minimum emissivity. Mainly important when X-rays are not present.
!
	DO I=MDST,DEND
	  IF(ETA(I) .LT. 1.0E-280_LDP)THEN
	    ETA(I)=1.0E-280_LDP
	    ETA_NOSCAT(I)=1.0E-280_LDP
	  END IF
	END DO
!
! The continuum source function is defined by:
!                                              S= ZETA + THETA.J
	ZETA=0.0_LDP; THETA=0.0_LDP
	DO I=MDST,DEND
	  ZETA(I)=ETA(I)/CHI(I)
	  THETA(I)=CHI_SCAT(I)/CHI(I)
	END DO
	IF(.NOT. LST_DEPTH_ONLY)THEN
	  CALL MPI_ALLREDUCE(ZETA,TA,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	    CALL MPI_ERR_CHECK(IERR,'MPI_ALLREDUCE -- ZETA -- COMP_OPAC')
	  ZETA(1:ND)=TA(1:ND)
	  CALL MPI_ALLREDUCE(THETA,TA,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	    CALL MPI_ERR_CHECK(IERR,'MPI_ALLREDUCE -- THETA -- COMP_OPAC')
	  THETA(1:ND)=TA(1:ND)
	END IF
!
! Store TOTAL continuum line emissivity and opacity.
!
	IF(LST_DEPTH_ONLY)THEN
	  ETA_CONT(ND)=ETA(ND); CHI_CONT(ND)=CHI(ND)
	ELSE
	  CALL MPI_ALLREDUCE(ETA,ETA_CONT,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	    CALL MPI_ERR_CHECK(IERR,'MPI_ALLREDUCE -- ETA -- COMP_OPAC')
	  CALL MPI_ALLREDUCE(CHI,CHI_CONT,ND,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	    CALL MPI_ERR_CHECK(IERR,'MPI_ALLREDUCE -- CONT -- COMP_OPAC')
	END IF
!
	RETURN
	END
