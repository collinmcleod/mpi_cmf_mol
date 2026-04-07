!
! Routine designed to convert the SE%STEQ and SE%BA structures into the
! full STEQ and BA matrices that are used to find the corrections.
! The conversion is done one depth at a time. At each depth, the DIAGONAL
! matrix should be passed first.
!
	SUBROUTINE GENERATE_FULL_MATRIX_V3(C_MAT,STEQ_VEC,POPS,
	1                REPLACE,ZERO_STEQ,
	1                NT,ND,NION,NUM_BNDS,BAND_INDX,DIAG_INDX,DEPTH_INDX,
	1                FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	USE CONTROL_VARIABLE_MOD, ONLY : LTE_MODEL, USE_ELEC_HEAT_BAL, INCL_MOL_RXN, REPLACE_ALL
	USE MOL_RXN_MOD
	IMPLICIT NONE
!
! Altered 25-Aug-2022 : Collin McLeod changed equation replacement again. Molecular terms
!                       were not being analytically eliminated in a consistent way
!
! Altered 23-Feb-2022 : Collin McLeod changed equation replacement for molecules + C/O. Now includes
!                       molecular terms in ion equation where appropriate
!
! Lines 378-386 commented out by Collin McLeod to avoid replacing molecular equations with incorrect
! ionization equations
!
! Altered 11-Jan-2022 : Collin McLeod added prints for BA_ASCI before the matrix is scaled by populations
!
! Altered 10-Oct-2021 : Changed treatment of number conservation to treat molecules (edited by Collin McLeod)
!
! Altered 13-Mar-2014 : Issues with crude electron-energy balance equation when
!                          non-thermal ionization was included.
! Altered 30-Jan-2002 : Changed to V2
!                          REPLACE and ZERO_STEQ now passed in call.
! Altered 10-Sep-2001 : Using DIAG_INDX as logical variable in an IF statement
! Created 05-Apr-2001
!
!******************************************************************************
!
	INTEGER NT
	INTEGER ND
	INTEGER NION
	INTEGER NUM_BNDS
	INTEGER BAND_INDX
	INTEGER DIAG_INDX
	INTEGER DEPTH_INDX
!
        REAL(KIND=LDP) POPS(NT,ND)
        REAL(KIND=LDP) C_MAT(NT,NT)
	REAL(KIND=LDP) STEQ_VEC(NT)
	LOGICAL ZERO_STEQ(NT)		!Which vectors to be zeroed.
	LOGICAL REPLACE(NION)		!Which equations to be replaced
	LOGICAL FIRST_MATRIX		!Fist call to GENERATE_FULL_MATRIX.
	LOGICAL LAST_MATRIX		!Last call to   "       "    "
!
! When true we replace those equations indicated by the passed logical vector REPLACE.
! We always use the passed logical vector REPLACE if not updating the diagonal band.
!
	LOGICAL USE_PASSED_REP
!
! Local matrices and vectors to facilitate construction of C_MAT and STEQ_VEC.
!
! We use ?_ION to create the ionization balance equations for each
! ionization stage.
!
	REAL(KIND=LDP) C_ION(NION,NT)
	REAL(KIND=LDP) STEQ_ION(NION)
!
! We use ?_NC for the Number conservation equation for each species.
! We don't replace them directly into C_MAT and STEQ_VEC for ease of
! programming.
!
	REAL(KIND=LDP) C_NC(NUM_SPECIES,NT)
	REAL(KIND=LDP) STEQ_NC(NUM_SPECIES)
!
	REAL(KIND=LDP) G_SUM(NT)
	REAL(KIND=LDP) EDGE_SUM(NT)
	REAL(KIND=LDP) SUM,T1
	INTEGER NS
!
! FAC is used as a scale factor to determine at what depth the ground-state
! equilibrium equation is replaced by the ionization equation.
!
!	REAL(KIND=LDP), PARAMETER :: FAC=1.0D+05
!	REAL(KIND=LDP), PARAMETER :: FAC=1.0D+02
	REAL(KIND=LDP), SAVE ::  FAC=1.0E+02_LDP
!
	LOGICAL DIAG_BAND
!
! For consistency with the old version of CMFGEN, we only replace the
! ground-state equations over a continuous set of depths. REPLACE
! is used to indicate whether the current depth is to be replaced
! (as determined from the DIAGONAL band). Replace now passed.
!
	INTEGER, SAVE, ALLOCATABLE ::  REP_CNT(:)
!
	INTEGER, SAVE :: LST_DEPTH_INDX=0
!
	INTEGER I,J,K,L,JJ,II,LL,N,I1,I2
	INTEGER EQ
	INTEGER ID
	INTEGER ISPEC,ISPEC1
!
	INTEGER ERROR_LU,WARNING_LU
	INTEGER LUER,LUWARN,LUOUT
	EXTERNAL ERROR_LU,WARNING_LU
!
! Some variables to handle the number conservation for molecules
!
	INTEGER OXYSPEC,CARBSPEC
	INTEGER, ALLOCATABLE :: MOL_SPECS(:)
	LOGICAL IGNORE_MOL_LEV
	LOGICAL FILE_PRES
	CHARACTER(LEN=80) FILENAME
!
	INTEGER, PARAMETER :: NINDX=5
	INTEGER, SAVE :: WRITE_INDX(NINDX)	
!
! To save typing.
!
	LUER=ERROR_LU()
	LUWARN=WARNING_LU()
        K=DEPTH_INDX
	DIAG_BAND=.FALSE.
	IF(DIAG_INDX .EQ. BAND_INDX)DIAG_BAND=.TRUE.
	IF(DEPTH_INDX .NE. LST_DEPTH_INDX)THEN
	  LST_DEPTH_INDX=DEPTH_INDX
	  IF(.NOT. DIAG_BAND)THEN
	    WRITE(LUER,*)'Error in GENERATE_FULL_MATRIX'
	    WRITE(LUER,*)'Diagonal band must be passed first at each depth'
	    WRITE(LUER,*)'DEPTH_INDX=',DEPTH_INDX
	    STOP
	  END IF
	END IF
!
	INQUIRE(FILE='BA_REPLACEMENT_FACTOR',EXIST=FILE_PRES)
	IF(FILE_PRES)THEN
	  CALL GET_LU(I,'BA_REP')
	  OPEN(UNIT=I,FILE='BA_REPLACEMENT_FACTOR',STATUS='OLD',ACTION='READ')
	    READ(I,*)FAC
	  CLOSE(UNIT=I)
	END IF
!
	C_MAT(:,:)=0.0_LDP
	C_ION(:,:)=0.0_LDP
	C_NC(:,:)=0.0_LDP
!
	IF(DIAG_BAND)THEN
	  STEQ_VEC(:)=0.0_LDP
	  STEQ_ION(:)=0.0_LDP
	  STEQ_NC(:)=0.0_LDP
!
!	  IF(FIRST_MATRIX)THEN
!            I=SE(4)%LNK_TO_IV(1101)
!	     WRITE(99,*)I
!            WRITE(99,*)SE(1)%BA(1,I,DIAG_INDX,1),SE(1)%BA(2,I,DIAG_INDX,1)
!            WRITE(99,*)SE(3)%BA(1,I,DIAG_INDX,1),SE(3)%BA(2,I,DIAG_INDX,1)
!            WRITE(99,*)SE(4)%BA(1,I,DIAG_INDX,1),SE(4)%BA(2,I,DIAG_INDX,1)
!	  END IF
!
	END IF
!
	IF( .NOT. ALLOCATED(REP_CNT))THEN
	  ALLOCATE (REP_CNT(NION)); REP_CNT(:)=0
	END IF
!
! Map the small BA array onto the full BA array
!
! NB: C_ION(1,:) refers to to the ionization/recombination equation
!                for ion 1 (e.g. CI in the carbon sequence). It is
!                dN(CI)/dt. Since the equation (i.e. BA(I,:,:,:) with
!                I > ATMD(ID)%NXzV refers to dN/dt for the recombining
!                level (i.e. C2) we need a - sign when we evaluate
!                C_ION and STEQ_ION.
!
	DO ISPEC=1,NUM_SPECIES
	   DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO J=1,SE(ID)%N_IV
		 JJ=SE(ID)%LNK_TO_F(J)
		 DO I=1,SE(ID)%N_SE-1
		    EQ=SE(ID)%EQ_IN_BA(I)
		    C_MAT(EQ,JJ)=C_MAT(EQ,JJ)+SE(ID)%BA(I,J,BAND_INDX,K)
		 END DO
		 DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
		    C_ION(ID,JJ)=C_ION(ID,JJ)-SE(ID)%BA(I,J,BAND_INDX,K)
		 END DO
		 C_NC(ISPEC,JJ)=C_NC(ISPEC,JJ)+SE(ID)%BA(SE(ID)%N_SE,J,BAND_INDX,K)
	      END DO
	   END DO
	END DO
!
! Add additional terms to C_NC and STEQ_NC for molecular terms which fall outside the species in question
!


!__________________________
	IF (INCL_MOL_RXN) THEN
	   DO ISPEC=1,NUM_SPECIES
	      I=SPECIES_BEG_ID(ISPEC)
	      II=SPECIES_END_ID(ISPEC)-1
	      IF (IS_MOLECULE(ISPEC) .OR. ANY(SPECIESNAME_LIST .EQ. SPECIES(ISPEC) ) ) THEN
		 DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
		    IF (ID .NE. 0) THEN
		       DO J=1,N_MOL_LEVS
			  JJ = LNK_MOL_TO_FULL(J)
			  IF (JJ .LT. ATM(I)%EQXzV .OR. JJ .GT. ATM(II)%EQXzV + ATM(II)%NXzV) THEN
			     DO L=1,N_MOL_LEVS
				LL = LNK_MOL_TO_FULL(L)
				IF (LL .GE. ATM(ID)%EQxZV .AND. LL .LT. ATM(ID)%EQXzV + ATM(ID)%NXzV) THEN
				   C_NC(ISPEC,JJ)=C_NC(ISPEC,JJ)+MTOT_BA(L,J,K)
!				   C_ION(ID,JJ) = C_ION(ID,JJ) - MTOT_BA(L,J,K)
!				   WRITE(LUER,*) 'Adding molecular term in generate_full_matrix: ' ************************************************
!				   WRITE(LUER,*) 'ISPEC, SPECIES: ',ISPEC, ' ', SPECIES(ISPEC)
!				   WRITE(LUER,*) 'JJ, L, J, K: ', JJ, ' ', L, ' ', J, ' ', K
!				   WRITE(LUER,*) 'Added above to C_NC'
				END IF
			     END DO
			  END IF
		       END DO
		    END IF
		 END DO
	      END IF
	   END DO
	END IF
!________________________________
!
! Now do STEQ
!
	IF(DIAG_BAND)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,SE(ID)%N_SE-1
	        EQ=SE(ID)%EQ_IN_BA(I)
                STEQ_VEC(EQ)=STEQ_VEC(EQ)+SE(ID)%STEQ(I,K)
	      END DO
	      DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
	        STEQ_ION(ID)=STEQ_ION(ID)-SE(ID)%STEQ(I,K)
	      END DO
	      STEQ_NC(ISPEC)=STEQ_NC(ISPEC)+SE(ID)%STEQ(SE(ID)%N_SE,K)
	    END DO
	  END DO
	  STEQ_VEC(NT-1)=STEQ_ED(K)
	  STEQ_VEC(NT)=STEQ_T(K)
	END IF
!
! Update the ionization equations for when X-rays are included.
! The following is photoionizations/recombinations which change z by 2.
! Only other process allowed are /\z=1. The corrections to the ionization
! equations follow from simple algebraic manipulations.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-2
            I=SE(ID)%XRAY_EQ
	    IF(I .NE. 0)THEN
	      DO J=1,SE(ID)%N_IV
	        JJ=SE(ID)%LNK_TO_F(J)
	        C_ION(ID+1,JJ)=C_ION(ID+1,JJ)-SE(ID)%BA(I,J,BAND_INDX,K)
	      END DO
	      IF(DIAG_BAND)STEQ_ION(ID+1)=STEQ_ION(ID+1)-SE(ID)%STEQ(I,K)
	    END IF
	  END DO
	END DO
!
! Update the ionization equation to account for molecular reactions
! Note that, unlike the photoionization terms, molecular terms are with respect to 
! the ion being destroyed, so we do not need a - sign
!
!_____________________________
!
	IF (INCL_MOL_RXN) THEN
	   DO ISPEC=1,NUM_SPECIES
	      IF (IS_MOLECULE(ISPEC) .OR. ANY(SPECIESNAME_LIST .EQ. SPECIES(ISPEC) )) THEN
		 DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1 !,SPECIES_BEG_ID(ISPEC) !SPECIES_END_ID(ISPEC)-1
		    IF (ID .NE. 0) THEN
		       DO I=1,N_MOL_LEVS
			  II = LNK_MOL_TO_FULL(I)
			  DO J=1,N_MOL_LEVS+2
			     JJ = LNK_MOL_TO_FULL(J)
!       We only consider levels in this particular ion
			     ! Should it be the opposite? Only include levels not in this ion? added .NOT.
			     IF ( (II .GE. ATM(ID)%EQXzV) .AND. (II .LT. (ATM(ID)%EQXzV+ATM(ID)%NXzV)) ) THEN
!			     IF ( .NOT. (II .GE. ATM(SPECIES_BEG_ID(ISPEC))%EQXzV .AND. II .LE. EQ_SPECIES(ISPEC) ) ) THEN
				C_ION(ID,JJ) = C_ION(ID,JJ) + MTOT_BA(I,J,K)
!				C_ION(ID,JJ) = C_ION(ID,JJ) + MTOT_BA(I,J,K)
!			     WRITE(LUER,*) 'Updating ionization equation with molecular terms: '
!			     WRITE(LUER,*) 'ID, JJ, I, J, K: ',ID, JJ, I, J, K
!			     WRITE(LUER,*) '---------------------------------'
			     END IF
			  END DO
			  IF ( (II .GE. ATM(ID)%EQXzV) .AND. (II .LT. (ATM(ID)%EQXzV+ATM(ID)%NXzV)) .AND. DIAG_BAND) THEN
			     STEQ_ION(ID) = STEQ_ION(ID) + MTOT_STEQ(I,K)
!			     STEQ_ION(ID) = STEQ_ION(ID) + MTOT_STEQ(I,K)
!			     IF (SPECIES(SPECIES_LNK(ID)) .EQ. 'OXY' .AND. K .EQ. 51) THEN
!				WRITE(*,*) 'IN GENERATE_FULL_MATRIX, adding STEQ to oxygen: '
!				WRITE(*,*) 'ION: ', ID
!				WRITE(*,*) 'STEQ FROM MOL RXNS: ',MTOT_STEQ(I,K)
!			     END IF
!			  WRITE(LUER,*) 'Updating STEQ for ion with molecular terms: '
!			  WRITE(LUER,*) 'ID, I, K: ', ID, I, K
!			  WRITE(LUER,*) '*******************************************'
			  END IF
		       END DO
		    END IF
		 END DO
!END DO
	      ELSE
		 CYCLE
	      END IF
	   END DO
	END IF
!_______________________________
!
! Allow for advection terms in the ionization equations.
!
	IF(DIAG_BAND)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      STEQ_ION(ID)=STEQ_ION(ID)+SE(ID)%STEQ_ADV(K)
	    END DO
	  END DO
	END IF
	CALL SAVE_STEQ_ION(STEQ_ION,DEPTH_INDX,NION,DST,DEND)
!
! Allow for advection terms in the ionization equations. Since we are using linear derivatives,
! the terms, at each depth, are identical. We only need to worry about the sign of the terms.
!

	IF(BA_ADV_TERM(DIAG_INDX,K) .NE. 0.0_LDP)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      IF(SE(ID)%STRT_ADV_ID(K) .EQ. SPECIES_BEG_ID(ISPEC))THEN
	        DO L=SE(ID)%STRT_ADV_ID(K),ID
	          DO J=1,ATM(L)%NXzV
	            JJ=SE(L)%LNK_TO_F(J)
	            C_ION(ID,JJ)=C_ION(ID,JJ)-BA_ADV_TERM(BAND_INDX,K)
	          END DO
	        END DO
	      ELSE
	        DO L=ID+1,SE(ID)%END_ADV_ID(K)
	          DO J=1,ATM(L)%NXzV
	            JJ=SE(L)%LNK_TO_F(J)
	            C_ION(ID,JJ)=C_ION(ID,JJ)+BA_ADV_TERM(BAND_INDX,K)
	          END DO
	        END DO
!
!Including DxZV for last ionization state
!
	        L=SE(ID)%END_ADV_ID(K)
	        J=ATM(L)%NXzV
	        JJ=SE(L)%LNK_TO_F(J)+1
	        C_ION(ID,JJ)=C_ION(ID,JJ)+BA_ADV_TERM(BAND_INDX,K)
	      END IF
	    END DO
	  END DO
	END IF
!
! We now replace any equations that need replacing.
!
! Insert the number conservation equation for each species.
!
! Significant changes inserted here by Collin McLeod: molecules do not get their own number conservation equations
! carbon monoxide population is used to update both carbon and oxygen number conservation, etc.
!
	OXYSPEC=-1
	CARBSPEC=-1
	ALLOCATE(MOL_SPECS(N_SPECIES))
	MOL_SPECS = 0
	DO J=1,NUM_SPECIES
	   DO I=1,N_SPECIES
	      IF (TRIM(ADJUSTL(SPECIES(J))) .EQ. TRIM(ADJUSTL(SPECIESNAME_LIST(I)))) MOL_SPECS(I)=J
	      
	      IF (PRODUCT(MOL_SPECS) .GT. 1) EXIT
	   END DO
	END DO
!	WRITE(*,*) 'MOL_SPECS IN GENERATE_FULL_MATRIX_V3:'
!	WRITE(*,*) 'Number of species: ',N_SPECIES
!	DO J=1,N_SPECIES
!	   WRITE(*,*) 'SPECIESNAME(',J,')= ',SPECIESNAME_LIST(J)
!	   WRITE(*,*) 'MOL_SPECS(',J,')= ',MOL_SPECS(J)
!	   CALL FLUSH(6)
!	END DO
	IF (ANY(MOL_SPECS .LT. 0)) THEN
	   WRITE(LUER,*)'Error in GENERATE_FULL_MATRIX:'
	   WRITE(LUER,*)'Molecules included but requisite atoms not located'
	   WRITE(LUER,*)'MOL atom IDs: ',MOL_SPECS
	   STOP
	END IF
!
	DO ISPEC=1,NUM_SPECIES
	   IF(SPECIES_PRES(ISPEC))THEN
	      IF (IS_MOLECULE(ISPEC)) THEN
!_____________________________
		 DO II=1,LEN(SPECIESNAME_LIST)
		    IF (INDEX(MOL_ATOMS(ISPEC,1),TRIM(ADJUSTL(SPECIESNAME_LIST(II)))) .NE. 0) I1 = II
		    IF (INDEX(MOL_ATOMS(ISPEC,2),TRIM(ADJUSTL(SPECIESNAME_LIST(II)))) .NE. 0) I2 = II
		 END DO
!
		 I = EQ_SPECIES(MOL_SPECS(I1))
		 IF (MOL_ATOMS(ISPEC,2) .NE. 'NONE') J = EQ_SPECIES(MOL_SPECS(I2))

		 C_MAT(I,:) = C_MAT(I,:) + MOL_ATM_NUMS(ISPEC,1)*C_NC(ISPEC,:)
		 IF (MOL_ATOMS(ISPEC,2) .NE. 'NONE') C_MAT(J,:) = C_MAT(J,:) + MOL_ATM_NUMS(ISPEC,2)*C_NC(ISPEC,:)
!_____________________________
		
!		IF (SPECIES(ISPEC) .EQ. 'COMOL') THEN
!		   C_MAT(I,:) = C_MAT(I,:) + C_NC(ISPEC,:)
!		   C_MAT(J,:) = C_MAT(J,:) + C_NC(ISPEC,:)
!		ELSE IF (SPECIES(ISPEC) .EQ. 'C2MOL') THEN
!		   C_MAT(J,:) = C_MAT(J,:) + 2*C_NC(ISPEC,:)
!		ELSE IF (SPECIES(ISPEC) .EQ. 'O2MOL') THEN
!		   C_MAT(I,:) = C_MAT(I,:) + 2*C_NC(ISPEC,:)
!		ELSE IF (SPECIES(ISPEC) .EQ. 'CO2MOL') THEN
!		   C_MAT(I,:) = C_MAT(I,:) + 2*C_NC(ISPEC,:)
!		   C_MAT(J,:) = C_MAT(J,:) + C_NC(ISPEC,:)
!		ELSE IF (SPECIES(ISPEC) .EQ. 'C2OMOL') THEN
!		   C_MAT(I,:) = C_MAT(I,:) + C_NC(ISPEC,:)
!		   C_MAT(J,:) = C_MAT(J,:) + 2*C_NC(ISPEC,:)
!		END IF
	      ELSE
		 I=EQ_SPECIES(ISPEC)
		 C_MAT(I,:)=C_MAT(I,:) + C_NC(ISPEC,:)
	      END IF
	   END IF
	END DO
!       
	DO ISPEC=1,NUM_SPECIES
	   IF(DIAG_BAND .AND. SPECIES_PRES(ISPEC))THEN
	      IF (IS_MOLECULE(ISPEC)) THEN
!___________________________________
		 DO II=1,LEN(SPECIESNAME_LIST)
		    IF (INDEX(MOL_ATOMS(ISPEC,1),TRIM(ADJUSTL(SPECIESNAME_LIST(II)))) .NE. 0) I1 = II
		    IF (INDEX(MOL_ATOMS(ISPEC,2),TRIM(ADJUSTL(SPECIESNAME_LIST(II)))) .NE. 0) I2 = II
		 END DO
		 I = EQ_SPECIES(MOL_SPECS(I1))
		 IF (MOL_ATOMS(ISPEC,2) .NE. 'NONE') J = EQ_SPECIES(MOL_SPECS(I2))

		 STEQ_VEC(I) = STEQ_VEC(I) + MOL_ATM_NUMS(ISPEC,1)*STEQ_NC(ISPEC)
		 IF (MOL_ATOMS(ISPEC,2) .NE. 'NONE') STEQ_VEC(J) = STEQ_VEC(J) + MOL_ATM_NUMS(ISPEC,2)*STEQ_NC(ISPEC)
!___________________________________
!		I = EQ_SPECIES(OXYSPEC)
!		J = EQ_SPECIES(CARBSPEC)
!		IF (SPECIES(ISPEC) .EQ. 'COMOL') THEN
!		   STEQ_VEC(I) = STEQ_VEC(I) + STEQ_NC(ISPEC)
!		   STEQ_VEC(J) = STEQ_VEC(J) + STEQ_NC(ISPEC)
!		ELSE IF (SPECIES(ISPEC) .EQ. 'C2MOL') THEN
!		   STEQ_VEC(J) = STEQ_VEC(J) + 2*STEQ_NC(ISPEC)
!		ELSE IF (SPECIES(ISPEC) .EQ. 'O2MOL') THEN
!		   STEQ_VEC(I) = STEQ_VEC(I) + 2*STEQ_NC(ISPEC)
!		ELSE IF (SPECIES(ISPEC) .EQ. 'CO2MOL') THEN
!		   STEQ_VEC(I) = STEQ_VEC(I) + 2*STEQ_NC(ISPEC)
!		   STEQ_VEC(J) = STEQ_VEC(J) + STEQ_NC(ISPEC)
!		ELSE IF (SPECIES(ISPEC) .EQ. 'C2OMOL') THEN
!		   STEQ_VEC(I) = STEQ_VEC(I) + STEQ_NC(ISPEC)
!		   STEQ_VEC(J) = STEQ_VEC(J) + 2*STEQ_NC(ISPEC)
!		END IF
	      ELSE
		 I=EQ_SPECIES(ISPEC)
		 STEQ_VEC(I)=STEQ_VEC(I) + STEQ_NC(ISPEC)
	      END IF
	   END IF
	END DO
!
! Charge conservation and radiative equilibrium equations.
!
	C_MAT(NT-1,:)=BA_ED(:,BAND_INDX,DEPTH_INDX)
	IF(USE_ELEC_HEAT_BAL)THEN
	  IF(DEPTH_INDX .EQ. 1)WRITE(6,*)'Using EHB for T equation'
	  C_MAT(NT,:)  =BA_T_EHB(:,BAND_INDX,DEPTH_INDX)
	  IF(DIAG_BAND)STEQ_VEC(NT)=STEQ_T_EHB(DEPTH_INDX)
	ELSE
	  C_MAT(NT,:)  =BA_T(:,BAND_INDX,DEPTH_INDX)
	END IF
!
! Crude section to create ELECTRON cooling equation.
! Must be done after all equations are done, but before GS equation is replaced.
!
! The following was omitted 13-Mar-2014. Problem with non-thermal ionization. Needs modifcation.
!
	IF(DEPTH_INDX .LE. -12)THEN
          DO ISPEC=1,NUM_SPECIES
	    IF(SPECIES_PRES(ISPEC))THEN
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        NS=ATM(ID)%NXzV
                G_SUM(1:NS)=0.0_LDP
                EDGE_SUM(1:NS)=0.0_LDP
                DO J=1,ATM(ID)%NXzV_F
                  I=ATM(ID)%F_TO_S_XzV(J)
                  G_SUM(I)=G_SUM(I)+ATM(ID)%GXzV_F(J)
                  EDGE_SUM(I)=EDGE_SUM(I)+ATM(ID)%EDGEXzV_F(J)*ATM(ID)%GXzV_F(J)
                END DO
!
	        IF(DIAG_BAND)THEN
	          SUM=0.0_LDP
                  DO I=1,ATM(ID)%NXzV
                    T1=EDGE_SUM(I)/G_SUM(I)
	            EQ=SE(ID)%EQ_IN_BA(I)
	            SUM=SUM+5.27296E-03_LDP*T1*STEQ_VEC(EQ)
	          END DO
	          STEQ_VEC(NT)=STEQ_VEC(NT)+SUM
	        END IF
!
	        DO J=1,NT
	          SUM=0.0_LDP
                  DO I=1,ATM(ID)%NXzV
                    T1=EDGE_SUM(I)/G_SUM(I)
	            EQ=SE(ID)%EQ_IN_BA(I)
	            SUM=SUM+5.27296E-03_LDP*T1*C_MAT(EQ,J)
	          END DO
	          C_MAT(NT,J)=C_MAT(NT,J)+SUM
	        END DO
	      END DO
	    END IF
	  END DO
	END IF
!
! Check if we wish to use the ionization conservation equation.
! We only use the diagonal band to do this. We use REPLACE
! array so that we only use the conservation equation over a
! sequence of optical depths.
!
! Altered Jan 12 by Collin McLeod to avoid replacing any molecular g.s. equations
! Note that, since molecules are treated as one level, the g.s. and ion equations 
! should be identical (and thus replacement is unnecessary)
!
	IF(DIAG_BAND .AND. DEPTH_INDX .EQ. 1)REP_CNT(:)=0
        IF(.NOT. USE_PASSED_REP .AND. DIAG_BAND)THEN
	   DO ID=1,NION
	      IF (REPLACE_ALL) THEN
		 REPLACE(ID)=.TRUE.
	      ELSE
		 REPLACE(ID)=.FALSE.
	      END IF
	    IF (INCL_MOL_RXN) THEN
	       IF( ATM(ID)%XzV_PRES .AND. .NOT. IS_MOLECULE(SPECIES_LNK(ID)) ) THEN !.NOT. ANY(SPECIESNAME_LIST .EQ. SPECIES(SPECIES_LNK(ID)) ) ) THEN
!
!       Don't replace equations which have molecular reaction terms--can cause an inconsistency
!       IF (.NOT. (SPECIES_LNK(ID) .EQ. CARBSPEC .OR. SPECIES_LNK(ID) .EQ. OXYSPEC) ) THEN
		  EQ=ATM(ID)%EQXzV
		  N=ATM(ID)%NXzV
		  IF( ABS(C_MAT(EQ,EQ))*ATM(ID)%XzV(1,K) .GT.
	1	       FAC*ABS(C_MAT(EQ,EQ+N))*ATM(ID)%DXzV(K) )THEN
!
!       Commented out 23-Aug-2022 to test running without equation replacement
!
		     REPLACE(ID)=.TRUE.
		     REP_CNT(ID)=REP_CNT(ID)+1
		  END IF
	       END IF
	    ELSE
	       IF( ATM(ID)%XzV_PRES )  THEN
!       
		  EQ=ATM(ID)%EQXzV
		  N=ATM(ID)%NXzV
		  IF( ABS(C_MAT(EQ,EQ))*ATM(ID)%XzV(1,K) .GT.
	1	       FAC*ABS(C_MAT(EQ,EQ+N))*ATM(ID)%DXzV(K) )THEN
!
!       Commented out 23-Aug-2022 to test running without equation replacement
!       
		     REPLACE(ID)=.TRUE.
		     REP_CNT(ID)=REP_CNT(ID)+1
		  END IF
	       END IF
	    END IF
	 END DO
	END IF
	
	IF(DIAG_BAND)CALL SAVE_REPLACE(REPLACE,ION_ID,DEPTH_INDX,NION,DST,DEND)	
!
! In all cases, we replace the ground state equation.
! We need the check on XzV_PRES since we initially set
! all the REPLACE to true.
!
	DO ID=1,NION
          IF(REPLACE(ID) .AND. ATM(ID)%XzV_PRES)THEN
	    C_MAT(ATM(ID)%EQXzV,:)=C_ION(ID,:)
	    IF(DIAG_BAND)STEQ_VEC(ATM(ID)%EQXzV)=STEQ_ION(ID)
	  END IF
	END DO
!
	IF(DIAG_BAND .AND. DEPTH_INDX .EQ. ND)THEN
	  WRITE(LUWARN,'(/,/,1X,A,/)')' Equation selection in generate_full_matrix_v3.f'
	  DO ID=1,NION
	    IF(REP_CNT(ID) .GT. 0)THEN
              WRITE(LUWARN,'(1X,A,T9,A,I3,A)')
	1         TRIM(ION_ID(ID)),' g.s. eq. replaced by ionization eq. at ',
	1         REP_CNT(ID),' depths.'
	      END IF
	  END DO
	  WRITE(LUWARN,'(A)')' '
	END IF
!
!	IF(K .EQ. 1 .AND. DIAG_BAND)THEN
!	  CALL WR2D(STEQ_VEC,NT,1,'STEQ_VEC_D1',94)
!	  OPEN(UNIT=96,FILE='BA_ASCI_D1',STATUS='UNKNOWN')
!	    CALL WR2D(C_MAT,NT,NT,'C_MAT_VEC_D1',96)
!	  CLOSE(UNIT=96)
!	END IF
!	IF(K .EQ. 61 .AND. DIAG_BAND)THEN
!	  CALL WR2D(STEQ_VEC,NT,1,'STEQ_VEC_D10',94)
!	  OPEN(UNIT=96,FILE='BA_ASCI_D61',STATUS='UNKNOWN',POSITION='APPEND')
!	    CALL WR2D(C_MAT,NT,NT,'C_MAT_VEC_D10',96)
!	  CLOSE(UNIT=96)
!	END IF
!
! Fix any populations.
!
	IF(USE_PASSED_REP)THEN
	  DO I=1,NT
	     IF(ZERO_STEQ(I))STEQ_VEC(I)=0.0_LDP
	  END DO
	ELSE
	  CALL FIXPOP_IN_BA_V3(C_MAT,STEQ_VEC,ZERO_STEQ,
	1       NT,ND,NION,DIAG_BAND,K,
	1       FIRST_MATRIX,LAST_MATRIX)
	END IF
!
	IF(LTE_MODEL)THEN
	  CALL ADJUST_CMAT_TO_LTE(C_MAT,STEQ_VEC,DIAG_BAND,DEPTH_INDX,NT)
	END IF
!
! Print the ba matrices before scaling (for comparison with molecular reactions
!
	CALL GET_LU(LUOUT,'In generate_full_matrix_v3')
!
	IF(K .EQ. 31 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D31',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 30 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D30',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 32 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D32',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 33 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D33',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 34 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D34',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 29 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D29',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 28 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D28',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 41 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D41',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 51 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D51',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 5 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D5',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 2 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D2',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 60 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D60',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 57 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D57',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 63 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D63',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 55 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D55',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
!
! Extra writes for checks related to Luc's models
!
	IF(K .EQ. 91 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	   OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D91',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 92 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D92',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF	
	IF(K .EQ. 135 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D135',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 136 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D136',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 137 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D137',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 138 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D138',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 140 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D140',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
!
	IF(K .EQ. 183 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D183',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 184 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D184',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 185 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D185',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 186 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D186',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 187 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D187',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 188 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D188',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 189 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D189',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 190 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D190',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 184 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D184',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
 	IF(K .EQ. 177 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D177',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 88 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D88',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
 	IF(K .EQ. 196 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D196',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
 	IF(K .EQ. 309 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D309',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
 	IF(K .EQ. 111 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D111',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
 	IF(K .EQ. 201 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_UNSCALED_D201',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF

!
! Scale the BA matrix so that we solve for the fractional corrections to
! the populations. This seems to yield better solutions for large matrices,
! especially when we subsequently condition the matrices.
!
	L=DEPTH_INDX-DIAG_INDX+BAND_INDX
	DO J=1,NT
	  DO I=1,NT
	    C_MAT(I,J)=C_MAT(I,J)*POPS(J,L)
	  END DO
	END DO
!
	IF(K .EQ. 31 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_D31',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
!
	IF(K .EQ. 41 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_D41',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 51 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_D51',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 5 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_D5',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 60 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_D60',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 57 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_D57',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 63 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_D63',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
	IF(K .EQ. 55 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=LUOUT,FILE='BA_ASCI_N_D55',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',LUOUT)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',LUOUT)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',LUOUT)
	  CLOSE(UNIT=LUOUT)
	END IF
!
!
!	WRITE(LUER,*)'NT= ',NT
!
	RETURN
	END
