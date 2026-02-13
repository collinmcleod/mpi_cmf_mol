	SUBROUTINE READ_ARNAUD_ION_DATA(ND)
	USE SET_KIND_MODULE
	USE STEQ_DATA_MOD               !, :: ONLY SE%XRAY_EQ
	USE MOD_CMFGEN
	USE MOD_NON_THERM
	IMPLICIT NONE
	INTEGER ND
!
! Altered : 12-Jan-2025 : Unimportant access error to array fixed.
!                         Bug fix -- issues for double ionization when X-rays not included. 
! Altered : 26-Apr-2019 : better error handling and reporting.
!
	INTEGER I,J,K
	INTEGER IT1,IT2
	INTEGER IT,ISPEC
	INTEGER ID,ID_ION
	INTEGER RD_AT_NO
	INTEGER IOS
	INTEGER NM_POS
	CHARACTER(LEN=120) TMP_NAME
	CHARACTER(LEN=120) STRING
	CHARACTER(LEN=40) TMP_ION_NAME
!
	INTEGER, PARAMETER :: LU_ER=6
	INTEGER, PARAMETER :: LU_IN=7
	INTEGER, PARAMETER :: LU_OUT=10
	INTEGER LU_WARN
!
	REAL(KIND=LDP), PARAMETER :: Hz_TO_eV=4.1356691_LDP
	REAL(KIND=LDP) T1,T2
!
	LOGICAL MATCH_FOUND
	INTEGER COUNT_OCCUR,WARNING_LU
	EXTERNAL COUNT_OCCUR,WARNING_LU
!
	LU_WARN=WARNING_LU()
	OPEN(UNIT=LU_OUT,FILE='NON_THERM_ION_SUM',STATUS='UNKNOWN')
!	CALL SET_LINE_BUFFERING(LU_OUT)
!
! File containing the parameters for the ionization potential
! and direct-ionization cross section used in Arnaud & Rothenflug 1985
! These parameters are then used in
!
	IF(MYPE .EQ. 0)WRITE(LU_ER,'(/,1X,A)')'Opening arnaud_rothenflug.dat'
	OPEN(UNIT=LU_IN,FILE='arnaud_rothenflug.dat',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*)'Unble to open arnaud_rothenflug.dat in READ_ION'
	    WRITE(LU_ER,*)'IOS=',IOS
	    STOP
	  END IF
	  IF(MYPE .EQ. 0)WRITE(6,*)'Opened arnaud_rothenflug.dat'; FLUSH(UNIT=6)
!
! All names must begin in the same column.
!
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(LU_IN,'(A)',END=100)STRING
	  END DO
	  NM_POS=INDEX(STRING,'  Ion')+2
	  BACKSPACE(LU_IN)
!
	  IT=0
	  THD(1:MAX_NUM_THD)%PRES=.FALSE.
	  DO WHILE(1 .EQ. 1)
!
! Skip blank lines and comments.
!
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	      READ(LU_IN,'(A)',END=100)STRING
	    END DO
!
! Find species and ionization state match.
!
	    READ(STRING,*,IOSTAT=IOS)RD_AT_NO,IT2
	    IF(IOS .NE. 0)THEN
	      WRITE(LU_ER,*)'Error reading in Arnaud atomic No.in READ_ARNAUD_ION_DATA'
	      WRITE(LU_ER,*)TRIM(STRING)
	      STOP
	    END IF
	    MATCH_FOUND=.FALSE.
	    DO ISPEC=1,NUM_SPECIES
	      IF( NINT( AT_NO(ISPEC) ) .EQ. RD_AT_NO .AND. SPECIES_BEG_ID(ISPEC) .NE. 0)THEN
	        IT1=RD_AT_NO+1-IT2
	        WRITE(LU_OUT,*)RD_AT_NO,AT_NO(ISPEC),IT1,SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	          IF( IT1 .EQ. NINT(ATM(ID)%ZXzV) )THEN
	             MATCH_FOUND=.TRUE.
	             IT=IT+1
	             IF(IT .GT. MAX_NUM_THD)THEN
	               WRITE(LU_ER,*)'Error reading in Arnaud ionization data in READ_ARNAUD_ION_DATA'
	               WRITE(LU_ER,*)'IT exceeds MAX_NUM_THD: MAX_NUM_THD=',MAX_NUM_THD
	               WRITE(LU_ER,*)TRIM(STRING)
	               STOP
	             END IF
	             THD(IT)%LNK_TO_SPECIES=ISPEC
	             THD(IT)%LNK_TO_ION=ID
	             THD(IT)%NTAB=0
	             IF(COUNT_OCCUR(STRING,'.') .GT. 1)THEN
	               READ(STRING,*) THD(IT)%ZION,THD(IT)%N_ION_EL,THD(IT)%PQN,THD(IT)%ANG,THD(IT)%ION_POT, &
	                   THD(IT)%A_COL,THD(IT)%B_COL,THD(IT)%C_COL,THD(IT)%D_COL
	             ELSE
	               READ(STRING,*) THD(IT)%ZION,THD(IT)%N_ION_EL,THD(IT)%PQN,THD(IT)%ANG,THD(IT)%ION_POT, THD(IT)%NTAB
	               ALLOCATE(THD(IT)%XTAB(THD(IT)%NTAB),STAT=IOS)
	               IF(IOS .EQ. 0)ALLOCATE(THD(IT)%YTAB(THD(IT)%NTAB),STAT=IOS)
	               IF(IOS .NE. 0)THEN
	                 WRITE(6,*)'Unable to allocate XTAB and YTAB in READ_ARNAUD'
	                 STOP
	               END IF
	               READ(LU_IN,*)THD(IT)%XTAB
	               READ(LU_IN,*)THD(IT)%YTAB
	               THD(IT)%XTAB=THD(IT)%XTAB*THD(IT)%ION_POT
	             END IF
	             THD(IT)%PRES=.TRUE.
!
! Find level in atom. We assume ground term, unless the lower term (indicated by a *) is specified.
!
	             THD(IT)%ATOM_STATES=0
	             THD(IT)%N_STATES=0
	             TMP_NAME=ADJUSTL(STRING(NM_POS:))
	             IF(TMP_NAME(1:1) .EQ. '*')THEN
	               TMP_NAME=TMP_NAME(2:)
	             ELSE
	               J=INDEX(ATM(ID)%XzVLEVNAME_F(1),'[')-1
	               IF(J .GT. 0)THEN
	                 TMP_NAME=ATM(ID)%XZVLEVNAME_F(1)(1:J)
	               ELSE
	                 TMP_NAME=TRIM(ATM(ID)%XZVLEVNAME_F(1))
	               END IF
	             END IF
	             DO I=1,MIN(30,ATM(ID)%NXzV_F)
	               K=INDEX(ATM(ID)%XzVLEVNAME_F(I),'[')-1
	               IF(K .LT. 0)K=LEN_TRIM(ATM(ID)%XzVLEVNAME_F(I))
	                IF(INDEX(TMP_NAME,ATM(ID)%XZVLEVNAME_F(I)(1:K)) .NE. 0)THEN
	                  THD(IT)%N_STATES=THD(IT)%N_STATES+1
	                  THD(IT)%ATOM_STATES(THD(IT)%N_STATES)=I
	               END IF
	             END DO
	             IF(THD(IT)%N_STATES .EQ. 0)THEN
	               WRITE(LU_ER,*)'Error in READ_ARNAUD_ION_DATA'
	               WRITE(LU_ER,*)'Lower level not found'
	               WRITE(LU_ER,*)TRIM(STRING)
	               STOP
	             END IF
!
! Find level(s) in ion.
!
	             TMP_NAME=ADJUSTL(STRING(NM_POS:))
	             IF(MYPE .EQ. 0)THEN
	               IF(ID .EQ. 1)THEN
	                 WRITE(LU_ER,*)'Check WARNINGS to see what Arnaud cross-sectons have been read in'
	                 WRITE(LU_WARN,'(//,A3,A7,A28,A23,3X,A,5X,A20,/)')  &
	                                   'ID','Ion','Arn. final state',' GS level', &
                                           'Xray Eqn.','Final State'
	               END IF
	             END IF
	             TMP_ION_NAME=' '
	             IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)THEN
	               THD(IT)%SUM_GION=ATM(ID)%GIONXzV_F
	               THD(IT)%N_ION_ROUTES=1
	               THD(IT)%ION_LEV(1)=SE(ID)%XRAY_EQ
	               TMP_ION_NAME='Final ion'
	             ELSE
	               THD(IT)%SUM_GION=0.0_LDP
	               THD(IT)%N_ION_ROUTES=0
	               ID_ION=ID+1
	               DO I=1,ATM(ID_ION)%NXzV_F
	                 J=INDEX(ATM(ID_ION)%XZVLEVNAME_F(I),'[')-1
	                 IF(J .LT. 0)J=LEN_TRIM(ATM(ID_ION)%XZVLEVNAME_F(I))
	                 IF(INDEX(TRIM(TMP_NAME),ATM(ID_ION)%XZVLEVNAME_F(I)(1:J)) .NE. 0)THEN
	                   THD(IT)%N_ION_ROUTES = THD(IT)%N_ION_ROUTES + 1
	                   THD(IT)%ION_LEV(THD(IT)%N_ION_ROUTES)=I
	                   THD(IT)%SUM_GION=THD(IT)%SUM_GION+ATM(ID_ION)%GXzV_F(I)
	                   TMP_ION_NAME=ATM(ID_ION)%XZVLEVNAME_F(I)
	                 END IF
	               END DO
	               IF(INDEX(TMP_NAME,'&') .NE. 0 .AND. THD(IT)%N_ION_ROUTES .EQ.0)THEN
	                 THD(IT)%SUM_GION=ATM(ID+1)%GIONXzV_F
	                 TMP_ION_NAME=ATM(ID+1)%XZVLEVNAME_F(1)
	                 THD(IT)%N_ION_ROUTES=1
	                 IF(SE(ID)%XRAY_EQ .EQ. 0)THEN
	                    THD(IT)%ION_LEV(THD(IT)%N_ION_ROUTES)=ATM(ID)%NXzV+1
	                 ELSE
	                   THD(IT)%ION_LEV(THD(IT)%N_ION_ROUTES)=SE(ID)%XRAY_EQ
	                 END IF
	               END IF
	             END IF
	             IF(MYPE .EQ. 0)THEN
	               DO J=1,THD(IT)%N_ION_ROUTES
	                 WRITE(LU_WARN,'(I3,A7,3X,A30,3X,A20,7X,I5,5X,A20)')           &
	                           ID,TRIM(ION_ID(ID)),TRIM(TMP_NAME),               &
	                           TRIM(ATM(ID)%XZVLEVNAME_F(1)),                    &
	                           THD(IT)%ION_LEV(J),TRIM(TMP_ION_NAME)
	                END DO
	             END IF
!
	             IF(THD(IT)%N_ION_ROUTES .EQ. 0)THEN
	               IF(MYPE .EQ. 0)THEN
	                 WRITE(6,*)'Error in RD_ARNAUD_ION_DATA -- unmatched ion level name'
	                 WRITE(6,*)TRIM(STRING)
	                 WRITE(6,*)'Extracted ion level name is: ',TRIM(TMP_NAME)
	               END IF
	               T1=ATM(ID_ION)%EDGEXzV_F(1)-ATM(ID_ION)%EDGEXzV_F(ATM(ID_ION)%NXzV_F)
	               T1=ATM(ID)%EDGEXzV_F(1)+T1
	               T2=ATM(ID)%EDGEXzV_F(1)+ATM(ID_ION)%EDGEXzV_F(1)
	               IF(MYPE .EQ. 0)WRITE(6,*)THD(IT)%ION_POT,Hz_TO_eV*T1,Hz_TO_eV*T2
	               IF(THD(IT)%ION_POT .GT. Hz_TO_eV*T1 .AND. THD(IT)%ION_POT .LT. Hz_TO_eV*T2)THEN
	                 IF(MYPE .EQ. 0)WRITE(6,*)'Continuing execution as level not included'
	                 THD(IT)%NTAB=0
	                 IT=IT-1
	               ELSE
	                 STOP
	               END IF
	             END IF
!
	             EXIT
	          END IF
	        END DO
	        EXIT
              END IF
	    END DO
	    IF(.NOT. MATCH_FOUND .AND. COUNT_OCCUR(STRING,'.') .LT. 2)THEN
	      READ(LU_IN,'(A)')STRING
	      READ(LU_IN,'(A)')STRING
	    END IF
	  END DO
100	CONTINUE
	NUM_THD=IT
	CLOSE(UNIT=LU_IN)
!
	IF(MYPE .EQ. 0)THEN
	  CALL GET_LU(J,'Summary of cross sections')
	  OPEN(UNIT=J,STATUS='UNKNOWN',FILE='ARNAUD_SUMMARY')
	  DO IT=1,NUM_THD
	    WRITE(J,*) THD(IT)%PRES,THD(IT)%ZION,THD(IT)%N_ION_EL,THD(IT)%PQN,THD(IT)%ANG,THD(IT)%ION_POT, &
	         THD(IT)%A_COL,THD(IT)%B_COL,THD(IT)%C_COL,THD(IT)%D_COL
	  END DO	
	  CLOSE(J)
	END IF
!
	DO IT=1,NUM_THD
	  IF(THD(IT)%PRES)THEN
!
	    I=THD(IT)%N_STATES
	    WRITE(STRING,'(I3,A,10(I2,A))')I,'(',(THD(IT)%ATOM_STATES(J),',',J=1,I)
	    STRING(LEN_TRIM(STRING):)=')'
!
	    I=THD(IT)%N_ION_ROUTES
            WRITE(TMP_NAME,'(I3,A,10(I3,A))')I,'(',(THD(IT)%ION_LEV(J),',',J=1,I)
	    TMP_NAME(LEN_TRIM(TMP_NAME):)=')'
	    I=1
	    DO WHILE(TMP_NAME(I+1:I+1) .NE. ')')
	      IF(TMP_NAME(I:I+1) .EQ. '   ')TMP_NAME(I:)=TMP_NAME(I+1:)
	      I=I+1
	    END DO
!
	    WRITE(LU_OUT,'(I4,2X,I2,A5,I4,2X,A6,3X,A,T60,A)') &
	      IT,THD(IT)%LNK_TO_SPECIES,SPECIES_ABR(THD(IT)%LNK_TO_SPECIES), &
	      THD(IT)%LNK_TO_ION,ION_ID(THD(IT)%LNK_TO_ION), &
	      TRIM(STRING),TRIM(TMP_NAME)
	  END IF
	END DO
!
	CLOSE(LU_OUT)
	RETURN
	END
