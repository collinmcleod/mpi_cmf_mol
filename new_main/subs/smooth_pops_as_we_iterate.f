	SUBROUTINE SMOOTH_POPS_AS_WE_ITERATE(POPS,STEQ_VALS,ND,NT)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	USE MPI
	IMPLICIT NONE
!
! Altered 11-Jun-2025 - Substantial changes made to improve smootinh in MPI.
!                          Changes made on galah (15-May-2025) for running models of 10 Lac.
! Altered 05-May-2025 - Substantial changes
!                     - POPS only directedly updated by root process.
! Altered 19-May-2019 _ Still under devlopment.
!
	INTEGER ND
	INTEGER NT
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) STEQ_VALS(NT,ND)
	INTEGER, ALLOCATABLE ::  DONE(:)
!
	INTEGER, SAVE :: LUOUT=0
	INTEGER I,J,K,L
	INTEGER ISPEC
	INTEGER ID,LOOP
	REAL(KIND=LDP) T0,T1,T2,T3
	REAL(KIND=LDP) RAT23,RAT12,RAT01
	REAL(KIND=LDP) SUM_STEQ
	LOGICAL DID_FULL_L
!!
	IF(LUOUT .EQ. 0 .AND. MYPE .EQ. 0)THEN
	  CALL GET_LU(LUOUT,'In SMOOTH_POP_AS_WE_ITERATE')
	  OPEN(UNIT=LUOUT,FILE='SMOOTH_POP_CORRECTIONS',STATUS='UNKNOWN',ACTION='WRITE')
	END IF
!
	IF(MYPE .EQ. 0)THEN
	  WRITE(LUOUT,*)'Maximum correction is [sign corrected]',-MINVAL(STEQ_VALS)
	  WRITE(LUOUT,*)'Minimum correction is [sign corrected]',-MAXVAL(STEQ_VALS)
	  WRITE(LUOUT,'(5X,6A,3(7X,A),5(3X,A))')
	1          '   L',' ISP','  ID','   I','    J',' DJ','RAT23','RAT12','RAT01',
	1          '  OPOP(L)','  New Val','OPOP(L+1)','OPOP(L+2)','OPOP(L+3)'
!
	 ALLOCATE(DONE(NT)); DONE=0
!
	  DID_FULL_L=.FALSE.
	  DO L=ND-5,1,-1
	    SUM_STEQ=SUM(STEQ_VALS(:,L))
	    IF(SUM_STEQ .EQ. 0.0_LDP)THEN
	      WRITE(6,*)'SUM_STEQ zero for depth ',L
	      IF(.NOT. DID_FULL_L)THEN
	        POPS(1:NT-2,L)=POPS(1:NT-2,L+1)
	        DID_FULL_L=.TRUE.
	        WRITE(6,*)'Replaced depth (SUM_STEQ) ',L
	      END IF
	    ELSE
	      DID_FULL_L=.FALSE.
	      DO ISPEC=1,NUM_SPECIES
	        DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	          K=0; IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)K=1
	          DO I=1,ATM(ID)%NXzV+K
	            J=ATM(ID)%EQXzV+I-1
	            IF(DONE(J) .GE. 2)THEN
	            ELSE IF(STEQ_VALS(J,L) .GE. 0.998_LDP .OR. STEQ_VALS(J,L) .LT. -1.0E+05_LDP)THEN
	              T3=LOG10(POPS(J,L+3))-LOG10(POP_SPECIES(L+3,ISPEC))
	              T2=LOG10(POPS(J,L+2))-LOG10(POP_SPECIES(L+2,ISPEC))
	              T1=LOG10(POPS(J,L+1))-LOG10(POP_SPECIES(L+1,ISPEC))
	              T0=LOG10(POPS(J,L))-LOG10(POP_SPECIES(L,ISPEC))
	              RAT23=T2-T3
	              RAT12=T1-T2
	              RAT01=T0-T1
	              IF( ABS(RAT23-RAT12) .LT. 1.0_LDP)THEN
	                IF( ABS(RAT01-RAT12) .GT. 5.0_LDP .AND. ABS(RAT23) .LT.  3.0_LDP)THEN
	                  T3=(10**RAT12)*(POPS(J,L+1)/POP_SPECIES(L+1,ISPEC))*POP_SPECIES(L,ISPEC)
	                  WRITE(LUOUT,'(A3,2X,4I4,I5,2X,L1,8ES12.3)')' CM',L,ISPEC,ID,I,J,DONE(J),
	1                   RAT23,RAT12,RAT01,POPS(J,L),T3,POPS(J,L+1:L+3)
	                  POPS(J,L)=T3
	                  DONE(J)=DONE(J)+1
	                ELSE
	                  T3=(10**RAT12)*(POPS(J,L+1)/POP_SPECIES(L+1,ISPEC))*POP_SPECIES(L,ISPEC)
	                  WRITE(LUOUT,'(A3,2X,4I4,I5,2X,I1,8ES12.3)')'NCM',L,ISPEC,ID,I,J,DONE(J),
	1                     RAT23,RAT12,RAT01,POPS(J,L),T3,POPS(J,L+1:L+3)
	                END IF
	              END IF
	            END IF
	          END DO
	        END DO
	      END DO
	    END IF
	  END DO
!
	  WRITE(6,*)'See SMOOTH_POP_CORRECTIONS for a listing of changes'
	  DEALLOCATE(DONE)
	END IF
	K=NT*ND
	CALL MPI_BCAST(POPS,K,MPI_DOUBLE_PRECISION,IZERO,MPI_COMM_WORLD,IERR)
!
	RETURN
	END
