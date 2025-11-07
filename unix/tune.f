C
C Routine to replace VMS TUNE routine for collecting tuning statistics for
C a program section.
C
C Usage:
C          CALL TUNE(1,'Unique ID')
C              .......
C              Section of code.
C              .......
C          CALL TUNE(2,'Same Unique ID')
C
C To print results:
C
C           CALL TUNE(3,' ')
C
C Output from print:
C
C           CPUTOT() : CPU time accumulated in each IDENT
C           WALLTOT() : Total wall time bewteen each IDENT.
C
C NB: Total time for each code section is accumulated on successive calls.
C i.e. TUNE(1,'Unique ID') does not initialize the counters.
C
	MODULE TUNE_DATA
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: MAX_IDS=100
!
        REAL(KIND=LDP), SAVE :: T0,T1,dT,OVERHEAD
        REAL(KIND=LDP), SAVE :: ST_CPU(MAX_IDS)
	REAL(KIND=LDP), SAVE :: CPUTOT(MAX_IDS)
        REAL(KIND=LDP), SAVE :: WALLTOT(MAX_IDS)
	INTEGER, SAVE :: STACK(MAX_IDS)
        CHARACTER(LEN=20), SAVE ::  IDLIST(MAX_IDS)
        CHARACTER(LEN=20), SAVE :: LAST_1_CALL=' '
        CHARACTER(LEN=20), SAVE :: LAST_2_CALL=' '
!
	REAL(KIND=LDP), SAVE :: RR0
	INTEGER*8, SAVE :: IEND_WALL
	INTEGER*8, SAVE :: IC0,IR0,IM0,IT1
	INTEGER*8, SAVE :: IST_WALL(MAX_IDS)
	INTEGER, SAVE :: NUM_IDS
        INTEGER, SAVE :: NUM_STACK
!
	INTEGER LUOUT
	INTEGER CURRENT_ID
	INTEGER ACTIVE_ID
        INTEGER I,LU,TERM_OUT
	INTEGER, SAVE :: CNT=0
	EXTERNAL TERM_OUT
!
	LOGICAL, PARAMETER :: OUT_DIAGNOSTICS=.FALSE. 
!
	LOGICAL, SAVE :: ONLY_PROC_ZERO=.TRUE.
	LOGICAL, SAVE :: DO_CALL_DIAGNOSTICS=.FALSE.
	LOGICAL, SAVE :: FIRST_TOO_MANY=.TRUE.
	LOGICAL, SAVE :: FIRST_UNMATCHED=.TRUE.
	LOGICAL, SAVE :: FIRSTTIME=.TRUE.
	CHARACTER(LEN=30), SAVE :: FILE_NAME
!
	END MODULE TUNE_DATA
!
	SUBROUTINE INIT_TUNE_ROUTINE(SET_ALL_PROCESSES,SET_CALL_DIAGNOSTICS)
	USE SET_KIND_MODULE
	USE TUNE_DATA
	LOGICAL SET_ALL_PROCESSES
	LOGICAL SET_CALL_DIAGNOSTICS
!
	IF(SET_ALL_PROCESSES)ONLY_PROC_ZERO=.FALSE.
	IF(SET_CALL_DIAGNOSTICS)DO_CALL_DIAGNOSTICS=.TRUE.
!
	RETURN
	END
!
        SUBROUTINE TUNE(LRUN,IDENT)
	USE SET_KIND_MODULE
	USE TUNE_DATA
	IMPLICIT NONE
!
! Altered 07-Nov-2025 : Default is to output tuning information only for processor zero.
!                         By calling INIT_TUNE_ROUTINE you can change the default to write 
!                         output for all processors. The same routine can also output call 
!                         information for process 0 for debugging help.
! Altered 01-Sep-2025 : Replace ETIME by call to CPU_TIME.
! Altered 18-Feb-2013 : MAX_IDS increased. STACK introduced.
!                       Routine should now be much more efficient, with less instructions per call.
! Altered 08-Mar-2010 : Change variable for system clock to 8 bytes.
!                          This prevents loss of elapsed time due to clock rollover.
! Altered 11-Nov-2000 : Call to F90 SYSTEM_CLOCK routine implemented.
!                       Wall time now returened as in original VMS routine.
!                       Counters now initialized if LRUN=0 is passed.
!
	INTEGER LRUN
	CHARACTER*(*) IDENT
!
	IF(MYPE .NE. 0 .AND. ONLY_PROC_ZERO)RETURN
!
	IF(DO_CALL_DIAGNOSTICS)THEN
	  IF(MYPE .EQ. 0)THEN
	    IF(LRUN .EQ. 1 .AND. IDENT .NE. LAST_2_CALL)THEN
	      IF(CNT .NE. 0)WRITE(6,'(I4,2X,A,I10)')LRUN,LAST_2_CALL,CNT
	      WRITE(6,'(I4,2X,A,I10)')LRUN,TRIM(IDENT)
	      FLUSH(6)
	      CNT=0
	    ELSE IF(LRUN .EQ. 2)THEN
	      IF(IDENT .NE. LAST_2_CALL)THEN
	         WRITE(6,'(I4,2X,A,I10)')LRUN,TRIM(IDENT); FLUSH(6)
	      ELSE
	         CNT=CNT+1
	      END IF
	    END IF
	    RETURN
	  END IF
	END IF
!
	IF(FIRSTTIME)THEN
          FIRSTTIME=.FALSE.
          DO  I=1,MAX_IDS
            ST_CPU(I)=0.0_LDP
            IST_WALL(I)=0.0_LDP
            CPUTOT(I)=0.0_LDP
            WALLTOT(I)=0.0_LDP
	    IDLIST(I)=' '
	    STACK(I)=1
          END DO
	  ACTIVE_ID=0
	  NUM_IDS=0
	  NUM_STACK=0
	  CALL SYSTEM_CLOCK(IC0,IR0,IM0);    RR0=IR0
          CALL CPU_TIME(T0)
	  CALL CPU_TIME(T1)
          OVERHEAD=2.0_LDP*(T1-T0)
	  WRITE(FILE_NAME,'(I3.3)')MYPE; FILE_NAME='TIMING_'//FILE_NAME
	  CALL GET_LU(LUOUT,'TIMING file in TUNE')
	  OPEN(UNIT=LUOUT,STATUS='REPLACE',FILE=FILE_NAME)
	  WRITE(LUOUT,*)' '
	  WRITE(LUOUT,*)'Overhead is ',OVERHEAD
	  WRITE(LUOUT,*)'   Count rate for wall clock is',IR0
	  WRITE(LUOUT,*)'Maximum count for wall clock is',IM0
	  WRITE(LUOUT,*)' '
	  CLOSE(LUOUT)
        ENDIF
!
! If LRUN =1, we are beginning the TIME bracket. Therefore we find the
! correct storage location first.
!
	IF (LRUN .EQ. 1) THEN
	  DO I=NUM_IDS,1,-1
            IF (IDENT .EQ. IDLIST(I))THEN
	      NUM_STACK=NUM_STACK+1
	      STACK(NUM_STACK)=I
	      CALL SYSTEM_CLOCK(IST_WALL(I))
	      CALL CPU_TIME(ST_CPU(I))
	      RETURN	
	    END IF
	  END DO
	  NUM_IDS=NUM_IDS+1
	  IF(NUM_IDS .LE. MAX_IDS)THEN
	    I=NUM_IDS
	    IDLIST(I)=IDENT
	    NUM_STACK=NUM_STACK+1
	    STACK(NUM_STACK)=I
	    CALL SYSTEM_CLOCK(IST_WALL(I))
	    CALL CPU_TIME(ST_CPU(I))
	    RETURN	
	  END IF
	  IF(FIRST_TOO_MANY)THEN
	    WRITE (LUOUT,'(A)')' ***** TOO MANY TUNING POINTS '
	    WRITE (LUOUT,'(A)')' Current TUNE points follow:'
	    DO I=1,NUM_IDS
	      WRITE(6,'(A)')TRIM(IDLIST(I))
	    END DO
	    FIRST_TOO_MANY=.FALSE.
	  END IF
	  RETURN
!
	ELSE IF (LRUN .EQ. 2) THEN
!
! If LRUN=2, we are ending the TIME bracket. Therefore we call the timing
! routine first.
!
! If TUNE has called been called correctly, then STACK should always be set correctly.
!
          CALL CPU_TIME(T0)
	  CALL SYSTEM_CLOCK(IEND_WALL)
	  ACTIVE_ID=STACK(NUM_STACK)
	  IF (IDENT .EQ. IDLIST(ACTIVE_ID))THEN
	    CURRENT_ID=ACTIVE_ID
	    NUM_STACK=NUM_STACK-1
	  ELSE
	    LU=TERM_OUT()
	    WRITE(LU,*)'Error in TUNE: STACK not correct'
	    WRITE(LU,*)'Make sure inner TUNE section is fully contained in outer TUNE section'
            WRITE(LU,*)LRUN,TRIM(IDENT)
	    WRITE(LU,*)' '
	    WRITE(LU,*)'A printout of the STACK follows'
	    WRITE(LU,*)' '
	    DO I=1,NUM_STACK
              WRITE(LU,'(I7,I7,T20,A)')I,STACK(I),TRIM(IDLIST(STACK(I)))
	    END DO
	    STOP
	    CURRENT_ID=0
            DO I=1,MAX_IDS
	      IF (IDENT.EQ.IDLIST(I))THEN
	        CURRENT_ID=I
	        EXIT
	      END IF
	    END DO
	  END IF
!
	  IF(CURRENT_ID .NE. 0)THEN
	    I=CURRENT_ID
	    dT=T0-ST_CPU(I)-OVERHEAD
	    CPUTOT(I)=CPUTOT(I)+MAX(0.0_LDP,dT)
	    IT1=IEND_WALL-IST_WALL(I)
	    IF(IT1 .LT. 0)IT1=IT1+IM0
	    WALLTOT(I)=WALLTOT(I)+IT1/RR0
	  ELSE IF(FIRST_UNMATCHED)THEN
	    FIRST_UNMATCHED=.FALSE.
	    LU=TERM_OUT()
	    WRITE(LU,*)' Warning ***** UNMATCHED TUNING POINT: ',TRIM(IDENT)
	  END IF
	  RETURN
C
	ELSE IF (LRUN .EQ. 3) THEN
	  CALL GET_LU(LUOUT,'TIMING file in TUNE(3)')
	  OPEN(UNIT=LUOUT,STATUS='OLD',ACTION='WRITE',POSITION='APPEND',FILE=FILE_NAME)
	  WRITE(LUOUT,'(8X,''Identifier'',11x,''Elapsed'',11x,''  CPU'')')
	  WRITE(LUOUT,'(29x,''  Time '',11x,''  Time'')')
	  DO I=1,MAX_IDS
	    IF (IDLIST(I).EQ.' ') EXIT
	    WRITE(LUOUT,'(1X,a24,f15.6,2X,f15.6)')
	1   IDLIST(I),WALLTOT(I),CPUTOT(I)
          END DO
	  CLOSE(LUOUT)
	ELSE IF(LRUN .EQ. 0) THEN
          DO  I=1,MAX_IDS
            ST_CPU(I)=0.0_LDP
            IST_WALL(I)=0.0_LDP
            CPUTOT(I)=0.0_LDP
            WALLTOT(I)=0.0_LDP
	    IDLIST(I)=' '
	    STACK(I)=1
          END DO
	  ACTIVE_ID=0
	  NUM_IDS=0
	  NUM_STACK=0
	ELSE
	  LU=TERM_OUT()
	  WRITE (LU,'(A)')' ***** ILLEGAL VALUE OF LRUN IN CALL TO TUNE '
	  WRITE(LU,*)' LRUN=',LRUN
	  STOP
	ENDIF

	RETURN
	END
