!
! Routine to write a 2D Matrix out as a "MATRIX" . A maximum of ten
! numbers are written across the page.
!
	SUBROUTINE WR2D_GATH_MPI_V1(A,N,DST,DEND,ND,MES,SYMB,WR_INDEX,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 10-Feb-2002: SYMB installed as pass option
!                      Index I & J now written out.
!                      Changed to V2
! Altered 13-Dec-1989 - Implicit none installed. I index written out.
!
	INTEGER N,ND,LU
	INTEGER DST,DEND
	REAL(KIND=LDP) A(N,DST:DEND)
	LOGICAL WR_INDEX
	CHARACTER(LEN=*) MES
	CHARACTER(LEN=*) SYMB
!
! Local variables
!
	REAL(KIND=LDP), ALLOCATABLE ::  MAT(:,:)
        INTEGER, ALLOCATABLE, SAVE :: SCAT_DISP(:)
        INTEGER, ALLOCATABLE, SAVE :: SCAT_SIZE(:)
!
	REAL(KIND=LDP) T1
	INTEGER MS,MF,ML,I,J
	CHARACTER(LEN=80) FORM
	INTEGER, PARAMETER :: IZERO=0
!
	INTEGER K
	INTEGER IERR
	include 'mpif.h'
!
	IF(MYPE .EQ. 0)ALLOCATE(MAT(N,ND))
	IF(. NOT. ALLOCATED(SCAT_DISP))THEN
	  ALLOCATE(SCAT_SIZE(0:NTHREAD-1))
	  ALLOCATE(SCAT_DISP(0:NTHREAD-1))
	END IF
!
	CALL SET_IDISP_ISEND(SCAT_DISP,SCAT_SIZE,ND,NTHREAD)
	SCAT_DISP=SCAT_DISP*N
	SCAT_SIZE=SCAT_SIZE*N
!	
	CALL MPI_GATHERV(A,SCAT_SIZE(MYPE),MPI_DOUBLE_PRECISION,
	1      MAT,SCAT_SIZE,SCAT_DISP,MPI_DOUBLE_PRECISION,IZERO,
	1      MPI_COMM_WORLD,IERR)
!
	IF(MYPE .NE. 0)RETURN
	WRITE(LU,'(/,1X,A)')MES
!
! If WR_INDEX is not set, we simply write the array values.
!
	IF(.NOT. WR_INDEX)THEN
	  MS=1
	  DO ML=0,ND-1,10
	    MF=ML+10
	    IF(MF .GT. ND)MF=ND
	    WRITE(LU,'()')
	    DO I=1,N
	      WRITE(LU,'(1X,10ES12.4)')(MAT(I,J),J=MS,MF)
	    END DO
	    MS=MS+10
	  END DO
	  DEALLOCATE(MAT)
	  RETURN
	END IF
!
! Since WR_INDEX must be set, we need to determine the corect format to use.
!
! We Insert a 0 after the I so that FORM is always exactly 8 characters long
! before we enter the M section.
!
	IF(N .LT. 100)THEN
	   FORM='(1X,I02,'
	ELSE IF(N .LT. 1000)THEN
	   FORM='(1X,I03,'
	ELSE IF(N .LT. 10000)THEN
	   FORM='(1X,I04,'
	ELSE
	   T1=N
	   T1=LOG10(T1)+1
	   FORM='(1X,I'
	   WRITE(FORM(5:6),'(I2.2)')INT(T1)
	   FORM(8:8)=','
	END IF
!
	IF(ND .LT. 10)THEN
	  FORM=TRIM(FORM)//'''('',I1,'')'','
	ELSE IF(ND .LT. 100)THEN
	  FORM=TRIM(FORM)//'''('',I2,'')'','
	ELSE IF(ND.LT. 1000)THEN
	  FORM=TRIM(FORM)//'''('',I3,'')'','
	ELSE IF(ND .LT. 10000)THEN
	  FORM=TRIM(FORM)//'''('',I4,'')'','
	ELSE
	   T1=ND
	   T1=LOG10(T1)+1
	   FORM=TRIM(FORM)//'''('',I'
	   WRITE(FORM(14:15),'(I2.2)')INT(T1)
	   FORM=TRIM(FORM)//','')'','
	END IF
	IF(SYMB(1:1) .NE. ' ')THEN
	  FORM=TRIM(FORM)//''''//SYMB(1:1)//''',1X,10ES12.4)'
	ELSE
	  FORM=TRIM(FORM)//'1X,10ES12.4)'
	END IF
!
	MS=1
	DO ML=0,ND-1,10
	  MF=ML+10
	  IF(MF .GT. ND)MF=ND
	  IF(ML .EQ. 0 .OR. N .GT. 1)WRITE(LU,'()')
	  DO I=1,N
	    WRITE(LU,FORM)I,MS,(MAT(I,J),J=MS,MF)
	  END DO
	  MS=MS+10
	END DO
	DEALLOCATE(MAT)
!
	RETURN
	END
