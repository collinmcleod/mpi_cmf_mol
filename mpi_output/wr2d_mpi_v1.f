!
! Routine to write a 2D Matrix out as a "MATRIX" . A maximum of ten
! numbers are written across the page.
!
	SUBROUTINE WR2D_MPI_V1(A,N,DST,DEND,ND,MES,SYMB,WR_INDEX,LU)
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
	REAL(KIND=LDP), ALLOCATABLE ::  MAT(:,:)
	LOGICAL WR_INDEX
	CHARACTER*(*) MES
	CHARACTER*(*) SYMB
!
	REAL(KIND=LDP) T1
	INTEGER MS,MF,ML,I,J
	CHARACTER(LEN=80) FORM
	INTEGER, PARAMETER :: IZERO=0
!
	LOGICAL, SAVE :: FIRST=.TRUE.
	INTEGER, SAVE :: NUM_DEPTHS_PER_THREAD
	INTEGER IERR
	INTEGER NTRANS
	INTEGER ND_BIG
!
	IF(DST .EQ. 1)THEN
	  ALLOCATE(MAT(N,ND))
	ELSE
	  ALLOCATE(MAT(1,1))
	END IF
	WRITE(6,*)'Calling gather',DST; FLUSH(UNIT=6)
	I=N*ND; CALL GATHER_VEC_TO_ROOT_MPI_V1(A,MAT,I)
	WRITE(6,*)'Called gather',DST; FLUSH(UNIT=6)
!
	IF(DST .NE. 1)THEN
	  DEALLOCATE(MAT)
	  RETURN
	END IF
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
