!
! F90, Computer independent routined to return a 20 character string
! containg a formatted date and time.

	SUBROUTINE DATE_TIME(D_T_STRING)
	IMPLICIT NONE
	CHARACTER*(*) D_T_STRING

	CHARACTER*3 MONTHS(12)
	DATA MONTHS/'Jan','Feb','Mar','Apr','May','Jun',&
	           &'Jul','Aug','Sep','Oct','Nov','Dec'/
	CHARACTER*8 CDATE
	CHARACTER*10 CTIME
	INTEGER ELEMENTS(8)

	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU

	!D_T_STRING=' '
	!IF(LEN(D_T_STRING) .LT. 20)THEN
	!  LUER=ERROR_LU()
	!  WRITE(LUER,*)'String in DATE_TIME must be at least 20 characters'
	!  STOP
	!END IF

	CALL DATE_AND_TIME(DATE=CDATE,TIME=CTIME,VALUES=ELEMENTS)

	D_T_STRING(1:3)=CDATE(7:8)//'-'
	D_T_STRING(4:6)=MONTHS(ELEMENTS(2))
	D_T_STRING(7:12)='-'//CDATE(1:4)//' '
	D_T_STRING(13:14)=CTIME(1:2)
	D_T_STRING(15:18)=':'//CTIME(3:4)//':'
	D_T_STRING(19:20)=CTIME(5:6)

	RETURN
	END subroutine date_time
