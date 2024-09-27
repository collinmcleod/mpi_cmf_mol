!
! Module for vectors required when computing the continuum radiation field,
! with blanketing. Some are simply scratch vectors. Others
! are required in CMFGEN_SUB and in the variation routines.
! Vectors required for dJ are in VAR_RAD_MOD.
!
	MODULE LUMINOSITY_MOD
	USE SET_KIND_MODULE
!
	REAL(KIND=LDP), ALLOCATABLE :: dE_WORK_LUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: DEP_RAD_EQ_LUM(:)        !ND - Integrated departure from radiative equilibrium
	REAL(KIND=LDP), ALLOCATABLE :: DIELUM(:)      		!ND - Dielectronic line emission luminosity.
	REAL(KIND=LDP), ALLOCATABLE :: DJDt_LUM(:)      	!ND - DJDT correction to integrated flux.
	REAL(KIND=LDP), ALLOCATABLE :: LLUMST(:)      		!ND - Line luminosity.
	REAL(KIND=LDP), ALLOCATABLE :: MECH_LUM(:)     		!ND - Mechanical luminosity
	REAL(KIND=LDP), ALLOCATABLE :: RAD_DECAY_LUM(:)      	!ND - Luminosity as a function of depth
	REAL(KIND=LDP), ALLOCATABLE :: RLUMST(:)      		!ND - 
	REAL(KIND=LDP), ALLOCATABLE :: SHOCK_POWER_LUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: XRAY_LUM_0p1(:)
	REAL(KIND=LDP), ALLOCATABLE :: XRAY_LUM_1keV(:)
	REAL(KIND=LDP), ALLOCATABLE :: XRAY_LUM_TOT(:)

	END MODULE LUMINOSITY_MOD
!
	SUBROUTINE SET_LUMINOSITY_MOD(ND)
	USE SET_KIND_MODULE
	USE LUMINOSITY_MOD
	IMPLICIT NONE
	INTEGER ND
!
	INTEGER IOS
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	IOS=0
        IF(IOS .EQ. 0)ALLOCATE ( dE_WORK_LUM(ND),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( DEP_RAD_EQ_LUM(ND),STAT=IOS )    	!Depature from radiative equilibrium.
	IF(IOS .EQ. 0)ALLOCATE ( DIELUM(ND),STAT=IOS )    		!Dielectronic line emission luminosity.
	IF(IOS .EQ. 0)ALLOCATE ( DJDt_LUM(ND),STAT=IOS )  	  	!DJDt correction to integrated flux.
	IF(IOS .EQ. 0)ALLOCATE ( LLUMST(ND),STAT=IOS )    		!Line luminosity.
	IF(IOS .EQ. 0)ALLOCATE ( MECH_LUM(ND),STAT=IOS )		!Mechanical luminosity
        IF(IOS .EQ. 0)ALLOCATE ( RAD_DECAY_LUM(ND) )
	IF(IOS .EQ. 0)ALLOCATE ( RLUMST(ND),STAT=IOS )			!Luminosity as a function of depth
        IF(IOS .EQ. 0)ALLOCATE ( SHOCK_POWER_LUM(ND) )
        IF(IOS .EQ. 0)ALLOCATE ( XRAY_LUM_0p1(ND) )
        IF(IOS .EQ. 0)ALLOCATE ( XRAY_LUM_1keV(ND) )
        IF(IOS .EQ. 0)ALLOCATE ( XRAY_LUM_TOT(ND) )
!
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error allocating memory in SET_LUMINOSITY_MOD'
	  WRITE(LUER,*)'STAT=',IOS
	  STOP
	END IF
!
	RETURN
	END
