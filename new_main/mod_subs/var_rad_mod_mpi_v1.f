!
! This moduls declares vectors and arrays required to compute dJ in the continuum section of
! the code. Storage is allocated by a call to SET_VAR_RAD_MOD_V2.
!
	MODULE VAR_RAD_MOD_MPI_V1
	USE SET_KIND_MODULE
!
	REAL(KIND=LDP), ALLOCATABLE :: DIFFW(:)         !NT - Variation of diffusion approx. at inner boundary.
!
	REAL(KIND=LDP), ALLOCATABLE :: VK(:,:)    	!ND,ND - Coef. matrix of %CHI vector
	REAL(KIND=LDP), ALLOCATABLE :: FC(:,:)    	!ND,ND - Coef. of %EMIS vector in angular equ.
	REAL(KIND=LDP), ALLOCATABLE :: F2DA(:,:)    	!ND,ND - Coef. of %CHi in angular equ.
	REAL(KIND=LDP), ALLOCATABLE :: FA(:)    		!ND
!
	REAL(KIND=LDP), ALLOCATABLE :: F2DAEXT(:,:)    	!NDMAX,NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: FCEXT(:,:)    	!NDMAX,NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: FAEXT(:)    	!NDMAX-
!
	REAL(KIND=LDP), ALLOCATABLE :: dJ_DIF_d_T_EXT(:)          !NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: dJ_DIF_d_dTdR_EXT(:)       !NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: dJ_DIF_d_T(:)              !NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: dJ_DIF_d_dTdR(:)           !NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: RHS_dHdCHI(:,:)            !NDMAX,ND -
	REAL(KIND=LDP), ALLOCATABLE :: dRSQH_DIF_d_T(:)           !NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: dRSQH_DIF_d_dTdR(:)        !NDMAX -
!
	END MODULE VAR_RAD_MOD_MPI_V1
!
! 
! Subroutine to allocate vectors and arrays.
!
	SUBROUTINE SET_VAR_RAD_MOD_MPI_V1(DST,DEND,ND,
	1                NDEXT,NT,NUM_BNDS,NM,MAX_SIM,NM_KI,ACCURATE,ALLOCATE_TX)
	USE SET_KIND_MODULE
	USE VAR_RAD_MOD_MPI_V1
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER DST,DEND
	INTEGER NDEXT
	INTEGER NT
	INTEGER NUM_BNDS
	INTEGER NM
	INTEGER NM_KI
	INTEGER MAX_SIM
	LOGICAL ACCURATE
	LOGICAL ALLOCATE_TX
!
	INTEGER IOS
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	IOS=0
	IF(IOS .EQ. 0)ALLOCATE( DIFFW(NT),STAT=IOS)       !Diffusion variation
	IF(IOS .EQ. 0)ALLOCATE( VK(ND,ND),STAT=IOS)       !Coef. matrix of %CHI vector
	IF(IOS .EQ. 0)ALLOCATE( FC(ND,ND),STAT=IOS)       !Coef. of %EMIS vector in angular equ.
	IF(IOS .EQ. 0)ALLOCATE( F2DA(ND,ND),STAT=IOS)     !Coef. of %CHi in angular equ.
	IF(IOS .EQ. 0)ALLOCATE( FA(ND),STAT=IOS)          !
!
	IF(ACCURATE)THEN
	  IF(IOS .EQ. 0)ALLOCATE( F2DAEXT(NDEXT,NDEXT),STAT=IOS)  	!These arrays don't need to be
	  IF(IOS .EQ. 0)ALLOCATE( FCEXT(NDEXT,NDEXT),STAT=IOS)    	!NDEXT,NDEXT - contiguous as for PERTJD.
	  IF(IOS .EQ. 0)ALLOCATE( FAEXT(NDEXT), STAT=IOS)
	END IF
!
! Variation arrays
!
	IF(IOS .EQ. 0)ALLOCATE ( dJ_DIF_d_T_EXT(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dJ_DIF_d_dTdR_EXT(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dJ_DIF_d_T(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dJ_DIF_d_dTdR(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RHS_dHdCHI(NDEXT,ND),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dRSQH_DIF_d_T(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dRSQH_DIF_d_dTdR(NDEXT),STAT=IOS )
!
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in SET_VAR_RAD_MOD_V2'
	  WRITE(LUER,*)'Unable to allocate required memory'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
	END
