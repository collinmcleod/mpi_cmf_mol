
! Data module for FG_J_CMF. Data placed in this module is automatically
! saved between subroutine calls..
!
	MODULE FG_J_CMF_MOD_MPI_V1
	USE SET_KIND_MODULE
	USE MPI
	IMPLICIT NONE
!
! The *_STORE routines are used to store the radiation field, as computed
! on the previous call to FG_J_CMF.
! The *_PREV routines are used to store the radiation field, as computed
! for the previous frequency.
!
! The *_PREV routines are updated from the *_STORE routines when NEW_FREQ
! is .TRUE.
!
!***************************************************************************
!***************************************************************************
!
	TYPE RAY_DATA
!
	  REAL(KIND=LDP), ALLOCATABLE :: R_RAY(:)
	  REAL(KIND=LDP), ALLOCATABLE :: Z(:)
	  REAL(KIND=LDP), ALLOCATABLE :: GAM(:)
	  REAL(KIND=LDP), ALLOCATABLE :: DTAU(:)
	  REAL(KIND=LDP), ALLOCATABLE :: AV(:)
	  REAL(KIND=LDP), ALLOCATABLE :: CV(:)
	  INTEGER, ALLOCATABLE :: REXT_PNT(:)
!
! Variables specific to short characteristic approach of solving the
! transfer equation.
!
	  REAL(KIND=LDP), ALLOCATABLE :: I_P_PREV(:)
	  REAL(KIND=LDP), ALLOCATABLE :: I_P_STORE(:)
	  REAL(KIND=LDP), ALLOCATABLE :: I_M_PREV(:)
	  REAL(KIND=LDP), ALLOCATABLE :: I_M_STORE(:)
	  REAL(KIND=LDP), ALLOCATABLE :: dGAMdR(:)
!
	  REAL(KIND=LDP), ALLOCATABLE :: A0(:)
	  REAL(KIND=LDP), ALLOCATABLE :: A1(:)
	  REAL(KIND=LDP), ALLOCATABLE :: A2(:)
	  REAL(KIND=LDP), ALLOCATABLE :: A3(:)
	  REAL(KIND=LDP), ALLOCATABLE :: A4(:)
!
	  REAL(KIND=LDP), ALLOCATABLE :: I_P(:)
	  REAL(KIND=LDP), ALLOCATABLE :: I_M(:)
!
	  REAL(KIND=LDP), ALLOCATABLE :: JQW(:)
	  REAL(KIND=LDP), ALLOCATABLE :: HQW(:)
	  REAL(KIND=LDP), ALLOCATABLE :: KQW(:)
	  REAL(KIND=LDP), ALLOCATABLE :: NQW(:)
	  REAL(KIND=LDP), ALLOCATABLE :: HMIDQW(:)
	  REAL(KIND=LDP), ALLOCATABLE :: NMIDQW(:)
!
	  REAL(KIND=LDP) IB_AQW(8)
	  REAL(KIND=LDP) OB_AQW(8)
!
! Dimension ND. These are used to select the correct location on the
! fine ray grids.
!
	  INTEGER, ALLOCATABLE :: J_PNT(:)
	  INTEGER, ALLOCATABLE :: H_PNT(:)
	  INTEGER, ALLOCATABLE :: NI_RAY
!
	END TYPE RAY_DATA
	TYPE(RAY_DATA), ALLOCATABLE :: RAY(:)
!
! PAR_MOM will be used to store all summed moments for each process.
! We use the pointer to point to a specifc location in PAR_MOM. See below.
!
	REAL(KIND=LDP), TARGET, ALLOCATABLE :: PAR_MOM(:)
	REAL(KIND=LDP), ALLOCATABLE :: MOM_STORE(:)
	REAL(KIND=LDP), POINTER :: PAR_JNU(:)
	REAL(KIND=LDP), POINTER :: PAR_HNU(:)
	REAL(KIND=LDP), POINTER :: PAR_KNU(:)
	REAL(KIND=LDP), POINTER :: PAR_NNU(:)
!
	REAL(KIND=LDP), POINTER :: PAR_IB_VEC(:)
	REAL(KIND=LDP), POINTER :: PAR_OB_VEC(:)
!
! Dimensiond ND_EXT,4
!
	REAL(KIND=LDP), ALLOCATABLE :: V_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_COEF(:,:)
!
!
! Dimensioned ND+ND_ADD=ND_EXT
!
	REAL(KIND=LDP), ALLOCATABLE :: R_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_R_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: Z_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_ETA_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_CHI_EXT(:)
!
! Dimensioned NRAY_MAX
!
	REAL(KIND=LDP), ALLOCATABLE :: dCHIdR(:)
	REAL(KIND=LDP), ALLOCATABLE :: Q(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_RAY(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: dCHIdR_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: SOURCE_RAY(:)
!
!
	REAL(KIND=LDP), ALLOCATABLE :: OLDCHI(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLDCHI_STORE(:)
!
!***************************************************************************
!***************************************************************************
!
! Variables specific to short characteristic approach of solving the
! transfer equation.
!
! As these variableis change from freqency to frequency, and are not need
! when the routine is called successively for the same frequencty, we do
! not need to store these for every ray.
!
	REAL(KIND=LDP), ALLOCATABLE :: EE(:)
	REAL(KIND=LDP), ALLOCATABLE :: E0(:)
	REAL(KIND=LDP), ALLOCATABLE :: E1(:)
	REAL(KIND=LDP), ALLOCATABLE :: E2(:)
	REAL(KIND=LDP), ALLOCATABLE :: E3(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: SOURCE_PRIME(:)
	REAL(KIND=LDP), ALLOCATABLE :: S(:)
	REAL(KIND=LDP), ALLOCATABLE :: dS(:)
!
	REAL(KIND=LDP) PREVIOUS_FREQ
	REAL(KIND=LDP) VDOP_FRAC_SAV
!
	INTEGER NRAY_MAX
	INTEGER ND_EXT
	INTEGER ND_ADD
	INTEGER NP_MAX
	INTEGER IDMIN,IDMAX
!
	CHARACTER(LEN=10) OLD_SOLUTION_OPTIONS
	CHARACTER(LEN=10) SOLUTION_METHOD
	LOGICAL INSERT
	LOGICAL FIRST_TIME
	LOGICAL NEW_R_GRID
	LOGICAL WRITE_IP
	LOGICAL DIF_OR_ZF
!
	DATA VDOP_FRAC_SAV/-1001.0_LDP/ 	!Absurd value.
	DATA FIRST_TIME/.TRUE./
	DATA PREVIOUS_FREQ/0.0_LDP/
!
! MPI related data
!
	INTEGER IERR
	INTEGER ERRORCODE
	INTEGER NUM_RAYS_PER_THREAD
!
	SAVE
	END MODULE FG_J_CMF_MOD_MPI_V1
!
!
! 
!
! Routine to compute the Eddington F, G and RSQN_ON_RSQJ Eddington factors
! for a single frequency. The transfer is done in the comoving-frame with
! using the INTEGRAL equation approach.
!
! The INTEGRAL equation approach is more stable. Because monotonic cubic
! interpolation is used for the Source function in the INTEGRAL equation
! approach, the intensities are guaranteed positive (or zero).
!
! The moments J, H, K, and N are computed by the routine. 
!
! NB:
!	F = K / J
!	G=  N / H
!
! Routine also returns I+, so that observers flux can be computed. Note because
! the outer boundary can be thick, I+ is actually I+ - I- and hence is a flux
! like variable (This I+ = 2v at outer boundary).
!
! The radiation field at the previous (bluer) frequency is stored locally.
! Coupling to this frequency is controlled by the variable INIT. If INIT is
! true, the atmosphere is treated with the assumption that V=0.
!
! INIT must be TRUE on the very first call.
!
!
	SUBROUTINE FG_J_CMF_MPI_V1(ETA,CHI,ESEC,V,SIGMA,R,P,JNU,FEDD,
	1                  RETURNED_IN_HBC,RETURNED_OUT_HBC,IPLUS_P,
	1                  FREQ,dLOG_NU,INNER_BND_METH,DBB,IC,
	1                  VDOP_VEC,VDOP_FRAC,REXT_FAC,
	1                  METHOD,SOLUTION_OPTIONS,TRAPFORJ,
	1                  THK,INIT,NEW_FREQ,NC,NP,ND)
	USE SET_KIND_MODULE
	USE MOD_RAY_MOM_STORE
	USE FG_J_CMF_MOD_MPI_V1
	IMPLICIT NONE
!
! Altered 03-Jul-2025 : Fixed bug for NI=2 when using thin outer boundary condition.
! Altered 12-Jun-2025 : Converted from MPI_DOUBLE_PRECISION to MY_MPI_DP
! Altered 18-May-2025 : Improved error checking for outer boundary condition (10-Jun-2025).
! Altered 17-Mar-2025 : Bug fix when checking whether NI=1.
!                       Error messages fixed to point to correct routines.
! Altered 16-Feb-2025 : Fixed bug with computation of NUM_RAYS_PER_THREAD.
!
	INTEGER NC,NP,ND
	REAL(KIND=LDP) ETA(ND),CHI(ND),ESEC(ND)
	REAL(KIND=LDP) V(ND),SIGMA(ND),R(ND),P(NP)
!
! NB: J,H,K,N refer to the first 4 moments of the radiation field.
!
	REAL(KIND=LDP) JNU(ND)
	REAL(KIND=LDP) FEDD(ND)
	REAL(KIND=LDP) IPLUS_P(NP)
!
! VDOP_VEC(I) is the minimum DOPPLER width for all species at depth I. It will include
! both a turbulent, and and thermal contribution for the ionization species with the
! highest mass. VDOP_FRAC is used to set the minimum velocity step size along a ray.
!
	REAL(KIND=LDP) OB_VEC(8)
	REAL(KIND=LDP) IB_VEC(8)
	REAL(KIND=LDP) VDOP_VEC(ND)
	REAL(KIND=LDP) VDOP_FRAC
	REAL(KIND=LDP) REXT_FAC		!Factor to scale RMAX by if thick atmosphere.
	REAL(KIND=LDP) RETURNED_IN_HBC
	REAL(KIND=LDP) RETURNED_OUT_HBC
!
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) IC
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) dLOG_NU
!
	CHARACTER*(*) SOLUTION_OPTIONS
	CHARACTER(LEN=6) METHOD
	CHARACTER(LEN=6) N_TYPE
	CHARACTER(LEN=*) INNER_BND_METH
!
! Use "Thick" boundary condition. at outer boundary. Only noted when INIT
! is true. All subsequent frequencies will use the same boundary condition
! independent of the passed value (Until INIT is set to TRUE again).
!
	LOGICAL TRAPFORJ
	LOGICAL THK
!
! First frequency -- no frequency coupling.

	LOGICAL INIT
!
! Upon leaving this routine the radiation field along each ray is stored. This
! will provide the blue wing information necessary for the next frequency.
! This routine may, however, be used in an iterative loop. In this case the
! "blue wing" information should remain unaltered between calls.
! NEW_FREQ indicates that a new_frequency is being passed, and hence the "blue
! wing" information should be updated.
!
	LOGICAL NEW_FREQ	
!
	INTEGER, PARAMETER :: ND_ADD_MAX=24
	INTEGER, PARAMETER :: NC_PNT_SRCE=2
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
!
	INTEGER, SAVE :: FREQ_CNT
	INTEGER ACCESS_F
        INTEGER LU_IP
	INTEGER IOS,REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC
!
! The following arrays do not need to be stored, and hence can be created
! each time.
!
	REAL(KIND=LDP) CV_BOUND(NP)		!Outer boundary V
	REAL(KIND=LDP) I_M_IN_BND(NP)		!Inner boundary
	REAL(KIND=LDP) IBOUND(NP)		!Incident intensity on outer boundary.
!
	INTEGER N_ERR_MAX,FG_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL(KIND=LDP) FG_ERR_ON_FREQ
	INTEGER FG_ERR_TYPE
	COMMON /FG_J_CMF_ERR/FG_ERR_ON_FREQ(N_ERR_MAX),FG_ERR_TYPE(N_ERR_MAX),FG_ERR_CNT
	LOGICAL NEG_AV_VALUE
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL(KIND=LDP), PARAMETER :: ONE=1
	INTEGER, PARAMETER :: NINS=4
!
	LOGICAL, PARAMETER :: LFALSE=.FALSE.
	LOGICAL, PARAMETER :: LTRUE=.TRUE.
!
	INTEGER NI_SMALL
	INTEGER I,J,K,LS
	INTEGER IPROC
	INTEGER NI
	INTEGER NP_TMP
!
	REAL(KIND=LDP) DBC
	REAL(KIND=LDP) I_CORE
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) DELTA_Z
	REAL(KIND=LDP) ALPHA
	REAL(KIND=LDP) ESEC_POW
	REAL(KIND=LDP) BETA
	REAL(KIND=LDP) VINF
	REAL(KIND=LDP) RMAX,DEL_R_FAC
	REAL(KIND=LDP) MU,dZ,PSQ
	REAL(KIND=LDP) DEL_R
!
	LOGICAL, PARAMETER :: DEFINE_AT_MID_POINTS=.FALSE.    !TRUE.
!
! Include MPI definitions.
!
!	INCLUDE 'mpif.h'
!
! Funtion to return LS as a function of MYPE an I (=IPROC loop counter).
! For example, with 4 process, we have
!
! Pocessor   0, 1, 2, 3, 0, 1, 2, 3, 0,  1,  2,  3,
! LS         1, 2, 3, 4, 8, 7, 6, 5, 9, 10, 11, 12, etc
!
	INTEGER GET_LS
	GET_LS(MYPE,NTHREAD,I)=MOD(I,2)*(MYPE+(I-1)*NTHREAD+1)+MOD(I+1,2)*(I*NTHREAD-MYPE)
!
!
	CALL TUNE(1,'FG_FULL')
!
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NTHREAD,IERR)
	NUM_RAYS_PER_THREAD=(NP+NTHREAD-1)/NTHREAD
!
	IF(FIRST_TIME .AND. MYPE .EQ. 0)THEN
	  WRITE(6,'(3(2X,A3,I5))')'NC=',NC,'ND=',ND,'NP=',NP
	  WRITE(6,*)'R(1), R(ND)=',R(1),R(ND)
	  WRITE(6,*)'P(1), P(NP)=',P(1),P(NP)
	  WRITE(6,*)'Using trapazoidal rules for quadrature weights?',TRAPFORJ
	END IF
!
! Allocate data for moments which will be used to construct the Eddington factors.
!
	LUER=ERROR_LU()
	IF(.NOT. ALLOCATED(JNU_STORE))THEN
	  ND_STORE=ND
	  ALLOCATE (JNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (KNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (R_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_REL_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RMID_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EXT_RMID_STORE(ND+1),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating JNU_STORE block in FG_J_CMF_MPI_V1: Status=',IOS
	    CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	    STOP	
          END IF
	  ALLOCATE(RAY(NP))
	END IF
!
! This routine is messy because well can extra points to each ray to increase the accuracy of
! the calculation. The addition of extra points is ray depndent.
!
	IF(INIT)THEN
	  R_STORE(1:ND)=R(1:ND)
	  GAM_REL_STORE(1:ND)=1.0_LDP/SQRT(1.0_LDP-(V(1:ND)/2.99792458E+10_LDP)**2)
	  DO I=1,ND-1
	    RMID_STORE(I)=0.5_LDP*(R(I)+R(I+1))
	    EXT_RMID_STORE(I+1)=0.5_LDP*(R(I)+R(I+1))
	  END DO
	  EXT_RMID_STORE(1)=R(1); EXT_RMID_STORE(ND+1)=R(ND)
	END IF
!
! Check to see whether we have a new R grid, or the solution options have
! changed. This can only happen when INIT is TRUE.
!
	NEW_R_GRID=.FALSE.
	IF(INIT .AND. .NOT. FIRST_TIME)THEN
	  DO I=1,ND
	    IF(R(I) .NE. R_EXT(ND_ADD+I))THEN
	      NEW_R_GRID=.TRUE.
	      J=ERROR_LU()
	      IF(MYPE .EQ. 0)WRITE(J,*)'Warning: Updating RGRID in FG_J_CMF_MPI_V1'
	      EXIT
	    END IF
	  END DO
	  IF(VDOP_FRAC .NE. VDOP_FRAC_SAV
	1       .OR.  OLD_SOLUTION_OPTIONS .NE. SOLUTION_OPTIONS)NEW_R_GRID=.TRUE.
	END IF
!
	IF(NEW_R_GRID .OR. (INIT .AND. FIRST_TIME))THEN
	  IF(MYPE .EQ. 0)WRITE(6,*)'About to compute angular quadratue weights'; FLUSH(UNIT=6)
	  CALL SET_ANG_QW_MPI(R,P,NC,ND,NP,TRAPFORJ)
	END IF
	IF(FIRST_TIME .AND. MYPE .EQ. 0)THEN
	  WRITE(6,*)'Ray assignments to thread 0 for ray integrations'
	  WRITE(6,'(4A10)')'LS','MYPE','NTHREAD','Loop'
	END IF
!
	DO IPROC=1,NUM_RAYS_PER_THREAD
	  LS=GET_LS(MYPE,NTHREAD,IPROC)
	  IF(FIRST_TIME .AND. MYPE .EQ. 0)THEN
	    WRITE(6,'(4I10)')LS,MYPE,NTHREAD,IPROC; FLUSH(UNIT=6)
	  END IF
	  IF(LS .GT. NP)EXIT
!
! This association was done to save time, and could probably beimproved.
!
	  RAY(LS)%OB_AQW(1)=RAY(LS)%JQW(1)
	  RAY(LS)%OB_AQW(2)=RAY(LS)%HQW(1)
	  RAY(LS)%OB_AQW(3)=RAY(LS)%KQW(1)
	  RAY(LS)%OB_AQW(4)=RAY(LS)%NQW(1)
	  RAY(LS)%OB_AQW(5)=RAY(LS)%JQW(1)
	  RAY(LS)%OB_AQW(6)=RAY(LS)%HQW(1)
	  RAY(LS)%OB_AQW(7)=RAY(LS)%KQW(1)
	  RAY(LS)%OB_AQW(8)=RAY(LS)%NQW(1)
	    
	  RAY(LS)%IB_AQW(1)=RAY(LS)%JQW(ND)
	  RAY(LS)%IB_AQW(2)=RAY(LS)%HQW(ND)
	  RAY(LS)%IB_AQW(3)=RAY(LS)%KQW(ND)
	  RAY(LS)%IB_AQW(4)=RAY(LS)%NQW(ND)
	  RAY(LS)%IB_AQW(5)=RAY(LS)%JQW(ND)
	  RAY(LS)%IB_AQW(6)=RAY(LS)%HQW(ND)
	  RAY(LS)%IB_AQW(7)=RAY(LS)%KQW(ND)
	  RAY(LS)%IB_AQW(8)=RAY(LS)%NQW(ND)
	END DO
!	
! Deallocate all allocated rays if we are using a diferent solution technique.
! This option will only be used when testing, since in CMFGEN we will always use
! the same atmospheric structure.
!
	IF( ALLOCATED(R_EXT) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R_EXT )
	  DEALLOCATE ( LOG_R_EXT )
	  DEALLOCATE ( V_EXT )
	  DEALLOCATE ( Z_EXT )
	  DEALLOCATE ( SIGMA_EXT )
	  DEALLOCATE ( ETA_EXT )
	  DEALLOCATE ( CHI_EXT )
	  DEALLOCATE ( LOG_ETA_EXT )
	  DEALLOCATE ( LOG_CHI_EXT )
!
	  DEALLOCATE ( EE )
	  DEALLOCATE ( E0 )
	  DEALLOCATE ( E1 )
	  DEALLOCATE ( E2 )
	  DEALLOCATE ( E3 )
!
	  DEALLOCATE ( SOURCE_PRIME )
	  DEALLOCATE ( S )
	  DEALLOCATE ( dS )
	  DEALLOCATE ( MOM_STORE )
	  DEALLOCATE ( PAR_MOM )
!
	  DEALLOCATE ( V_COEF )
	  DEALLOCATE ( SIGMA_COEF )
	  DEALLOCATE ( ETA_COEF )
	  DEALLOCATE ( CHI_COEF )
	  DEALLOCATE ( V_RAY )
	  DEALLOCATE ( SIGMA_RAY )
	  DEALLOCATE ( ETA_RAY )
	  DEALLOCATE ( CHI_RAY )
	  DEALLOCATE ( dCHIdR_RAY )
	  DEALLOCATE ( SOURCE_RAY )
	  DEALLOCATE ( dCHIdR )
	  DEALLOCATE ( Q )
!
	  DO IPROC=1,NUM_RAYS_PER_THREAD
	    LS=GET_LS(MYPE,NTHREAD,IPROC)
	    IF(LS .GT. NP)EXIT
!
! These arrays are only used when using the INTEGERAL apporach.
!
	    DEALLOCATE ( RAY(LS)%Z )
	    DEALLOCATE ( RAY(LS)%GAM )
	    DEALLOCATE ( RAY(LS)%DTAU )
	    DEALLOCATE ( RAY(LS)%REXT_PNT )
	    DEALLOCATE ( RAY(LS)%R_RAY )
	    DEALLOCATE ( RAY(LS)%NI_RAY )
	    DEALLOCATE ( RAY(LS)%J_PNT )
	    DEALLOCATE ( RAY(LS)%H_PNT )
!
	    DEALLOCATE ( RAY(LS)%I_P )
	    DEALLOCATE ( RAY(LS)%I_M )
!
	    DEALLOCATE ( RAY(LS)%I_P_PREV )
	    DEALLOCATE ( RAY(LS)%I_P_STORE )
	    DEALLOCATE ( RAY(LS)%I_M_PREV )
	    DEALLOCATE ( RAY(LS)%I_M_STORE )
	    DEALLOCATE ( RAY(LS)%dGAMdR )
	    DEALLOCATE ( RAY(LS)%AV )
	    DEALLOCATE ( RAY(LS)%CV )
!
	    DEALLOCATE ( RAY(LS)%A0 )
	    DEALLOCATE ( RAY(LS)%A1 )
	    DEALLOCATE ( RAY(LS)%A2 )
	    DEALLOCATE ( RAY(LS)%A3 )
	    DEALLOCATE ( RAY(LS)%A4 )
!
	  END DO
	  IF(MYPE .EQ. 0)WRITE(6,*)'Done deallocations'; FLUSH(UNIT=6)
	END IF
	VDOP_FRAC_SAV=VDOP_FRAC
!
! 
!
! Set up the revised grid to improve computational accuracy. Unles we are
! carrying out tests, these need only be constructed once.
!
! Note: MPI_ABORT should be sued when stoping a process due to an error.
! It stops all processors.
!
	IF(FIRST_TIME .OR. .NOT. ALLOCATED(R_EXT) )THEN
!
	  OLD_SOLUTION_OPTIONS=SOLUTION_OPTIONS
	  IF(SOLUTION_OPTIONS .EQ. 'INT/INS')THEN
	    SOLUTION_METHOD='INTEGRAL'
	    INSERT=.TRUE.
	  ELSE IF(SOLUTION_OPTIONS(1:3) .EQ. 'INT')THEN
	    SOLUTION_METHOD='INTEGRAL'
	    INSERT=.FALSE.
	  ELSE
	    IF(MYPE .EQ. 0)THEN
	      J=ERROR_LU()
	      WRITE(J,*)'Error in FG_J_CMF_MPI_V1: Invalid solution option'
	      WRITE(J,*)SOLUTION_OPTIONS
	      CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	      STOP
	    END IF
	  END IF
!
	  IF(INNER_BND_METH .NE. 'DIFFUSION' .AND.
	1     INNER_BND_METH .NE. 'ZERO_FLUX' .AND.
	1     INNER_BND_METH .NE. 'PNT_SRCE' .AND.
	1     INNER_BND_METH .NE. 'SCHUSTER')THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in FG_J_CMF_MPI_V1: Invalid inner boundary condition'
	    WRITE(J,*)INNER_BND_METH
	    CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	    STOP
	  END IF
!
          ND_ADD=0
          IF(THK)ND_ADD=ND_ADD_MAX
          ND_EXT=ND+ND_ADD
!
	  ALLOCATE ( R_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_R_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( V_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( Z_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_ETA_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_CHI_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating R_EXT block in FG_J_CMF_MPI_V1: Status=',IOS
	    CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	    STOP	
          END IF
!
! Compute the extended R grid, excluding inserted points.
!
	  DO I=1,ND
	    R_EXT(ND_ADD+I)=R(I)
	  END DO
	  IF(THK)THEN
	    IF(REXT_FAC .GT. 1.0_LDP .AND. REXT_FAC .LT. 10.0_LDP)THEN
	       RMAX=REXT_FAC*R(1)
	    ELSE IF(V(ND) .LT. 10.0_LDP .AND. R(1)/R(ND) .GE. 9.99_LDP)THEN
	      RMAX=10.0_LDP*R(1)		!Stellar wind
	    ELSE IF(V(ND) .GT. 10.0_LDP .OR. V(1) .GT. 2.0E+04_LDP)THEN
	      RMAX=1.5_LDP*R(1)		!SN model
	    ELSE
	      RMAX=MIN(10.0_LDP,SQRT(R(1)/R(ND)))*R(1)
	    END IF
	    ALPHA=R(1)+(R(1)-R(2))
	    DEL_R_FAC=EXP( LOG(RMAX/ALPHA)/(ND_ADD-3) )
	    R_EXT(1)=RMAX
	    R_EXT(4)=RMAX/DEL_R_FAC
	    R_EXT(2)=R_EXT(1)-0.1_LDP*(R_EXT(1)-R_EXT(4))
	    R_EXT(3)=R_EXT(1)-0.4_LDP*(R_EXT(1)-R_EXT(4))
	    DO I=5,ND_ADD-1
	      R_EXT(I)=R_EXT(I-1)/DEL_R_FAC
	    END DO
	    R_EXT(ND_ADD)=ALPHA
	  END IF
!
! Compute VEXT and R_EXT. We assume a BETA velocity law at large R.
!
	  V_EXT(ND_ADD+1:ND_EXT)=V(1:ND)
	  SIGMA_EXT(ND_ADD+1:ND_EXT)=SIGMA(1:ND)
	  IF(THK)THEN
	    BETA=(SIGMA(1)+1.0_LDP)*(R(1)/R(ND)-1.0_LDP)
            VINF=V(1)/(1-R(ND)/R(1))**BETA
	    DO I=1,ND_ADD
	      V_EXT(I)=VINF*(1.0_LDP-R_EXT(ND_EXT)/R_EXT(I))**BETA
	      SIGMA_EXT(I)=BETA/(R_EXT(I)/R_EXT(ND_EXT)-1.0_LDP)-1.0_LDP
	    END DO
	    J=ERROR_LU()
	    IF(MYPE .EQ. 0)THEN
	      WRITE(J,'(A)')' '
	      WRITE(J,*)'Using thick boundary condition in FG_J_CMF_MPI_V1'
	      WRITE(J,'(2(A,ES16.8,3X))')' R(1)=',R(1),'RMAX=',RMAX
	      WRITE(J,'(2(A,ES16.8,3X))')' V(1)=',V(1),'VMAX=',V_EXT(1)
	    END IF
	  ELSE
	    IF(MYPE .EQ. 0)THEN
	      WRITE(J,'(A)')' '
	      WRITE(J,*)'Using thin boundary condition in FG_J_CMF_MPI_V1'
	    END IF
	  END IF
!
! Define zone used to extrapolate opacities and emissivities.
!
	  IF(ND_ADD .NE. 0)THEN
	    IDMIN=1
	    T1=R_EXT(1)/R(1)
	    IF(T1 .LT. R(1)/R(ND/3))THEN
	      T1=MIN(3.0_LDP,T1)
	      DO I=1,ND
	        IF(R(1)/R(I) .GT. T1)THEN
	          IDMAX=MAX(4,I)
	          EXIT
	        END IF
	      END DO
	      IDMAX=MIN(IDMAX,ND/6)
	    ELSE
	      IDMAX=MIN(10,ND/6)
	    END IF
	    IF(MYPE .EQ. 0)THEN
	      WRITE(6,'(A,I4,A)')' Using depth 1 and',IDMAX,' to extrapolate opacities'
	      WRITE(6,*)'Done EXT allocations';FLUSH(UNIT=6)
	    END IF
	  END IF
!
! Comnpute the extend Z grid. Then work out the maximim number of points per ray. 
! This will allow us to allocate the required memory.
!
	  NRAY_MAX=0
	  DO IPROC=1,NUM_RAYS_PER_THREAD
	    LS=GET_LS(MYPE,NTHREAD,IPROC)
	    IF(LS .GT. NP)EXIT
	    NI_SMALL=ND_EXT-(LS-NC-1)
	    IF(LS .LE. NC+1)NI_SMALL=ND_EXT
!
	    DO I=1,NI_SMALL
	      IF(R_EXT(I) .EQ. P(LS))THEN
	        Z_EXT(I)=0.0_LDP
	      ELSE
	        Z_EXT(I)=SQRT(R_EXT(I)*R_EXT(I)-P(LS)*P(LS))
	      END IF
	    END DO
!
	    K=1
	    T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:NI_SMALL-ND_ADD))
	    DO I=1,NI_SMALL-1
	      T1=(Z_EXT(I)*V_EXT(I)/R_EXT(I)-Z_EXT(I+1)*V_EXT(I+1)/R_EXT(I+1))/T2
	      IF(.NOT. THK .AND. NI_SMALL .EQ. 2)T1=MAX(2.01_LDP,T1)
	      IF(T1 .GT. 1.0_LDP)K=K+INT(T1)
	      K=K+1
	    END DO
            NRAY_MAX=MAX(NRAY_MAX,K)
!
	    ALLOCATE ( RAY(LS)%Z(K),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE ( RAY(LS)%REXT_PNT(K),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE ( RAY(LS)%R_RAY(K),STAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error allocating Z block in FG_J_CMF_MPI_V1: Status=',IOS
	      CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	      STOP	
            END IF
!
! 
!
! Define the grid along each ray. NB: REXT_PNT is used to indicate the
! interpolation interval on the original (EXT) grid.
!
	    K=1
	    RAY(LS)%Z(1)=Z_EXT(1)
	    RAY(LS)%REXT_PNT(1)=1
	    T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:NI_SMALL-ND_ADD))
	    DO I=1,NI_SMALL-1
	      T1=(Z_EXT(I)*V_EXT(I)/R_EXT(I)-Z_EXT(I+1)*V_EXT(I+1)/R_EXT(I+1))/T2
	      IF(.NOT. THK .AND. NI_SMALL .EQ. 2)T1=MAX(2.01_LDP,T1)
	      IF(T1 .GT. 1.0_LDP)THEN
	        DELTA_Z=(Z_EXT(I+1)-Z_EXT(I))/(INT(T1)+1)
	        DO J=1,INT(T1)
	          K=K+1
	          RAY(LS)%Z(K)=RAY(LS)%Z(K-1)+DELTA_Z
	          RAY(LS)%REXT_PNT(K)=I
	        END DO
	      END IF
	      K=K+1
	      RAY(LS)%Z(K)=Z_EXT(I+1)
	      RAY(LS)%REXT_PNT(K)=I
	    END DO
	    RAY(LS)%NI_RAY=K
!
	    PSQ=P(LS)*P(LS)
	    RAY(LS)%R_RAY(1)=R_EXT(1)
	    DO I=2,RAY(LS)%NI_RAY
	      RAY(LS)%R_RAY(I)=SQRT(RAY(LS)%Z(I)*RAY(LS)%Z(I)+PSQ)
	    END DO
	    RAY(LS)%R_RAY(RAY(LS)%NI_RAY)=R_EXT(NI_SMALL)
!
! 
!
! J_PNT and H_PNT are used to indicate the postion of J(I) amd H(I) along the ray
! so that J and H can be computed on our regular radius grid.
!
	    ALLOCATE ( RAY(LS)%J_PNT(ND),STAT=IOS); RAY(LS)%J_PNT(:)=0
	    IF(IOS .EQ. 0)ALLOCATE ( RAY(LS)%H_PNT(ND),STAT=IOS); RAY(LS)%H_PNT(:)=0
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error allocating J_PNT block in FG_J_CMF_MPI_V1: Status=',IOS
	      CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	      STOP	
            END IF
!
! LS=NP with only 1 depth point must be treated separately.
!
	    IF(RAY(LS)%NI_RAY .EQ. 1)THEN
	      RAY(LS)%J_PNT=1
	    ELSE
	      NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	      K=1
	      DO I=1,NI_SMALL
	        DO WHILE(RAY(LS)%J_PNT(I) .EQ. 0)
	          IF(R(I) .LE. RAY(LS)%R_RAY(K) .AND. R(I) .GE. RAY(LS)%R_RAY(K+1))THEN
	            IF( (RAY(LS)%R_RAY(K)-R(I)) .LT. (R(I)-RAY(LS)%R_RAY(K+1)) )THEN
	              RAY(LS)%J_PNT(I)=K
	            ELSE
	              RAY(LS)%J_PNT(I)=K+1
	            END IF
	          ELSE
	            K=K+1
	          END IF
	        END DO
	      END DO
	    END IF
!
	    DO I=1,NI_SMALL
	      K=RAY(LS)%J_PNT(I)
	      IF(RAY(LS)%J_PNT(I) .LT. 1 .OR. RAY(LS)%J_PNT(I) .GT. RAY(LS)%NI_RAY)THEN
	        WRITE(LUER,*)'Error setting J_PNT in FG_J_CMF_MPI_V1 -- invalid values: MYPE=',MYPE
	        WRITE(LUER,*)'Depth=',I,'Ray=',LS,'J_PNT value=',RAY(LS)%J_PNT(I)
	        WRITE(LUER,*)'NP=',NP,'R(1)=',R(1),'R_RAY(1,LS)=',RAY(LS)%R_RAY(1)
	        WRITE(LUER,*)'NI_SMALL=',NI_SMALL
	        FLUSH(LUER)
	        CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	        STOP
	      ELSE IF( ABS(RAY(LS)%R_RAY(K)-R(I))/R(I) .GT. 1.0E-12_LDP)THEN
	        WRITE(LUER,*)'Error setting J_PNT in FG_JCMF_MPI_V1 -- invalid values: MYPE=',MYPE
	        WRITE(LUER,*)'Fractional difference is',ABS(RAY(LS)%R_RAY(K)-R(I))/R(I)
	        WRITE(LUER,*)'Depth=',I,'Ray=',LS,'J_PNT value=',RAY(LS)%J_PNT(I)
	        WRITE(LUER,*)'RAY(LS)%R_RAY(K),R(I)=',RAY(LS)%R_RAY(K),R(I)
	        WRITE(LUER,*)'NI_SMALL=',NI_SMALL,'NC,ND=',NC,ND
	        FLUSH(LUER)
	        CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	        STOP
	      END IF
	    END DO
!
! CV is defined on the nodes. Need to interpolate to the midpoint of R.
!
	    K=1
	    DO I=1,NI_SMALL-1
	      DO WHILE(RAY(LS)%H_PNT(I) .EQ. 0)
	        T1=0.5_LDP*(R(I)+R(I+1))
	        IF(T1 .LE. RAY(LS)%R_RAY(K) .AND. T1 .GE. RAY(LS)%R_RAY(K+1))THEN
	          RAY(LS)%H_PNT(I)=K
	        ELSE
	          K=K+1
	        END IF
	      END DO
	    END DO
	  END DO
!
! 
!
	  ALLOCATE ( V_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_COEF(ND_EXT,4),STAT=IOS )
!
!***************************************************************************
!***************************************************************************
!
	  K=NRAY_MAX
	  IF(IOS .EQ. 0)ALLOCATE ( V_RAY(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_RAY(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_RAY(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_RAY(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( dCHIdR_RAY(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SOURCE_RAY(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( dCHIdR(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( Q(K),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating V_RAY etc in FG_J_CMF_MPI_V1: Status=',IOS
	    CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	    STOP	
          END IF
!
	  K=NRAY_MAX
	  IF(IOS .EQ. 0)ALLOCATE ( EE(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( E0(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( E1(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( E2(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( E3(K),STAT=IOS )

	  IF(IOS .EQ. 0)ALLOCATE (SOURCE_PRIME(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE (S(K),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE (dS(K),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating INTEGRAL block in FG_J_CMF_MPI_V1: Status=',IOS
	    CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	    STOP	
          END IF
!
! Set up stodarge and pointers for computation of the moments, and partial moments
! at the outer boudaries.
!
	  IOS=0
	  IF(IOS .EQ. 0)ALLOCATE ( PAR_MOM(4*ND+16),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( MOM_STORE(4*ND+16),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error -- unable to allocate PAR_MOM or MOM_STORE: Error=',IOS
	    CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	    CALL MPI_FINALIZE(IERR)
	    STOP
	  END IF
	  PAR_JNU=>PAR_MOM(1:ND)
	  PAR_HNU=>PAR_MOM(ND+1:2*ND)
	  PAR_KNU=>PAR_MOM(2*ND+1:3*ND)
	  PAR_NNU=>PAR_MOM(3*ND+1:4*ND)
	  PAR_IB_VEC=>PAR_MOM(4*ND+1:4*ND+8)
	  PAR_OB_VEC=>PAR_MOM(4*ND+9:4*ND+16)
!
!***************************************************************************************
!
! Allocation of ray variables. The use of te funtion GET_LS ensures we always allocate
! the same processor to the same rays.
!
	  DO IPROC=1,NUM_RAYS_PER_THREAD
	    LS=GET_LS(MYPE,NTHREAD,IPROC)
	    IF(LS .GT. NP)EXIT
	    NI=RAY(LS)%NI_RAY 
	    ALLOCATE( RAY(LS)%I_P(NI),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%I_M(NI),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%I_P_PREV(NI),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%I_P_STORE(NI),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%I_M_PREV(NI),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%I_M_STORE(NI),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%dGAMdR(NI),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%A0(NI),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%A1(NI),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%A2(NI),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%A3(NI),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE( RAY(LS)%A4(NI),STAT=IOS )
!
	    K=RAY(LS)%NI_RAY
	    IF(IOS .EQ. 0)ALLOCATE ( RAY(LS)%AV(K),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( RAY(LS)%CV(K),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( RAY(LS)%DTAU(K),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( RAY(LS)%GAM(K),STAT=IOS )
!
	  END DO
!
	  FIRST_TIME=.FALSE.
	ELSE
	  IF(OLD_SOLUTION_OPTIONS .NE. SOLUTION_OPTIONS)THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in FG_J_CMF_MPI_V1'
	    WRITE(J,*)'Can''t switch SOLUTION_OPTIONS while runing code'
	    WRITE(J,*)'New setting:',SOLUTION_OPTIONS
	    WRITE(J,*)'Old setting:',OLD_SOLUTION_OPTIONS
	    CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	    STOP
	  END IF
	END IF
!
!
!
	NEG_AV_VALUE=.FALSE.
!
! Perform initializations.
!
	IF(INIT)THEN
!
! Insert extra points into radius grid. Not all points will be used along a ray.
! We only insert additional points in the interval between Z(NI)=0 and Z(NI-1),
! and between Z(NI-1) and  Z(NI-2).
!
	  FG_ERR_ON_FREQ(:)=0.0_LDP
	  FG_ERR_TYPE(:)=0
	  FG _ERR_CNT=0
!
	  LOG_R_EXT(1:ND_EXT)=LOG(R_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(V_COEF,V_EXT,LOG_R_EXT,ND_EXT)
	  CALL MON_INT_FUNS_V2(SIGMA_COEF,SIGMA_EXT,LOG_R_EXT,ND_EXT)
!
	  DO I=1,ND
	    IF(SIGMA(I) .LT. -1.0_LDP)THEN
	      WRITE(LUER,*)'Warnining Error in FG_J_CMF_MPI_V1 - SIGMA .LT. -1.0D0'
	      WRITE(LUER,*)I,SIGMA(I)
	      EXIT
	    END IF
	  END DO
!
	  DO IPROC=1,NUM_RAYS_PER_THREAD
	    LS=GET_LS(MYPE,NTHREAD,IPROC)
	    IF(LS .GT. NP)EXIT
	    NI=RAY(LS)%NI_RAY  
	    RAY(LS)%I_P_PREV(:)=0.0_LDP
	    RAY(LS)%I_M_PREV(:)=0.0_LDP
!
! The check on SIGMA is to allow for the possibility of a non-monotonic velcoity law in the
! region where V is small, and hece where SIGAMA is unimportant.
!
	    DO I=1,RAY(LS)%NI_RAY
	      K=RAY(LS)%REXT_PNT(I)
	      T1=LOG(RAY(LS)%R_RAY(I)/R_EXT(K))
	      V_RAY(I)=((V_COEF(K,1)*T1+V_COEF(K,2))*T1+V_COEF(K,3))*T1+V_COEF(K,4)
	      SIGMA_RAY(I)=((SIGMA_COEF(K,1)*T1+SIGMA_COEF(K,2))*T1+SIGMA_COEF(K,3))*T1+SIGMA_COEF(K,4)
	      IF(SIGMA_RAY(I) .LE. -1.0_LDP)SIGMA_RAY(I)=-0.999_LDP
	    END DO
!
! Compute GAMMA. This section is straight from the subroutine GAMMA, except
! That _EXT has been added to V, SIGMA, and R.
!
! We assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
!  	    (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	    NI=RAY(LS)%NI_RAY
	    DO I=1,RAY(LS)%NI_RAY
	      MU=RAY(LS)%Z(I)/RAY(LS)%R_RAY(I)
	      T1=3.33564E-06_LDP*V_RAY(I)/RAY(LS)%R_RAY(I)
	      RAY(LS)%GAM(I)=T1*( 1.0_LDP+SIGMA_RAY(I)*(MU**2) )
	    END DO
!
	    DO I=1,RAY(LS)%NI_RAY
	      MU=RAY(LS)%Z(I)/RAY(LS)%R_RAY(I)
	      T1=3.33564E-06_LDP*V_RAY(I)/RAY(LS)%R_RAY(I)
	      RAY(LS)%dGAMdR(I)=RAY(LS)%GAM(I)*SIGMA_RAY(I)/RAY(LS)%R_RAY(I)
	      J=MAX(I-1,1); K=MIN(NI,I+1)
	      RAY(LS)%dGAMdR(I)=RAY(LS)%dGAMdR(I)+T1*(MU**2)*
	1           (SIGMA_RAY(J)-SIGMA_RAY(K))/(RAY(LS)%R_RAY(J)-RAY(LS)%R_RAY(K))
	    END DO
!
	  END DO		!LS Loop
!
!
	ELSE IF(NEW_FREQ)THEN
	  DO IPROC=1,NUM_RAYS_PER_THREAD
	    LS=GET_LS(MYPE,NTHREAD,IPROC)
	    IF(LS .GT. NP)EXIT
	    NI=RAY(LS)%NI_RAY 
	    RAY(LS)%I_P_PREV(:)=RAY(LS)%I_P_STORE(:)
	    RAY(LS)%I_M_PREV(:)=RAY(LS)%I_M_STORE(:)
	  END DO 
	  IF(FREQ .GE. PREVIOUS_FREQ)THEN
	    WRITE(ERROR_LU(),*)'Error in FG_J_CMF_MPI_V1'
	    WRITE(ERROR_LU(),*)'Frequencies must be monotonically decreasng'
	    CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	    STOP
	  END IF
	ELSE IF(FREQ .NE. PREVIOUS_FREQ)THEN
	   WRITE(ERROR_LU(),*)'Error in FG_J_CMF_MPI_V1'
	   WRITE(ERROR_LU(),*) 'Frequencies must not change if NEW_FREQ=.FALSE.'
	   CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	   STOP
	END IF
	PREVIOUS_FREQ=FREQ
!
! 
!
! Compute CHI_EXT, and ETA_EXT. CHI_EXT could be saved, as it doesn't change
! during the iteration procedure. ETA does however change (since it depends
! of J) and thus ETA_EXT must be re-computed on each entry.
!
! The first checks whether we may have negative line opacities due
! to stimulated emission. In such a case we simply assume an 1/r^2
! extrapolation.
!
! We also interpolate in ESEC, since ESEC (in the absence of negative
! absorption) provides a lower bound to the opacity. NB: When CHI is much
! larger then ESEC its variation with r dominates, and it is possible to
! extrapolate CHI below ESEC.
!
	CALL TUNE(1,'FG_CHI_BEG')
	IF(ND_ADD .NE. 0)THEN
	  IF(CHI(IDMIN) .LE. ESEC(IDMIN) .OR. CHI(IDMAX) .LE. ESEC(IDMAX))THEN
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.0_LDP)ESEC_POW=2.0_LDP
	    DO I=1,ND_ADD
	      CHI_EXT(I)=CHI(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      dCHIdR(I)= -ESEC_POW*CHI_EXT(I)/R_EXT(I)
	    END DO
	  ELSE
	    ALPHA=LOG( (CHI(IDMAX)-ESEC(IDMAX)) / (CHI(IDMIN)-ESEC(IDMIN)) )
	1          /LOG(R(IDMIN)/R(IDMAX))
	    IF(ALPHA .LT. 2.0_LDP)ALPHA=2.0_LDP
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.0_LDP)ESEC_POW=2.0_LDP
   	    DO I=1,ND_ADD
	      T1=(CHI(IDMIN)-ESEC(IDMIN))*(R(IDMIN)/R_EXT(I))**ALPHA
	      T2=ESEC(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      CHI_EXT(I)=T1+T2
	      dCHIdR(I)=(-ALPHA*T1-ESEC_POW*T2)/R_EXT(I)
	    END DO
	  END IF
	  DO I=ND_ADD+1,ND_EXT
	    CHI_EXT(I)=CHI(I-ND_ADD)
	  END DO
!
! We limit alpha to 3.5 to avoid excess envelope emission. If alpha were
! 3 we would get a logarithmic flux divergence as we increase the volume.
!
	  ALPHA=LOG(ETA(IDMAX)/ETA(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	  IF(ALPHA .LT. 3.5_LDP)ALPHA=3.5_LDP
	  DO I=1,ND_ADD
	    ETA_EXT(I)=ETA(IDMIN)*(R(IDMIN)/R_EXT(I))**ALPHA
	    IF(ETA_EXT(I) .LE. 1.0E-280_LDP)ETA_EXT(I)=1.0E-280_LDP
	  END DO
	  DO I=ND_ADD+1,ND_EXT
	    ETA_EXT(I)=ETA(I-ND_ADD)
	  END DO
	ELSE
	  CHI_EXT(1:ND)=CHI(1:ND)		!NB: In this can ND=ND_EXT
	  ETA_EXT(1:ND)=ETA(1:ND)
	END IF
!
	LOG_CHI_EXT(1:ND_EXT)=LOG(CHI_EXT(1:ND_EXT))
	CALL MON_INT_FUNS_V2(CHI_COEF,LOG_CHI_EXT,LOG_R_EXT,ND_EXT)
	LOG_ETA_EXT(1:ND_EXT)=LOG(ETA_EXT(1:ND_EXT))
	CALL MON_INT_FUNS_V2(ETA_COEF,LOG_ETA_EXT,LOG_R_EXT,ND_EXT)
	CALL TUNE(2,'FG_CHI_BEG')
!
! 
!************************************************************************
!************************************************************************
!
	CALL TUNE(1,'FG_MAIN_LOOP')
!
	PAR_JNU=0.0_LDP
	PAR_HNU=0.0_LDP
	PAR_KNU=0.0_LDP
	PAR_NNU=0.0_LDP
	PAR_OB_VEC=0.0_LDP
	PAR_IB_VEC=0.0_LDP
	IPLUS_P=0.0_LDP
!
! Enter loop to perform integration along each ray.
!
	DO IPROC=1,(NP-1)/NTHREAD+1
	  LS=GET_LS(MYPE,NTHREAD,IPROC)
	  IF(LS .GT. NP)EXIT
	  NI=RAY(LS)%NI_RAY
!
	  IF(RAY(LS)%NI_RAY .EQ. 1)THEN
	    RAY(LS)%I_P(1)=0.0_LDP
	    RAY(LS)%I_M(1)=0.0_LDP
	    NI=1
	    GOTO 1000
	  END IF
!
! NB: dCHIdR = LOG(CHI)/LOG(R) * CHI/R
!
	  DO I=1,RAY(LS)%NI_RAY
	    K=RAY(LS)%REXT_PNT(I)
	    IF(RAY(LS)%R_RAY(I) .EQ. R_EXT(K))THEN
	      CHI_RAY(I)=CHI_EXT(K)
	      ETA_RAY(I)=ETA_EXT(K)
	      dCHIdR_RAY(I)=CHI_COEF(K,3)*CHI_RAY(I)/RAY(LS)%R_RAY(I)
	    ELSE
	      T1=LOG(RAY(LS)%R_RAY(I)/R_EXT(K))
	      T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	      CHI_RAY(I)=EXP(T2)
	      T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	      ETA_RAY(I)=EXP(T2)
	      T2=(3.0_LDP*CHI_COEF(K,1)*T1+2.0_LDP*CHI_COEF(K,2))*T1+CHI_COEF(K,3)
	            dCHIdR_RAY(I)=T2*CHI_RAY(I)/RAY(LS)%R_RAY(I)
	    END IF
	  END DO
!
! By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum
! calculation for the first frequency.
!
	  IF(INIT)THEN
	    Q(1:NI)=0.0_LDP
	    SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
	  ELSE
	    Q(1:NI)=RAY(LS)%GAM(1:NI)/dLOG_NU
	    CHI_RAY(1:NI)=CHI_RAY(1:NI)+Q(1:NI)
	    dCHIdR_RAY(1:NI)=dCHIdR_RAY(1:NI)+RAY(LS)%dGAMdR(1:NI)/dLOG_NU
	    Q(1:NI)=Q(1:NI)/CHI_RAY(1:NI)
	    SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
	  END IF
!
!
!
	  IF(NEW_FREQ)THEN
!
! Compute the optical depth increments. This code is from TAU, and NORDTAU. We
! check that the Euler-Mauclarin correction is not too large. This is mainly
! done to prevent negative optical depths. The check is still necessary when
! we re using monotonic interpolation, since the monotonic interpolation only
! applies to CHI. CHI_RAY contains an additional term due to the frequency
! derivative. As non relativistic, DTAU is the same for both directions.
!
	    IF(METHOD .EQ. 'ZERO')THEN
	      DO I=1,NI-1
	        dZ=RAY(LS)%Z(I)-RAY(LS)%Z(I+1)
	        RAY(LS)%DTAU(I)=0.5_LDP*(CHI_RAY(I)+CHI_RAY(I+1))*dZ
	      END DO
	    ELSE
	      DO I=1,NI-1
	        dZ=RAY(LS)%Z(I)-RAY(LS)%Z(I+1)
	        RAY(LS)%DTAU(I)=0.5_LDP*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*RAY(LS)%Z(I+1)/RAY(LS)%R_RAY(I+1) -
	1            dCHIdR_RAY(I)*RAY(LS)%Z(I)/RAY(LS)%R_RAY(I) )/6.0_LDP )
	        IF( CHI_RAY(I) .LT. CHI_RAY(I+1) )THEN
	          RAY(LS)%DTAU(I)=MAX(CHI_RAY(I)*dZ,RAY(LS)%DTAU(I))
	          RAY(LS)%DTAU(I)=MIN(CHI_RAY(I+1)*dZ,RAY(LS)%DTAU(I))
	        ELSE
	          RAY(LS)%DTAU(I)=MIN(CHI_RAY(I)*dZ,RAY(LS)%DTAU(I))
	          RAY(LS)%DTAU(I)=MAX(CHI_RAY(I+1)*dZ,RAY(LS)%DTAU(I))
	        END IF
	      END DO
	    END IF
!
! Compute the functions used to evaluate the weighted integral over the
! polynomial fit to the source function.
!
! NB:  En(I)= EXP(-DTAU) {Integral[0 to DTAU] t^n EXP(t) dt }/ DTAU^n
!
	    DO I=1,NI-1
	      T1=RAY(LS)%DTAU(I)
	      EE(I)=0.0_LDP
	      IF(T1 .LT. 700.0_LDP)EE(I)=EXP(-T1)
	      IF(T1 .GE. 40.0_LDP)THEN
	      ELSE IF(T1 .GT. 0.5_LDP)THEN
	        E0(I)=1.0_LDP-EE(I)
	        E1(I)=1.0_LDP-E0(I)/T1
	        E2(I)=1.0_LDP-2.0_LDP*E1(I)/T1
	        E3(I)=1.0_LDP-3.0_LDP*E2(I)/T1
	      ELSE IF(T1 .GT. 0.1_LDP)THEN
	        E3(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1             (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	1             (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP*
	1             (1.0_LDP-T1/10.0_LDP*(1.0_LDP-T1/11.0_LDP*
	1             (1.0_LDP-T1/12.0_LDP*(1.0_LDP-T1/13.0_LDP)))))))) )
	        E2(I)=T1*( 1.0_LDP-E3(I) )/3.0_LDP
	        E1(I)=T1*( 1.0_LDP-E2(I) )/2.0_LDP
	        E0(I)=T1*( 1.0_LDP-E1(I) )
	      ELSE
	        E3(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1             (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	1             (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP) ))))
	        E2(I)=T1*( 1.0_LDP-E3(I) )/3.0_LDP
	        E1(I)=T1*( 1.0_LDP-E2(I) )/2.0_LDP
	        E0(I)=T1*( 1.0_LDP-E1(I) )
	      END IF
	     END DO
!
	     DO I=1,NI-1
	      T1=RAY(LS)%DTAU(I)
	      IF(T1 .GE. 40.0_LDP)THEN
	        RAY(LS)%A0(I)=EE(I)
	        RAY(LS)%A1(I)=(6.0_LDP-12.0_LDP/T1)/T1/T1
	        RAY(LS)%A2(I)=1.0_LDP-RAY(LS)%A1(I)
	        RAY(LS)%A3(I)=(2.0_LDP-6.0_LDP/T1)/T1
	        RAY(LS)%A4(I)=(4.0_LDP-6.0_LDP/T1)/T1-1.0_LDP
	      ELSE
	        RAY(LS)%A0(I)=EE(I)
	        RAY(LS)%A1(I)=E0(I)-3.0_LDP*E2(I)+2.0_LDP*E3(I)
	        RAY(LS)%A2(I)=3.0_LDP*E2(I)-2.0_LDP*E3(I)
	        RAY(LS)%A3(I)=RAY(LS)%DTAU(I)*(E1(I)-2.0_LDP*E2(I)+E3(I))
	        RAY(LS)%A4(I)=RAY(LS)%DTAU(I)*(E3(I)-E2(I))
	      END IF
	    END DO
	  END IF
!
!
! ******************* INWARD DIRECTED RAYS *********************************
!
! Compute the Source function for inward directed rays, and find the
! monotonic interpolating polynomial.
!
	  SOURCE_PRIME(1:NI)=SOURCE_RAY(1:NI)+Q(1:NI)*RAY(LS)%I_M_PREV(1:NI)
	  DO I=1,NI-1
	    S(I)=(SOURCE_PRIME(I+1)-SOURCE_PRIME(I))/RAY(LS)%DTAU(I)
	  END DO
!
! Now compute the derivatives node I.
!
	  dS(1)=S(1) +(S(1)-S(2))*RAY(LS)%DTAU(1)/(RAY(LS)%DTAU(1)+RAY(LS)%DTAU(2))
	  DO I=2,NI-1
	    dS(I)=(S(I-1)*RAY(LS)%DTAU(I)+S(I)*RAY(LS)%DTAU(I-1))/
	1                     (RAY(LS)%DTAU(I-1)+RAY(LS)%DTAU(I))
	  END DO
	  dS(NI)=S(NI-1)+(S(NI-1)-S(NI-2))*RAY(LS)%DTAU(NI-1)/
	1                       (RAY(LS)%DTAU(NI-2)+RAY(LS)%DTAU(NI-1))
!
! Adjust first derivatives so that function is monotonic in each interval.
!
	  dS(1)=( SIGN(ONE,S(1))+SIGN(ONE,dS(1)) )*
	1                      MIN(ABS(S(1)),0.5_LDP*ABS(dS(1)))
	  DO I=2,NI-1
	    dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1             MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(dS(I)))
	  END DO
	  dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*
	1               MIN(ABS(S(NI-1)),0.5_LDP*ABS(dS(NI)))
!
          RAY(LS)%I_M(1)=0.0_LDP
	  DO I=1,NI-1
	    RAY(LS)%I_M(I+1)=RAY(LS)%I_M(I)*RAY(LS)%A0(I)+ (
	1             SOURCE_PRIME(I)*RAY(LS)%A1(I)
	1        +    SOURCE_PRIME(I+1)*RAY(LS)%A2(I)
	1        +    dS(I)*RAY(LS)%A3(I) 
	1        +    dS(I+1)*RAY(LS)%A4(I) )
	  END DO
!
	  DO I=1,NI-1
	    IF(RAY(LS)%I_M(I+1) .LT. 0)THEN
	      WRITE(6,*)'Error for I_M, LS=',LS,ND,NI,FREQ
	      WRITE(6,'(2X,A,13(8X,A))')'I','IMP1','IMPR','DTAU','   Q',
	1                     '  SI','SIP1',' dSI','dSIP',
	1                     '  A0','  A1','  A2','  A3','  A4'
	      WRITE(6,'(I3,13ES12.4)')I,RAY(LS)%I_M(I+1),RAY(LS)%I_M_PREV(I+1),RAY(LS)%DTAU(I),Q(I),
	1                   SOURCE_PRIME(I),SOURCE_PRIME(I+1),dS(I),dS(I+1),
	1                   RAY(LS)%A0(I),RAY(LS)%A1(I),RAY(LS)%A2(I),RAY(LS)%A3(I),RAY(LS)%A4(I)
	    END IF
	  END DO
!
! ******************* OUTWARD DIRECTED RAYS *********************************
!
! Compute the Source function for outward directed rays, and find the
! monotonic interpolating polynomial.
!
	  SOURCE_PRIME(1:NI)=SOURCE_RAY(1:NI)+Q(1:NI)*RAY(LS)%I_P_PREV(1:NI)
	  DO I=1,NI-1
	    S(I)=(SOURCE_PRIME(I+1)-SOURCE_PRIME(I))/RAY(LS)%DTAU(I)
	  END DO
!
! Now compute the derivatives at node I.
!
	  dS(1)=S(1) +(S(1)-S(2))*RAY(LS)%DTAU(1)/(RAY(LS)%DTAU(1)+RAY(LS)%DTAU(2))
	  DO I=2,NI-1
	    dS(I)=(S(I-1)*RAY(LS)%DTAU(I)+S(I)*RAY(LS)%DTAU(I-1))/
	1                  (RAY(LS)%DTAU(I-1)+RAY(LS)%DTAU(I))
	  END DO
	  dS(NI)=S(NI-1)+(S(NI-1)-S(NI-2))*RAY(LS)%DTAU(NI-1)/
	1                  (RAY(LS)%DTAU(NI-2)+RAY(LS)%DTAU(NI-1))
!
! Adjust the first derivatives so that function is monotonic in each interval.
!
	  dS(1)=( SIGN(ONE,S(1))+SIGN(ONE,dS(1)) )*MIN(ABS(S(1)),0.5_LDP*ABS(dS(1)))
	  DO I=2,NI-1
	    dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(dS(I)))
	  END DO
	  dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*MIN(ABS(S(NI-1)),0.5_LDP*ABS(dS(NI)))
!
	  IF(INNER_BND_METH .EQ. 'DIFFUSION' .AND. LS .LE. NC)THEN
	    I_CORE=( ETA_RAY(NI)+ DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND) )/CHI_RAY(NI)
	  ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	    I_CORE=RAY(LS)%I_M(NI)
	  ELSE IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	    I_CORE=RAY(LS)%I_M(NI)
	    IF(LS .LE. NC_PNT_SRCE)I_CORE=IC+RAY(LS)%I_M(NI)
	  ELSE
	    I_CORE=IC
	  END IF
	  IF(LS .LE. NC)THEN
	    RAY(LS)%I_P(NI)=I_CORE
	  ELSE
	    RAY(LS)%I_P(NI)=RAY(LS)%I_M(NI)
	  END IF
	  DO I=NI-1,1,-1
	    RAY(LS)%I_P(I)=RAY(LS)%I_P(I+1)*RAY(LS)%A0(I)+ (
	1             SOURCE_PRIME(I+1)*RAY(LS)%A1(I)
	1        +    SOURCE_PRIME(I)*RAY(LS)%A2(I)
	1        -    dS(I+1)*RAY(LS)%A3(I)
	1        -    dS(I)*RAY(LS)%A4(I) )
	  END DO
!
	  DO I=1,NI
	    IF(RAY(LS)%I_P(I) .LT. 0)THEN
	      WRITE(6,*)'Error for RAY(LS)%I_P, LS=',LS,ND,NI
	      WRITE(6,'(2X,A,11(7X,A))')'I','  IM',' IMP','DTAU','   Q','   S','  dS',
	1                     '  A0','  A1','  A2','  A3','  A4'
	      WRITE(6,'(I3,11ES12.4)')I,RAY(LS)%I_P(I),RAY(LS)%I_P_PREV(I),RAY(LS)%DTAU(I),
	1                   Q(I),SOURCE_PRIME(I),dS(I),
	1                   RAY(LS)%A0(I),RAY(LS)%A1(I),RAY(LS)%A2(I),RAY(LS)%A3(I),RAY(LS)%A3(I)
	    END IF
	  END DO
C
C Note that V=AV(1)-IBOUND.
C
1000	  CONTINUE
	  K=RAY(LS)%J_PNT(1)
	  CV_BOUND(LS)=0.5_LDP*(RAY(LS)%I_P(K)-RAY(LS)%I_M(K))
	  I_M_IN_BND(LS)=RAY(LS)%I_M(NI)
!
! Compute the mean intensity like variable U at each grid point. Used to compute J.
!
	  DO I=1,RAY(LS)%NI_RAY
	    RAY(LS)%AV(I)=0.5_LDP*( RAY(LS)%I_P(I)+RAY(LS)%I_M(I) )
	  END DO
!
! We define CV on the gregular grid. We perform the interpolation onto
! the original R grid as we compute H.
!
	  DO I=1,RAY(LS)%NI_RAY
	    RAY(LS)%CV(I)= 0.5_LDP*( RAY(LS)%I_P(I)-RAY(LS)%I_M(I) )
	  END DO
!
	  RAY(LS)%I_P_STORE(:)=RAY(LS)%I_P(:)
	  RAY(LS)%I_M_STORE(:)=RAY(LS)%I_M(:)
!
! J and K are always evaluated on the nodes.
!
	  DO I=1,MIN(ND,ND-(LS-NC-1))
	    K=RAY(LS)%J_PNT(I)
	    PAR_JNU(I)=PAR_JNU(I)+RAY(LS)%JQW(I)*RAY(LS)%AV(K)
	    PAR_KNU(I)=PAR_KNU(I)+RAY(LS)%KQW(I)*RAY(LS)%AV(K)
	  END DO
!
	  IF(DEFINE_AT_MID_POINTS)THEN
            DO I=1,MIN(ND,ND-(LS-NC-1))-1
              K=RAY(LS)%H_PNT(I)
              T2=0.5_LDP*(R(I)+R(I+1))
              T1=(T2-RAY(LS)%R_RAY(K))/(RAY(LS)%R_RAY(K+1)-RAY(LS)%R_RAY(K))
              PAR_HNU(I)=PAR_HNU(I)+RAY(LS)%HMIDQW(I)*((1.0_LDP-T1)*RAY(LS)%CV(K)+T1*RAY(LS)%CV(K+1))
              PAR_NNU(I)=PAR_NNU(I)+RAY(LS)%NMIDQW(I)*((1.0_LDP-T1)*RAY(LS)%CV(K)+T1*RAY(LS)%CV(K+1))
            END DO
	    HN_DEF_ON_NODES=.FALSE.
	  ELSE
	    DO I=1,MIN(ND,ND-(LS-NC-1))
	      K=RAY(LS)%J_PNT(I)
	      PAR_HNU(I)=PAR_HNU(I)+RAY(LS)%HQW(I)*RAY(LS)%CV(K)
	      PAR_NNU(I)=PAR_NNU(I)+RAY(LS)%NQW(I)*RAY(LS)%CV(K)
	    END DO
	    HN_DEF_ON_NODES=.TRUE.
	  END IF
	  IPLUS_P(LS)=2.0_LDP*CV_BOUND(LS)
!
! This procedure works for the INTEGRAL or DIFFERENCE approach.
!
	  K=RAY(LS)%J_PNT(1)
	  T1=RAY(LS)%AV(K)+CV_BOUND(LS)
	  T2=0.0_LDP; T2=MAX(T2,RAY(LS)%AV(K)-CV_BOUND(LS))
	  PAR_OB_VEC(1:4)=PAR_OB_VEC(1:4)+T1*RAY(LS)%OB_AQW(1:4)
	  PAR_OB_VEC(5:8)=PAR_OB_VEC(5:8)+T2*RAY(LS)%OB_AQW(5:8)
!
	  IF(LS .LE. NC+1)THEN
	    T1=RAY(LS)%AV(K)-0.5_LDP*I_M_IN_BND(LS)
	    T2=I_M_IN_BND(LS)
	    PAR_IB_VEC(1:4)=PAR_IB_VEC(1:4)+T1*RAY(LS)%IB_AQW(1:4)
	    PAR_IB_VEC(5:8)=PAR_IB_VEC(5:8)+T2*RAY(LS)%IB_AQW(5:8)
	  END IF
!
	END DO  	!Processor loop
	CALL TUNE(2,'FG_MAIN_LOOP')
!
!
!***************************************************************************
!***************************************************************************
!
	CALL TUNE(1,'REDUCTIONS')
!
! Zero boundary conditions.
!
	HBC=0.0_LDP			!H/J at model outer boundary.
	NBC=0.0_LDP			!N/J at model outer boundary.
	IN_HBC=0.0_LDP
!
! Initialize intensity matrices.
!
	JNU_STORE(:)=0.0_LDP			!1:ND
	HNU_STORE(:)=0.0_LDP
	KNU_STORE(:)=0.0_LDP
	NNU_STORE(:)=0.0_LDP
!
! Sum the partial moments, sending the result to all processes. MPI_ALLREDUCE altomatically causes
! all process to sync at this locations.
!
	I=4*ND+16
	CALL MPI_ALLREDUCE(PAR_MOM,MOM_STORE,I,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,ierr)
	JNU_STORE(:)=MOM_STORE(1:ND)
	HNU_STORE(:)=MOM_STORE(ND+1:2*ND)
	KNU_STORE(:)=MOM_STORE(2*ND+1:3*ND)
	NNU_STORE(:)=MOM_STORE(3*ND+1:4*ND)
	IF(HN_DEF_ON_NODES)THEN
	  HNU_AT_OB=HNU_STORE(1); NNU_AT_OB=NNU_STORE(1)
	ELSE
	END IF
!
	CALL MPI_ALLREDUCE(MPI_IN_PLACE,IPLUS_P,NP,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,ierr)
!
	K=4*ND
	JPLUS_IB=MOM_STORE(K+1);  HPLUS_IB=MOM_STORE(K+2);  KPLUS_IB=MOM_STORE(K+3);    NPLUS_IB=MOM_STORE(K+4)
	JMIN_IB=MOM_STORE(K+5);   HMIN_IB=MOM_STORE(K+6);   KMIN_IB=MOM_STORE(K+7);     NMIN_IB=MOM_STORE(K+8)
	JPLUS_OB=MOM_STORE(K+9);  HPLUS_OB=MOM_STORE(K+10); KPLUS_OB=MOM_STORE(K+11);   NPLUS_OB=MOM_STORE(K+12)
	JMIN_OB=MOM_STORE(K+13);  HMIN_OB=MOM_STORE(K+14);  KMIN_OB=MOM_STORE(K+15);   NMIN_OB=MOM_STORE(K+16)
!
	IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
!	  HNU_AT_IB=DBB/R(ND)/CHI_RAY(NI)/3.0_LDP
!	  NNU_AT_IB=DBB/R(ND)/CHI_RAY(NI)/5.0_LDP
	  HNU_AT_IB=DBB/R(ND)/CHI(ND)/3.0_LDP
	  NNU_AT_IB=DBB/R(ND)/CHI(ND)/5.0_LDP
	ELSE IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	  T1=0.0_LDP; T2=0.0_LDP; HNU_AT_IB=0.0_LDP; NNU_AT_IB=0.0_LDP
	  DO IPROC=1,(NC_PNT_SRCE-1)/NTHREAD+1
	    LS=GET_LS(MYPE,NTHREAD,IPROC)
	    IF(LS .GT. NC_PNT_SRCE)EXIT
	    T1=T1+0.5_LDP*IC*RAY(LS)%HQW(ND)
	    T2=T2+0.5_LDP*IC*RAY(LS)%NQW(ND)
	  END DO
	  CALL MPI_ALLREDUCE(T1,HNU_AT_IB,IONE,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,ierr)
	  CALL MPI_ALLREDUCE(T2,NNU_AT_IB,IONE,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,ierr)
	ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	  HNU_AT_IB=0.0_LDP
	  NNU_AT_IB=0.0_LDP
	ELSE
	  T1=0.0_LDP; T2=0.0_LDP; HNU_AT_IB=0.0_LDP; NNU_AT_IB=0.0_LDP
	  DO IPROC=1,NC/NTHREAD+1
	    LS=GET_LS(MYPE,NTHREAD,IPROC)
	    IF(LS .GT. NC+1)EXIT
	    T1=T1+0.5_LDP*(IC-I_M_IN_BND(LS))*RAY(LS)%HQW(ND)
	    T2=T2+0.5_LDP*(IC-I_M_IN_BND(LS))*RAY(LS)%NQW(ND)
	  END DO
	  CALL MPI_ALLREDUCE(T1,HNU_AT_IB,IONE,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,ierr)
	  CALL MPI_ALLREDUCE(T2,NNU_AT_IB,IONE,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,ierr)
	END IF
!
! Get boundary conditons for moment calculations.
!
	T1=0.0_LDP; IN_HBC=0.0_LDP
	DO IPROC=1,NC/NTHREAD+1
	  LS=GET_LS(MYPE,NTHREAD,IPROC)
	  IF(LS .GT. NC+1)EXIT
	  T1=T1 + RAY(LS)%HQW(ND)*I_M_IN_BND(LS)
	END DO
	CALL MPI_ALLREDUCE(T1,IN_HBC,IONE,MY_MPI_DP,MPI_SUM,MPI_COMM_WORLD,ierr)
!
	CALL TUNE(2,'REDUCTIONS')
!
!
!
! Compute the boundary Eddington factors.
!
	HBC=HNU_AT_OB/JNU_STORE(1)
	NBC=NNU_AT_OB/JNU_STORE(1)
	IN_HBC=IN_HBC/(2.0_LDP*JNU_STORE(ND)-IC)
	RETURNED_OUT_HBC=HBC
	RETURNED_IN_HBC=IN_HBC
!
! Compute the Eddington factor, F. This required in CMFGEN to
! compute the K moment (from J computed by MOM_J_CMF).
!
	DO I=1,ND
	  JNU(I)=JNU_STORE(I)
	  FEDD(I)=KNU_STORE(I)/JNU_STORE(I)
	END DO
!
! Store frequencies at which errors occurred, and give an indication of the
! error.
!
! The was a leftover from difference approach, and should not be nbeeded here.

	DO J=1,N_ERR_MAX
	  IF(.NOT. NEG_AV_VALUE)EXIT
	  IF(FG_ERR_ON_FREQ(J) .EQ. FREQ)THEN
	    FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+1
	    NEG_AV_VALUE=.FALSE.
	  ELSE IF(J .EQ. FG_ERR_CNT+1)THEN
	    FG_ERR_CNT=J
	    FG_ERR_TYPE(J)=1
	    FG_ERR_ON_FREQ(J)=FREQ
	    NEG_AV_VALUE=.FALSE.
	  END IF
	END DO
!
	IF(MYPE .EQ. 0)THEN
	  DO I=1,ND
	    IF(JNU_STORE(I) .LT. 0)THEN
	      OPEN(UNIT=7,FILE='FG_J_CMF_MPI_V1_ERRORS',STATUS='UNKNOWN')
	        WRITE(7,*)'FREQ=',FREQ
	        WRITE(7,'(3X,A,5X,4(5X,A,5X))')'I','JNU_STORE','ETA','CHI','ESEC'
	        DO J=1,ND
	          WRITE(7,'(X,I5,4ES16.6)')J,JNU_STORE(J),ETA(J),CHI(J),ESEC(J)
	        END DO
	      CLOSE(UNIT=7)
	      J=ERROR_LU()
	      WRITE(J,*)'Error on FG_J_CMF_MPI_V1 --- negative mean intensities.'
	      WRITE(J,*)'Check out file FG_J_CMF_MPI_V1_ERRORS for aditional information.'
	      WRITE(J,*)'Halting code execution.'
	      CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
	      STOP
	    END IF
	  END DO
	END IF
!
! This will only output IP_FG_DATA if a file IP_FG_DATA exists.
! This is only for debugging purposes.
!
	IF(MYPE .EQ. 0)THEN
          ACCESS_F=5
	  LU_IP=56
          IF(INIT)THEN
	    INQUIRE(FILE='IP_FG_DATA',EXIST=WRITE_IP)
	    IF(WRITE_IP)THEN
	      CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
              I=WORD_SIZE*(NP+1)/UNIT_SIZE
              CALL WRITE_DIRECT_INFO_V3(NP,I,'20-Aug-2000','IP_FG_DATA',LU_IP)
              OPEN(UNIT=LU_IP,FILE='IP_FG_DATA',FORM='UNFORMATTED',
	1         ACCESS='DIRECT',STATUS='UNKNOWN',RECL=I,IOSTAT=IOS)
	      FREQ_CNT=0
              WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
              WRITE(LU_IP,REC=ACCESS_F)(P(I),I=1,NP)
	    END IF
	  END IF
          IF(WRITE_IP)THEN
	    T1=0.0_LDP
            IF(FREQ_CNT .NE. 0)READ(LU_IP,REC=ACCESS_F+FREQ_CNT)(IBOUND(LS),LS=1,NP),T1
	    IF(T1 .NE. FREQ)FREQ_CNT=FREQ_CNT+1
            WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
            WRITE(LU_IP,REC=ACCESS_F+FREQ_CNT)(CV_BOUND(LS),LS=1,NP),FREQ
          END IF
	END IF
	CALL TUNE(2,'FG_FULL')
!
	RETURN
	END SUBROUTINE FG_J_CMF_MPI_V1
!
! Declares arrays required for the angular quadrature weights, and
! computes their vaule. Routine can be called more than once, although
! the size of the declared arrays must not change.
!
	SUBROUTINE SET_ANG_QW_MPI(R,P,NC,ND,NP,TRAPFORJ)
	USE SET_KIND_MODULE
	USE FG_J_CMF_MOD_MPI_V1
	IMPLICIT NONE
!
! Created 11-Aug-2024
!
	INTEGER NC,ND,NP
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) P(NP)
	LOGICAL TRAPFORJ
!
! External function calls.
!
	EXTERNAL JWEIGHT_V2,HWEIGHT_V2,KWEIGHT_V2,NWEIGHT_V2
	EXTERNAL JTRPWGT_V2,HTRPWGT_V2,KTRPWGT_V2,NTRPWGT_V2
	EXTERNAL ERROR_LU
	INTEGER ERROR_LU
!
	REAL(KIND=LDP), ALLOCATABLE :: AQW(:,:)
!
	INTEGER LU_ER
	INTEGER LS
	INTEGER IOS
	INTEGER IPROC
	LOGICAL MID
	LOGICAL AQW_ALLOCATED
!
	INTEGER GET_LS,I
	GET_LS(MYPE,NTHREAD,I)=MOD(I,2)*(MYPE+(I-1)*NTHREAD+1)+MOD(I+1,2)*(I*NTHREAD-MYPE)
!
	AQW_ALLOCATED=.FALSE.
!
! Allocate wory arrays.
!
	ALLOCATE(AQW(ND,NP),STAT=IOS)
!
	DO IPROC=1,NUM_RAYS_PER_THREAD
	  IF(IPROC .GT. 1)THEN
	  ELSE IF(TRAPFORJ)THEN
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,JTRPWGT_V2)
	  ELSE
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,JWEIGHT_V2)
	  END IF
	  LS=GET_LS(MYPE,NTHREAD,IPROC)
	  IF(LS .GT. NP)EXIT
	  IF( ALLOCATED(RAY(LS)%JQW) )AQW_ALLOCATED=.TRUE.
	  IF(.NOT. AQW_ALLOCATED)ALLOCATE(RAY(LS)%JQW(ND))
	  RAY(LS)%JQW(1:ND)=AQW(1:ND,LS)
	END DO
!
	DO IPROC=1,NUM_RAYS_PER_THREAD
	  IF(IPROC .GT. 1)THEN
	  ELSE IF(TRAPFORJ)THEN
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,HTRPWGT_V2)
	  ELSE
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,HWEIGHT_V2)
	  END IF
	  LS=GET_LS(MYPE,NTHREAD,IPROC)
	  IF(LS .GT. NP)EXIT
	  IF(.NOT. AQW_ALLOCATED)ALLOCATE(RAY(LS)%HQW(ND))
	  RAY(LS)%HQW(1:ND)=AQW(1:ND,LS)
	END DO
!
	DO IPROC=1,NUM_RAYS_PER_THREAD
	  IF(IPROC .GT. 1)THEN
	  ELSE IF(TRAPFORJ)THEN
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,KTRPWGT_V2)
	  ELSE
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,KWEIGHT_V2)
	  END IF
	  LS=GET_LS(MYPE,NTHREAD,IPROC)
	  IF(LS .GT. NP)EXIT
	  IF(.NOT. AQW_ALLOCATED)ALLOCATE(RAY(LS)%KQW(ND))
	  RAY(LS)%KQW(1:ND)=AQW(1:ND,LS)
	END DO
!
	DO IPROC=1,NUM_RAYS_PER_THREAD
	  IF(IPROC .GT. 1)THEN
	  ELSE IF(TRAPFORJ)THEN
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,NTRPWGT_V2)
	  ELSE
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,NWEIGHT_V2)
	  END IF
	  LS=GET_LS(MYPE,NTHREAD,IPROC)
	  IF(LS .GT. NP)EXIT
	  IF(.NOT. AQW_ALLOCATED)ALLOCATE(RAY(LS)%NQW(ND))
	  RAY(LS)%NQW(1:ND)=AQW(1:ND,LS)
	END DO
!
!
	DO IPROC=1,NUM_RAYS_PER_THREAD
	  MID=.TRUE.
	  IF(IPROC .GT. 1)THEN
	  ELSE IF(TRAPFORJ)THEN
	    CALL GENANGQW_V2(AQW,R,P,NC,ND,NP,HTRPWGT_V2,MID)
	  ELSE
	    CALL GENANGQW_V2(AQW,R,P,NC,ND,NP,HWEIGHT_V2,MID)
	  END IF
	  LS=GET_LS(MYPE,NTHREAD,IPROC)
	  IF(LS .GT. NP)EXIT
	  IF(.NOT. AQW_ALLOCATED)ALLOCATE(RAY(LS)%HMIDQW(ND))
	  RAY(LS)%HMIDQW(1:ND)=AQW(1:ND,LS)
	END DO
!
	DO IPROC=1,NUM_RAYS_PER_THREAD
	  MID=.TRUE.
	  IF(IPROC .GT. 1)THEN
	  ELSE IF(TRAPFORJ)THEN
	    CALL GENANGQW_V2(AQW,R,P,NC,ND,NP,HTRPWGT_V2,MID)
	  ELSE
	    CALL GENANGQW_V2(AQW,R,P,NC,ND,NP,HWEIGHT_V2,MID)
	  END IF
	  LS=GET_LS(MYPE,NTHREAD,IPROC)
	  IF(LS .GT. NP)EXIT
	  IF(.NOT. AQW_ALLOCATED)ALLOCATE(RAY(LS)%NMIDQW(ND))
	  RAY(LS)%NMIDQW(1:ND)=AQW(1:ND,LS)
	END DO
!
	DEALLOCATE(AQW)
!
	RETURN
	END SUBROUTINE SET_ANG_QW_MPI
