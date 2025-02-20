	MODULE EXT_REL_GRID_MPI_V1
	  USE SET_KIND_MODULE
	  REAL(KIND=LDP), ALLOCATABLE :: R_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: LOG_R_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: Z_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: V_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: VDOP_VEC_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: SIGMA_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: ETA_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: CHI_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: LOG_ETA_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: LOG_CHI_EXT(:)
	  REAL(KIND=LDP), ALLOCATABLE :: TMP_VEC(:)
!
	  REAL(KIND=LDP), ALLOCATABLE :: CHI_COEF(:,:)
	  REAL(KIND=LDP), ALLOCATABLE :: ETA_COEF(:,:)
!
	  REAL(KIND=LDP), ALLOCATABLE :: CHI_RAY(:)
	  REAL(KIND=LDP), ALLOCATABLE :: ETA_RAY(:)
!
! PAR_MOM will be used to store all summed moments for each process.
! We use the pointer to point to a specifc location in PAR_MOM. See
! below.
!
	REAL(KIND=LDP), TARGET, ALLOCATABLE :: PAR_MOM(:)
	REAL(KIND=LDP), POINTER :: PAR_JNU(:)
	REAL(KIND=LDP), POINTER :: PAR_HNU(:)
	REAL(KIND=LDP), POINTER :: PAR_KNU(:)
	REAL(KIND=LDP), POINTER :: PAR_NNU(:)
!
	REAL(KIND=LDP), POINTER :: PAR_IB_VEC(:)
	REAL(KIND=LDP), POINTER :: PAR_OB_VEC(:)
!
	  INTEGER ND_EXT
	  INTEGER ND_ADD
	  INTEGER NP_MAX
	  INTEGER IDMIN,IDMAX
	  INTEGER, PARAMETER :: ND_ADD_MAX=24
	  LOGICAL, SAVE ::  FIRST_TIME=.TRUE.
!
	END MODULE EXT_REL_GRID_MPI_V1
!
	SUBROUTINE CMF_FORMAL_REL_MPI_V1
	1           (ETA,CHI,ESEC,V,SIGMA,R,P,
	1            JNU,FEDD,RET_HNU_AT_IB,RET_HNU_AT_OB,IPLUS,
	1            FREQ,dLOG_NU,B_PLANCK,DBB,
	1            INNER_BND_METH,THICK_OB,
	1            VDOP_VEC,VDOP_FRAC,REXT_FAC,
	1            DUST_SCAT_OPAC,G_HEN_GREEN,USE_HEN_GREEN,
	1            METHOD,INITIALIZE,NEW_FREQ,
	1            NC,NP,DST,DEND,ND)
	USE SET_KIND_MODULE
        USE EXT_REL_GRID_MPI_V1
	USE MOD_SPACE_GRID_MPI_V1
	USE MOD_RAY_MOM_STORE
	USE MPI
	IMPLICIT NONE
!
! Created 19-Jan-2025 -- Under development.
!
! Altered 26-Mar-2020: Increased N_STORE by a factor of 5.
! Altered 19-Nov-2016: CHI_COEF, ETA_COEF, and VDOP_VEC_EXT were not being deallocated.
!                        Only reallocate vectors if ND_EXT changes.
! Altered 04-Apr-2013: Fixed bug. NP_LIMIT(not NP) is limit when THICK_OB=FALSE.
! Altered 08-Jan-2012: Changed to V4 -- added REXT_FAC to call.
! Altered 14-May-2009: Altered value of IDMAX (for ETA and CHI interpolaton).
!                        IDMAX & IDMIN now stored in FG_J_CMF_MOD_V11.
        INTEGER NC
        INTEGER DST,DEND,ND
        INTEGER NP
!
	REAL(KIND=LDP) R(ND)
        REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) P(NP)
	REAL(KIND=LDP) ETA(ND)	
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) VDOP_VEC(ND)
	REAL(KIND=LDP) VDOP_FRAC
	REAL(KIND=LDP) DELV_FRAC_FG
	REAL(KIND=LDP) REXT_FAC
!
	REAL(KIND=LDP) RET_HNU_AT_IB
	REAL(KIND=LDP) RET_HNU_AT_OB
!
! NB: J,H,K,N refer to the first 4 moments of the radiation field.
!
        REAL(KIND=LDP) JNU(ND)
	REAL(KIND=LDP) FEDD(ND)
        REAL(KIND=LDP) IPLUS(NP)
!
	REAL(KIND=LDP) B_PLANCK
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) dLOG_NU
!
	CHARACTER*6 METHOD
!
	CHARACTER(LEN=*) INNER_BND_METH
!
! Use "Thick" boundary condition. at outer boundary. Only noted when INITIALIZE
! is true. All subsequent frequencies will use the same boundary condition
! independent of the passed value (Until INITIALIZE is set to TRUE again).
!
	LOGICAL THICK_OB
!
	REAL(KIND=LDP) DUST_SCAT_OPAC(ND)
	REAL(KIND=LDP) G_HEN_GREEN
	LOGICAL USE_HEN_GREEN
!
! First frequency -- no frequency coupling.
!
	LOGICAL INITIALIZE
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
! Local variables.
!
	INTEGER, PARAMETER :: IONE=1
	LOGICAL, PARAMETER :: LFALSE=.FALSE.
	LOGICAL, PARAMETER :: LTRUE=.TRUE.
!
	REAL(KIND=LDP) ETA_SCAT(ND)
	REAL(KIND=LDP) DBC
	REAL(KIND=LDP) I_CORE
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) dBdTAU
	REAL(KIND=LDP) ALPHA
	REAL(KIND=LDP) ESEC_POW
	REAL(KIND=LDP) BETA
	REAL(KIND=LDP) VINF
	REAL(KIND=LDP) RMAX
	REAL(KIND=LDP) DEL_R_FAC
	REAL(KIND=LDP) NU_ON_dNU
	REAL(KIND=LDP) MU_AT_RMAX
	REAL(KIND=LDP) HQW_AT_RMAX
	CHARACTER(LEN=20) BOUNDARY
!
	INTEGER NDM1
	INTEGER I,J,IERR
	INTEGER K,IP,ID,IOS
	INTEGER NP_LIMIT
	INTEGER NRAY
!
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) SPEED_OF_LIGHT
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU, SPEED_OF_LIGHT
	LOGICAL NEW_R_GRID
	LOGICAL REALLOCATE_GRID
!
	INTEGER NUM_RAYS_PER_THREAD
	INTEGER GET_IP
	INTEGER IPROC
	GET_IP(MYPE,NTHREAD,I)=MOD(I,2)*(MYPE+(I-1)*NTHREAD+1)+MOD(I+1,2)*(I*NTHREAD-MYPE)
!
!
!
	IF(INITIALIZE)THEN
	  NU_ON_dNU=0.0_LDP
	ELSE
	  NU_ON_dNU=1.0_LDP/dLOG_NU
	END IF
	NUM_RAYS_PER_THREAD=(NP+NTHREAD-1)/NTHREAD
!
! Allocate data for moments which will be used to construct the Eddington
! factors.
!
	IF(.NOT. ALLOCATED(JNU_STORE))THEN
	  ND_STORE=ND
	  ALLOCATE (JNU_STORE(ND))
	  ALLOCATE (HNU_STORE(ND))
	  ALLOCATE (KNU_STORE(ND))
	  ALLOCATE (NNU_STORE(ND))
	  ALLOCATE (R_STORE(ND))
	  ALLOCATE (GAM_REL_STORE(ND))
	  ALLOCATE (RMID_STORE(ND))
	  ALLOCATE (EXT_RMID_STORE(ND+1))
	END IF
!
	IF(INITIALIZE)THEN
	  R_STORE(1:ND)=R(1:ND)
	  GAM_REL_STORE(1:ND)=1.0_LDP/SQRT(1.0_LDP-(V(1:ND)/2.99792458E+05_LDP)**2)
	  DO I=1,ND-1
	    RMID_STORE(I)=0.5_LDP*(R(I)+R(I+1))
	    EXT_RMID_STORE(I+1)=0.5_LDP*(R(I)+R(I+1))
	  END DO
	  EXT_RMID_STORE(1)=R(1); EXT_RMID_STORE(ND+1)=R(ND)
	END IF
!
! Check to see whether we have a new R grid, or the solution options have
! change. This can only happen when INIT is TRUE.
!
	LUER=ERROR_LU()
	NEW_R_GRID=.FALSE.
	REALLOCATE_GRID=.FALSE.
	IF(INITIALIZE .AND. .NOT. FIRST_TIME)THEN
	  IF(ND_EXT .NE. ND+ND_ADD)THEN
	    NEW_R_GRID=.TRUE.
	    REALLOCATE_GRID=.TRUE.
	  ELSE
	    DO I=1,ND
	      IF(R(I) .NE. R_EXT(ND_ADD+I))THEN
	        NEW_R_GRID=.TRUE.
	        WRITE(LUER,*)'Updating RGRID in CMF_FORMAL_REL_V4'
	        EXIT
	      END IF
	    END DO
	  END IF
	ELSE IF(INITIALIZE)THEN
	  NEW_R_GRID=.TRUE.
	END IF
!
! Deallocate all allocated rays if we are using a diferent solution technique.
! This option will only be used when testing, since in CMFGEN we will always use
! the same atmospheric structure.
!
	IF(ALLOCATED(R_EXT) .AND. REALLOCATE_GRID)THEN
	  DEALLOCATE ( R_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( LOG_R_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( V_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( Z_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( SIGMA_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( ETA_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( CHI_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( LOG_ETA_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( LOG_CHI_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( TMP_VEC, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( VDOP_VEC_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( CHI_COEF, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( ETA_COEF, STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error deallocating R_EXT etc in CMF_FORMA_REL_V4'
	    STOP
	  END IF
	END IF
!
! Set up the revised grid to improve computational accuracy. Unles we are
! carrying out tests, these need only be constructed once.
!
	IF(FIRST_TIME .OR. .NOT. ALLOCATED(R_EXT) )THEN
!
          ND_ADD=0
          IF(THICK_OB)ND_ADD=ND_ADD_MAX
          ND_EXT=ND+ND_ADD
!
	  ALLOCATE ( R_EXT(ND_EXT),STAT=IOS );          R_EXT=0.0_LDP
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_R_EXT(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( V_EXT(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( VDOP_VEC_EXT(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( Z_EXT(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_EXT(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_EXT(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_EXT(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_ETA_EXT(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_CHI_EXT(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( TMP_VEC(ND_EXT),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error allocating memore in CMF_FORMAL_REL_V4'
	    STOP
	  END IF
	END IF
!
	IF(NEW_R_GRID)THEN
!
! Compute the extended R grid, excluding inserted points.
!
	  DO I=1,ND
	    R_EXT(ND_ADD+I)=R(I)
	  END DO
	  IF(THICK_OB)THEN
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
	    DEL_R_FAC=EXP( LOG(RMAX/ALPHA)/(ND_ADD-4) )
	    R_EXT(1)=RMAX
	    R_EXT(5)=RMAX/DEL_R_FAC
	    R_EXT(2)=R_EXT(1)-0.001_LDP*(R_EXT(1)-R_EXT(5))
	    R_EXT(3)=R_EXT(1)-0.1_LDP*(R_EXT(1)-R_EXT(5))
	    R_EXT(4)=R_EXT(1)-0.4_LDP*(R_EXT(1)-R_EXT(5))
	    DO I=5,ND_ADD-1
	      R_EXT(I)=R_EXT(I-1)/DEL_R_FAC
	    END DO
	    R_EXT(ND_ADD)=ALPHA
!
	  END IF
	  C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
!
! Set up stodarge and pointers for computation of the moments, and
! partial moments
! at the outer boudaries.
!
	  IF(FIRST_TIME)THEN
            IOS=0
	    ALLOCATE ( PAR_MOM(4*ND+16),STAT=IOS )
	    ALLOCATE ( PAR_IB_VEC(8),STAT=IOS )
	    ALLOCATE ( PAR_OB_VEC(8),STAT=IOS )
            PAR_JNU=>PAR_MOM(1:ND)
            PAR_HNU=>PAR_MOM(ND+1:2*ND)
            PAR_KNU=>PAR_MOM(2*ND+1:3*ND)
            PAR_NNU=>PAR_MOM(3*ND+1:4*ND)
            PAR_IB_VEC=>PAR_MOM(4*ND+1:4*ND+8)
            PAR_OB_VEC=>PAR_MOM(4*ND+9:4*ND+16)
	    IF(MYPE .EQ. 0)THEN
	      WRITE(6,*)'PAR_MOM allocated in CMF_REL...'
	      FLUSH(UNIT=6)
	    END IF
	  END IF
!
!
! Compute VEXT and R_EXT. We assume a BETA velocity law at large R.
!
	  V_EXT(ND_ADD+1:ND_EXT)=V(1:ND)
	  VDOP_VEC_EXT(ND_ADD+1:ND_EXT)=VDOP_VEC(1:ND)
	  SIGMA_EXT(ND_ADD+1:ND_EXT)=SIGMA(1:ND)
	  IF(THICK_OB)THEN
	    BETA=(SIGMA(1)+1.0_LDP)*(R(1)/R(ND)-1.0_LDP)
            VINF=V(1)/(1-R(ND)/R(1))**BETA
	    DO I=1,ND_ADD
	      V_EXT(I)=VINF*(1.0_LDP-R_EXT(ND_EXT)/R_EXT(I))**BETA
	      SIGMA_EXT(I)=BETA/(R_EXT(I)/R_EXT(ND_EXT)-1.0_LDP)-1.0_LDP
	    END DO
	    VDOP_VEC_EXT(1:ND_ADD)=VDOP_VEC(1)
	    WRITE(LUER,*)'   Using thick boundary condition in CMF_FORM_REL'
	    WRITE(LUER,'(7X,2(A,ES16.8,3X))')' R(1)=',R(1),'RMAX=',RMAX
	    WRITE(LUER,'(7X,2(A,ES16.8,3X))')' V(1)=',V(1),'VMAX=',V_EXT(1)
	  END IF
	  LOG_R_EXT(1:ND_EXT)=LOG(R_EXT(1:ND_EXT))
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
	      IDMAX=ND/6
	    END IF
	    WRITE(6,'(4X,A,I4,A)')'Using depth 1 and',IDMAX,' to extrapolate opacities'
	  END IF
!
	END IF
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
	IF(ND_ADD .NE. 0)THEN
	  IF(CHI(IDMIN) .LE. ESEC(IDMIN) .OR. CHI(IDMAX) .LE. ESEC(IDMAX))THEN
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.0_LDP)ESEC_POW=2.0_LDP
	    DO I=1,ND_ADD
	      CHI_EXT(I)=CHI(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
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
! This must be after the call to DEFINE_GRID so that RAY_POINTS_INSERTED is defined.
!
	IF(FIRST_TIME .OR. NEW_R_GRID)THEN
	  CALL DEFINE_GRID_MPI_V1(R_EXT,V_EXT,VDOP_VEC_EXT,VDOP_FRAC,ND_EXT,R,P,ND,NC,NP)
	  J=0
	  OPEN(UNIT=7,FILE='MU_VALUE_CHK',STATUS='UNKNOWN')
	  IF(MYPE .EQ. 0)THEN
	    WRITE(7,'(A)')' '
	    WRITE(7,'(A)')' Comparison of MU(cmf) and MU(obs) at outer boundary (CMF_FORMAL_REL_V4)'
	    WRITE(7,'(A)')' The  first MU(obs) is the transformed value of MU(cmf)'
	    WRITE(7,'(A)')' The second MU(obs) is simply computed from P(ip) and RMAX'
	    WRITE(7,'(A)')' '
	    WRITE(7,'(2X,A,4(6X,A,2X))')'IP',' MU(cmf)',' HQW(cmf)',' MU(obs)',' MU(obs)'
	    FLUSH(UNIT=6)
	  END IF
	  T1=0.0_LDP
	  DO IPROC=1,NUM_RAYS_PER_THREAD
	    IP=GET_IP(MYPE,NTHREAD,IPROC)
	    IF(IP .GT. NP)EXIT
	    J=MAX(J,RAY(IP)%NZ)
	    MU_AT_RMAX=RAY(IP)%MU_P(RAY(IP)%LNK(1))
	    HQW_AT_RMAX=2.0_LDP*ray(ip)%HQW_P(1)
!	    WRITE(7,'(I4,4ES16.6)')IP,MU_AT_RMAX,HQW_AT_RMAX,
!	1                 (MU_AT_RMAX+V(1)/C_KMS)/(1.0D0+MU_AT_RMAX*V(1)/C_KMS),
!	1                 SQRT( (R(1)-P(IP))*(R(1)+P(IP)) )/R(1)
	  END DO
	  IF(MYPE .EQ. 0)CLOSE(UNIT=7)
	  IF( ALLOCATED(CHI_RAY) ) DEALLOCATE (CHI_RAY)
	  IF( ALLOCATED(ETA_RAY) ) DEALLOCATE (ETA_RAY)
	  CALL MPI_ALLREDUCE(MPI_IN_PLACE,J,IONE,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
	  ALLOCATE (CHI_RAY(J),ETA_RAY(J))
!
	  IF(USE_HEN_GREEN)THEN
	    WRITE(6,*)'Calling CMF_REL_DUST_QW',J,NP;FLUSH(UNIT=6)
	    CALL CMF_REL_DUST_QW(G_HEN_GREEN,ND,NP)
	  END IF
!
          DO IPROC=1,NUM_RAYS_PER_THREAD
	    IP=GET_IP(MYPE,NTHREAD,IPROC)
	    IF(IP .GT. NP)EXIT
	    IF(ALLOCATED(RAY(IP)%ETA_M))DEALLOCATE(RAY(IP)%ETA_M,RAY(IP)%ETA_P)
            NRAY=RAY(IP)%NZ
            ALLOCATE (RAY(IP)%ETA_M(NRAY))
            ALLOCATE (RAY(IP)%ETA_P(NRAY))
	  END DO
	  IF(MYPE .EQ. 0)WRITE(6,*)'Allocated RAY(IP)%ETA_M'; FLUSH(UNIT=6)
	END IF
!
! For SN we can have a hollow core. In this case we need to store the
! inward radiation field for use in calculating the outward radiation field.
! Because of the expansion, the comoving frequencies do not match, hence
! the need for storage of results at earlier frequencies.
!
	IF(INNER_BND_METH .EQ. 'HOLLOW')THEN
	  IF(.NOT. ALLOCATED(FREQ_STORE))THEN
	    N_STORE=2.0_LDP*V(ND)                   !/5.0D0
	    WRITE(6,*)'N_STORE=',N_STORE
	    WRITE(6,*)'VDOP_FRAC=',VDOP_FRAC
	    WRITE(6,*)'MIN(VDOP)=',MINVAL(VDOP_VEC)
	    ALLOCATE (FREQ_STORE(0:N_STORE-1))
	    DO IP=1,NC
	      ALLOCATE (RAY(IP)%I_IN_BND_STORE(0:N_STORE-1))
	    END DO
	    IF(NEW_R_GRID)THEN
	      BETA=V(ND)/2.99792458E+05_LDP
	      DO IP=1,NC
	        T1=SQRT( (R(ND)-P(IP))*(R(ND)+P(IP)) )/R(ND)
	        T2=1.0_LDP-BETA*(T1+BETA)/(1.0_LDP+BETA*T1)
	        RAY(IP)%FREQ_CONV_FAC=1.0_LDP/GAM_REL_STORE(ND)/GAM_REL_STORE(ND)/T2/(1.0_LDP-BETA*T1)
	        WRITE(6,'(5ES14.4)')RAY(IP)%FREQ_CONV_FAC,GAM_REL_STORE(ND),T1,T2,BETA
	      END DO
	    END IF
	  END IF
!
	  IF(INITIALIZE)THEN
	    CUR_LOC=-1
	    FREQ_STORE=0.0_LDP
	    DO IP=1,NC
	      RAY(IP)%I_IN_BND_STORE=0.0_LDP
	    END DO
	  END IF
	END IF
!
	IF(RAY_POINTS_INSERTED)THEN
	  LOG_CHI_EXT(1:ND_EXT)=LOG(CHI_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(CHI_COEF,LOG_CHI_EXT,LOG_R_EXT,ND_EXT)
	  LOG_ETA_EXT(1:ND_EXT)=LOG(ETA_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(ETA_COEF,LOG_ETA_EXT,LOG_R_EXT,ND_EXT)
	END IF
!
	IF(INITIALIZE .AND. NEW_FREQ)THEN
          DO IPROC=1,NUM_RAYS_PER_THREAD
	    IP=GET_IP(MYPE,NTHREAD,IPROC)
	    IF(IP .GT. NP)EXIT
	    RAY(IP)%I_P=0.0_LDP; RAY(IP)%I_M=0.0_LDP
	    RAY(IP)%I_P_PREV=0.0_LDP; RAY(IP)%I_M_PREV=0.0_LDP
	    RAY(IP)%I_P_SAVE=0.0_LDP; RAY(IP)%I_M_SAVE=0.0_LDP
	  END DO	
	  HNU_AT_OB_PREV=0.0_LDP; NNU_AT_OB_PREV=0.0_LDP
	  HNU_AT_IB_PREV=0.0_LDP; NNU_AT_IB_PREV=0.0_LDP
	ELSE IF(NEW_FREQ)THEN
          DO IPROC=1,NUM_RAYS_PER_THREAD
	    IP=GET_IP(MYPE,NTHREAD,IPROC)
	    IF(IP .GT. NP)EXIT
	    RAY(IP)%I_P=0.0_LDP; RAY(IP)%I_M=0.0_LDP
	    RAY(IP)%I_P_PREV=RAY(IP)%I_P_SAVE
	    RAY(IP)%I_M_PREV=RAY(IP)%I_M_SAVE
	    RAY(IP)%I_P_SAVE=0.0_LDP;   RAY(IP)%I_M_SAVE=0.0_LDP
	  END DO	
	  HNU_AT_OB_PREV=HNU_AT_OB; NNU_AT_OB_PREV=NNU_AT_OB
	  HNU_AT_IB_PREV=HNU_AT_IB; NNU_AT_IB_PREV=NNU_AT_IB
	END IF
!
! If no points have been inserted on the rays, CHI and ETA are the same for
! all rays.
!
	IF(.NOT. RAY_POINTS_INSERTED)THEN
	  CHI_RAY(1:ND_EXT)=CHI_EXT(1:ND_EXT)
	  ETA_RAY(1:ND_EXT)=ETA_EXT(1:ND_EXT)
	END IF
!
	PAR_MOM=0.0_LDP	
!
! If using the HOLLOW core option, we need to determine location to
! store inner boundary intensity. We do it here, since the storage
! location hase the same pointer for all rays.
!
	IF(INNER_BND_METH .eq. 'HOLLOW')THEN
          IF(CUR_LOC .EQ. -1)THEN
            CUR_LOC=0
            FREQ_STORE(CUR_LOC)=FREQ
          ELSE IF(FREQ_STORE(CUR_LOC) .NE. FREQ)THEN
            CUR_LOC=MOD(CUR_LOC+1,N_STORE)
            FREQ_STORE(CUR_LOC)=FREQ
	  END IF
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
! Determine radiative transfer along each p-ray
!
	NP_LIMIT=NP-1
	IF(THICK_OB)NP_LIMIT=NP
	dBdTAU=DBB/CHI(ND)			!dB/dTAU
	IPLUS=0.0_LDP
!
	IF(RAY_POINTS_INSERTED)THEN
          DO IPROC=1,NUM_RAYS_PER_THREAD
	    IP=GET_IP(MYPE,NTHREAD,IPROC)
	    IF(IP .GT. NP_LIMIT)EXIT
!
	    NRAY=RAY(IP)%NZ
	    IF(RAY_POINTS_INSERTED)THEN
              K=1
              DO I=1,RAY(IP)%NZ
100	        CONTINUE
                IF( RAY(IP)%R_RAY(I) .EQ. R_EXT(K))THEN
                  CHI_RAY(I)=CHI_EXT(K)
                  ETA_RAY(I)=ETA_EXT(K)
	        ELSE IF(RAY(IP)%R_RAY(I) .GT. R_EXT(K+1))THEN
                  T1=LOG(RAY(IP)%R_RAY(I)/R_EXT(K))
                  T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
                  CHI_RAY(I)=EXP(T2)
                  T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
                  ETA_RAY(I)=EXP(T2)
	        ELSE
	          K=K+1
                  GOTO 100
	        END IF
              END DO
	    END IF
!
! Solve using Relativistic Formal Integral
!
	    IF(USE_HEN_GREEN)THEN
	      CALL ADD_ETA_DUST(DUST_SCAT_OPAC,ETA_RAY,RAY(IP)%ETA_M,RAY(IP)%ETA_P,IP,NRAY,R,ND,NP)
              CALL SOLVE_CMF_FORMAL_V3(CHI_RAY,RAY(IP)%ETA_M,RAY(IP)%ETA_P,
	1               IP,FREQ,NU_ON_dNU,INNER_BND_METH,b_planck,dBdTAU,NRAY,NP,NC)
	    ELSE
              CALL SOLVE_CMF_FORMAL_MPI_V1(CHI_RAY,ETA_RAY,IP,FREQ,NU_ON_dNU,INNER_BND_METH,b_planck,dBdTAU,NRAY,NP,NC)
	    END IF
	  END DO
!
	ELSE
!
	  CALL TUNE(1,'FG_SOLVE')
          DO IPROC=1,NUM_RAYS_PER_THREAD
	    IP=GET_IP(MYPE,NTHREAD,IPROC)
	    IF(IP .GT. NP_LIMIT)EXIT
	    NRAY=RAY(IP)%NZ
	    IF(USE_HEN_GREEN)THEN
	      CALL ADD_ETA_DUST(DUST_SCAT_OPAC,ETA_RAY,RAY(IP)%ETA_M,RAY(IP)%ETA_P,IP,NRAY,R,ND,NP)
              CALL SOLVE_CMF_FORMAL_V3(CHI_RAY,RAY(IP)%ETA_M,RAY(IP)%ETA_P,
	1               IP,FREQ,NU_ON_dNU,INNER_BND_METH,b_planck,dBdTAU,NRAY,NP,NC)
            ELSE
	      CALL SOLVE_CMF_FORMAL_MPI_V1(CHI_RAY,ETA_RAY,IP,FREQ,NU_ON_dNU,INNER_BND_METH,b_planck,dBdTAU,NRAY,NP,NC)
	    END IF
	  END DO
	  CALL TUNE(2,'FG_SOLVE')
	END IF
!
! Integrate over p to get J and K.
!
	CALL TUNE(1,'JVAL')
        DO IPROC=1,NUM_RAYS_PER_THREAD
	  IP=GET_IP(MYPE,NTHREAD,IPROC)
	  IF(IP .GT. NP)EXIT
	  DO ID=1,MIN(ND,NP-IP+1)
	    T1=RAY(IP)%I_P(RAY(IP)%LNK(ID))
	    T2=RAY(IP)%I_M(RAY(IP)%LNK(ID))
            PAR_Jnu(ID)=PAR_Jnu(ID)+T1*ray(ip)%Jqw_p(ID)+T2*ray(ip)%Jqw_m(ID)
            PAR_Hnu(ID)=PAR_Hnu(ID)+T1*ray(ip)%Hqw_p(ID)+T2*ray(ip)%Hqw_m(ID)
            PAR_Knu(ID)=PAR_Knu(ID)+T1*ray(ip)%Kqw_p(ID)+T2*ray(ip)%Kqw_m(ID)
            PAR_Nnu(ID)=PAR_Nnu(ID)+T1*ray(ip)%Nqw_p(ID)+T2*ray(ip)%Nqw_m(ID)
	  END DO
	END DO
	CALL TUNE(2,'JVAL')
!
	IF(USE_HEN_GREEN)THEN
	  ETA_SCAT=0.0_LDP
          DO IPROC=1,NUM_RAYS_PER_THREAD
	    IP=GET_IP(MYPE,NTHREAD,IPROC)
	    IF(IP .GT. NP_LIMIT)EXIT
	    DO ID=1,MIN(ND,NP-IP+1)
	      T1=RAY(IP)%ETA_P(RAY(IP)%LNK(ID))-ETA_RAY(RAY(IP)%LNK(ID))
	      T2=RAY(IP)%ETA_M(RAY(IP)%LNK(ID))-ETA_RAY(RAY(IP)%LNK(ID))
              ETA_SCAT(ID)=ETA_SCAT(ID)+T1*ray(ip)%Jqw_p(ID)+T2*ray(Ip)%Jqw_m(ID)
	    END DO
	  END DO
	END IF
!
	ID=ND
         DO IPROC=1,NUM_RAYS_PER_THREAD
	   IP=GET_IP(MYPE,NTHREAD,IPROC)
	   IF(IP .GT. NC+1)EXIT
	   T1=RAY(IP)%I_P(RAY(IP)%LNK(ID))
	   T2=RAY(IP)%I_M(RAY(IP)%LNK(ID))
	   PAR_IB_VEC(1)=PAR_IB_VEC(1)+T1*ray(ip)%Jqw_p(ND)
	   PAR_IB_VEC(2)=PAR_IB_VEC(2)+T1*ray(ip)%Hqw_p(ND)
	   PAR_IB_VEC(3)=PAR_IB_VEC(3)+T1*ray(ip)%Kqw_p(ND)
	   PAR_IB_VEC(4)=PAR_IB_VEC(4)+T1*ray(ip)%Nqw_p(ND)
	   PAR_IB_VEC(5)=PAR_IB_VEC(5) +T2*ray(ip)%Jqw_m(ND)
	   PAR_IB_VEC(6)=PAR_IB_VEC(6) -T2*ray(ip)%Hqw_m(ND)    !- to make +ve
	   PAR_IB_VEC(7)=PAR_IB_VEC(7) +T2*ray(ip)%Kqw_m(ND)
	   PAR_IB_VEC(8)=PAR_IB_VEC(8) -T2*ray(ip)%Nqw_m(ND)
	END DO
!
! Evaluate half moments at outer boundary, and store intensity at the outer boundary.
!
	ID=1
        DO IPROC=1,NUM_RAYS_PER_THREAD
	  IP=GET_IP(MYPE,NTHREAD,IPROC)
	  IF(IP .GT. NP)EXIT                  !?MP_LIMIT
	  T1=RAY(IP)%I_P(RAY(IP)%LNK(ID))
	  T2=RAY(IP)%I_M(RAY(IP)%LNK(ID))
	  PAR_OB_VEC(1) = PAR_OB_VEC(1)  + T1*ray(ip)%Jqw_p(1)
	  PAR_OB_VEC(2) = PAR_OB_VEC(2) + T1*ray(ip)%Hqw_p(1)
	  PAR_OB_VEC(3) = PAR_OB_VEC(3) + T1*ray(ip)%Kqw_p(1)
	  PAR_OB_VEC(4) = PAR_OB_VEC(4) + T1*ray(ip)%Nqw_p(1)
	  PAR_OB_VEC(5) = PAR_OB_VEC(5) + T2*ray(ip)%Jqw_m(1)
	  PAR_OB_VEC(6) = PAR_OB_VEC(6) - T2*ray(ip)%Hqw_m(1)    !- to make +ve
	  PAR_OB_VEC(7) = PAR_OB_VEC(7) + T2*ray(ip)%Kqw_m(1)
	  PAR_OB_VEC(8) = PAR_OB_VEC(8) - T2*ray(ip)%Nqw_m(1)
	  IPLUS(IP) = T1-T2
	END DO
!
	I=4*ND+16
	CALL MPI_ALLREDUCE(MPI_IN_PLACE,PAR_MOM,I,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
        JNU_STORE(:)=PAR_MOM(1:ND)
        HNU_STORE(:)=PAR_MOM(ND+1:2*ND)
        KNU_STORE(:)=PAR_MOM(2*ND+1:3*ND)
        NNU_STORE(:)=PAR_MOM(3*ND+1:4*ND)
!
        K=4*ND
        JPLUS_IB=PAR_MOM(K+1);  HPLUS_IB=PAR_MOM(K+2);  KPLUS_IB=PAR_MOM(K+3);    NPLUS_IB=PAR_MOM(K+4)
        JMIN_IB=PAR_MOM(K+5);   HMIN_IB=PAR_MOM(K+6);   KMIN_IB=PAR_MOM(K+7);     NMIN_IB=PAR_MOM(K+8)
        JPLUS_OB=PAR_MOM(K+9);  HPLUS_OB=PAR_MOM(K+10); KPLUS_OB=PAR_MOM(K+11);   NPLUS_OB=PAR_MOM(K+12)
        JMIN_OB=PAR_MOM(K+13);  HMIN_OB=PAR_MOM(K+14);  KMIN_OB=PAR_MOM(K+15);    NMIN_OB=PAR_MOM(K+16)
!
! Save intensity for integration at next frequency.
!
	CALL TUNE(1,'FGP_SAVE')
        DO IPROC=1,NUM_RAYS_PER_THREAD
	  IP=GET_IP(MYPE,NTHREAD,IPROC)
	  IF(IP .GT. NP)EXIT                   !?MP_LIMIT
	  DO ID=1,RAY(IP)%NZ
	    RAY(IP)%I_P_SAVE(ID)=RAY(IP)%I_P(ID)
	    RAY(IP)%I_M_SAVE(ID)=RAY(IP)%I_M(ID)
	  END DO
	END DO
	CALL TUNE(2,'FGP_SAVE')
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
	JNU=JNU_STORE
	FEDD=KNU_STORE/JNU_STORE
	HN_DEF_ON_NODES=.TRUE.
	HNU_AT_OB=Hnu_store(1);  NNU_AT_OB=Nnu_store(1)
	HNU_AT_IB=Hnu_store(ND); NNU_AT_IB=Nnu_store(ND)
!
	HBC=HNU_AT_OB/JNU_STORE(1)
	NBC=NNU_AT_OB/JNU_STORE(1)
	IF(HBC .LT. 0.0_LDP)HBC=0.0_LDP
	IF(NBC .LT. 0.0_LDP)NBC=0.0_LDP
!
	RET_HNU_AT_OB=HNU_AT_OB
	RET_HNU_AT_IB=HNU_AT_IB
!
	T1=1.0E-200_LDP
!
	FIRST_TIME=.FALSE.
	RETURN
	END
