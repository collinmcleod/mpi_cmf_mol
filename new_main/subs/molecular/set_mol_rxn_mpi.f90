!
! Subroutine to compute important parameters related to the molecular reactions.  These are
  ! AI_AR_MOL
  ! dlnAI_AR_MOL_dlnT
  ! COOL_MOL
! These arrays are used to constrain LTE properties and compute cooling rates (so as to keep energy conserved)
! The energy exchange may not be perfect due to the fudge we are using--assuming carbon monoxide emits in LTE
! 
! Unlike charge exchange reactions, the molecular reactions do not have easily computed reaction rates
! In fact, many of the included reactions either have reverse reactions included explicitly as separate reactions
! or their reverse reaction rates are negligible (for instance, some reactions have three products and thus have
! a negligible cross-section for the reverse rates)
!
  SUBROUTINE SET_MOL_RXN(ND)
    USE MOL_RXN_MOD
    USE MOD_CMFGEN
    USE SET_KIND_MODULE
    IMPLICIT NONE
!    
! Edited 22-Sep-2021: Changed things again--because I am not computing reverse reaction rates (they are either 
! included explicitly or assumed to be negligible) this routine only computes overall cooling rates
!
! Edited 22-Sep-2021: Significantly altered the overall structure of the subroutine--The primary loop is now over 
!                     the reaction list.  Unlike set_chg_exch_v4, this routine now needs to be called only once
!                     in the main program run, and will loop over all the reactions inside this subroutine
!                     This routine DOES NOT need to be called for each included ion
! Created 22-Sep-2021: based on set_chg_exch_v4, which handles the same issues for the charge exchange reactions 
!
    REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
    COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
    INTEGER ND
!
! Atomic Variables
!
    INTEGER ID
    INTEGER I_S,I_F
    INTEGER N_S,N_F
    INTEGER EQSPEC
    INTEGER EQHYD
!
    REAL(KIND=LDP) GION
!
! Local variables
!
    REAL(KIND=LDP) G_MOL_VEC(ND)
    REAL(KIND=LDP) DG_MOL_VEC(ND)
    REAL(KIND=LDP) T_VEC(ND)
    REAL(KIND=LDP) T1
!
    INTEGER I,J,K,L
    INTEGER IOS
    INTEGER,PARAMETER :: LU=119
    CHARACTER(LEN=30) LOC_NAME
    LOGICAL LEVEL_SET
    LOGICAL, SAVE :: FIRST_TIME
    REAL(KIND=LDP), ALLOCATABLE LHS_ENERGY(:,:),RHS_ENERGY(:,:),ENERGY_CHANGE(:,:)
    DATA FIRST_TIME/.TRUE./
!
!
    IF (.NOT. DO_MOL_RXNS) RETURN
!
    IF (FIRST_TIME) THEN
       WRITE(LU,*) ' '
       WRITE(LU,*) 'Check of reverse reaction factors for molecular reactions, at depth 1'
       WRITE(LU,*) 'The last AI/AR ratio is the conversion factor for a given reaction'
       WRITE(LU,*)' '
       WRITE(LU,'(1X,A,2X,A,2X,A,6X,A,9X,A,8X,A)')'Reaction','Col','ID','AI/AR','G_EXP','GION'
       FIRST_TIME=.FALSE.
    END IF
!
    IF (INITIALIZE_ARRAYS) THEN
       AI_AR_MOL(:,:)=1.0_LDP
       DLNAI_AR_MOL_DLNT(:,:)=0.0_LDP
       COOL_MOL(:,:)=0.0_LDP
       INITIALIZE_ARRAYS=.FALSE.
       WRITE(LU,*)' '
       WRITE(LU,'(1X,A,ES14.4)')'Temperature (10^4 K): ',T(1)
       WRITE(LU,*)' '
    END IF
!
    ALLOCATE(LHS_ENERGY(N_MOL_RXN,ND),STAT=IOS)
    IF (IOS .EQ. 0) ALLOCATE(RHS_ENERGY(N_MOL_RXN,ND),STAT=IOS)
    IF (IOS .NE. 0) THEN
       WRITE(LUER1,*)'Error in SET_MOL_RXN'
       WRITE(LUER1,*)'Unable to allocate memory for local variables'
       STOP
    END IF
!
! Now determine the energy change of various reactions, using the ionization energies. Bond energies may need to 
! be included later.
! I am assuming that any net change of energy due to a molecular reaction will be transferred to the electron 
! thermal pool, which may or may not be a good assumption
!

!       DO J=1,N_MOL_RXN
!          DO ID=1,NUM_IONS-1
!             LEVEL_SET = .FALSE.
!             DO K=1,4
                
                
    DO ID=1,NUM_IONS-1
       IF(ATM(ID)%XzV_PRES) THEN
       DO J=1,N_MOL_RXN
          LEVEL_SET=.FALSE.
          G_MOL_VEC(1:ND) = 0.0_LDP !Initialize
          DG_MOL_VEC(1:ND) = 0.0_LDP
          DO K=1,4
             IF (ID .EQ. ID_ION_MOL(J,K)) THEN
                I_F = 
                DO I_F=1,N_F
                   I_S=F_TO_S(I_F)
                   IF(I_S .EQ. LEV_IN_ION_MOL(J,K)) THEN
                      LEVEL_SET=.TRUE.
                      IF (MOL_Z_CHG(J,K) .GT. 1) THEN
                         T_VEC(1:ND)=HDKT*(EDGE_F(I_F)-EDGE_F(1))/T(1:ND)
                      ELSE
                         T_VEC(1:ND)=HDKT*EDGE_F(I_F)/T(1:ND)
                      END IF
                      DO L=1,ND
                         T1 = EXP(T_VEC(L))
                         G_MOL_VEC(L)=G_MOL_VEC(L)+G_F(I_F)*T1
                         DG_MOL_VEC(L)=DG_MOL_VEC(L)-G_F(I_F)*T_VEC(L)*T1
                      END DO
                   END IF
                END DO
!
                IF (LEVEL_SET) THEN

! Because of differing conventions, I may have to calculate the cooling rate from my reactions in a different
! way from how the charge exchange reactions are done... I will have to come back to this

!                   DO L=1,ND
!                      COOL_MOL(J,L)=COOL_MOL(J,L) + (DG_MOL_VEC(L)/G_MOL_VEC(L))*T(L)/HDKT
                      
                   IF (MOL_RXN_IS_FINAL_ION(J,K)) THEN
                      G_MOL_VEC(1:ND) = G_MOL_VEC(1:ND)/GION
                      DG_MOL_VEC(1:ND) = DG_MOL_VEC(1:ND)/GION
                   END IF
!
                   IF (K .EQ. 1 .OR. K .EQ. 2) THEN
                      AI_AR_MOL(J,1:ND)=AI_AR_MOL(J,1:ND)*G_MOL_VEC(1:ND)
                      WRITE(LU,'(1X,3I6,3ES14.4)')J,K,ID_ION_MOL(J,K),AI_AR_MOL(J,1),G_MOL_VEC(1),GION
                      DLNAI_AR_MOL_DLNT(J,1:ND)=DLNAI_AR_MOL_DLNT(J,1:ND)+DG_MOL_VEC(1:ND)/G_MOL_VEC(1:ND)           
!                   
                   ELSE
                      AI_AR_MOL(J,1:ND)=AI_AR_MOL(J,1:ND)/G_MOL_VEC(1:ND)
                      WRITE(LU,'(1X,3I6,3ES14.4)')J,K,ID_ION_MOL(J,K),AI_AR_MOL(J,1),G_MOL_VEC(1),GION
                      DLNAI_AR_MOL_DLNT(J,1:ND)=DLNAI_AR_MOL_DLNT(J,1:ND)-DG_MOL_VEC(1:ND)/G_MOL_VEC(1:ND)
                   END IF
                ELSE
                   WRITE(LUER1,*)'*********************************'
                   WRITE(LUER1,*)'Error in SET_MOL_RXN'
                   WRITE(LUER1,*)'Species matches, but level name does not'
                   WRITE(LUER1,*)'Molecular reaction number ',J
                   WRITE(LUER1,*)'Species number ',K
                   WRITE(LUER1,*)'Check for consistency with level names'
                   WRITE(LUER1,*)'*********************************'
                END IF  !Level set
             END IF     !Species match
          END DO        !Over K
       END DO           !Over J
    END IF
 END DO
!
       RETURN
     END SUBROUTINE SET_MOL_RXN
