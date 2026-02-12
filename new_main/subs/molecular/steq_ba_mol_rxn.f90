!
  ! Routine to modify the statistical equilibrium (STEQ) and variation of statistical equilibrium (BA) equations
  ! for the contributions from molecular reactions.  This routine must be preceded by calls of RD_MOL_RXN, 
  ! SET_MOL_RXN_LEV, and calls of SET_MOL_RXN for all included species
  !
  SUBROUTINE STEQ_BA_MOL_RXN(POPS,T,ED,NT,ND,DIAG_INDX,UPDATE_BA)
    USE MOL_RXN_MOD
    USE STEQ_DATA_MOD
!    USE MOD_CMFGEN
    IMPLICIT NONE
!
! Edited 18-Feb-2022: Moved most of the important output to print_mol_info.f90, which now takes care of
!                     split levels and photoionization/recombination rates
! Edited 6-Dec-2021: added more printing routines to look at rates for each species
! Edited 27-Sep-2021: Finished treatment of ELECTRON and NOATOM reactions, changed the way routine updates electron
!                     population
! Edited 23-Sep-2021: Redesigned to fit philosophy for molecular reactions, reverse reaction rates are not computed
! reverse reaction rates are included explicitly as separate reactions or assumed to be negligible
! Created 22-Sep-2021 by Collin McLeod, based on existing routine STEQ_BA_CHG_EXCH_V3
! See that routine for additional information on variables
!
! External/argument variables
!
    INTEGER NT
    INTEGER ND
    INTEGER DIAG_INDX
!
    REAL(8) POPS(NT,ND)
    REAL(8) T(ND)
    REAL(8) ED(ND)
    LOGICAL UPDATE_BA
!
    INTEGER LU_ER
    INTEGER ERROR_LU
    EXTERNAL ERROR_LU
!
! Local variables
!
    INTEGER I,J,K,L          
    INTEGER E1,E2,E3,E4
    INTEGER L1,L2,L3,L4
    INTEGER L11,L12,L13,L14
    INTEGER L21,L22,L23,L24
    INTEGER L31,L32,L33,L34
    INTEGER L41,L42,L43,L44
    INTEGER NIV
    INTEGER II
    INTEGER ID,ID1,ID2,ID3,ID4
    REAL(8) RATE_K,SCALED_K
    REAL(8), ALLOCATABLE :: RATE_KS(:,:)
    REAL(8) RATE_K_ION
    REAL(8) ALPHA_1,ALPHA_2,ALPHA_3,ALPHA_4
    REAL(8) DRATEDT
    REAL(8) DLN_K_DLNT,DLN_K_ION_DLNT
    REAL(8) FRD_R,REV_R
    REAL(8) T1,TVAL,TREAL1,TREAL2,TREAL3,TVAL4
    REAL(8) TEN4
!
    ! Some added variables useful for when I want to print out lots of information
!
    LOGICAL :: PRINT_RXN_INFO = .TRUE.
    LOGICAL :: PRINT_DEBUG_INFO = .TRUE.
    INTEGER :: MOL_FILE_NUM = 69
    CHARACTER(30) :: MOL_OUTPUT_FILE = 'MOLECULAR_REACTION_OUTPUT'
    INTEGER :: IOS,LDIFF,CI1,CI2
    CHARACTER(30) :: MOL_OUT_FILE_UNSORT = 'MOLECULAR_REACTIONS_OUT'
    CHARACTER(30) :: MOL_OUT_FILE_SORT = 'MOLECULAR_REACTIONS_SORTED'
    REAL(8), ALLOCATABLE :: MOL_RXN_RATES_ALL(:,:), MOL_RXN_RATES_SORTED(:,:,:)
!
    CHARACTER(30) :: MOL_RXNS_SPECIES_FILE='MOL_RXNS_BY_SPECIES'
!    REAL(8), ALLOCATABLE :: TOTAL_CHANGE_BY_SPECIES(:,:)
!    INTEGER, ALLOCATABLE :: TOP_RXNS_BY_SPECIES(:,:,:,:)
!    CHARACTER(20), ALLOCATABLE :: REACTANTS(:), PRODUCTS(:), 
    CHARACTER(20), ALLOCATABLE :: SPECIESIDLIST(:)
    CHARACTER(20) :: SPECIESID1
    INTEGER :: SPECIESNUM
    INTEGER, ALLOCATABLE :: MOL_LEVID_LIST(:), MOL_ATMID_LIST(:)
!
    LOGICAL :: WRITE_OUT_POPS=.TRUE.
    LOGICAL :: IS_OPEN

    CHARACTER(30) :: BA_CHK_FILE='MOL_RXN_BA_CHK'
!
    INTEGER, ALLOCATABLE :: LEVS_ACCOUNTED_FOR(:)
    INTEGER :: NLEVS_ACCOUNTED
!
    IF (.NOT. DO_MOL_RXNS) RETURN
!
    TEN4 = 10.0**4 !A useful factor to convert temperature units
!
    LU_ER=ERROR_LU()
!
    IF (PRINT_RXN_INFO) THEN
       ALLOCATE(MOL_RXN_RATES_ALL(ND,N_MOL_RXN+1),STAT=IOS)
       IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_RATES_SORTED(ND,12,2),STAT=IOS)
       IF (IOS .EQ. 0) THEN
          MOL_RXN_RATES_ALL = 0.0
          MOL_RXN_RATES_SORTED = 0.0
       ELSE
          WRITE(LU_ER,*)'Error in STEQ_BA_MOL_RXN: '
          WRITE(LU_ER,*)'Unable to allocate arrays for checking rates'
       END IF
    END IF
!
    ALLOCATE(RATE_KS(ND,N_MOL_RXN),STAT=IOS)
    IF (IOS .EQ. 0) THEN
       RATE_KS = 0.0
    ELSE
       WRITE(LU_ER,*)'Error in STEQ_BA_MOL_RXN: '
       WRITE(LU_ER,*)'Unable to allocate arrays for rate constants'
    END IF
!
    IF (.NOT. ALLOCATED(LVL_IN_MOL_RXNS)) ALLOCATE(LVL_IN_MOL_RXNS(NT),STAT=IOS)
    IF (IOS .NE. 0) THEN
       WRITE(LU_ER,*)'Error in STEQ_BA_MOL_RXN: '
       WRITE(LU_ER,*)'Unable to allocate LVL_IN_MOL_RXNS'
    END IF
    LVL_IN_MOL_RXNS(:)=.FALSE.
!
    IF(.NOT. ALLOCATED(MTOT_BA)) THEN
       ALLOCATE(MTOT_BA(N_MOL_LEVS,N_MOL_LEVS+2,ND),STAT=IOS)
       IF (IOS .NE. 0) WRITE(LU_ER,*) 'Error in steq_ba_mol_rxn: unable to allocate mtot_ba'
    END IF
    IF(.NOT. ALLOCATED(MTOT_STEQ)) THEN
       ALLOCATE(MTOT_STEQ(N_MOL_LEVS,ND),STAT=IOS)
       IF (IOS .NE. 0) WRITE(LU_ER,*) 'Error in steq_ba_mol_rxn: unable to allocate mtot_steq'
    END IF
    MTOT_STEQ(:,:)=0.0D0
    MTOT_BA(:,:,:)=0.0D0
    IF (.NOT. ALLOCATED(LEVS_ACCOUNTED_FOR)) THEN
       ALLOCATE(LEVS_ACCOUNTED_FOR(N_MOL_LEVS),STAT=IOS)
    END IF
    LEVS_ACCOUNTED_FOR(:)=0
    NLEVS_ACCOUNTED=0
    IF (.NOT. ALLOCATED(LNK_MOL_TO_FULL)) THEN
       ALLOCATE(LNK_MOL_TO_FULL(N_MOL_LEVS+2),STAT=IOS)
       IF (IOS .NE. 0) WRITE(LU_ER,*) 'Error in steq_ba_mol_rxn: unable to allocate LNK_MOL_TO_FULL'
    END IF
    LNK_MOL_TO_FULL(:)=0
    IF (.NOT. ALLOCATED(LNK_FULL_TO_MOL)) THEN
       ALLOCATE(LNK_FULL_TO_MOL(NT),STAT=IOS)
       IF (IOS .NE. 0) WRITE(LU_ER,*) 'Error in steq_ba_mol_rxn: unable to allocate LNK_FULL_TO_MOL'
    END IF
    LNK_FULL_TO_MOL(:)=0
    !
! MTOT_BA stores all the molecular reaction contributions so that they can later be included
! in the ion equation for each ion (in the big matrix)
!
    DO J=1,N_MOL_RXN
       IF (MOL_REACTION_AVAILABLE(J)) THEN
!
! Set proper level associations
!
          L1=LEV_IN_POPS_MOL(J,1)
          L2=LEV_IN_POPS_MOL(J,2)
          L3=LEV_IN_POPS_MOL(J,3)
          L4=LEV_IN_POPS_MOL(J,4)
!
! Set links for values in MTOT arrays
!
          IF (LNK_MOL_TO_FULL(N_MOL_LEVS+1) .EQ. 0) LNK_MOL_TO_FULL(N_MOL_LEVS+1)=NT-1
          IF (LNK_FULL_TO_MOL(NT-1) .EQ. 0) LNK_FULL_TO_MOL(NT-1)=N_MOL_LEVS+1
          IF (LNK_MOL_TO_FULL(N_MOL_LEVS+2) .EQ. 0) LNK_MOL_TO_FULL(N_MOL_LEVS+2)=NT
          IF (LNK_FULL_TO_MOL(NT) .EQ. 0) LNK_FULL_TO_MOL(NT)=N_MOL_LEVS+2
          IF (.NOT. ANY(LEVS_ACCOUNTED_FOR .EQ. L1)) THEN
             NLEVS_ACCOUNTED = NLEVS_ACCOUNTED + 1
             LEVS_ACCOUNTED_FOR(NLEVS_ACCOUNTED)=L1
             LNK_MOL_TO_FULL(NLEVS_ACCOUNTED)=L1
             LNK_FULL_TO_MOL(L1)=NLEVS_ACCOUNTED
          END IF
          IF (.NOT. ANY(LEVS_ACCOUNTED_FOR .EQ. L2)) THEN
             IF (L2 .GT. 0) THEN
                NLEVS_ACCOUNTED = NLEVS_ACCOUNTED + 1
                LEVS_ACCOUNTED_FOR(NLEVS_ACCOUNTED)=L2
                LNK_MOL_TO_FULL(NLEVS_ACCOUNTED)=L2
                LNK_FULL_TO_MOL(L2)=NLEVS_ACCOUNTED
             END IF
          END IF
          IF (.NOT. ANY(LEVS_ACCOUNTED_FOR .EQ. L3)) THEN
             NLEVS_ACCOUNTED = NLEVS_ACCOUNTED + 1
             LEVS_ACCOUNTED_FOR(NLEVS_ACCOUNTED)=L3
             LNK_MOL_TO_FULL(NLEVS_ACCOUNTED)=L3
             LNK_FULL_TO_MOL(L3)=NLEVS_ACCOUNTED
          END IF
          IF (.NOT. ANY(LEVS_ACCOUNTED_FOR .EQ. L4)) THEN
             IF (L4 .GT. 0) THEN
                NLEVS_ACCOUNTED = NLEVS_ACCOUNTED + 1
                LEVS_ACCOUNTED_FOR(NLEVS_ACCOUNTED)=L4
                LNK_MOL_TO_FULL(NLEVS_ACCOUNTED)=L4
                LNK_FULL_TO_MOL(L4)=NLEVS_ACCOUNTED
             END IF
          END IF
!
! Now I will update my logical array to indicate that certain levels in pops are involved in 
! the molecular reactions
!
          IF (L1 .GT. 0) LVL_IN_MOL_RXNS(L1)=.TRUE.
          IF (L2 .GT. 0) LVL_IN_MOL_RXNS(L2)=.TRUE.
          IF (L3 .GT. 0) LVL_IN_MOL_RXNS(L3)=.TRUE.
          IF (L4 .GT. 0) LVL_IN_MOL_RXNS(L4)=.TRUE.
!
! I set the levid variables which will be used in later write statements
!

          IF (.NOT. ALLOCATED(MOL_LEVID_LIST)) THEN
             ALLOCATE(MOL_LEVID_LIST(14),STAT=IOS)
             IF (IOS .EQ. 0) ALLOCATE(SPECIESIDLIST(14),STAT=IOS)
             IF (IOS .EQ. 0) ALLOCATE(MOL_ATMID_LIST(14),STAT=IOS)
             MOL_LEVID_LIST = 0

             SPECIESIDLIST(1)='CI'
             SPECIESIDLIST(2)='C2'
             SPECIESIDLIST(3)='OI'
             SPECIESIDLIST(4)='O2'
             SPECIESIDLIST(5)='COMI'
             SPECIESIDLIST(6)='COM2'
             SPECIESIDLIST(7)='C2MI'
             SPECIESIDLIST(8)='C2M2'
             SPECIESIDLIST(9)='O2MI'
             SPECIESIDLIST(10)='O2M2'
             SPECIESIDLIST(11)='CO2MI'
             SPECIESIDLIST(12)='CO2M2'
             SPECIESIDLIST(13)='C2OMI'
             SPECIESIDLIST(14)='C2OM2'
          END IF
!
             DO K=1,4
                DO I=1,14
                   IF (MOL_LEVID_LIST(I) .NE. 0) CYCLE
                   IF(TRIM(ADJUSTL(MOL_RXN_SPECIESID(J,K))) .EQ. TRIM(ADJUSTL(SPECIESIDLIST(I)))) THEN
                      IF (K .EQ. 1) THEN
                         MOL_LEVID_LIST(I) = L1
                      ELSE IF (K .EQ. 2) THEN
                         MOL_LEVID_LIST(I) = L2
                      ELSE IF (K .EQ. 3) THEN
                         MOL_LEVID_LIST(I) = L3
                      ELSE IF (K .EQ. 4) THEN
                         MOL_LEVID_LIST(I) = L4
                      END IF
                      MOL_ATMID_LIST(I) = ID_ION_MOL(J,K)
                   END IF
                END DO
             END DO

!
! I have introduced a convention whereby the level ID is -1 if the reactant/product is an electron,
! and -2 if the product is a photon (in which case it is assumed to immediately deposit its energy in the electron
! thermal pool)
!
! By my conventions, an electron (if involved) can only appear in position 2, and a photon (or NOATOM) will only
! appear in position 4.  These checks will be added at the end
!
          IF (L1 .GT. 0 .AND. L2 .GT. 0 .AND. L3 .GT. 0 .AND. L4 .GT. 0) THEN 
             !Check whether electrons/photons involved, if not then go ahead
!
! Refer to the ion states of the involved species
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,1)) THEN
                ID1=ID_ION_MOL(J,1)
             ELSE
                ID1=ID_ION_MOL(J,1)-1
             END IF
             L11=SE(ID1)%LNK_TO_IV(L1)
             L12=SE(ID1)%LNK_TO_IV(L2)
             L13=SE(ID1)%LNK_TO_IV(L3)
             L14=SE(ID1)%LNK_TO_IV(L4)
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,2)) THEN
                ID2=ID_ION_MOL(J,2)
             ELSE
                ID2=ID_ION_MOL(J,2)-1
             END IF
             L21=SE(ID2)%LNK_TO_IV(L1)
             L22=SE(ID2)%LNK_TO_IV(L2)
             L23=SE(ID2)%LNK_TO_IV(L3)
             L24=SE(ID2)%LNK_TO_IV(L4)
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,3)) THEN
                ID3=ID_ION_MOL(J,3)
             ELSE
                ID3=ID_ION_MOL(J,3)-1
             END IF
             L31=SE(ID3)%LNK_TO_IV(L1)
             L32=SE(ID3)%LNK_TO_IV(L2)
             L33=SE(ID3)%LNK_TO_IV(L3)
             L34=SE(ID3)%LNK_TO_IV(L4)
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,4)) THEN
                ID4=ID_ION_MOL(J,4)
             ELSE
                ID4=ID_ION_MOL(J,4)-1
             END IF
             L41=SE(ID4)%LNK_TO_IV(L1)
             L42=SE(ID4)%LNK_TO_IV(L2)
             L43=SE(ID4)%LNK_TO_IV(L3)
             L44=SE(ID4)%LNK_TO_IV(L4)
!
! Now associate the relevant statistical equilibrium equations, using the links from STEQ_DATA_MOD
!
             IF (MOL_RXN_IS_FINAL_ION(J,1)) THEN
                E1=SE(ID1)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,1))
             ELSE
                E1 = LEV_IN_ION_MOL(J,1)
             END IF
!
             IF (MOL_RXN_IS_FINAL_ION(J,2)) THEN
                E2=SE(ID2)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,2))
             ELSE
                E2 = LEV_IN_ION_MOL(J,2)
             END IF
!
             IF (MOL_RXN_IS_FINAL_ION(J,3)) THEN
                E3=SE(ID3)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,3))
             ELSE
                E3 = LEV_IN_ION_MOL(J,3)
             END IF
!
             IF (MOL_RXN_IS_FINAL_ION(J,4)) THEN
                E4=SE(ID4)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,4))
             ELSE
                E4 = LEV_IN_ION_MOL(J,4)
             END IF
!
             IF (E1*E2*E3*E4 .EQ. 0) THEN
                WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: '
                WRITE(LU_ER,*) 'At least one of E1...E4 is zero--check ion levels'
                WRITE(LU_ER,*) E1,E2,E3,E4
                STOP
             END IF
!
             II = DIAG_INDX
!
             DO L=1,ND
                TVAL=MAX(T(L),MOL_RXN_TLO(J))
                TVAL=MIN(TVAL,MOL_RXN_THI(J))
                TVAL4 = TEN4*TVAL
!
! Calculate rate constants and their variation with temperature
!
                IF((MOL_RXN_TYPE(J) .EQ. 1) .OR. (MOL_RXN_TYPE(J) .EQ. 2)) THEN
                   TREAL1 = EXP(-MOL_RXN_COEF(J,3)/(TVAL4))
                   TREAL2 = MOL_RXN_COEF(J,1)*((TVAL4/300.)**MOL_RXN_COEF(J,2))
                   RATE_K = TREAL1*TREAL2
                   DLN_K_DLNT = MOL_RXN_COEF(J,2) + (MOL_RXN_COEF(J,3)/(TVAL4))
                ELSE
                   LU_ER=ERROR_LU()
                   WRITE(LU_ER,*)'Error in STEQ_BA_MOL_RXN: '
                   WRITE(LU_ER,*)'Invalid reaction type passed,'
                   WRITE(LU_ER,*)'Reaction type: ',MOL_RXN_TYPE(J)
                   STOP
                END IF
!                
                RATE_KS(L,J) = RATE_K
                SCALED_K = RXN_SCALE_FAC*MOL_RXN_RATE_SCALE(J)*RATE_K
                IF (MOL_RXN_SPECIESID(J,3) .EQ. 'COMI' .OR. MOL_RXN_SPECIESID(J,4) .EQ. 'COMI') THEN
                   SCALED_K = SCALED_K*CO_SCALE_FAC
                END IF
!
                IF((TVAL .GT. MOL_RXN_THI(J)) .OR. (TVAL .LT. MOL_RXN_TLO(J))) DLN_K_DLNT = 0.0D0
!
                FRD_R=SCALED_K*POPS(L1,L)*POPS(L2,L) 
!
! The forward reaction rate, in my case the only rate considered. Note that we account for multiple product states
! rate is number of reactions /cm**3 /s
!
! Update statistical equilibrium
!
                SE(ID1)%STEQ(E1,L)=SE(ID1)%STEQ(E1,L) - FRD_R
                SE(ID2)%STEQ(E2,L)=SE(ID2)%STEQ(E2,L) - FRD_R
!
                SE(ID3)%STEQ(E3,L)=SE(ID3)%STEQ(E3,L) + FRD_R
                SE(ID4)%STEQ(E4,L)=SE(ID4)%STEQ(E4,L) + FRD_R                
!
                MTOT_STEQ(LNK_FULL_TO_MOL(L1),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L1),L) - FRD_R
                MTOT_STEQ(LNK_FULL_TO_MOL(L2),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L2),L) - FRD_R
                MTOT_STEQ(LNK_FULL_TO_MOL(L3),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L3),L) + FRD_R
                MTOT_STEQ(LNK_FULL_TO_MOL(L4),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L4),L) + FRD_R
!
                DRATEDT = (FRD_R*DLN_K_DLNT)/TVAL
!
                IF (UPDATE_BA) THEN
                   ALPHA_1 = SCALED_K*POPS(L1,L)
                   ALPHA_2 = SCALED_K*POPS(L2,L)
                   ALPHA_3 = SCALED_K*POPS(L3,L)
                   ALPHA_4 = SCALED_K*POPS(L4,L)
!
                   IF (L1 .EQ. L2) THEN
!
! When the reactants are the same, the change to the molecule/atom will be twice the reaction rate
! There's an additional factor of two when considering dR1/dR1 since the rate is propto R1^2
!                      
                      ID=ID1
                      NIV=SE(ID)%N_IV
                      IF(L11 .NE. 0) SE(ID)%BA(E1,L11,II,L)=SE(ID)%BA(E1,L11,II,L) - 4*ALPHA_1
                      SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) - 2*DRATEDT
! 
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) - &
4*ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) - &
2*DRATEDT
!
                      ID=ID3
                      NIV=SE(ID)%N_IV
                      IF(L31 .NE. 0) SE(ID)%BA(E3,L31,II,L)=SE(ID)%BA(E3,L31,II,L) + 2*ALPHA_1
                      SE(ID)%BA(E3,NIV,II,L)=SE(ID)%BA(E3,NIV,II,L) + DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) + &
2*ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) + &
DRATEDT
!
                      ID=ID4
                      NIV=SE(ID)%N_IV
                      IF(L41 .NE. 0) SE(ID)%BA(E4,L41,II,L)=SE(ID)%BA(E4,L41,II,L) + 2*ALPHA_1
                      SE(ID)%BA(E4,NIV,II,L)=SE(ID)%BA(E4,NIV,II,L) + DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(L1),L) + &
2*ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(NT),L) + &
DRATEDT
!

                   ELSE IF (L3 .EQ. L4) THEN
!
                      ID=ID1
                      NIV=SE(ID)%N_IV
                      IF(L11 .NE. 0) SE(ID)%BA(E1,L11,II,L)=SE(ID)%BA(E1,L11,II,L) - ALPHA_2
                      IF(L12 .NE. 0) SE(ID)%BA(E1,L12,II,L)=SE(ID)%BA(E1,L12,II,L) - ALPHA_1
                      SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) - DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) - &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L2),L) - &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) - &
DRATEDT
!
                      ID=ID2
                      NIV=SE(ID)%N_IV
                      IF(L21 .NE. 0) SE(ID)%BA(E2,L21,II,L)=SE(ID)%BA(E2,L21,II,L) - ALPHA_2
                      IF(L22 .NE. 0) SE(ID)%BA(E2,L22,II,L)=SE(ID)%BA(E2,L22,II,L) - ALPHA_1
                      SE(ID)%BA(E2,NIV,II,L)=SE(ID)%BA(E2,NIV,II,L) - DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L1),L) - &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L2),L) - &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(NT),L) - &
DRATEDT
!
                      ID=ID3
                      NIV=SE(ID)%N_IV
                      IF(L31 .NE. 0) SE(ID)%BA(E3,L31,II,L)=SE(ID)%BA(E3,L31,II,L) + 2*ALPHA_2
                      IF(L32 .NE. 0) SE(ID)%BA(E3,L32,II,L)=SE(ID)%BA(E3,L32,II,L) + 2*ALPHA_1
                      SE(ID)%BA(E3,NIV,II,L)=SE(ID)%BA(E3,NIV,II,L) + 2*DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) + &
2*ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L2),L) + &
2*ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) + &
2*DRATEDT
!
                   ELSE
!
                      ID=ID1
                      NIV=SE(ID)%N_IV
                      IF(L11 .NE. 0) SE(ID)%BA(E1,L11,II,L)=SE(ID)%BA(E1,L11,II,L) - ALPHA_2
                      IF(L12 .NE. 0) SE(ID)%BA(E1,L12,II,L)=SE(ID)%BA(E1,L12,II,L) - ALPHA_1
                      SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) - DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) - &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L2),L) - &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) - &
DRATEDT
!
                      ID=ID2
                      NIV=SE(ID)%N_IV
                      IF(L21 .NE. 0) SE(ID)%BA(E2,L21,II,L)=SE(ID)%BA(E2,L21,II,L) - ALPHA_2
                      IF(L22 .NE. 0) SE(ID)%BA(E2,L22,II,L)=SE(ID)%BA(E2,L22,II,L) - ALPHA_1
                      SE(ID)%BA(E2,NIV,II,L)=SE(ID)%BA(E2,NIV,II,L) - DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L1),L) - &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L2),L) - &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(NT),L) - &
DRATEDT
!
                      ID=ID3
                      NIV=SE(ID)%N_IV
                      IF(L31 .NE. 0) SE(ID)%BA(E3,L31,II,L)=SE(ID)%BA(E3,L31,II,L) + ALPHA_2
                      IF(L32 .NE. 0) SE(ID)%BA(E3,L32,II,L)=SE(ID)%BA(E3,L32,II,L) + ALPHA_1
                      SE(ID)%BA(E3,NIV,II,L)=SE(ID)%BA(E3,NIV,II,L) + DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) + &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L2),L) + &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) + &
DRATEDT
!
                      ID=ID4
                      NIV=SE(ID)%N_IV
                      IF(L41 .NE. 0) SE(ID)%BA(E4,L41,II,L)=SE(ID)%BA(E4,L41,II,L) + ALPHA_2
                      IF(L42 .NE. 0) SE(ID)%BA(E4,L42,II,L)=SE(ID)%BA(E4,L42,II,L) + ALPHA_1
                      SE(ID)%BA(E4,NIV,II,L)=SE(ID)%BA(E4,NIV,II,L) + DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(L1),L) + &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(L2),L) + &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(NT),L) + &
DRATEDT
!
                   END IF
!
                !    Add some print statements for debugging, use variable PRINT_RXN_INFO to toggle
        IF (PRINT_DEBUG_INFO) THEN
           IF (L .EQ. 1 .AND. J .EQ. 1) THEN
              OPEN(UNIT=MOL_FILE_NUM,FILE=MOL_OUTPUT_FILE,STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
              WRITE(MOL_FILE_NUM,*)'***************************************************************'
              WRITE(MOL_FILE_NUM,*)'Info for molecular reactions during evaluation, final iteration'
              WRITE(MOL_FILE_NUM,*)'First line: reaction'
              WRITE(MOL_FILE_NUM,*)'Second line: Depth, Temperature (10^4 K), Species 1 population, Species 2 population'
              WRITE(MOL_FILE_NUM,*)'Third line: Reaction coefficients Alpha, Beta, Gamma'
              WRITE(MOL_FILE_NUM,*)'Fourth line: Rate constant, total rate, DrateDT'
              WRITE(MOL_FILE_NUM,*)'Fifth line: IDs of levels involved, in order. IDs are skipped for electrons or photons'
              WRITE(MOL_FILE_NUM,*)'Sixth line: alphas 1-4 where applicable'
              WRITE(MOL_FILE_NUM,*)'**************************************************************'
              WRITE(MOL_FILE_NUM,*)
              IF (IOS .NE. 0) THEN
                 WRITE(LUER1,*)'Error opening file to print molecular reaction debugging info'
              END IF
           END IF
           WRITE(MOL_FILE_NUM,'(A16,I2,A1,2X,A8,1X,A1,1X,A8,1X,A2,1X,A8,1X,A1,1X,A8)')'Reaction Number ',&
                J,':',MOL_RXN_SPECIESID(J,1),'+',MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4)
           WRITE(MOL_FILE_NUM,'(I3,A1,1X,F7.5,A1,1X,ES15.5,A1,1X,ES15.5)')L,',',T(L),',',POPS(L1,L),',',&
                POPS(L2,L)
           WRITE(MOL_FILE_NUM,'(ES15.5,A1,1X,F10.3,A1,1X,F15.5)')MOL_RXN_COEF(J,1),',',MOL_RXN_COEF(J,2),&
                ',',MOL_RXN_COEF(J,3)
           WRITE(MOL_FILE_NUM,'(ES15.5,A1,1X,ES15.5,A1,1X,ES15.5)')RATE_K,',',FRD_R,',',DRATEDT
           WRITE(MOL_FILE_NUM,'(I4,A1,1X,I4,A1,1X,I4,A1,1X,I4)')L1,',',L2,',',L3,',',L4
           WRITE(MOL_FILE_NUM,'(ES11.4,A1,1X,ES11.4,A1,1X,ES11.4,A1,1X,ES11.4)') ALPHA_1,',',ALPHA_2,&
                ',',ALPHA_3,',',ALPHA_4
           WRITE(MOL_FILE_NUM,'(A,4I5)') 'ION IDS: ',ID1,ID2,ID3,ID4
           WRITE(MOL_FILE_NUM,'(A,4I5)') 'LEV_IN_IONS: ',LEV_IN_ION_MOL(J,:)
           WRITE(MOL_FILE_NUM,'(A,4I6)') 'E Values: ',E1,E2,E3,E4
           WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L1 Links: ', L11,L12,L13,L14
           WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L2 Links: ', L21,L22,L23,L24
           WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L3 Links: ', L31,L32,L33,L34
           WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L4 Links: ', L41,L42,L43,L44
           IF ( J .EQ. 1 .OR. J .EQ. 3) THEN
              WRITE(MOL_FILE_NUM,*) 'ION_LEV_TO_EQ_PNT (ID2): '
              WRITE(MOL_FILE_NUM,*) 'SIZE: ', SIZE(SE(ID2)%ION_LEV_TO_EQ_PNT)
              WRITE(MOL_FILE_NUM,*) SE(ID2)%ION_LEV_TO_EQ_PNT
           END IF
           WRITE(MOL_FILE_NUM,*)
        END IF
!
                END IF
!
                IF (PRINT_RXN_INFO) THEN
                   MOL_RXN_RATES_ALL(L,J) = FRD_R
                   IF (ABS(FRD_R) .GT. ABS(MOL_RXN_RATES_SORTED(L,12,1))) THEN
                     CHECKLOOP: DO CI1=1,12
                        IF (ABS(FRD_R) .GT. ABS(MOL_RXN_RATES_SORTED(L,CI1,1))) THEN
                           CI2 = CI1
                           EXIT CHECKLOOP
                        END IF
                     END DO CHECKLOOP
                     IF (CI2 .EQ. 12) THEN
                        MOL_RXN_RATES_SORTED(L,12,1) = FRD_R
                        MOL_RXN_RATES_SORTED(L,12,2) = J
                     ELSE
                        MOL_RXN_RATES_SORTED(L,CI2+1:12,:) = MOL_RXN_RATES_SORTED(L,CI2:11,:)
                        MOL_RXN_RATES_SORTED(L,CI2,1) = FRD_R
                        MOL_RXN_RATES_SORTED(L,CI2,2) = J
                     END IF
                   END IF
                END IF
!
             END DO
!
          ELSE IF(L2 .EQ. -1) THEN !Indicates a reaction involving an electron as a REACTANT
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,1)) THEN
                ID1=ID_ION_MOL(J,1)
             ELSE
                ID1=ID_ION_MOL(J,1)-1
             END IF
             L11=SE(ID1)%LNK_TO_IV(L1)
             L12=NT-1
             L13=SE(ID1)%LNK_TO_IV(L3)
             L14=SE(ID1)%LNK_TO_IV(L4)
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,3)) THEN
                ID3=ID_ION_MOL(J,3)
             ELSE
                ID3=ID_ION_MOL(J,3)-1
             END IF
             L31=SE(ID3)%LNK_TO_IV(L1)
             L32=NT-1
             L33=SE(ID3)%LNK_TO_IV(L3)
             L34=SE(ID3)%LNK_TO_IV(L4)
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,4)) THEN
                ID4=ID_ION_MOL(J,4)
             ELSE
                ID4=ID_ION_MOL(J,4)-1
             END IF
             L41=SE(ID4)%LNK_TO_IV(L1)
             L42=NT-1
             L43=SE(ID4)%LNK_TO_IV(L3)
             L44=SE(ID4)%LNK_TO_IV(L4)
!
! Now associate the relevant statistical equilibrium equations, using the links from STEQ_DATA_MOD
!
             IF (MOL_RXN_IS_FINAL_ION(J,1)) THEN
                E1=SE(ID1)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,1))
             ELSE
                E1 = LEV_IN_ION_MOL(J,1)
             END IF
!
             IF (MOL_RXN_IS_FINAL_ION(J,3)) THEN
                E3=SE(ID3)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,3))
             ELSE
                E3 = LEV_IN_ION_MOL(J,3)
             END IF
!
             IF (MOL_RXN_IS_FINAL_ION(J,4)) THEN
                E4=SE(ID4)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,4))
             ELSE
                E4 = LEV_IN_ION_MOL(J,4)
             END IF
!
             IF (E1*E3*E4 .EQ. 0) THEN
                LU_ER=ERROR_LU()
                WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: '
                WRITE(LU_ER,*) 'At least one of E1...E4 is zero--check ion levels'
                WRITE(LU_ER,*) 'This is an electron reaction: E2 is excluded'
                WRITE(LU_ER,*) E1,E3,E4
                STOP
             END IF
!
             II = DIAG_INDX
!
             DO L=1,ND
                TVAL=MAX(T(L),MOL_RXN_TLO(J))
                TVAL=MIN(TVAL,MOL_RXN_THI(J))
                TVAL4 = TEN4*TVAL
!
! Calculate rate constants and their variation with temperature
!
                IF((MOL_RXN_TYPE(J) .EQ. 1) .OR. (MOL_RXN_TYPE(J) .EQ. 2)) THEN
                   TREAL1 = EXP(-MOL_RXN_COEF(J,3)/(TVAL4))
                   TREAL2 = (TVAL4/300.)**MOL_RXN_COEF(J,2)
                   RATE_K = MOL_RXN_COEF(J,1)*TREAL1*TREAL2
                   DLN_K_DLNT = MOL_RXN_COEF(J,2) + (MOL_RXN_COEF(J,3)/(TVAL4))
                ELSE
                   LU_ER=ERROR_LU()
                   WRITE(LU_ER,*)'Error in STEQ_BA_MOL_RXN: '
                   WRITE(LU_ER,*)'Invalid reaction type passed,'
                   WRITE(LU_ER,*)'Reaction type: ',MOL_RXN_TYPE(J)
                   STOP
                END IF
!
                SCALED_K = RXN_SCALE_FAC*MOL_RXN_RATE_SCALE(J)*RATE_K
                IF (MOL_RXN_SPECIESID(J,3) .EQ. 'COMI' .OR. MOL_RXN_SPECIESID(J,4) .EQ. 'COMI') THEN
                   SCALED_K = SCALED_K*CO_SCALE_FAC
                END IF
!
                IF((TVAL .GT. MOL_RXN_THI(J)) .OR. (TVAL .LT. MOL_RXN_TLO(J))) DLN_K_DLNT = 0.0D0
!
                FRD_R=SCALED_K*POPS(L1,L)*ED(L) 
!
! The forward reaction rate, in my case the only rate considered
!
! Update statistical equilibrium
! Electrons are updated separately -- they don't need to be adjusted here
!
                SE(ID1)%STEQ(E1,L)=SE(ID1)%STEQ(E1,L) - FRD_R
!
                SE(ID3)%STEQ(E3,L)=SE(ID3)%STEQ(E3,L) + FRD_R
                SE(ID4)%STEQ(E4,L)=SE(ID4)%STEQ(E4,L) + FRD_R                
!
                MTOT_STEQ(LNK_FULL_TO_MOL(L1),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L1),L) - FRD_R
                MTOT_STEQ(LNK_FULL_TO_MOL(L3),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L3),L) + FRD_R
                MTOT_STEQ(LNK_FULL_TO_MOL(L4),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L4),L) + FRD_R
!
                IF (UPDATE_BA) THEN
                   DRATEDT = (FRD_R*DLN_K_DLNT)/TVAL
                   ALPHA_1 = SCALED_K*POPS(L1,L)
                   ALPHA_2 = SCALED_K*ED(L)
                   ALPHA_3 = SCALED_K*POPS(L3,L)
                   ALPHA_4 = SCALED_K*POPS(L4,L)
!
                   IF (L3 .EQ. L4) THEN
!
                      ID=ID1
                      NIV=SE(ID)%N_IV
                      IF(L11 .NE. 0) SE(ID)%BA(E1,L11,II,L)=SE(ID)%BA(E1,L11,II,L) - ALPHA_2
                      IF(L12 .NE. 0) SE(ID)%BA(E1,L12,II,L)=SE(ID)%BA(E1,L12,II,L) - ALPHA_1
                      SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) - DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) = MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) - &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT-1),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT-1),L) - &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) = MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) - &
DRATEDT
!
                      ID=ID3
                      NIV=SE(ID)%N_IV
                      IF(L31 .NE. 0) SE(ID)%BA(E3,L31,II,L)=SE(ID)%BA(E3,L31,II,L) + 2*ALPHA_2
                      IF(L32 .NE. 0) SE(ID)%BA(E3,L32,II,L)=SE(ID)%BA(E3,L32,II,L) + 2*ALPHA_1
                      SE(ID)%BA(E3,NIV,II,L)=SE(ID)%BA(E3,NIV,II,L) + 2*DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) = MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) + &
2*ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT-1),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT-1),L) + &
2*ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) = MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) + &
2*DRATEDT
!
                   ELSE
!
                      ID=ID1
                      NIV=SE(ID)%N_IV
                      IF(L11 .NE. 0) SE(ID)%BA(E1,L11,II,L)=SE(ID)%BA(E1,L11,II,L) - ALPHA_2
                      IF(L12 .NE. 0) SE(ID)%BA(E1,L12,II,L)=SE(ID)%BA(E1,L12,II,L) - ALPHA_1
                      SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) - DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) = MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) - &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT-1),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT-1),L) - &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) = MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) - &
DRATEDT
!
                      ID=ID3
                      NIV=SE(ID)%N_IV
                      IF(L31 .NE. 0) SE(ID)%BA(E3,L31,II,L)=SE(ID)%BA(E3,L31,II,L) + ALPHA_2
                      IF(L32 .NE. 0) SE(ID)%BA(E3,L32,II,L)=SE(ID)%BA(E3,L32,II,L) + ALPHA_1
                      SE(ID)%BA(E3,NIV,II,L)=SE(ID)%BA(E3,NIV,II,L) + DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) = MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) + &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT-1),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT-1),L) + &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) = MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) + &
DRATEDT
!
                      ID=ID4
                      NIV=SE(ID)%N_IV
                      IF(L41 .NE. 0) SE(ID)%BA(E4,L41,II,L)=SE(ID)%BA(E4,L41,II,L) + ALPHA_2
                      IF(L42 .NE. 0) SE(ID)%BA(E4,L42,II,L)=SE(ID)%BA(E4,L42,II,L) + ALPHA_1
                      SE(ID)%BA(E4,NIV,II,L)=SE(ID)%BA(E4,NIV,II,L) + DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(L1),L) = MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(L1),L) + &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(NT-1),L)=MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(NT-1),L) + &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(NT),L) = MTOT_BA(LNK_FULL_TO_MOL(L4),LNK_FULL_TO_MOL(NT),L) + &
DRATEDT
!
                   END IF
                   
                !    Add some print statements for debugging, use variable PRINT_RXN_INFO to toggle
                   IF (PRINT_DEBUG_INFO) THEN
                   WRITE(MOL_FILE_NUM,'(A16,I2,A1,2X,A8,1X,A1,1X,A8,1X,A2,1X,A8,1X,A1,1X,A8)')'Reaction Number ',&
J,':',MOL_RXN_SPECIESID(J,1),'+',MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4)
                   WRITE(MOL_FILE_NUM,'(I3,A1,1X,F7.5,A1,1X,ES15.5,A1,1X,ES15.5)')L,',',T(L),',',POPS(L1,L),',',&
ED(L)
                   WRITE(MOL_FILE_NUM,'(ES15.5,A1,1X,F10.3,A1,1X,F15.5)')MOL_RXN_COEF(J,1),',',MOL_RXN_COEF(J,2),&
',',MOL_RXN_COEF(J,3)
                   WRITE(MOL_FILE_NUM,'(ES15.5,A1,1X,ES15.5,A1,1X,ES15.5)')RATE_K,',',FRD_R,',',DRATEDT
                   WRITE(MOL_FILE_NUM,'(I4,A1,1X,I4,A1,1X,I4,A1,1X,I4)')L1,',',L2,',',L3,',',L4
                   WRITE(MOL_FILE_NUM,'(ES11.4,A1,1X,ES11.4,A1,1X,ES11.4,A1,1X,ES11.4)') ALPHA_1,',',ALPHA_2,&
',',ALPHA_3,',',ALPHA_4
                   WRITE(MOL_FILE_NUM,'(A,4I5)') 'ION IDS: ',ID1,ID2,ID3,ID4
                   WRITE(MOL_FILE_NUM,'(A,4I5)') 'LEV_IN_IONS: ',LEV_IN_ION_MOL(J,:)
                   WRITE(MOL_FILE_NUM,'(A,4I6)') 'E Values: ',E1,-1,E3,E4
                   WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L1 Links: ', L11,L12,L13,L14
                   WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L3 Links: ', L31,L32,L33,L34
                   WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L4 Links: ', L41,L42,L43,L44
                   WRITE(MOL_FILE_NUM,*)
                   END IF
                END IF
!
                IF (PRINT_RXN_INFO) THEN
                   MOL_RXN_RATES_ALL(L,J) = FRD_R
                   IF (ABS(FRD_R) .GT. ABS(MOL_RXN_RATES_SORTED(L,12,1))) THEN
                     CHECKLOOP1: DO CI1=1,12
                        IF (ABS(FRD_R) .GT. ABS(MOL_RXN_RATES_SORTED(L,CI1,1))) THEN
                           CI2 = CI1
                           EXIT CHECKLOOP1
                        END IF
                     END DO CHECKLOOP1
                     IF (CI2 .EQ. 12) THEN
                        MOL_RXN_RATES_SORTED(L,12,1) = FRD_R
                        MOL_RXN_RATES_SORTED(L,12,2) = J
                     ELSE
                        MOL_RXN_RATES_SORTED(L,CI2+1:12,:) = MOL_RXN_RATES_SORTED(L,CI2:11,:)
                        MOL_RXN_RATES_SORTED(L,CI2,1) = FRD_R
                        MOL_RXN_RATES_SORTED(L,CI2,2) = J
                     END IF
                  END IF
               END IF
            END DO
!
          ELSE IF (L4 .EQ. -2) THEN !A reaction with NOATOM listed for the 4th species, these are radiative 
             !association reactions.  We simply ignore the 4th spot in the reaction data
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,1)) THEN
                ID1=ID_ION_MOL(J,1)
             ELSE
                ID1=ID_ION_MOL(J,1)-1
             END IF
             L11=SE(ID1)%LNK_TO_IV(L1)
             L12=SE(ID1)%LNK_TO_IV(L2)
             L13=SE(ID1)%LNK_TO_IV(L3)
!             L14=SE(ID1)%LNK_TO_IV(L4)
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,2)) THEN
                ID2=ID_ION_MOL(J,2)
             ELSE
                ID2=ID_ION_MOL(J,2)-1
             END IF
             L21=SE(ID2)%LNK_TO_IV(L1)
             L22=SE(ID2)%LNK_TO_IV(L2)
             L23=SE(ID2)%LNK_TO_IV(L3)
!             L24=SE(ID2)%LNK_TO_IV(L4)
!
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,3)) THEN
                ID3=ID_ION_MOL(J,3)
             ELSE
                ID3=ID_ION_MOL(J,3)-1
             END IF
             L31=SE(ID3)%LNK_TO_IV(L1)
             L32=SE(ID3)%LNK_TO_IV(L2)
             L33=SE(ID3)%LNK_TO_IV(L3)
!             L34=SE(ID3)%LNK_TO_IV(L4)
!
! Now associate the relevant statistical equilibrium equations, using the links from STEQ_DATA_MOD
!
             IF (MOL_RXN_IS_FINAL_ION(J,1)) THEN
                E1=SE(ID1)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,1))
             ELSE
                E1 = LEV_IN_ION_MOL(J,1)
             END IF
!
             IF (MOL_RXN_IS_FINAL_ION(J,3)) THEN
                E3=SE(ID3)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,3))
             ELSE
                E3 = LEV_IN_ION_MOL(J,3)
             END IF
!
             IF (MOL_RXN_IS_FINAL_ION(J,2)) THEN
                E2=SE(ID2)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,2))
             ELSE
                E2 = LEV_IN_ION_MOL(J,2)
             END IF
!
             IF (E1*E2*E3 .EQ. 0) THEN
                LU_ER=ERROR_LU()
                WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: '
                WRITE(LU_ER,*) 'At least one of E1...E3 is zero--check ion levels'
                WRITE(LU_ER,*) 'This is a radiative association reaction, there is no E4'
                WRITE(LU_ER,*) E1,E2,E3
                STOP
             END IF
!
             II = DIAG_INDX
!
             DO L=1,ND
                TVAL=MAX(T(L),MOL_RXN_TLO(J))
                TVAL=MIN(TVAL,MOL_RXN_THI(J))
                TVAL4 = TEN4*TVAL
!
! Calculate rate constants and their variation with temperature
!
                IF((MOL_RXN_TYPE(J) .EQ. 1) .OR. (MOL_RXN_TYPE(J) .EQ. 2)) THEN
                   TREAL1 = EXP(-MOL_RXN_COEF(J,3)/(TVAL4))
                   TREAL2 = MOL_RXN_COEF(J,1)*((TVAL4/300.)**MOL_RXN_COEF(J,2))
                   RATE_K = TREAL1*TREAL2
                   DLN_K_DLNT = MOL_RXN_COEF(J,2) + (MOL_RXN_COEF(J,3)/(TVAL4))
                ELSE
                   LU_ER=ERROR_LU()
                   WRITE(LU_ER,*)'Error in STEQ_BA_MOL_RXN: '
                   WRITE(LU_ER,*)'Invalid reaction type passed,'
                   WRITE(LU_ER,*)'Reaction type: ',MOL_RXN_TYPE(J)
                   STOP
                END IF
!
                SCALED_K = RXN_SCALE_FAC*MOL_RXN_RATE_SCALE(J)*RATE_K
                IF (MOL_RXN_SPECIESID(J,3) .EQ. 'COMI' .OR. MOL_RXN_SPECIESID(J,4) .EQ. 'COMI') THEN
                   SCALED_K = SCALED_K*CO_SCALE_FAC
                END IF
!
                IF((TVAL .GT. MOL_RXN_THI(J)) .OR. (TVAL .LT. MOL_RXN_TLO(J))) DLN_K_DLNT = 0.0D0
!
                FRD_R=SCALED_K*POPS(L1,L)*POPS(L2,L) 
!
! The forward reaction rate, in my case the only rate considered
!
! Update statistical equilibrium
!
                SE(ID1)%STEQ(E1,L)=SE(ID1)%STEQ(E1,L) - FRD_R
                SE(ID2)%STEQ(E2,L)=SE(ID2)%STEQ(E2,L) - FRD_R
!
                SE(ID3)%STEQ(E3,L)=SE(ID3)%STEQ(E3,L) + FRD_R
!
                MTOT_STEQ(LNK_FULL_TO_MOL(L1),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L1),L) - FRD_R
                MTOT_STEQ(LNK_FULL_TO_MOL(L2),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L2),L) - FRD_R
                MTOT_STEQ(LNK_FULL_TO_MOL(L3),L)=MTOT_STEQ(LNK_FULL_TO_MOL(L3),L) + FRD_R
!
                IF (UPDATE_BA) THEN
                   DRATEDT = (FRD_R*DLN_K_DLNT)/TVAL
                   ALPHA_1 = SCALED_K*POPS(L1,L)
                   ALPHA_2 = SCALED_K*POPS(L2,L)
                   ALPHA_3 = SCALED_K*POPS(L3,L)
!
                   IF (L1 .EQ. L2) THEN
!                      
                      ID=ID1
                      NIV=SE(ID)%N_IV
                      IF(L11 .NE. 0) SE(ID)%BA(E1,L11,II,L)=SE(ID)%BA(E1,L11,II,L) - 4*ALPHA_1
                      SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) - 2*DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) = MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) - &
4*ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) = MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) - &
2*DRATEDT
!
                      ID=ID3
                      NIV=SE(ID)%N_IV
                      IF(L31 .NE. 0) SE(ID)%BA(E3,L31,II,L)=SE(ID)%BA(E3,L31,II,L) + 2*ALPHA_1
                      SE(ID)%BA(E3,NIV,II,L)=SE(ID)%BA(E3,NIV,II,L) + DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) = MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) + &
2*ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) = MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) + &
DRATEDT
!
                   ELSE
!
                      ID=ID1
                      NIV=SE(ID)%N_IV
                      IF(L11 .NE. 0) SE(ID)%BA(E1,L11,II,L)=SE(ID)%BA(E1,L11,II,L) - ALPHA_2
                      IF(L12 .NE. 0) SE(ID)%BA(E1,L12,II,L)=SE(ID)%BA(E1,L12,II,L) - ALPHA_1
                      SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) - DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L1),L) - &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(L2),L) - &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L1),LNK_FULL_TO_MOL(NT),L) - &
DRATEDT
!
                      ID=ID2
                      NIV=SE(ID)%N_IV
                      IF(L21 .NE. 0) SE(ID)%BA(E2,L21,II,L)=SE(ID)%BA(E2,L21,II,L) - ALPHA_2
                      IF(L22 .NE. 0) SE(ID)%BA(E2,L22,II,L)=SE(ID)%BA(E2,L22,II,L) - ALPHA_1
                      SE(ID)%BA(E2,NIV,II,L)=SE(ID)%BA(E2,NIV,II,L) - DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L1),L) - &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(L2),L) - &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L2),LNK_FULL_TO_MOL(NT),L) - &
DRATEDT
!
                      ID=ID3
                      NIV=SE(ID)%N_IV
                      IF(L31 .NE. 0) SE(ID)%BA(E3,L31,II,L)=SE(ID)%BA(E3,L31,II,L) + ALPHA_2
                      IF(L32 .NE. 0) SE(ID)%BA(E3,L32,II,L)=SE(ID)%BA(E3,L32,II,L) + ALPHA_1
                      SE(ID)%BA(E3,NIV,II,L)=SE(ID)%BA(E3,NIV,II,L) + DRATEDT
!
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L1),L) + &
ALPHA_2
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L2),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(L2),L) + &
ALPHA_1
                      MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L)=MTOT_BA(LNK_FULL_TO_MOL(L3),LNK_FULL_TO_MOL(NT),L) + &
DRATEDT
!
                   END IF
!
                END IF
!
                IF (PRINT_RXN_INFO) THEN
                   MOL_RXN_RATES_ALL(L,J) = FRD_R
                   IF (ABS(FRD_R) .GT. ABS(MOL_RXN_RATES_SORTED(L,12,1))) THEN
                     CHECKLOOP2: DO CI1=1,12
                        IF (ABS(FRD_R) .GT. ABS(MOL_RXN_RATES_SORTED(L,CI1,1))) THEN
                           CI2 = CI1
                           EXIT CHECKLOOP2
                        END IF
                     END DO CHECKLOOP2
                     IF (CI2 .EQ. 12) THEN
                        MOL_RXN_RATES_SORTED(L,12,1) = FRD_R
                        MOL_RXN_RATES_SORTED(L,12,2) = J
                     ELSE
                        MOL_RXN_RATES_SORTED(L,CI2+1:12,:) = MOL_RXN_RATES_SORTED(L,CI2:11,:)
                        MOL_RXN_RATES_SORTED(L,CI2,1) = FRD_R
                        MOL_RXN_RATES_SORTED(L,CI2,2) = J
                     END IF
                   END IF
                !Add some print statements for debugging, use variable PRINT_RXN_INFO to toggle
                   IF (PRINT_DEBUG_INFO) THEN
                   WRITE(MOL_FILE_NUM,'(A16,I2,A1,2X,A8,1X,A1,1X,A8,1X,A2,1X,A8,1X,A1,1X,A8)')'Reaction Number ',&
J,':',MOL_RXN_SPECIESID(J,1),'+',MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4)
                   WRITE(MOL_FILE_NUM,'(I3,A1,1X,F7.5,A1,1X,ES15.5,A1,1X,ES15.5)')L,',',T(L),',',POPS(L1,L),',',&
POPS(L2,L)
                   WRITE(MOL_FILE_NUM,'(ES15.5,A1,1X,F10.3,A1,1X,F15.5)')MOL_RXN_COEF(J,1),',',MOL_RXN_COEF(J,2),&
',',MOL_RXN_COEF(J,3)
                   WRITE(MOL_FILE_NUM,'(ES15.5,A1,1X,ES15.5,A1,1X,ES15.5)')RATE_K,',',FRD_R,',',DRATEDT
                   WRITE(MOL_FILE_NUM,'(I4,A1,1X,I4,A1,1X,I4,A1,1X,I4)')L1,',',L2,',',L3,',',L4
                   WRITE(MOL_FILE_NUM,'(ES11.4,A1,1X,ES11.4,A1,1X,ES11.4)') ALPHA_1,',',ALPHA_2,&
',',ALPHA_3
                   WRITE(MOL_FILE_NUM,'(A,3I5)') 'ION IDS: ',ID1,ID2,ID3
                   WRITE(MOL_FILE_NUM,'(A,4I5)') 'LEV_IN_IONS: ',LEV_IN_ION_MOL(J,:)
                   WRITE(MOL_FILE_NUM,'(A,3I6)') 'E Values: ',E1,E2,E3
                   WRITE(MOL_FILE_NUM,'(A,5X,3I5)')'L1 Links: ', L11,L12,L13
                   WRITE(MOL_FILE_NUM,'(A,5X,3I5)')'L2 Links: ', L21,L22,L23
                   WRITE(MOL_FILE_NUM,'(A,5X,3I5)')'L3 Links: ', L31,L32,L33
                   WRITE(MOL_FILE_NUM,*)
!
                END IF
             END IF
!
          END DO
!
       END IF
!
    END IF
!
 END DO
!
 INITIALIZE_ARRAYS_MOL=.TRUE. !So that set_mol_rxn will reset the arrays when next called
!
    INQUIRE(UNIT=MOL_FILE_NUM,OPENED=IS_OPEN)
    IF (IS_OPEN) THEN
       CALL FLUSH(MOL_FILE_NUM)
       CLOSE(MOL_FILE_NUM)
    END IF
!
    IF (PRINT_RXN_INFO) THEN
       !Write a file with all the molecular reaction rates, sorted by depth       
       OPEN(UNIT=MOL_FILE_NUM+1,FILE=MOL_OUT_FILE_UNSORT,STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
       WRITE(MOL_FILE_NUM+1,*)'************************************************************************************'
       WRITE(MOL_FILE_NUM+1,*)'Info for molecular reactions, final iteration'
       WRITE(MOL_FILE_NUM+1,*)'Shown: total reaction rates (in cm^-3 s^-1), 10 depths at a time'
       WRITE(MOL_FILE_NUM+1,*)'************************************************************************************'
       WRITE(MOL_FILE_NUM+1,*)
       L=1
       WRITELOOP1: DO 
          IF (L .LE. ND-9) THEN
             WRITE(MOL_FILE_NUM+1,'(A7,I3,A1,I3,1X,A13,T30,10ES15.5)')'Depths ',L,'-',L+9,'T (10^4 K) = ',T(L:L+9)
             WRITE(MOL_FILE_NUM+1,*)
             DO J=1,N_MOL_RXN
                WRITE(MOL_FILE_NUM+1,'(A5,A1,1X,A5,A2,1X,A5,A1,1X,A5,A1,T30,10ES15.5)')MOL_RXN_SPECIESID(J,1),'+',&
MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4),':',MOL_RXN_RATES_ALL(L:L+9,J)
             END DO
             WRITE(MOL_FILE_NUM+1,*)
             L=L+10
          ELSE
             LDIFF=ND-L
             WRITE(MOL_FILE_NUM+1,'(A7,I3,A1,I3,1X,A13,T30,10ES15.5)')'Depths ',L,'-',ND,'T (10^4 K) = ',T(L:ND)
             DO J=1,N_MOL_RXN
                WRITE(MOL_FILE_NUM+1,'(A5,A1,1X,A5,A2,1X,A5,A1,1X,A5,A1,T30,10ES15.5)')MOL_RXN_SPECIESID(J,1),'+',&
MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4),':',MOL_RXN_RATES_ALL(L:ND,J)
             END DO
             EXIT WRITELOOP1
          END IF
       END DO WRITELOOP1
       CALL FLUSH(MOL_FILE_NUM+1)
       CLOSE(MOL_FILE_NUM+1)
       !Write a file which records the largest reaction rates at each depth, sorted by rate
       OPEN(UNIT=MOL_FILE_NUM+2,FILE=MOL_OUT_FILE_SORT,STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
       WRITE(MOL_FILE_NUM+2,*)'************************************************************************************'
       WRITE(MOL_FILE_NUM+2,*)'Info for molecular reactions, final iteration'
       WRITE(MOL_FILE_NUM+2,*)'Shown: total reaction rates (in cm^-3 s^-1), 12 largest reactions at each depth'
       WRITE(MOL_FILE_NUM+2,*)'************************************************************************************'
       WRITE(MOL_FILE_NUM+2,*)
       DO L=1,ND
          WRITE(MOL_FILE_NUM+2,'(A,I3,1X,A,ES15.5)')'Depth: ',L,'T (10^4 K): ',T(L)
          WRITE(MOL_FILE_NUM+2,*)
          DO CI1=1,12
             J=MOL_RXN_RATES_SORTED(L,CI1,2)
             WRITE(MOL_FILE_NUM+2,'(A5,A1,1X,A5,A2,1X,A5,A1,1X,A6,A1,5X,ES15.5)')MOL_RXN_SPECIESID(J,1),'+',&
MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4),':', MOL_RXN_RATES_SORTED(L,CI1,1)
          END DO
          WRITE(MOL_FILE_NUM+2,*)
       END DO
       CALL FLUSH(MOL_FILE_NUM+2)
       CLOSE(MOL_FILE_NUM+2)
    END IF
!
    !Another section for printing out information on a per-species basis
!
    IF (PRINT_RXN_INFO) THEN
!
!       ALLOCATE(TOP_RXNS_BY_SPECIES(ND,14,2,5),STAT=IOS)
!       IF (IOS .EQ. 0) ALLOCATE(TOTAL_CHANGE_BY_SPECIES(ND,14),STAT=IOS)
!       IF (IOS .EQ. 0) ALLOCATE(REACTANTS(2),STAT=IOS)
!       IF (IOS .EQ. 0) ALLOCATE(PRODUCTS(2),STAT=IOS)
!       IF (IOS .NE. 0) WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: Unable to allocate sorting arrays for printing info'
!!
!       TOP_RXNS_BY_SPECIES = N_MOL_RXN+1
!       TOTAL_CHANGE_BY_SPECIES = 0.0D0
!!
!       DO L=1,ND
!          DO J=1,N_MOL_RXN
!             REACTANTS(1) = MOL_RXN_SPECIESID(J,1)
!             REACTANTS(2) = MOL_RXN_SPECIESID(J,2)
!             PRODUCTS(1) = MOL_RXN_SPECIESID(J,3)
!             PRODUCTS(2) = MOL_RXN_SPECIESID(J,4)
!!
!             DO SPECIESNUM=1,12
!                SPECIESID1 = SPECIESIDLIST(SPECIESNUM)
!                IF (SPECIESID1 .EQ. REACTANTS(1) .OR. SPECIESID1 .EQ. REACTANTS(2)) THEN
!                   DO ID=1,5
!                      IF (ABS(MOL_RXN_RATES_ALL(L,J)) .GT. ABS(MOL_RXN_RATES_ALL(L,TOP_RXNS_BY_SPECIES(L,SPECIESNUM,2,ID)))) THEN
!                         IF (ID .LT. 5) TOP_RXNS_BY_SPECIES(L,SPECIESNUM,2,ID+1:5) = TOP_RXNS_BY_SPECIES(L,SPECIESNUM,2,ID:4)
!                         TOP_RXNS_BY_SPECIES(L,SPECIESNUM,2,ID) = J
!                         EXIT
!                      END IF
!                   END DO
!                   IF (SPECIESID1 .EQ. REACTANTS(1) .AND. SPECIESID1 .EQ. REACTANTS(2)) THEN
!                      TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM) = TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM)-2.0*MOL_RXN_RATES_ALL(L,J)
!                   ELSE
!                      TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM) = TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM)-1.0*MOL_RXN_RATES_ALL(L,J)
!                   END IF
!                END IF
!!
!                IF (SPECIESID1 .EQ. PRODUCTS(1) .OR. SPECIESID1 .EQ. PRODUCTS(2)) THEN
!                   DO ID=1,5
!                      IF (ABS(MOL_RXN_RATES_ALL(L,J)) .GT. ABS(MOL_RXN_RATES_ALL(L,TOP_RXNS_BY_SPECIES(L,SPECIESNUM,1,ID)))) THEN
!                         IF (ID .LT. 5) TOP_RXNS_BY_SPECIES(L,SPECIESNUM,1,ID+1:5) = TOP_RXNS_BY_SPECIES(L,SPECIESNUM,1,ID:4)
!                         TOP_RXNS_BY_SPECIES(L,SPECIESNUM,1,ID) = J
!                         EXIT
!                      END IF
!                   END DO
!                   IF (SPECIESID1 .EQ. PRODUCTS(1) .AND. SPECIESID1 .EQ. PRODUCTS(2)) THEN
!                      TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM) = TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM)+2.0*MOL_RXN_RATES_ALL(L,J)
!                   ELSE
!                      TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM) = TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM)+1.0*MOL_RXN_RATES_ALL(L,J)
!                   END IF
!                END IF
!!
!             END DO
!!
!          END DO
!       END DO
!
! I have moved code for writing the file to a separate subroutine: print_mol_info.f90
!
       CALL PRINT_MOL_INFO(SPECIESIDLIST,MOL_RXN_RATES_ALL,POPS,MOL_LEVID_LIST,MOL_ATMID_LIST,ND,NT)

!       OPEN(UNIT=MOL_FILE_NUM+3,FILE=MOL_RXNS_SPECIES_FILE,STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)       
!       ID1 = MOL_FILE_NUM+3
!       WRITE(ID1,*)'****************************************************************************************************'
!       WRITE(ID1,*)'Molecular reaction info for each species, final iteration'
!       WRITE(ID1,*)'For each species, shown are top destruction reactions (on the left) and top formation reactions (on the right)'
!       WRITE(ID1,*)'****************************************************************************************************'
!       WRITE(ID1,*)
!
!       DO L=1,ND
!          WRITE(ID1,*)'*************************************************************************************************'
!          WRITE(ID1,'(A,I3,2X,A,F10.5,2X,A,ES15.5)')'Depth:',L,'Temperature (10^4 K): ',T(L),'ED:',ED(L)
!          DO SPECIESNUM=1,12
!             WRITE(ID1,*)'-----------------------------------------------------------------------------'
!             WRITE(ID1,*)
!             WRITE(ID1,'(A,A)')TRIM(SPECIESIDLIST(SPECIESNUM)),':'
!             WRITE(ID1,'(ES15.5,5X,ES15.5)') ATM(MOL_ATMID_LIST(SPECIESNUM)%APRXzV(1,L),ATM(MOL_ATMID_LIST(SPECIESNUM)%ARRXzV(1,L)
!            WRITE(ID1,*)
!            DO J=1,5
!               ID2 = TOP_RXNS_BY_SPECIES(L,SPECIESNUM,2,J)
!               ID3 = TOP_RXNS_BY_SPECIES(L,SPECIESNUM,1,J)
!               IF (ID2 .LE. SIZE(MOL_RXN_SPECIESID,1) .AND. ID3 .LE. SIZE(MOL_RXN_SPECIESID,1)) THEN
!                  WRITE(ID1,'(A6,A1,1X,A6,A2,1X,A6,A1,1X,A6,A1,1X,ES15.5,T50,A6,A1,1X,A6,A2,1X,A6,A1,1X,A6,A1,1X,ES15.5)')&
! MOL_RXN_SPECIESID(ID2,1),'+',MOL_RXN_SPECIESID(ID2,2),'->',MOL_RXN_SPECIESID(ID2,3),'+',MOL_RXN_SPECIESID(ID2,4),&
!':', MOL_RXN_RATES_ALL(L,ID2),MOL_RXN_SPECIESID(ID3,1),'+',MOL_RXN_SPECIESID(ID3,2),'->',MOL_RXN_SPECIESID(ID3,3),&
!'+',MOL_RXN_SPECIESID(ID3,4),':', MOL_RXN_RATES_ALL(L,ID3)
!                ELSE IF (ID2 .LE. SIZE(MOL_RXN_SPECIESID,1)) THEN
!                   WRITE(ID1,'(A6,A1,1X,A6,A2,1X,A6,A1,1X,A6,A1,1X,ES15.5)')&
! MOL_RXN_SPECIESID(ID2,1),'+',MOL_RXN_SPECIESID(ID2,2),'->',MOL_RXN_SPECIESID(ID2,3),'+',MOL_RXN_SPECIESID(ID2,4),&
!':', MOL_RXN_RATES_ALL(L,ID2)
!                ELSE IF (ID3 .LE. SIZE(MOL_RXN_SPECIESID,1)) THEN
!                   WRITE(ID1,'(T50,A6,A1,1X,A6,A2,1X,A6,A1,1X,A6,A1,1X,ES15.5)')&
!MOL_RXN_SPECIESID(ID3,1),'+',MOL_RXN_SPECIESID(ID3,2),'->',MOL_RXN_SPECIESID(ID3,3),&
!'+',MOL_RXN_SPECIESID(ID3,4),':', MOL_RXN_RATES_ALL(L,ID3)
!                ELSE
!                   WRITE(ID1,*)
!                END IF
!             END DO
!             WRITE(ID1,*)
!             WRITE(ID1,'(A,A,A,2X,ES15.5,2X,A,2X,ES15.5)')'Total change to ',TRIM(SPECIESIDLIST(SPECIESNUM)),&
!'from molecular reactions:',TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM),&
!'Scaled by pop:',TOTAL_CHANGE_BY_SPECIES(L,SPECIESNUM)/POPS(MOL_LEVID_LIST(SPECIESNUM),L)
!             WRITE(ID1,*)'-----------------------------------------------------------------------------'
!          END DO
!          WRITE(ID1,*)'*************************************************************************************************'
!          WRITE(ID1,*)
!       END DO
!
!       CALL FLUSH(ID1)
!       CLOSE(ID1)
!
    END IF
!
    IF (WRITE_OUT_POPS) THEN
       OPEN(UNIT=11,FILE='CI_POPS_INFO',STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
       WRITE(11,*)'Neutral Carbon Populations'
       WRITE(11,*) 
       DO L=1,ND
          WRITE(11,'(ES17.7)') POPS(MOL_LEVID_LIST(1),L)
       END DO
       WRITE(11,*) "levid list: ",MOL_LEVID_LIST
       CLOSE(11)
!
       OPEN(UNIT=11,FILE='C2_POPS_INFO',STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
       WRITE(11,*)'Ionized Carbon Populations'
       WRITE(11,*) 
       DO L=1,ND
          WRITE(11,'(ES17.7)') POPS(MOL_LEVID_LIST(2),L)
       END DO
       CLOSE(11)
!
       OPEN(UNIT=11,FILE='OI_POPS_INFO',STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
       WRITE(11,*)'Neutral Oxygen Populations'
       WRITE(11,*) 
       DO L=1,ND
          WRITE(11,'(ES17.7)') POPS(MOL_LEVID_LIST(3),L)
       END DO
       CLOSE(11)
!
       OPEN(UNIT=11,FILE='O2_POPS_INFO',STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
       WRITE(11,*)'Ionized Oxygen Populations'
       WRITE(11,*) 
       DO L=1,ND
          WRITE(11,'(ES17.7)') POPS(MOL_LEVID_LIST(4),L)
       END DO
       CLOSE(11)
    END IF
!
    IF (UPDATE_BA .AND. PRINT_DEBUG_INFO) THEN
       OPEN(UNIT=12,FILE='MOL_RXN_BA_CHK_D5',STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
       OPEN(UNIT=13,FILE='MOL_RXN_BA_CHK_IND_D5',STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
       WRITE(12,*)'Jacobian check for molecular reactions, depth = 5'
       WRITE(12,*)
       WRITE(12,*)'*****************************************************************************'
       WRITE(12,*)
!
       WRITE(13,'(A,5X,I4,A1,I4,A1,I4,A1,I4)') 'CI_BA_SIZE: ',SIZE(SE(MOL_ATMID_LIST(1))%BA,1),',',&
SIZE(SE(MOL_ATMID_LIST(1))%BA,2),',',SIZE(SE(MOL_ATMID_LIST(1))%BA,3),',',SIZE(SE(MOL_ATMID_LIST(1))%BA,4)
       WRITE(13,'(A,5X,I4,A1,I4,A1,I4,A1,I4)') 'C2_BA_SIZE: ',SIZE(SE(MOL_ATMID_LIST(2))%BA,1),',',&
SIZE(SE(MOL_ATMID_LIST(2))%BA,2),',',SIZE(SE(MOL_ATMID_LIST(2))%BA,3),',',SIZE(SE(MOL_ATMID_LIST(2))%BA,4)
       WRITE(13,'(A,5X,I4,A1,I4,A1,I4,A1,I4)') 'OI_BA_SIZE: ',SIZE(SE(MOL_ATMID_LIST(3))%BA,1),',',&
SIZE(SE(MOL_ATMID_LIST(3))%BA,2),',',SIZE(SE(MOL_ATMID_LIST(3))%BA,3),',',SIZE(SE(MOL_ATMID_LIST(3))%BA,4)
       WRITE(13,'(A,5X,I4,A1,I4,A1,I4,A1,I4)') 'O2_BA_SIZE: ',SIZE(SE(MOL_ATMID_LIST(4))%BA,1),',',&
SIZE(SE(MOL_ATMID_LIST(4))%BA,2),',',SIZE(SE(MOL_ATMID_LIST(4))%BA,3),',',SIZE(SE(MOL_ATMID_LIST(4))%BA,4)
       WRITE(13,'(A,5X,I4,A1,I4,A1,I4,A1,I4)') 'COMI_BA_SIZE: ',SIZE(SE(MOL_ATMID_LIST(5))%BA,1),',',&
SIZE(SE(MOL_ATMID_LIST(5))%BA,2),',',SIZE(SE(MOL_ATMID_LIST(5))%BA,3),',',SIZE(SE(MOL_ATMID_LIST(5))%BA,4)
       WRITE(13,'(A,5X,I4,A1,I4,A1,I4,A1,I4)') 'C2MI_BA_SIZE: ',SIZE(SE(MOL_ATMID_LIST(7))%BA,1),',',&
SIZE(SE(MOL_ATMID_LIST(7))%BA,2),',',SIZE(SE(MOL_ATMID_LIST(7))%BA,3),',',SIZE(SE(MOL_ATMID_LIST(7))%BA,4)
       WRITE(13,'(A,5X,I4,A1,I4,A1,I4,A1,I4)') 'O2MI_BA_SIZE: ',SIZE(SE(MOL_ATMID_LIST(9))%BA,1),',',&
SIZE(SE(MOL_ATMID_LIST(9))%BA,2),',',SIZE(SE(MOL_ATMID_LIST(9))%BA,3),',',SIZE(SE(MOL_ATMID_LIST(9))%BA,4)
       WRITE(13,'(A,5X,I4,A1,I4,A1,I4,A1,I4)') 'CO2MI_BA_SIZE: ',SIZE(SE(MOL_ATMID_LIST(11))%BA,1),',',&
SIZE(SE(MOL_ATMID_LIST(11))%BA,2),',',SIZE(SE(MOL_ATMID_LIST(11))%BA,3),',',SIZE(SE(MOL_ATMID_LIST(11))%BA,4)
!
       WRITE(13,*)
       WRITE(13,*)
!
       L=1
       LOOP2: DO
          IF (L .GT. NT) EXIT LOOP2
          IF (L+9 .LE. NT) THEN
             WRITE(12,'(T12,I3,T24,I3,T36,I3,T48,I3,T60,I3,T72,I3,T84,I3,T96,I3,T108,I3,T120,I3)') L,L+1,L+2,L+3,L+4,L+5,&
L+6,L+7,L+8,L+9
             WRITE(13,'(T12,I3,T24,I3,T36,I3,T48,I3,T60,I3,T72,I3,T84,I3,T96,I3,T108,I3,T120,I3)') L,L+1,L+2,L+3,L+4,L+5,&
L+6,L+7,L+8,L+9
!          WRITE(12,*)
!
             DO K=1,4
                WRITE(13,'(I4,A1,5X,10ES12.3)') MOL_LEVID_LIST(K),'*',SE(MOL_ATMID_LIST(K))%BA(1,L:L+9,II,5)
             END DO
!
             K=5
             DO
                WRITE(13,'(I4,A1,5X,10ES12.3)') MOL_LEVID_LIST(K),'*',SE(MOL_ATMID_LIST(K))%BA(1,L:L+9,II,5)
                WRITE(13,'(I4,A1,5X,10ES12.3)') MOL_LEVID_LIST(K+1),'*',SE(MOL_ATMID_LIST(K))%BA(2,L:L+9,II,5)
                WRITE(13,'(I4,A2,5X,10ES12.3)') MOL_LEVID_LIST(K),'$*',SE(MOL_ATMID_LIST(K))%BA(2,L:L+9,II,5)
                K=K+2
                IF (K .GT. 12) EXIT
             END DO
!
          ELSE
             WRITE(13,'(T12,I3,T24,I3,T36,I3,T48,I3,T60,I3,T72,I3,T84,I3,T96,I3,T108,I3,T120,I3)') L,L+1,L+2,L+3,L+4,L+5,&
L+6,L+7,L+8,L+9
!          WRITE(12,*)
!
             DO K=1,4
                WRITE(13,'(I4,A1,5X,10ES12.3)') MOL_LEVID_LIST(K),'*',SE(MOL_ATMID_LIST(K))%BA(1,L:NT,II,5)
             END DO
!
             K=5
             DO
                WRITE(13,'(I4,A1,5X,10ES12.3)') MOL_LEVID_LIST(K),'*',SE(MOL_ATMID_LIST(K))%BA(1,L:NT,II,5)
                WRITE(13,'(I4,A1,5X,10ES12.3)') MOL_LEVID_LIST(K+1),'*',SE(MOL_ATMID_LIST(K))%BA(2,L:NT,II,5)
                WRITE(13,'(I4,A2,5X,10ES12.3)') MOL_LEVID_LIST(K),'$*',SE(MOL_ATMID_LIST(K))%BA(2,L:NT,II,5)
                K=K+2
                IF (K .GT. 12) EXIT
             END DO
          END IF
             
          WRITE(13,*)
          WRITE(13,*)
!
          L=L+10
       END DO LOOP2
!
       CLOSE(12)
       CLOSE(13)
!
    END IF
!
    RETURN
  END SUBROUTINE STEQ_BA_MOL_RXN
