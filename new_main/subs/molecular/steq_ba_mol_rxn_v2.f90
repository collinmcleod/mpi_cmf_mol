!
!
! Updated version of steq_ba_mol_rxn.f90, intended to clean up garbage in the routine and hopefully eliminate 
! any lingering bugs
!
! This routine modifies the individual BA and STEQ matrices for each species involved in molecular reactions.  
! Should be called after RD_MOL_RXN, SET_MOL_RXN_LEV, and VERIFY_MOL_RXN
!
  SUBROUTINE STEQ_BA_MOL_RXN_V2(POPS,T,ED,NT,ND,DIAG_INDX,UPDATE_BA)
    USE MOL_RXN_MOD
    USE STEQ_DATA_MOD
    USE MOD_NON_THERM
    IMPLICIT NONE
!
! Edited 21-Sep-2022: Added non-thermal (Compton) dissociation reactions
! Edited 16-Jun-2022: Added reaction rate scaling options
! Created 1-Apr-2022: Created by Collin McLeod to address issues with existing STEQ_BA_MOL_RXN file
!

    !Argument variables
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
! Output variables
!
    INTEGER LU_ER
    INTEGER ERROR_LU
    EXTERNAL ERROR_LU
!
! Local Variables
! 
    INTEGER I,J,K,L
    INTEGER E1,E2,E3,E4
    INTEGER L1,L2,L3,L4
    INTEGER L11,L12,L13,L14
    INTEGER L21,L22,L23,L24
    INTEGER L31,L32,L33,L34
    INTEGER L41,L42,L43,L44
    INTEGER LNK1,LNK2,LNK3,LNK4,LNKT
    INTEGER NIV, II
    INTEGER ID,ID1,ID2,ID3,ID4
!
    REAL(8) RATE_K,SCALED_K
    REAL(8), ALLOCATABLE :: RATE_KS(:,:)
    REAL(8) ALPHA_1,ALPHA_2
    REAL(8) DLNK_DLNT,DRATEDT
    REAL(8) TVAL,TVAL4,TEN4
    REAL(8) FRD_R
    REAL(8) T1,TR1,TR2,TR3,TR4
!
! Other variables I use (mainly for printing info)
!
    LOGICAL :: PRINT_DEBUG_INFO=.TRUE.
    INTEGER :: MOL_FILE_NUM=60
    CHARACTER(30) :: MOL_DEBUG_FILE='MOLECULAR_REACTION_OUTPUT'
    CHARACTER(30) :: MOL_OUT_FILE_SORT='MOLECULAR_REACTIONS_SORTED'
    INTEGER :: IOS,LDIFF,CI1,CI2
    LOGICAL :: IS_OPEN
    REAL(8), ALLOCATABLE :: MOL_RXN_RATES_ALL(:,:), MOL_RXN_RATES_SORTED(:,:,:)
!
!    CHARACTER(20), ALLOCATABLE :: SPECIESIDLIST(:)
    INTEGER, ALLOCATABLE :: MOL_LEVID_LIST(:), MOL_ATMID_LIST(:)
    CHARACTER(20) :: SPECIESID1
    INTEGER :: SPECIESNUM, NLEVS_ACCOUNTED,NSPECIES
    INTEGER, ALLOCATABLE :: LEVS_ACCOUNTED_FOR(:)
!
! Variables for Compton reactions (dissociation by non-thermal electrons)
!
    REAL(8), ALLOCATABLE :: INIT_CO_ENERGIES(:), INIT_CO_XCS(:)
    REAL(8), ALLOCATABLE :: INIT_CPO_ENERGIES(:), INIT_CPO_XCS(:)
    REAL(8), ALLOCATABLE :: INIT_COP_ENERGIES(:), INIT_COP_XCS(:)
    REAL(8), ALLOCATABLE :: INTERP_CO_XCS(:), INTERP_CPO_XCS(:), INTERP_COP_XCS(:)
    REAL(8) :: COMP_POP
    LOGICAL :: DO_COMP_READ
    CHARACTER(100) :: TSTR1
    INTEGER :: INIT_COMP_N
    REAL(8) :: TEN16,MELEC,EVTOERG
!
! Check whether molecular reactions are included
!
    IF (.NOT. DO_MOL_RXNS) RETURN
!
    TEN4=10.0**4 ! To convert temperature units
    TEN16=1.0D-16 ! For converting the tabulated cross-sections
    MELEC= 9.10938D-28 ! Mass of electron, g
    EVTOERG= 6.2415D11 ! Conversion factor between erg and eV
!
    LU_ER=ERROR_LU()
!
    ALLOCATE(MOL_RXN_RATES_ALL(ND,N_MOL_RXN+1),STAT=IOS)
    IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_RATES_SORTED(ND,N_SPECIES,2),STAT=IOS)
    IF(IOS .EQ. 0) THEN
       MOL_RXN_RATES_ALL(:,:)=0.0
       MOL_RXN_RATES_SORTED(:,:,:)=0.0
    ELSE
       WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: '
       WRITE(LU_ER,*) 'Unable to allocate reaction rate arrays'
    END IF
!
    ALLOCATE(RATE_KS(ND,N_MOL_RXN),STAT=IOS)
    IF(IOS .EQ. 0) THEN
       RATE_KS(:,:)=0.0
    ELSE
       WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: '
       WRITE(LU_ER,*) 'Unable to allocate array for rate constants'
    END IF
!
    IF (.NOT. ALLOCATED(LVL_IN_MOL_RXNS)) ALLOCATE(LVL_IN_MOL_RXNS(NT),STAT=IOS)
    IF(IOS .EQ. 0) THEN
       LVL_IN_MOL_RXNS(:)=.FALSE.
    ELSE
       WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: '
       WRITE(LU_ER,*) 'Unable to allocate LVL_IN_MOL_RXNS'
    END IF
!
    IF (.NOT. ALLOCATED(MTOT_BA)) ALLOCATE(MTOT_BA(N_MOL_LEVS,N_MOL_LEVS+2,ND),STAT=IOS)
    IF (.NOT. ALLOCATED(MTOT_STEQ)) ALLOCATE(MTOT_STEQ(N_MOL_LEVS,ND),STAT=IOS)
    IF(IOS .EQ. 0) THEN
       MTOT_BA(:,:,:)=0.0D0
       MTOT_STEQ(:,:)=0.0D0
    ELSE
       WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: '
       WRITE(LU_ER,*) 'Unable to allocate MTOT arrays'
    END IF
!
    IF (.NOT. ALLOCATED(LNK_MOL_TO_FULL)) ALLOCATE(LNK_MOL_TO_FULL(N_MOL_LEVS+2),STAT=IOS)
    IF (.NOT. ALLOCATED(LNK_FULL_TO_MOL)) ALLOCATE(LNK_FULL_TO_MOL(NT),STAT=IOS)
    IF (.NOT. ALLOCATED(LEVS_ACCOUNTED_FOR)) ALLOCATE(LEVS_ACCOUNTED_FOR(N_MOL_LEVS),STAT=IOS)
    IF (IOS .EQ. 0) THEN
       LNK_MOL_TO_FULL(:)=0
       LNK_FULL_TO_MOL(:)=0
       LEVS_ACCOUNTED_FOR(:)=0
       NLEVS_ACCOUNTED=0
    ELSE
       WRITE(LU_ER,*) 'Error STEQ_BA_MOL_RXN: '
       WRITE(LU_ER,*) 'Unable to allocate MOL level links'
    END IF
!
! Establish the list of species relevant to molecules
! This is now accomplished in set_mol_rxn_lev when reading in reactions
!
    ALLOCATE(MOL_LEVID_LIST(N_SPECIES),STAT=IOS)
    IF (IOS .EQ. 0) ALLOCATE(MOL_ATMID_LIST(N_SPECIES),STAT=IOS)
    IF (IOS .EQ. 0) THEN
       MOL_LEVID_LIST(:)=0
       MOL_ATMID_LIST(:)=0       
    END IF
!
! Set certain MOL to FULL links
!
!electrons
    LNK_MOL_TO_FULL(N_MOL_LEVS+1)=NT-1 
    LNK_FULL_TO_MOL(NT-1)=N_MOL_LEVS+1
!Temperature
    LNK_MOL_TO_FULL(N_MOL_LEVS+2)=NT
    LNK_FULL_TO_MOL(NT)=N_MOL_LEVS+2
!
! Now start the main loop over all reactions
!
    DO J=1,N_MOL_RXN
       IF (MOL_REACTION_AVAILABLE(J)) THEN
          !
          ! Establish level links for species in this reaction
          !
          L1=LEV_IN_POPS_MOL(J,1)
          L2=LEV_IN_POPS_MOL(J,2)
          L3=LEV_IN_POPS_MOL(J,3)
          L4=LEV_IN_POPS_MOL(J,4)
          !
          ! Set proper level links for swapping from MOL to FULL (necessary for MTOT arrays)
          ! Note that L2 and L4 can be negative, indicating an electron or a photon, respectively
          ! If that is the case, we want to ignore them for most associations
          !
          IF (.NOT. ANY(LEVS_ACCOUNTED_FOR .EQ. L1)) THEN
             NLEVS_ACCOUNTED=NLEVS_ACCOUNTED+1
             LEVS_ACCOUNTED_FOR(NLEVS_ACCOUNTED)=L1
             LNK_MOL_TO_FULL(NLEVS_ACCOUNTED)=L1
             LNK_FULL_TO_MOL(L1)=NLEVS_ACCOUNTED
          END IF
!
          IF (.NOT. ANY(LEVS_ACCOUNTED_FOR .EQ. L2) .AND. L2 .GT. 0) THEN
             NLEVS_ACCOUNTED=NLEVS_ACCOUNTED+1
             LEVS_ACCOUNTED_FOR(NLEVS_ACCOUNTED)=L2
             LNK_MOL_TO_FULL(NLEVS_ACCOUNTED)=L2
             LNK_FULL_TO_MOL(L2)=NLEVS_ACCOUNTED
          END IF
!
          IF (.NOT. ANY(LEVS_ACCOUNTED_FOR .EQ. L3)) THEN
             NLEVS_ACCOUNTED=NLEVS_ACCOUNTED+1
             LEVS_ACCOUNTED_FOR(NLEVS_ACCOUNTED)=L3
             LNK_MOL_TO_FULL(NLEVS_ACCOUNTED)=L3
             LNK_FULL_TO_MOL(L3)=NLEVS_ACCOUNTED
          END IF
!
          IF (.NOT. ANY(LEVS_ACCOUNTED_FOR .EQ. L4) .AND. L4 .GT. 0) THEN
             NLEVS_ACCOUNTED=NLEVS_ACCOUNTED+1
             LEVS_ACCOUNTED_FOR(NLEVS_ACCOUNTED)=L4
             LNK_MOL_TO_FULL(NLEVS_ACCOUNTED)=L4
             LNK_FULL_TO_MOL(L4)=NLEVS_ACCOUNTED
          END IF
!
          !
          ! Update my logical array to indicate that these levels are involved in molecular reactions
          !
          LVL_IN_MOL_RXNS(L1)=.TRUE.
          IF (L2 .GT. 0) LVL_IN_MOL_RXNS(L2)=.TRUE.
          LVL_IN_MOL_RXNS(L3)=.TRUE.
          IF (L4 .GT. 0) LVL_IN_MOL_RXNS(L4)=.TRUE.
!
          !
          ! Now I establish the list of atom and level IDS so I know which IDs correspond to 
          ! the species I am interested in
          !
          DO K=1,4 !Species in reaction
             DO I=1,N_SPECIES !Species I care about
                IF (MOL_LEVID_LIST(I) .NE. 0) CYCLE
                IF (TRIM(ADJUSTL(MOL_RXN_SPECIESID(J,K))) .EQ. TRIM(ADJUSTL(SPECIESID_LIST(I)))) THEN
                   IF (K .EQ. 1) THEN
                      MOL_LEVID_LIST(I)=L1
                   ELSE IF (K .EQ. 2) THEN
                      MOL_LEVID_LIST(I)=L2
                   ELSE IF (K .EQ. 3) THEN
                      MOL_LEVID_LIST(I)=L3
                   ELSE IF (K .EQ. 4) THEN
                      MOL_LEVID_LIST(I)=L4
                   END IF
                   MOL_ATMID_LIST(I)=ID_ION_MOL(J,K)
                END IF
             END DO
          END DO
!
          !
          ! Now establish the other relevant level links needed to update STEQ and BA properly
          ! We have to be careful about the treatment of electrons and photons
          ! These additional checks only have to be included for species 2 and 4
          !
          ! ION IDs
          !
!
          IF (MOL_RXN_IS_FINAL_ION(J,1)) THEN !First species
             ID1=ID_ION_MOL(J,1)-1
          ELSE
             ID1=ID_ION_MOL(J,1)
          END IF
!
          L11=SE(ID1)%LNK_TO_IV(L1)
!
          IF (L2 .GT. 0) THEN !I.e., if L2 is NOT an electron
             L12=SE(ID1)%LNK_TO_IV(L2)
          ELSE
             L12=NT-1 !Electron equation
          END IF
!
          L13=SE(ID1)%LNK_TO_IV(L3)
!
          IF (L4 .GT. 0) THEN !I.e., if L4 is NOT a photon
             L14=SE(ID1)%LNK_TO_IV(L4)
          ELSE
             L14=0
          END IF 
!
          IF (MOL_RXN_IS_FINAL_ION(J,2)) THEN !Second species
             ID2=ID_ION_MOL(J,2)-1
          ELSE
             ID2=ID_ION_MOL(J,2)
          END IF
!
          IF (L2 .GT. 0) THEN 
             L21=SE(ID2)%LNK_TO_IV(L1)
             L22=SE(ID2)%LNK_TO_IV(L2)
             L23=SE(ID2)%LNK_TO_IV(L3)
             IF (L4 .GT. 0) THEN
                L24=SE(ID2)%LNK_TO_IV(L4)
             ELSE
                L24=0
             END IF
          ELSE !If L2 is an electron, nothing gets updated explicitly
             L21=0
             L22=0
             L23=0
             L24=0
          END IF
!
          IF (MOL_RXN_IS_FINAL_ION(J,3)) THEN !Third species
             ID3=ID_ION_MOL(J,3)-1
          ELSE
             ID3=ID_ION_MOL(J,3)
          END IF
!
          L31=SE(ID3)%LNK_TO_IV(L1)
!
          IF (L2 .GT. 0) THEN !I.e., if L2 is NOT an electron
             L32=SE(ID3)%LNK_TO_IV(L2)
          ELSE
             L32=NT-1 !Electron equation
          END IF
!
          L33=SE(ID3)%LNK_TO_IV(L3)
!
          IF (L4 .GT. 0) THEN !I.e., if L4 is NOT a photon
             L34=SE(ID3)%LNK_TO_IV(L4)
          ELSE
             L34=0
          END IF 
!
          IF (MOL_RXN_IS_FINAL_ION(J,4)) THEN !Fourth species
             ID4=ID_ION_MOL(J,4)-1
          ELSE
             ID4=ID_ION_MOL(J,4)
          END IF
!
          IF (L4 .GT. 0) THEN 
             L41=SE(ID4)%LNK_TO_IV(L1)
             IF (L2 .GT. 0) THEN
                L42=SE(ID4)%LNK_TO_IV(L2)
             ELSE
                L42=NT-1 !Electron equation
             END IF
             L43=SE(ID4)%LNK_TO_IV(L3)
             L44=SE(ID4)%LNK_TO_IV(L4)
          ELSE !If L4 is a photon, nothing gets updated explicitly
             L41=0
             L42=0
             L43=0
             L44=0
          END IF
!         
          !
          ! Now we associate the levels with the appropriate statistical equilibrium equations
          ! The conventions for LEV_IN_ION are the same as for LEV_IN_POPS:  electrons have level -1, photons -2
          !
          IF (MOL_RXN_IS_FINAL_ION(J,1)) THEN !Then we need to use the ion links
             E1=SE(ID1)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,1))
          ELSE !Then the equation is just the level in ion
             E1=LEV_IN_ION_MOL(J,1)
          END IF
!
          IF (MOL_RXN_IS_FINAL_ION(J,2)) THEN !Then we need to use the ion links
             E2=SE(ID2)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,2))
          ELSE !Then the equation is just the level in ion
             E2=LEV_IN_ION_MOL(J,2) !Note that if this is an electron, then E2 will be -1
          END IF
!
          IF (MOL_RXN_IS_FINAL_ION(J,3)) THEN !Then we need to use the ion links
             E3=SE(ID3)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,3))
          ELSE !Then the equation is just the level in ion
             E3=LEV_IN_ION_MOL(J,3)
          END IF
!
          IF (MOL_RXN_IS_FINAL_ION(J,4)) THEN !Then we need to use the ion links
             E4=SE(ID4)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_MOL(J,4))
          ELSE !Then the equation is just the level in ion
             E4=LEV_IN_ION_MOL(J,4) !If a photon, then E4 will be -2
          END IF
!
          !Quick check for errors
          IF (E1*E2*E3*E4 .EQ. 0) THEN
             WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: '
             WRITE(LU_ER,*) 'At least one of E1...E4 is zero--check level associations'
             WRITE(LU_ER,*) 'E1-E4: ',E1,E2,E3,E4
          END IF
!
          II=DIAG_INDX
!
          DO L=1,ND
!
             TVAL=MAX(T(L),MOL_RXN_TLO(J))
             TVAL=MIN(TVAL,MOL_RXN_THI(J))
             TVAL4=TEN4*TVAL !TVAL4 is the temperature in Kelvin
             !
             ! Now calculate the rate constant and dlnk_dlnt
             !
             IF ( (MOL_RXN_TYPE(J) .EQ. 1) .OR. (MOL_RXN_TYPE(J) .EQ. 2) ) THEN
                TR1 = EXP(-MOL_RXN_COEF(J,3)/TVAL4) !e^(-gamma/T) Needs to be in Kelvin
                TR2 = (TVAL4/300.0D0)**MOL_RXN_COEF(J,2) !(T/300 K)^beta, T in Kelvin
                RATE_K = MOL_RXN_COEF(J,1)*TR1*TR2 ! the rate constant K
                DLNK_DLNT = MOL_RXN_COEF(J,2) + (MOL_RXN_COEF(J,3)/TVAL4) !dlnk/dlnt = beta + gamma/T (in Kelvin)
                ! If T is outside the bounds, the derivative is set to 0
                IF( (T(L) .GE. MOL_RXN_THI(J)) .OR. (T(L) .LE. MOL_RXN_TLO(J)) ) DLNK_DLNT=0.0D0
             ELSE IF (MOL_RXN_TYPE(J) .EQ. 3) THEN
                IF (.NOT. ALLOCATED(XKT)) THEN
                   WRITE(*,*) 'Compton reactions included but non-thermal electron spectrum not available'
                   STOP
                END IF
                IF (MOL_RXN_COEF(J,1) .LE. 1.5) THEN
                   IF (.NOT. ALLOCATED(INIT_CO_XCS)) THEN
                      DO_COMP_READ=.TRUE.
                      OPEN(UNIT=30,FILE='COMPTON_XC_CO',ACTION='READ')
                   END IF
                ELSE IF (MOL_RXN_COEF(J,1) .LE. 2.5) THEN
                   IF (.NOT. ALLOCATED(INIT_CPO_XCS)) THEN
                      DO_COMP_READ=.TRUE.
                      OPEN(UNIT=30,FILE='COMPTON_XC_CPO',ACTION='READ')
                   END IF
                ELSE IF (MOL_RXN_COEF(J,1) .LE. 3.5) THEN
                   IF (.NOT. ALLOCATED(INIT_COP_XCS)) THEN
                      DO_COMP_READ=.TRUE.
                      OPEN(UNIT=30,FILE='COMPTON_XC_COP',ACTION='READ')
                   END IF
                ELSE
                   WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN:'
                   WRITE(LU_ER,*) 'Unrecognized Compton reaction: rxn # ',MOL_RXN_COEF(J,1)
                   STOP
                END IF
                IF (DO_COMP_READ) THEN
                   DO
                      READ(30,'(A)') TSTR1
                      IF (INDEX(TSTR1,'!Number of') .NE. 0) THEN
                         READ(TSTR1,'(I3,A)') INIT_COMP_N
                         EXIT
                      END IF
                   END DO
                   IF (MOL_RXN_COEF(J,1) .LE. 1.5) THEN
                      ALLOCATE(INIT_CO_ENERGIES(INIT_COMP_N))
                      INIT_CO_ENERGIES(:)=0.0D0
                      ALLOCATE(INIT_CO_XCS(INIT_COMP_N))
                      INIT_CO_XCS(:)=0.0D0
                      ALLOCATE(INTERP_CO_XCS(NKT))
                      INTERP_CO_XCS(:)=0.0D0
                      READ(30,'(A)') TSTR1
                      DO I=1,INIT_COMP_N
                         READ(30,*) INIT_CO_ENERGIES(I),INIT_CO_XCS(I)
                   !                      READ(30,*) INIT_COMP_XCS(I)
                      END DO
                   ELSE IF (MOL_RXN_COEF(J,1) .LE. 2.5) THEN
                      ALLOCATE(INIT_CPO_ENERGIES(INIT_COMP_N))
                      INIT_CPO_ENERGIES(:)=0.0D0
                      ALLOCATE(INIT_CPO_XCS(INIT_COMP_N))
                      INIT_CPO_XCS(:)=0.0D0
                      ALLOCATE(INTERP_CPO_XCS(NKT))
                      INTERP_CPO_XCS(:)=0.0D0
                      READ(30,'(A)') TSTR1
                      DO I=1,INIT_COMP_N
                         READ(30,*) INIT_CPO_ENERGIES(I),INIT_CPO_XCS(I)
                   !                      READ(30,*) INIT_COMP_XCS(I)
                      END DO
                   ELSE IF (MOL_RXN_COEF(J,1) .LE. 3.5) THEN
                      ALLOCATE(INIT_COP_ENERGIES(INIT_COMP_N))
                      INIT_COP_ENERGIES(:)=0.0D0
                      ALLOCATE(INIT_COP_XCS(INIT_COMP_N))
                      INIT_COP_XCS(:)=0.0D0
                      ALLOCATE(INTERP_COP_XCS(NKT))
                      INTERP_COP_XCS(:)=0.0D0
                      READ(30,'(A)') TSTR1
                      DO I=1,INIT_COMP_N
                         READ(30,*) INIT_COP_ENERGIES(I),INIT_COP_XCS(I)
                   !                      READ(30,*) INIT_COMP_XCS(I)
                      END DO
                   END IF
                   CLOSE(30)
                   DO_COMP_READ=.FALSE.
                END IF
!                WRITE(LU_ER,*) 'STEQ_BA_MOL_RXN Version 28'
!                WRITE(LU_ER,*) 'RXN Identifier: ',MOL_RXN_COEF(J,1)
!                WRITE(LU_ER,*) 'Sigmas before interp:', INIT_COMP_XCS(1:10)
                !                   WRITE(LU_ER,*) 'Energies read in: ', INIT_COMP_ENERGIES
                !                   WRITE(LU_ER,*) 'Sigmas read in: ', INIT_COMP_XCS
                IF (MOL_RXN_COEF(J,1) .LE. 1.5) THEN
                   CALL MON_INTERP_FAST_V2(INTERP_CO_XCS,NKT,1,XKT,NKT,&
                     INIT_CO_XCS,INIT_COMP_N,INIT_CO_ENERGIES,INIT_COMP_N,'STEQ_BA_MOL_RXN_V2') !IDENT NEEDED
                ELSE IF (MOL_RXN_COEF(J,1) .LE. 2.5) THEN
                   CALL MON_INTERP_FAST_V2(INTERP_CPO_XCS,NKT,1,XKT,NKT,&
                     INIT_CPO_XCS,INIT_COMP_N,INIT_CPO_ENERGIES,INIT_COMP_N,'STEQ_BA_MOL_RXN_V2') !IDENT NEEDED
                ELSE IF (MOL_RXN_COEF(J,1) .LE. 3.5) THEN
                   CALL MON_INTERP_FAST_V2(INTERP_COP_XCS,NKT,1,XKT,NKT,&
                     INIT_COP_XCS,INIT_COMP_N,INIT_COP_ENERGIES,INIT_COMP_N,'STEQ_BA_MOL_RXN_V2') !IDENT NEEDED
                END IF
                !                   WRITE(LU_ER,*) 'Energies after interp: ', XKT
!                IF (MOL_RXN_COEF(J,1) .LE. 1.5) THEN
!                   WRITE(LU_ER,*) 'Sigmas after interp: ', INTERP_CO_XCS(10:20)
!                ELSE IF (MOL_RXN_COEF(J,1) .LE. 2.5) THEN
!                   WRITE(LU_ER,*) 'Sigmas after interp: ', INTERP_CPO_XCS(10:20)
!                ELSE IF (MOL_RXN_COEF(J,1) .LE. 3.5) THEN
!                   WRITE(LU_ER,*) 'Sigmas after interp: ', INTERP_COP_XCS(10:20)
!                END IF
!                WRITE(LU_ER,*) 'STEQ_BA_MOL_RXN Version 28'
                RATE_K=0.0D0
                TR1=EVTOERG/MELEC
                COMP_POP = 0.0D0
                DO I=1,NKT
                   IF (MOL_RXN_COEF(J,1) .LE. 1.5) THEN
                      RATE_K = RATE_K + YE(I,L)*(INTERP_CO_XCS(I))*DXKT(I)
!                      IF (L .EQ. 1 .AND. I .GE. 10 .AND. I .LE. 20) THEN
!                         WRITE(LU_ER,*) 'YE(I,L), SQRT(2*XKT(I)*TR1), INTERP_CO_XCS(I), DXKT(I):'
!                         WRITE(LU_ER,*) YE(I,L), SQRT(2.0D0*XKT(I)*TR1), INTERP_CO_XCS(I), DXKT(I)
!                      END IF
                   ELSE IF (MOL_RXN_COEF(J,1) .LE. 2.5) THEN
                      RATE_K = RATE_K + YE(I,L)*(INTERP_CPO_XCS(I))*DXKT(I)
                   ELSE IF (MOL_RXN_COEF(J,1) .LE. 3.5) THEN
                      RATE_K = RATE_K + YE(I,L)*(INTERP_COP_XCS(I))*DXKT(I)
                   END IF
                   COMP_POP = COMP_POP + (YE(I,L)/XKT(I))*DXKT(I)
                END DO
                COMP_POP=SQRT(MELEC/2.0D0)*COMP_POP !To reduce total number of operations
!
! YE is the degradation spectrum.  I want to evaluate the integral int( sigma*n_e*v*dv ). Apparently, the
! nt electron number density is given by YE/v*dv, so the v and 1/v cancel out, leaving the integral 
! sigma*YE*dv
!
                DLNK_DLNT=0.0D0
                !                RATE_K = RATE_K*TEN16
                !CALL RD_NON_THERM_ELEC_SPEC_V1(ND,LU_ER)             
             ELSE
                WRITE(LU_ER,*) 'Error in STEQ_BA_MOL_RXN: '
                WRITE(LU_ER,*) 'Reaction type not recognized--type is ',MOL_RXN_TYPE(J)
                WRITE(LU_ER,*) 'Allowed reaction types are 1, 2 and 3'
             END IF
!
             SCALED_K = RXN_SCALE_FAC*MOL_RXN_RATE_SCALE(J)*RATE_K
             IF (MOL_RXN_SPECIESID(J,3) .EQ. 'COMI' .OR. MOL_RXN_SPECIESID(J,4) .EQ. 'COMI') THEN
                SCALED_K = SCALED_K*CO_SCALE_FAC
             END IF
!
! Multiple Scalings applied to the rate constant:
! one to account for split product states (for physical consistency)
! one to reduce the rates of all reactions (for diagnostic purposes)
! one to reduce the rates of CO-producing reactions (for diagnostic purposes)
!
             RATE_KS(L,J) = SCALED_K ! Recorded for printing debugging info
             !
             ! FRD_R is the forward reaction rate (I am not computing reverse rates)
             ! If L2 is an electron, we use the electron density instead to calculate the rate
             ! unit is number of reactions /cm^3 /s
             !
             IF (L2 .GT. 0) THEN
                FRD_R = SCALED_K*POPS(L1,L)*POPS(L2,L) 
             ELSE IF (MOL_RXN_TYPE(J) .EQ. 3) THEN
                FRD_R = SCALED_K*POPS(L1,L) ! The electron population is accounted for in the evaluation of RATE_K
             ELSE
                FRD_R = SCALED_K*POPS(L1,L)*ED(L)
             END IF
             !
             ! Now update statistical equilibrium
             !
             SE(ID1)%STEQ(E1,L)=SE(ID1)%STEQ(E1,L) - FRD_R
             IF (L2 .GT. 0) SE(ID2)%STEQ(E2,L)=SE(ID2)%STEQ(E2,L) - FRD_R !No update is done for electrons
             SE(ID3)%STEQ(E3,L)=SE(ID3)%STEQ(E3,L) + FRD_R
             IF (L4 .GT. 0) SE(ID4)%STEQ(E4,L)=SE(ID4)%STEQ(E4,L) + FRD_R !No update is done for photons
             !
             ! We also record the contributions in the MTOT array
             !
             LNK1=LNK_FULL_TO_MOL(L1)
             IF (L2 .GT. 0) THEN
                LNK2=LNK_FULL_TO_MOL(L2)
             ELSE
                LNK2=LNK_FULL_TO_MOL(NT-1)
             END IF
             LNK3=LNK_FULL_TO_MOL(L3)
             IF (L4 .GT. 0) THEN
                LNK4=LNK_FULL_TO_MOL(L4)
             ELSE
                LNK4=0
             END IF
             LNKT = LNK_FULL_TO_MOL(NT)
!
             MTOT_STEQ(LNK1,L)=MTOT_STEQ(LNK1,L) - FRD_R
             IF (L2 .GT. 0) MTOT_STEQ(LNK2,L)=MTOT_STEQ(LNK2,L) - FRD_R
             MTOT_STEQ(LNK3,L)=MTOT_STEQ(LNK3,L) + FRD_R
             IF (L4 .GT. 0) MTOT_STEQ(LNK4,L)=MTOT_STEQ(LNK4,L) + FRD_R
!
             DRATEDT=(FRD_R*DLNK_DLNT)/TVAL !This accounts for the unit conversion (K to 10^4 K)
!
             IF (UPDATE_BA) THEN
!
                IF (MOL_RXN_TYPE(J) .NE. 3) THEN
                   ALPHA_2 = SCALED_K*POPS(L1,L) !The derivative w.r.t. Q2
                   IF (L2 .GT. 0) THEN
                      ALPHA_1 = SCALED_K*POPS(L2,L) !The derivative w.r.t. Q1
                   ELSE
                      ALPHA_1 = SCALED_K*ED(L) !Derivative if Q2 is an electron
                   END IF
                ELSE
!                   ALPHA_2 = (SCALED_K/COMP_POP)*POPS(L1,L)
                   ALPHA_2 = 0
                   IF (L2 .GT. 0) THEN
                      ALPHA_1 = (SCALED_K/COMP_POP)*POPS(L2,L)
                   ELSE
                      ALPHA_1 = SCALED_K
                   END IF
                END IF
!
                IF (L1 .EQ. L2) THEN 
                   !
                   ! This scenario requires special consideration due to the additional quadratic dependence
                   ! In this case, the rate is K*Q1^2, so derivatives w.r.t. Q1 are 2*K*Q1
                   ! In addition, the derivative of Q1 w.r.t. Q1 is actually 4*K*Q1 due to the fact that the 
                   ! rate of change of Q1 is actually double the reaction rate
                   ! Note also that this case precludes the possibility of Q2 being an electron
                   !
                   ID=ID1
                   NIV=SE(ID)%N_IV
                   IF (L11 .NE. 0) SE(ID)%BA(E1,L11,II,L)=SE(ID)%BA(E1,L11,II,L) - 4*ALPHA_2
                   SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) - 2*DRATEDT
                   !
                   ! Also record the contributions in the MTOT matrix
                   !
                   IF (L11 .NE. 0) MTOT_BA(LNK1,LNK1,L)=MTOT_BA(LNK1,LNK1,L)-4*ALPHA_2
                   MTOT_BA(LNK1,LNKT,L)=MTOT_BA(LNK1,LNKT,L)-2*DRATEDT
!
                   ID=ID3
                   NIV=SE(ID)%N_IV
                   IF (L31 .NE. 0) SE(ID)%BA(E3,L31,II,L)=SE(ID)%BA(E3,L31,II,L) + 2*ALPHA_2
                   SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) + DRATEDT
!
                   IF (L31 .NE. 0) MTOT_BA(LNK3,LNK1,L)=MTOT_BA(LNK3,LNK1,L)+2*ALPHA_2
                   MTOT_BA(LNK3,LNKT,L)=MTOT_BA(LNK3,LNKT,L)+DRATEDT
!
                   IF (L4 .GT. 0) THEN
                      ID=ID4
                      NIV=SE(ID)%N_IV
                      IF(L41 .NE. 0) SE(ID)%BA(E4,L41,II,L)=SE(ID)%BA(E4,L41,II,L)+2*ALPHA_2
                      SE(ID)%BA(E4,NIV,II,L)=SE(ID)%BA(E4,NIV,II,L)+DRATEDT
!
                      IF (L41 .NE. 0) MTOT_BA(LNK4,LNK1,L)=MTOT_BA(LNK4,LNK1,L)+2*ALPHA_2
                      MTOT_BA(LNK4,LNKT,L)=MTOT_BA(LNK4,LNKT,L)+DRATEDT
                   END IF
!
                ELSE
!
                   ID=ID1
                   NIV=SE(ID)%N_IV
                   IF (L11 .NE. 0) SE(ID)%BA(E1,L11,II,L)=SE(ID)%BA(E1,L11,II,L)-ALPHA_1
                   !
                   ! If Q2 is an electron, then L12 will be set to NT-1--the electron equation
                   ! Thus the update should work correctly with just this one line
                   !
                   IF (L12 .NE. 0) SE(ID)%BA(E1,L12,II,L)=SE(ID)%BA(E1,L12,II,L)-ALPHA_2
                   SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L)-DRATEDT
!
                   IF (L11 .NE. 0) MTOT_BA(LNK1,LNK1,L)=MTOT_BA(LNK1,LNK1,L)-ALPHA_1
                   IF (L12 .NE. 0) MTOT_BA(LNK1,LNK2,L)=MTOT_BA(LNK1,LNK2,L)-ALPHA_2
                   MTOT_BA(LNK1,LNKT,L)=MTOT_ba(LNK1,LNKT,L)-DRATEDT
!
                   IF (L2 .GT. 0) THEN !Do not update the electron equation explicitly
                      ID=ID2
                      NIV=SE(ID)%N_IV
                      IF (L21 .NE. 0) SE(ID)%BA(E2,L21,II,L)=SE(ID)%BA(E2,L21,II,L)-ALPHA_1
                      IF (L22 .NE. 0) SE(ID)%BA(E2,L22,II,L)=SE(ID)%BA(E2,L22,II,L)-ALPHA_2
                      SE(ID)%BA(E2,NIV,II,L)=SE(ID)%BA(E2,NIV,II,L)-DRATEDT
!
                      IF (L21 .NE. 0) MTOT_BA(LNK2,LNK1,L)=MTOT_BA(LNK2,LNK1,L)-ALPHA_1
                      IF (L22 .NE. 0) MTOT_BA(LNK2,LNK2,L)=MTOT_BA(LNK2,LNK2,L)-ALPHA_2
                      MTOT_BA(LNK2,LNKT,L)=MTOT_BA(LNK2,LNKT,L)-DRATEDT
                   END IF
!
                   ID=ID3
                   NIV=SE(ID)%N_IV
                   IF (L31 .NE. 0) SE(ID)%BA(E3,L31,II,L)=SE(ID)%BA(E3,L31,II,L)+ALPHA_1
                   IF (L32 .NE. 0) SE(ID)%BA(E3,L32,II,L)=SE(ID)%BA(E3,L32,II,L)+ALPHA_2
                   SE(ID)%BA(E3,NIV,II,L)=SE(ID)%BA(E3,NIV,II,L)+DRATEDT
!
                   IF (L31 .NE. 0) MTOT_BA(LNK3,LNK1,L)=MTOT_BA(LNK3,LNK1,L)+ALPHA_1
                   IF (L32 .NE. 0) MTOT_BA(LNK3,LNK2,L)=MTOT_BA(LNK3,LNK2,L)+ALPHA_2
                   MTOT_BA(LNK3,LNKT,L)=MTOT_BA(LNK3,LNKT,L)+DRATEDT
!
                   IF (L4 .GT. 0) THEN !Only do the updates if Q4 is an actual ion
                      ID=ID4
                      NIV=SE(ID)%N_IV
                      IF (L41 .NE. 0) SE(ID)%BA(E4,L41,II,L)=SE(ID)%BA(E4,L41,II,L)+ALPHA_1
                      IF (L42 .NE. 0) SE(ID)%BA(E4,L42,II,L)=SE(ID)%BA(E4,L42,II,L)+ALPHA_2
                      SE(ID)%BA(E4,NIV,II,L)=SE(ID)%BA(E4,NIV,II,L)+DRATEDT
!
                      IF (L41 .NE. 0) MTOT_BA(LNK4,LNK1,L)=MTOT_BA(LNK4,LNK1,L)+ALPHA_1
                      IF (L42 .NE. 0) MTOT_BA(LNK4,LNK2,L)=MTOT_BA(LNK4,LNK2,L)+ALPHA_2
                      MTOT_BA(LNK4,LNKT,L)=MTOT_BA(LNK4,LNKT,L)+DRATEDT
                   END IF
!
                END IF ! L1 == L2
!
             END IF ! UPDATE_BA
!
             IF (PRINT_DEBUG_INFO) THEN
                IF (L .EQ. 1 .AND. J .EQ. 1) THEN
                   OPEN(UNIT=MOL_FILE_NUM,FILE=MOL_DEBUG_FILE,STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
                   WRITE(MOL_FILE_NUM,*)'**************************************************************************'
                   WRITE(MOL_FILE_NUM,*)'Info for molecular reactions, final iteration'
                   WRITE(MOL_FILE_NUM,*)'First line: reaction string'
                   WRITE(MOL_FILE_NUM,*)'Second line: Depth, Temperature (10^4 K), Q1, Q2'
                   WRITE(MOL_FILE_NUM,*)'Third line: Alpha, Beta, Gamma'
                   WRITE(MOL_FILE_NUM,*)'Fourth line: Rate constant K, Forward rate, DRATEDT'
                   WRITE(MOL_FILE_NUM,*)'Fifth line: LEV_IN_POPS L1 through L4'
                   WRITE(MOL_FILE_NUM,*)'Sixth line: Alpha_1 (derivative w.r.t. Q1) and Alpha_2 (derivative w.r.t. Q2)'
                   WRITE(MOL_FILE_NUM,*)'Seventh line: Ion IDs 1-4'
                   WRITE(MOL_FILE_NUM,*)'Eighth line: LEV_IN_ION 1-4'
                   WRITE(MOL_FILE_NUM,*)'Ninth line: E values (equation numbers)'
                   WRITE(MOL_FILE_NUM,*)'Tenth line: L1 links'
                   WRITE(MOL_FILE_NUM,*)'Eleventh line: L2 links'
                   WRITE(MOL_FILE_NUM,*)'Twelfth line: L3 links'
                   WRITE(MOL_FILE_NUM,*)'Thirteenth line: L4 links'
                   WRITE(MOL_FILE_NUM,*)'**************************************************************************'
                   IF (IOS .NE. 0) WRITE(LU_ER,*) 'Error opening molecular debug info file'
                END IF
                WRITE(MOL_FILE_NUM,'(A16,I2,A1,2X,A8,1X,A1,1X,A8,1X,A2,1X,A8,1X,A1,1X,A8)')'Reaction Number ',&
J,':',MOL_RXN_SPECIESID(J,1),'+',MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4)
                IF (L1 .GE. 1 .AND. L2 .GE. 1) THEN
                   WRITE(MOL_FILE_NUM,'(I3,A1,1X,F7.5,A1,1X,ES15.5,A1,1X,ES15.5)')L,',',T(L),',',POPS(L1,L),',',POPS(L2,L)
                END IF
                WRITE(MOL_FILE_NUM,'(ES15.5,A1,1X,F10.3,A1,1X,F15.5)')MOL_RXN_COEF(J,1),',',MOL_RXN_COEF(J,2),',',&
MOL_RXN_COEF(J,3)
                WRITE(MOL_FILE_NUM,'(ES15.5,A1,1X,ES15.5,A1,1X,ES15.5)')RATE_K,',',FRD_R,',',DRATEDT
                WRITE(MOL_FILE_NUM,'(I4,A1,1X,I4,A1,1X,I4,A1,1X,I4)')L1,',',L2,',',L3,',',L4
                WRITE(MOL_FILE_NUM,'(A,ES11.4,A1,1X,ES11.4)')'Alphas: ',ALPHA_1,',',ALPHA_2
                WRITE(MOL_FILE_NUM,'(A,4I5)')'Ion IDs: ',ID1,ID2,ID3,ID4
                WRITE(MOL_FILE_NUM,'(A,4I5)')'Levels in ion: ',LEV_IN_ION_MOL(J,:)
                WRITE(MOL_FILE_NUM,'(A,4I6)')'E values: ',E1,E2,E3,E4
                WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L1 links: ',L11,L12,L13,L14
                WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L2 links: ',L21,L22,L23,L24
                WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L3 links: ',L31,L32,L33,L34
                WRITE(MOL_FILE_NUM,'(A,5X,4I5)')'L4 links: ',L41,L42,L43,L44
                WRITE(MOL_FILE_NUM,*)
             END IF
!
             !
             ! Record the reaction rates and sort the fastest for printing
             ! 
             MOL_RXN_RATES_ALL(L,J)=FRD_R
             IF (ABS(FRD_R) .GT. ABS(MOL_RXN_RATES_SORTED(L,12,1))) THEN
                CHECKLOOP: DO CI1=1,12
                   IF (ABS(FRD_R) .GT. ABS(MOL_RXN_RATES_SORTED(L,CI1,1))) THEN
                      CI2=CI1
                      EXIT CHECKLOOP
                   END IF
                END DO CHECKLOOP
                IF (CI2 .EQ. 12) THEN
                   MOL_RXN_RATES_SORTED(L,12,1)=FRD_R
                   MOL_RXN_RATES_SORTED(L,12,2)=J
                ELSE
                   MOL_RXN_RATES_SORTED(L,CI2+1:12,:) = MOL_RXN_RATES_SORTED(L,CI2:11,:)
                   MOL_RXN_RATES_SORTED(L,CI2,1)=FRD_R
                   MOL_RXN_RATES_SORTED(L,CI2,2)=J
                END IF
             END IF
!
          END DO ! Loop over depths
!
       END IF ! Reaction available
    END DO
!
    INITIALIZE_ARRAYS_MOL=.TRUE. !So that arrays will be reset on the next call to mol_rxn routines
!
    INQUIRE(UNIT=MOL_FILE_NUM,OPENED=IS_OPEN)
    IF (IS_OPEN) THEN
       CALL FLUSH(MOL_FILE_NUM)
       CLOSE(MOL_FILE_NUM)
    END IF
!
    !
    ! Write out a file with all the reaction rates at each depth
    !
    OPEN(UNIT=MOL_FILE_NUM,FILE=MOL_OUT_FILE_SORT,STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
    WRITE(MOL_FILE_NUM,*)'*****************************************************************************************'
    WRITE(MOL_FILE_NUM,*)'Info for molecular reactions, final iteration'
    WRITE(MOL_FILE_NUM,*)'Total reaction rates (cm^-3 s^-1), 10 depths at a time'
    WRITE(MOL_FILE_NUM,*)'*****************************************************************************************'
    WRITE(MOL_FILE_NUM,*)
    L=1
    WRITELOOP1: DO
       IF (L .LE. ND-9) THEN
          WRITE(MOL_FILE_NUM,'(A7,I3,A1,I3,1X,A13,T30,10ES15.5)')'Depths ',L,'-',L+9,'T (10^4 K) = ',T(L:L+9)
          WRITE(MOL_FILE_NUM,*)
          DO J=1,N_MOL_RXN
             WRITE(MOL_FILE_NUM,'(A5,A1,1X,A5,A2,1X,A5,A1,1X,A5,A1,T30,10ES15.5)')MOL_RXN_SPECIESID(J,1),'+',&
MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4),':',MOL_RXN_RATES_ALL(L:L+9,J)
          END DO
          WRITE(MOL_FILE_NUM,*)
          L=L+10
       ELSE
          LDIFF=ND-L
          WRITE(MOL_FILE_NUM,'(A7,I3,A1,I3,1X,A13,T30,10ES15.5)')'Depths ',L,'-',ND,'T (10^4 K) = ',T(L:ND)
          DO J=1,N_MOL_RXN
             WRITE(MOL_FILE_NUM,'(A5,A1,1X,A5,A2,1X,A5,A1,1X,A5,A1,T30,10ES15.5)')MOL_RXN_SPECIESID(J,1),'+',&
MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4),':',MOL_RXN_RATES_ALL(L:ND,J)
          END DO
          EXIT WRITELOOP1
       END IF
    END DO WRITELOOP1
    CALL FLUSH(MOL_FILE_NUM)
    CLOSE(MOL_FILE_NUM)
!
    WRITE(*,*) 'Calling print_mol_info'
    CALL FLUSH(6)
    CALL PRINT_MOL_INFO(MOL_RXN_RATES_ALL,POPS,MOL_LEVID_LIST,MOL_ATMID_LIST,ND,NT)
!
    RETURN
  END SUBROUTINE STEQ_BA_MOL_RXN_V2
