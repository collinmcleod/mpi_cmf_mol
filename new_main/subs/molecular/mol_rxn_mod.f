!
! Module containing variables for chemical reactions involving molecules. 
! Created by Collin McLeod, Aug 11 2021
! Compare to module for charge exchange reactions (chg_exch_mod_v3.f at time of creation of this file)

      MODULE MOL_RXN_MOD
      USE SET_KIND_MODULE
          
          INTEGER N_MOL_RXN_RD
          INTEGER N_MOL_RXN
          INTEGER LUER1

! Arrays for reaction information
! For now I'm using two types of reactions: standard chemical reactions (3 coefficients)
! and electron collision reactions (using functional cross-sections)
! See the data formats in MOL_RXNS_DATA
! All data is sourced from the UMIST database for Astrochemistry (udfa.ajmarkwick.net/index.php)
! Full citation to come

          INTEGER, ALLOCATABLE :: MOL_RXN_TYPE_RD(:)
          REAL(KIND=LDP), ALLOCATABLE :: MOL_RXN_COEF_RD(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: MOL_RXN_TLO_RD(:) !Temp ranges
          REAL(KIND=LDP), ALLOCATABLE :: MOL_RXN_THI_RD(:)
          LOGICAL, ALLOCATABLE :: MOL_RXN_INCL_RD(:) !Whether reaction is included
          CHARACTER(8), ALLOCATABLE :: MOL_RXN_SPECIESID_RD(:,:)
          CHARACTER(30), ALLOCATABLE :: MOL_RXN_LEV_NAME_RD(:,:)
          CHARACTER(30), ALLOCATABLE :: MOL_RXN_ALT_LEV_NAME_RD(:,:)
          LOGICAL, ALLOCATABLE :: MOL_RXN_IS_FINAL_ION_RD(:,:)

! All reactions will be assumed to occur to and from ground states only, for simplicity

!Arrays to be populated after accounting for split states

          INTEGER, ALLOCATABLE :: MOL_RXN_TYPE(:)
          REAL(KIND=LDP), ALLOCATABLE :: MOL_RXN_COEF(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: MOL_RXN_G(:,:) !degeneracy for upper and lower states (I think)
          REAL(KIND=LDP), ALLOCATABLE :: MOL_RXN_TLO(:) !Temp ranges
          REAL(KIND=LDP), ALLOCATABLE :: MOL_RXN_THI(:)
          LOGICAL, ALLOCATABLE :: MOL_RXN_INCL(:) !Whether reaction is included
          CHARACTER(8), ALLOCATABLE :: MOL_RXN_SPECIESID(:,:)
          CHARACTER(30), ALLOCATABLE :: MOL_RXN_LEV_NAME(:,:)
          CHARACTER(30), ALLOCATABLE :: MOL_RXN_ALT_LEV_NAME(:,:)
          LOGICAL, ALLOCATABLE :: MOL_RXN_IS_FINAL_ION(:,:) !Whether the given species is the final ionization state
                                                            !currently included in CMFGEN for that element/molecule
                                                            !Needs to be treated slightly differently

! Arrays for updating STEQ and BA arrays

          REAL(KIND=LDP), ALLOCATABLE :: MOL_Z_CHG(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: AI_AR_MOL(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: DLNAI_AR_MOL_DLNT(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: COOL_MOL(:,:)

          INTEGER, ALLOCATABLE :: ID_ION_MOL(:,:)
          INTEGER, ALLOCATABLE :: LEV_IN_POPS_MOL(:,:)
          INTEGER, ALLOCATABLE :: LEV_IN_ION_MOL(:,:)
          LOGICAL, ALLOCATABLE :: MOL_REACTION_AVAILABLE(:)

          REAL(KIND=LDP), ALLOCATABLE :: MOL_RXN_RATE_SCALE(:)

          INTEGER, PARAMETER :: N_COEF_MAX_MOL=3 !May need to be modified when including
                                                  !electron collision reactions

          LOGICAL INITIALIZE_ARRAYS_MOL
          LOGICAL DO_MOL_RXNS
          LOGICAL, ALLOCATABLE :: LVL_IN_MOL_RXNS(:)
!
          REAL(KIND=LDP), ALLOCATABLE :: MTOT_BA(:,:,:), MTOT_STEQ(:,:) !Needed to store the net contributions from molecular reactions
                                                                 !Used in generate_full_matrix to update the ion equation properly
          INTEGER :: N_MOL_LEVS
          INTEGER, ALLOCATABLE :: LNK_MOL_TO_FULL(:)
          INTEGER, ALLOCATABLE :: LNK_FULL_TO_MOL(:)
          REAL(KIND=LDP) :: RXN_SCALE_FAC
          REAL(KIND=LDP) :: CO_SCALE_FAC
          LOGICAL :: DO_NT
!

        END MODULE MOL_RXN_MOD
