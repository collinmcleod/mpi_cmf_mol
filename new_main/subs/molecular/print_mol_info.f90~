!
! A Program to write out molecular reaction information in a useful way
! migrated from steq_ba_mol_rxn.f90 due to increased coding requirements related to using 
! split states, as well as for a way to access photoionization/recombination information
!
  SUBROUTINE PRINT_MOL_INFO(MOL_RXN_RATES_ALL,POPS,MOL_LEVID_LIST,MOL_ATMID_LIST,ND,NT)
    USE MOL_RXN_MOD
    USE MOD_CMFGEN
    USE SET_KIND_MODULE
!
    IMPLICIT NONE
!
! Created 17-Feb-2022
!
    INTEGER NT
    INTEGER ND
    REAL(KIND=LDP) POPS(NT,ND)
!
    REAL(KIND=LDP) MOL_RXN_RATES_ALL(ND,N_MOL_RXN+1)
    INTEGER MOL_LEVID_LIST(N_SPECIES)
    INTEGER MOL_ATMID_LIST(N_SPECIES)
!
    INTEGER N_RXNS_UNSPLIT
    CHARACTER(20), ALLOCATABLE :: REACTANTS(:), PRODUCTS(:)
    REAL(KIND=LDP), ALLOCATABLE :: TOTAL_CHANGE_BY_SPECIES(:,:)
    INTEGER, ALLOCATABLE :: TOP_RXNS_BY_SPECIES(:,:,:,:)
    REAL(KIND=LDP), ALLOCATABLE :: NEW_MOL_RXN_RATES_ALL(:,:)
    INTEGER, ALLOCATABLE :: RXN_POINTERS(:), REVERSE_POINTERS(:)
    INTEGER :: I,J,K,L,IOS,RN1,RN2,ID1,ID2,ID3,ID4,ID1R,ID2R,ID3R,ID4R
    CHARACTER(20) :: OUTFILE_NAME
    REAL(KIND=LDP) :: TOT_PR,TOT_RR,TR1,TR2
    CHARACTER(200) :: RXN_STRING1,RXN_STRING2
    INTEGER, ALLOCATABLE :: SPECIES_RELEVANT_LEVELS(:,:)
    LOGICAL, ALLOCATABLE :: SPECIES_ISFINAL(:)
    CHARACTER(20) :: SPECIESID1
    LOGICAL :: MATCHED
!
    IF (.NOT. DO_MOL_RXNS) RETURN
!
    WRITE(*,*) 'Starting print_mol_info'
    CALL FLUSH(6)
!
    IF (.NOT. ALLOCATED(RXN_POINTERS)) ALLOCATE(RXN_POINTERS(N_MOL_RXN))
!
! Here I create some pointers so that reactions with split states can be printed all together
!
    J=1
    RXN_POINTERS(1)=1
    DO I=2,N_MOL_RXN
       IF ( (MOL_RXN_SPECIESID(I,1) .EQ. MOL_RXN_SPECIESID(I-1,1)) .AND. &
(MOL_RXN_SPECIESID(I,2) .EQ. MOL_RXN_SPECIESID(I-1,2)) .AND. &
(MOL_RXN_SPECIESID(I,3) .EQ. MOL_RXN_SPECIESID(I-1,3)) .AND. &
(MOL_RXN_SPECIESID(I,4) .EQ. MOL_RXN_SPECIESID(I-1,4)) )THEN
          RXN_POINTERS(I)=J
       ELSE
          J=J+1
          RXN_POINTERS(I)=J
       END IF
    END DO
!
    N_RXNS_UNSPLIT=J
    IF (.NOT. ALLOCATED(NEW_MOL_RXN_RATES_ALL)) ALLOCATE(NEW_MOL_RXN_RATES_ALL(ND,N_RXNS_UNSPLIT+1))
    NEW_MOL_RXN_RATES_ALL = 0.0_LDP
!
! Reverse pointers point from the new reaction list to the first matching reaction in the original list
! Needed for getting the rest of the reaction info out (like species IDs, etc.)
!
    IF (.NOT. ALLOCATED(REVERSE_POINTERS)) ALLOCATE(REVERSE_POINTERS(N_RXNS_UNSPLIT))
    J=1
    REVERSE_POINTERS(1)=1
    DO I=2,N_MOL_RXN
       IF (RXN_POINTERS(I) .EQ. RXN_POINTERS(I-1)) THEN
          CYCLE
       ELSE
          J=J+1
          REVERSE_POINTERS(J)=I
       END IF
    END DO
!
! Now group all the reaction rates which correspond to the same reaction w/ different split states
!
    DO K=1,ND
       DO I=1,N_MOL_RXN
          J = RXN_POINTERS(I)
          NEW_MOL_RXN_RATES_ALL(K,J) = NEW_MOL_RXN_RATES_ALL(K,J) + MOL_RXN_RATES_ALL(K,I)
       END DO
    END DO
!
    IF (.NOT. ALLOCATED(TOP_RXNS_BY_SPECIES)) ALLOCATE(TOP_RXNS_BY_SPECIES(ND,N_SPECIES,2,5),STAT=IOS)
    IF (IOS .EQ. 0 .AND. .NOT. ALLOCATED(TOTAL_CHANGE_BY_SPECIES)) ALLOCATE(TOTAL_CHANGE_BY_SPECIES(ND,N_SPECIES),STAT=IOS)
    IF (IOS .EQ. 0 .AND. .NOT. ALLOCATED(REACTANTS)) ALLOCATE(REACTANTS(2))
    IF (IOS .EQ. 0 .AND. .NOT. ALLOCATED(PRODUCTS)) ALLOCATE(PRODUCTS(2))
    IF (IOS .NE. 0) WRITE(*,*) 'Error allocating arrays in print_mol_info.f90'
!
    TOP_RXNS_BY_SPECIES = N_RXNS_UNSPLIT+1
    TOTAL_CHANGE_BY_SPECIES=0.0_LDP
!
! Now loop over all reactions and determine which are the most important creation/destruction reactions for each 
! species included in the model
!
    DO L=1,ND
       DO J=1,N_RXNS_UNSPLIT
          RN1 = REVERSE_POINTERS(J) !Needed to pull reaction info from module
          REACTANTS(1) = MOL_RXN_SPECIESID(RN1,1)
          REACTANTS(2) = MOL_RXN_SPECIESID(RN1,2)
          PRODUCTS(1) = MOL_RXN_SPECIESID(RN1,3)
          PRODUCTS(2) = MOL_RXN_SPECIESID(RN1,4)
          DO K=1,N_SPECIES
             SPECIESID1 = SPECIESID_LIST(K)
             IF (SPECIESID1 .EQ. REACTANTS(1) .OR. SPECIESID1 .EQ. REACTANTS(2)) THEN
                DO ID1=1,5
                   ID2 = TOP_RXNS_BY_SPECIES(L,K,2,ID1)
                   IF(ABS(NEW_MOL_RXN_RATES_ALL(L,J)) .GT. ABS(NEW_MOL_RXN_RATES_ALL(L,ID2))) THEN
                      IF (ID1 .LT. 5) TOP_RXNS_BY_SPECIES(L,K,2,ID1+1:5) = TOP_RXNS_BY_SPECIES(L,K,2,ID1:4)
                      TOP_RXNS_BY_SPECIES(L,K,2,ID1) = J
                      EXIT
                   END IF
                END DO
                IF (REACTANTS(1) .EQ. REACTANTS(2)) THEN
                   TOTAL_CHANGE_BY_SPECIES(L,K) = TOTAL_CHANGE_BY_SPECIES(L,K) - 2.0_LDP*NEW_MOL_RXN_RATES_ALL(L,J)
                ELSE
                   TOTAL_CHANGE_BY_SPECIES(L,K) = TOTAL_CHANGE_BY_SPECIES(L,K) - 1.0_LDP*NEW_MOL_RXN_RATES_ALL(L,J)
                END IF
             END IF
             IF (SPECIESID1 .EQ. PRODUCTS(1) .OR. SPECIESID1 .EQ. PRODUCTS(2)) THEN
                DO ID1=1,5
                   ID2 = TOP_RXNS_BY_SPECIES(L,K,1,ID1)
                   IF(ABS(NEW_MOL_RXN_RATES_ALL(L,J)) .GT. ABS(NEW_MOL_RXN_RATES_ALL(L,ID2))) THEN
                      IF (ID1 .LT. 5) TOP_RXNS_BY_SPECIES(L,K,1,ID1+1:5) = TOP_RXNS_BY_SPECIES(L,K,1,ID1:4)
                      TOP_RXNS_BY_SPECIES(L,K,1,ID1) = J
                      EXIT
                   END IF
                END DO
                IF (PRODUCTS(1) .EQ. PRODUCTS(2)) THEN
                   TOTAL_CHANGE_BY_SPECIES(L,K) = TOTAL_CHANGE_BY_SPECIES(L,K) + 2.0_LDP*NEW_MOL_RXN_RATES_ALL(L,J)
                ELSE
                   TOTAL_CHANGE_BY_SPECIES(L,K) = TOTAL_CHANGE_BY_SPECIES(L,K) + 1.0_LDP*NEW_MOL_RXN_RATES_ALL(L,J)
                END IF
             END IF
          END DO
       END DO
    END DO
!
! Now I need to determine certain properties for each of my considered species. I will do so by looping over the 
! list of molecular reactions
!
    IF (.NOT. ALLOCATED(SPECIES_RELEVANT_LEVELS)) ALLOCATE(SPECIES_RELEVANT_LEVELS(N_SPECIES,2))
    SPECIES_RELEVANT_LEVELS=1
    IF (.NOT. ALLOCATED(SPECIES_ISFINAL)) ALLOCATE(SPECIES_ISFINAL(N_SPECIES))
    SPECIES_ISFINAL=.TRUE.
!
    LOOPOVERSP: DO L=1,N_SPECIES
       LOOPOVERRXNS: DO J=1,N_MOL_RXN
          MATCHED=.FALSE.
          LOOPINRXN: DO I=1,4
             IF (MOL_RXN_SPECIESID(J,I) .EQ. SPECIESID_LIST(L)) THEN
                ID1=L
                ID2=I
                MATCHED=.TRUE.
                EXIT LOOPINRXN
             END IF
          END DO LOOPINRXN
          IF (MATCHED) THEN
             SPECIES_ISFINAL(L) = MOL_RXN_IS_FINAL_ION(J,ID2)
             IF (.NOT. MOL_RXN_IS_FINAL_ION(J,ID2)) THEN
                IF (LEV_IN_ION_MOL(J,ID2) .GT. SPECIES_RELEVANT_LEVELS(L,2)) THEN
                   SPECIES_RELEVANT_LEVELS(L,2) = LEV_IN_ION_MOL(J,ID2)
                END IF
             END IF
          END IF
       END DO LOOPOVERRXNS
    END DO LOOPOVERSP

!    DO J=1,N_MOL_RXN
!       DO I=1,4
!          DO L=1,14
!             IF (MOL_RXN_SPECIESID(J,I) .EQ. SPECIESIDLIST(L)) THEN
!                ID1=L
!                EXIT
!             END IF
!          END DO
!          IF (LEV_IN_ION_MOL(J,I) .GT. SPECIES_RELEVANT_LEVELS(ID1,2)) THEN
!             SPECIES_RELEVANT_LEVELS(ID1,2) = LEV_IN_ION_MOL(J,I)
!          END IF
!          SPECIES_ISFINAL(ID1) = MOL_RXN_IS_FINAL_ION(J,I)
!       END DO
!    END DO
!
    OUTFILE_NAME = 'MOL_RXNS_BY_SPECIES'
    ID1=69
    OPEN(UNIT=ID1,FILE=TRIM(OUTFILE_NAME),STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
    IF (IOS .NE. 0) WRITE(*,*) 'Error in print_mol_info: Unable to open file '//OUTFILE_NAME
!    
    WRITE(ID1,*) '***********************************************************************************************'
    WRITE(ID1,*) 'Molecular reaction info, final iteration'
    WRITE(ID1,*) 'Shown are fastest destruction reactions (left) and formation reactions (right) for each species'
    WRITE(ID1,*) 'Info shown is summed over all relevant (split) levels. For info on reactions to/from individual'
    WRITE(ID1,*) 'levels, see the file MOLECULAR_REACTIONS_OUT'
    WRITE(ID1,*) '***********************************************************************************************'
    WRITE(ID1,*)
!
    DO L=1,ND
       WRITE(ID1,*) '***********************************************************************************************'
       WRITE(ID1,'(A,I3,2X,A,F10.5,2X,A,ES15.5)') 'Depth:',L,'Temperature (10^4 K): ',T(L),'ED:',ED(L)
       DO K=1,N_SPECIES
          IF (MOL_LEVID_LIST(K) .NE. 0) THEN
             TR1=0
             DO I=1,N_LVLS_SPECIES(K)
                TR1=TR1+POPS(MOL_LEVID_LIST(K)+I-1,L)
             END DO
             WRITE(ID1,*) '-----------------------------------------------------------------------------------------------'
             WRITE(ID1,*)
             WRITE(ID1,'(A,A,ES15.5,A,A)') TRIM(SPECIESID_LIST(K)),'(POP=',TR1,')',':'
             WRITE(ID1,*)
             DO J=1,5
                ID2=TOP_RXNS_BY_SPECIES(L,K,2,J)
                ID3=TOP_RXNS_BY_SPECIES(L,K,1,J)
                !
                IF (ID2 .LE. N_RXNS_UNSPLIT) THEN
                   ID2R=REVERSE_POINTERS(ID2)
                   WRITE(RXN_STRING1,'(A6,A1,1X,A6,A2,1X,A6,A1,1X,A6,A1,1X,ES15.5)') MOL_RXN_SPECIESID(ID2R,1),'+',&
                        MOL_RXN_SPECIESID(ID2R,2),'->',MOL_RXN_SPECIESID(ID2R,3),'+',MOL_RXN_SPECIESID(ID2R,4),':',&
                        NEW_MOL_RXN_RATES_ALL(L,ID2)
                ELSE
                   RXN_STRING1 = ' '
                END IF
                IF (ID3 .LE. N_RXNS_UNSPLIT) THEN
                   ID3R=REVERSE_POINTERS(ID3)
                   WRITE(RXN_STRING2,'(A6,A1,1X,A6,A2,1X,A6,A1,1X,A6,A1,1X,ES15.5)') MOL_RXN_SPECIESID(ID3R,1),'+',&
                        MOL_RXN_SPECIESID(ID3R,2),'->',MOL_RXN_SPECIESID(ID3R,3),'+',MOL_RXN_SPECIESID(ID3R,4),':',&
                        NEW_MOL_RXN_RATES_ALL(L,ID3)
                ELSE
                   RXN_STRING2 = ' '
                END IF
!
                WRITE(ID1,'(A,T50,A,2X,A)') TRIM(RXN_STRING1),'|',TRIM(RXN_STRING2)
!
             END DO
!         
             WRITE(ID1,*)
             WRITE(ID1,'(A,A,A,2X,ES15.5,2X,A,2X,ES15.5)') 'Total change to ',TRIM(SPECIESID_LIST(K)),&
                  ' from molecular reactions:',TOTAL_CHANGE_BY_SPECIES(L,K),'Scaled by pop:',&
                  TOTAL_CHANGE_BY_SPECIES(L,K)/POPS(MOL_LEVID_LIST(K),L)
             WRITE(ID1,*)
!
             IF (.NOT. SPECIES_ISFINAL(K)) THEN
                TOT_PR=0.0_LDP
                TOT_RR=0.0_LDP
                WRITE(ID1,*) 'Photoionization/Recombination rates for this species: '
                WRITE(ID1,*)
                DO I=1,ATM(MOL_ATMID_LIST(K))%NXzV
                   TOT_PR = TOT_PR + ATM(MOL_ATMID_LIST(K))%APRXzV(I,L)
                   TOT_RR = TOT_RR + ATM(MOL_ATMID_LIST(K))%ARRXzV(I,L)
                END DO
                WRITE(ID1,'(A,1X,ES15.5,T50,A,1X,ES15.5)') 'PR: ',TOT_PR,'RR: ',TOT_RR
             END IF
!
             WRITE(ID1,*)'-----------------------------------------------------------------------------------------------'
          END IF
       END DO
       WRITE(ID1,*)'***********************************************************************************************'
       WRITE(ID1,*)
    END DO
    !
    CALL FLUSH(ID1)
    CLOSE(ID1)
!
! Add rates to the atom variables (for printing in PRRR)
!
    DO I=1,N_SPECIES
       IF (MOL_ATMID_LIST(I) .NE. 0) THEN
          SPECIESID1=SPECIESID_LIST(I)
          ID1=MOL_ATMID_LIST(I)
          ATM(ID1)%MOL_DRXzV(:)=0.0_LDP
          ATM(ID1)%MOL_CRXzV(:)=0.0_LDP
          DO J=1,N_MOL_RXN
             DO K=1,4
                IF (MOL_RXN_SPECIESID(J,K) .EQ. SPECIESID1) THEN
                   IF (.NOT. MOL_RXN_IS_FINAL_ION(J,K)) THEN
                      DO L=1,ND
                         IF (K .EQ. 1 .OR. K .EQ. 2) THEN
                            ATM(ID1)%MOL_DRXzV(L)=ATM(ID1)%MOL_DRXzV(L)+MOL_RXN_RATES_ALL(L,J)
                         ELSE IF (K .EQ. 3 .OR. K .EQ. 4) THEN
                            ATM(ID1)%MOL_CRXzV(L)=ATM(ID1)%MOL_CRXzV(L)+MOL_RXN_RATES_ALL(L,J)
                         END IF
                      END DO
                   END IF
                END IF
             END DO
          END DO
       END IF
    END DO         
!
  END SUBROUTINE PRINT_MOL_INFO
