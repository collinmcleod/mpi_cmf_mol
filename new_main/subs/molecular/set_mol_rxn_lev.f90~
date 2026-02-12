!
! Subroutine to link SPECIES identifiers in molecular reactions with the main program variables in CMFGEN
!
! Fewer rules than charge exchange reaction files
!
! Created 26-Aug-2021 by Collin McLeod
! Edited 10-Jan-2025 by Collin McLeod - altered flow so that excluded reactions are actually fully excluded
! Edited 26-Aug-2021 by Collin McLeod
! Edited 22-Sep-2021: Changed treatment of ELECTRON and NOATOM occurrences in reactions

      SUBROUTINE SET_MOL_RXN_LEV(ND,LUOUT)
      USE MOD_CMFGEN
      USE MOL_RXN_MOD
      USE SET_KIND_MODULE
      USE CONTROL_VARIABLE_MOD, ONLY : INCL_MOL_COMP
      IMPLICIT NONE

!
! Variables -- see set_chg_lev_id_v4.f
      INTEGER ND,LUOUT
!
      INTEGER, PARAMETER :: IZERO=0
!
      INTEGER ID
      INTEGER, ALLOCATABLE :: ID_POINTER(:,:)
      INTEGER, ALLOCATABLE :: LEV_CNT(:,:)
      INTEGER, ALLOCATABLE :: MOL_ID(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: G_SUM(:,:)
!
      INTEGER, ALLOCATABLE :: TMP_ID_ION(:,:)
      INTEGER, ALLOCATABLE :: TMP_LEV_IN_POPS(:,:)
      INTEGER, ALLOCATABLE :: TMP_LEV_IN_ION(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: TMP_G(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: TMP_Z(:,:)
      INTEGER, ALLOCATABLE :: USED_IS_LIST(:)
!
      REAL(KIND=LDP) T1,T2
      INTEGER NOUT
      INTEGER NIN
      INTEGER I,J,K,L
      INTEGER L1,L2,L3,L4,ML
      INTEGER ICNT
      INTEGER ICOUNT
      INTEGER IPOS
      INTEGER I_S,I_F
      INTEGER LST
      INTEGER II
      INTEGER IP
      INTEGER IOS
      INTEGER N_RXN_OMITTED
      INTEGER V1,V2,J1,J2
      CHARACTER(30) LOC_NAME,NAME1,NAME2
      CHARACTER(132) STRING1,STRING2
      LOGICAL LEVEL_SET
!
      INTEGER STR_LEN
      CHARACTER(10),ALLOCATABLE :: TMP_SP_ID_LIST(:)
      CHARACTER(10),ALLOCATABLE :: TMP_SP_NAME_LIST(:)
!
      IF(.NOT. DO_MOL_RXNS) RETURN
!
      ALLOCATE(LEV_CNT(N_MOL_RXN_RD,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(G_SUM(N_MOL_RXN_RD,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(ID_POINTER(N_MOL_RXN_RD,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(USED_IS_LIST(450),STAT=IOS)
      IF (IOS .NE. 0) THEN
         WRITE(LUER1,*)'Error in SET_MOL_RXN_LEV'
         WRITE(LUER1,*)'Unable to allocate arrays'
         WRITE(LUER1,*)'STAT= ',IOS
         STOP
      END IF
!
!     SPECIESID_LIST(1) = 'CI'
!     SPECIESID_LIST(2) = 'C2'
!     SPECIESID_LIST(3) = 'OI'
!     SPECIESID_LIST(4) = 'O2'
!     SPECIESID_LIST(5) = 'COMI'
!     SPECIESID_LIST(6) = 'COM2'
!     SPECIESID_LIST(7) = 'C2MI'
!     SPECIESID_LIST(8) = 'C2M2'
!     SPECIESID_LIST(9) = 'O2MI'
!     SPECIESID_LIST(10) = 'O2M2'
!     SPECIESID_LIST(11) = 'CO2MI'
!     SPECIESID_LIST(12) = 'CO2M2'
!
      !Initialize everything to zero
      LEV_CNT(:,:)=0
      G_SUM(:,:)=0.0_LDP
      ID_POINTER(:,:)=0
      USED_IS_LIST(:)=0
!
      !Allocate space for arrays of important information
!      
      IF (.NOT. ALLOCATED(MOL_RXN_IS_FINAL_ION_RD))ALLOCATE(MOL_RXN_IS_FINAL_ION_RD(N_MOL_RXN_RD,4),STAT=IOS)
      IF (IOS .NE. 0) THEN
         WRITE(LUER1,*)'Error in SET_MOL_RXN_LEV'
         WRITE(LUER1,*)'Unable to allocate memory for is_final_ion'
         STOP
      END IF
           
      !Count the levels and account for super levels

      DO ID=1,NUM_IONS-1
         IF(ATM(ID)%XzV_PRES) THEN
            DO J=1,N_MOL_RXN_RD
               DO K=1,4
                  IF(TRIM(ADJUSTL(MOL_RXN_SPECIESID_RD(J,K))) .EQ. TRIM(ADJUSTL(ION_ID(ID)))) THEN !Species ID match
                     ID_POINTER(J,K)=ID
                     MOL_RXN_IS_FINAL_ION_RD(J,K) = .FALSE.

! Now check against all (non-super) levels

                     I_S=0
                     L=1
                     USED_IS_LIST(:)=0
                     DO I_F=1,ATM(ID)%NXzV_F
                        LOC_NAME=ATM(ID)%XzVLEVNAME_F(I_F)
                        IF(INDEX(MOL_RXN_LEV_NAME_RD(J,K),'|') .EQ. 0) THEN
                           IPOS=INDEX(LOC_NAME,'|')
                           IF(IPOS .NE. 0) LOC_NAME=LOC_NAME(1:IPOS-1)
                        END IF
!*
                        IF(INDEX(MOL_RXN_LEV_NAME_RD(J,K),'[') .EQ. 0) THEN      !repeat above condition for 
                           IPOS=INDEX(LOC_NAME,'[')                              !different delimiter (used in CO)
                           IF(IPOS .NE. 0) LOC_NAME=LOC_NAME(1:IPOS-1)           !I have changed the naming convention for CO-make sure to check for "|" first
                        ENDIF
!*
                        IF (INDEX(MOL_RXN_LEV_NAME_RD(J,K),'Rovibration') .NE. 0 .OR. INDEX(MOL_RXN_LEV_NAME_RD(J,K),'X1Sigma') &
                             .NE. 0) THEN !get the v and J for rovibrational levels--should work for CO and SiO
                           CALL VJ_FROM_NAME(ATM(ID)%XzVLEVNAME_F(I_F),V1,J1)
                        ELSE
                           J1=0
                        END IF
!
                        IF(LOC_NAME .EQ. MOL_RXN_LEV_NAME_RD(J,K)) THEN
                           IF(I_S .EQ. 0) THEN
                              I_S=ATM(ID)%F_TO_S_XZv(I_F)
                              USED_IS_LIST(L)=I_S
                              L=L+1
                              LEV_CNT(J,K)=1
                           ELSE IF(ANY(USED_IS_LIST .EQ. ATM(ID)%F_TO_S_XZv(I_F))) THEN
                              CONTINUE !Same super level, no change needed
!                          ELSE IF (J1 .NE. 0) THEN
!                             CONTINUE
                           ELSE
                              USED_IS_LIST(L)=ATM(ID)%F_TO_S_XzV(I_F)
                              L=L+1
                              LEV_CNT(J,K)=LEV_CNT(J,K)+1 !Different super level from previous
                           END IF
                           G_SUM(J,K)=G_SUM(J,K)+ATM(ID)%GXzV_F(I_F)
                        END IF
                     END DO
                  ELSE IF (INDEX(MOL_RXN_SPECIESID_RD(J,K),'M2') .NE. 0) THEN !Is this a +1 ion
                        STR_LEN = LEN_TRIM(MOL_RXN_SPECIESID_RD(J,K))
                        IF (MOL_RXN_SPECIESID_RD(J,K)(STR_LEN:STR_LEN) .EQ. '2') THEN !This means it's a +1 ion
                        !Takes two steps because I can't rule out this being CO2M2, for instance
                           IF (MOL_RXN_SPECIESID_RD(J,K)(1:STR_LEN-1)//'I' .EQ. ION_ID(ID)) THEN !Match species
                              ID_POINTER(J,K)=ID
                              MOL_RXN_IS_FINAL_ION_RD(J,K)=.TRUE.
                              LEV_CNT(J,K)=1
                              G_SUM(J,K)=ATM(ID)%GIONXzV_F
                           END IF !Matches previous ion
                        END IF !Is a +1 ion
                     END IF !SPECIESID has a 2 in it
                  
               END DO !Over reaction variables
            END DO !Over number of molecular reactions
         END IF
      END DO !Over all ion IDs

!Add a loop to check for reactions which directly involve electrons or protons (H2, which is not handled correctly otherwise)

      DO J=1,N_MOL_RXN_RD
         DO K=1,4
            IF ((MOL_RXN_SPECIESID_RD(J,K) .EQ. 'ELEC') .AND. (MOL_RXN_LEV_NAME_RD(J,K) .EQ. 'ELECTRON')) THEN
               ID_POINTER(J,K)=0 !Placeholder value which indicates electrons
               MOL_RXN_IS_FINAL_ION_RD(J,K) = .FALSE.
               LEV_CNT(J,K)=1 !Free electron
               G_SUM(J,K)=1 !Free electron
            ELSE IF (MOL_RXN_SPECIESID_RD(J,K) .EQ. 'H2' .AND. (MOL_RXN_LEV_NAME_RD(J,K) .EQ. 'ION') ) THEN
               ID_POINTER(J,K)=1
               MOL_RXN_IS_FINAL_ION_RD(J,K) = .TRUE.
               LEV_CNT(J,K)=1
               G_SUM(J,K)=ATM(1)%GIONXzV_F
            END IF
         END DO
      END DO


! A loop which checks for reactions with only a single product (mostly radiative association reactions)

      DO J=1,N_MOL_RXN_RD
         DO K=1,4
            IF((MOL_RXN_SPECIESID_RD(J,K) .EQ. 'NOATOM') .AND. (MOL_RXN_LEV_NAME_RD(J,K) .EQ. 'NOLEVEL')) THEN
               ID_POINTER(J,K)=-1 !Placeholder value to denote nothing
               MOL_RXN_IS_FINAL_ION_RD(J,K) = .FALSE.
               LEV_CNT(J,K)=1 !Free photon
               G_SUM(J,K)=1 !Free photon
            END IF
         END DO
      END DO
!
      !Count total number of included reactions (depending on which ions are included in main program run)
!

      WRITE(LUER1,*) "INCL_MOL_COMP: ",INCL_MOL_COMP
      N_MOL_RXN=0
      N_RXN_OMITTED=0
      DO J=1,N_MOL_RXN_RD
         K=LEV_CNT(J,1)*LEV_CNT(J,2)*LEV_CNT(J,3)*LEV_CNT(J,4)
         IF (K .EQ. 0) THEN
            N_RXN_OMITTED = N_RXN_OMITTED+1
            MOL_RXN_INCL_RD(J)=.FALSE.
         ELSE IF ( (MOL_RXN_TYPE_RD(J) .EQ. 3) ) THEN
            IF (DO_NT .AND. INCL_MOL_COMP) THEN
               N_MOL_RXN = N_MOL_RXN + K
            ELSE
               N_RXN_OMITTED = N_RXN_OMITTED+1
               MOL_RXN_INCL_RD(J)=.FALSE.
               LEV_CNT(J,:) = 0
            END IF
         ELSE
            N_MOL_RXN = N_MOL_RXN+K
         END IF
      END DO
!
! Report how many reactions are included/excluded
      WRITE(LUER1,*)'Number of molecular reactions read in is ',N_MOL_RXN_RD
      WRITE(LUER1,*)'Number of molecular reactions omitted is ',N_RXN_OMITTED
      WRITE(LUER1,*)'Total molecular reactions (including split states) is ',N_MOL_RXN
      CALL FLUSH(LUER1)
!
      ALLOCATE(MOL_RXN_TYPE(N_MOL_RXN),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_RATE_SCALE(N_MOL_RXN),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_COEF(N_MOL_RXN,N_COEF_MAX_MOL),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_G(N_MOL_RXN,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_TLO(N_MOL_RXN),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_THI(N_MOL_RXN),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_INCL(N_MOL_RXN),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_SPECIESID(N_MOL_RXN,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_LEV_NAME(N_MOL_RXN,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_Z_CHG(N_MOL_RXN,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_REACTION_AVAILABLE(N_MOL_RXN),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(ID_ION_MOL(N_MOL_RXN,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(LEV_IN_POPS_MOL(N_MOL_RXN,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(LEV_IN_ION_MOL(N_MOL_RXN,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(AI_AR_MOL(N_MOL_RXN,ND),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(DLNAI_AR_MOL_DLNT(N_MOL_RXN,ND),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(COOL_MOL(N_MOL_RXN,ND),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_ID(N_MOL_RXN,ND),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_IS_FINAL_ION(N_MOL_RXN,4),STAT=IOS)
!
      IF (IOS .NE. 0)THEN
         WRITE(LUER1,*)'Error in SET_MOL_RXN_LEV'
         WRITE(LUER1,*)'Unable to allocate main memory'
         WRITE(LUER1,*)'STAT= ',IOS
         STOP
      END IF
!
      WRITE(LUER1,*)'Allocated main memory for molecular reactions'
!
      !Initialize things to zero
      ID_ION_MOL(:,:)=0
      LEV_IN_POPS_MOL(:,:)=0
      LEV_IN_ION_MOL(:,:)=0
      MOL_Z_CHG(:,:)=0.0_LDP
      AI_AR_MOL(:,:)=0.0_LDP
      DLNAI_AR_MOL_DLNT(:,:)=0.0_LDP
      COOL_MOL(:,:)=0.0_LDP
      INITIALIZE_ARRAYS_MOL=.TRUE.
!
! Allocate temporary arrays to be used in this routine
!
      IOS=0
      I=MAXVAL(LEV_CNT)
      WRITE(LUER1,*)'Maximum level count is ',I
      IF(IOS .EQ. 0) ALLOCATE(TMP_ID_ION(I,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(TMP_LEV_IN_POPS(I,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(TMP_LEV_IN_ION(I,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(TMP_G(I,4),STAT=IOS)
      IF(IOS .EQ. 0) ALLOCATE(TMP_Z(I,4),STAT=IOS)
      IF (IOS .NE. 0)THEN
         WRITE(LUER1,*)'Error in SET_MOL_RXN_LEV'
         WRITE(LUER1,*)'Unable to allocate memory for temporary arrays; STAT=',IOS
         STOP
      END IF
!
!Determine ionization stages and levels for the involved species, as well as whether the species is 
!actually present in any of the molecular reactions
!
      LST=0
      DO J=1,N_MOL_RXN_RD
         WRITE(LUOUT,'(/,A,I3,4(2X,A))')'Operating on molecular reaction J=',J,MOL_RXN_SPECIESID_RD(J,1:4)
         WRITE(LUOUT,'(A,3X,4I5)')'Super levels associated with each species:',&
              LEV_CNT(J,1),LEV_CNT(J,2),LEV_CNT(J,3),LEV_CNT(J,4)
         IF (LEV_CNT(J,1)*LEV_CNT(J,2)*LEV_CNT(J,3)*LEV_CNT(J,4) .NE. 0) THEN !Check whether the reaction is included
            TMP_G = 0.0_LDP
            DO K=1,4
               ID=ID_POINTER(J,K)
               IF (MOL_RXN_IS_FINAL_ION_RD(J,K)) THEN
                  TMP_Z(1,K)=ATM(ID)%ZXzV
                  TMP_ID_ION(1,K)=ID+1
                  TMP_LEV_IN_POPS(1,K)=ATM(ID)%EQXzV+ATM(ID)%NXzV
                  TMP_LEV_IN_ION(1,K)=1
                  TMP_G(1,K)=ATM(ID)%GIONXzV_F
               ELSE IF (MOL_RXN_SPECIESID_RD(J,K) .EQ. 'ELEC') THEN
                  TMP_Z(1,K)=-1
                  TMP_ID_ION(1,K)=0
                  TMP_LEV_IN_POPS(1,K)=-1 !Placeholder values which will not actually be used in the calculations
                  TMP_LEV_IN_ION(1,K)=-1
                  TMP_G(1,K)=1
               ELSE IF (MOL_RXN_SPECIESID_RD(J,K) .EQ. 'NOATOM') THEN
                  TMP_Z(1,K)=0
                  TMP_ID_ION(1,K)=-1
                  TMP_LEV_IN_POPS(1,K)=-2 !Placeholder values that I really hope I've accounted for correctly
                  TMP_LEV_IN_ION(1,K)=-2
                  TMP_G(1,K)=1
               ELSE
                  ICNT=0
                  DO I_F=1,ATM(ID)%NXzV_F
                     LOC_NAME=ATM(ID)%XzVLEVNAME_F(I_F)
                     IF(INDEX(MOL_RXN_LEV_NAME_RD(J,K),'[') .EQ. 0) THEN
                        IPOS=INDEX(LOC_NAME,'[')
                        IF (IPOS .NE. 0) LOC_NAME=LOC_NAME(1:IPOS-1)
                     END IF
!*
                     IF(INDEX(MOL_RXN_LEV_NAME_RD(J,K),'|') .EQ. 0) THEN
                        IPOS=INDEX(LOC_NAME,'|')
                        IF (IPOS .NE. 0) LOC_NAME=LOC_NAME(1:IPOS-1)
                     END IF                     
!*
                     IF (LOC_NAME .EQ. MOL_RXN_LEV_NAME_RD(J,K)) THEN
                        I_S=ATM(ID)%F_TO_S_XzV(I_F)
                        ML=0
                        DO I=1,ICNT
                           IF (I_S .EQ. TMP_LEV_IN_ION(I,K)) THEN
                              ML=I
                              EXIT
                           END IF
                        END DO
                        IF (ML .EQ. 0) THEN
                           ICNT=ICNT+1
                           ML=ICNT
                        END IF
                        !
                        TMP_ID_ION(ML,K)=ID
                        TMP_LEV_IN_ION(ML,K)=I_S
                        TMP_LEV_IN_POPS(ML,K)=ATM(ID)%EQXzV+I_S-1
                        TMP_Z(ML,K)=ATM(ID)%ZXzV-1.0_LDP
                        TMP_G(ML,K)=TMP_G(ML,K)+ATM(ID)%GXzV_F(I_F)
                     END IF !Matching level names
                  END DO !Loop over levels in full ion
               END IF
            END DO !Loop over K (species in given reaction)
!
            L=LST
            DO L4=1,LEV_CNT(J,4)
               DO L3=1,LEV_CNT(J,3)
                  DO L2=1,LEV_CNT(J,2)
                     DO L1=1,LEV_CNT(J,1)
                        L=L+1
                        DO K=1,4
                           IF(K .EQ. 1) ML=L1
                           IF(K .EQ. 2) ML=L2
                           IF(K .EQ. 3) ML=L3
                           IF(K .EQ. 4) ML=L4
                           ID_ION_MOL(L,K)=TMP_ID_ION(ML,K)
                           LEV_IN_ION_MOL(L,K)=TMP_LEV_IN_ION(ML,K)
                           LEV_IN_POPS_MOL(L,K)=TMP_LEV_IN_POPS(ML,K)
                           MOL_Z_CHG(L,K)=TMP_Z(ML,K)
                           MOL_ID(L,K)=J
                           MOL_RXN_G(L,K)=TMP_G(ML,K)
                        END DO !Loop over species in reaction
                     END DO !Over levels in species 1
                  END DO !Over levels in species 2
               END DO !Over levels in species 3
            END DO !Over levels in species 4
!
! Now save the reaction information from what was initially read in
!
            DO K=LST+1,L
               MOL_RXN_TYPE(K)=MOL_RXN_TYPE_RD(J)
               MOL_RXN_TLO(K)=MOL_RXN_TLO_RD(J)
               MOL_RXN_THI(K)=MOL_RXN_THI_RD(J)
               MOL_RXN_INCL(K)=MOL_RXN_INCL_RD(J)
               MOL_RXN_SPECIESID(K,1:4)=MOL_RXN_SPECIESID_RD(J,1:4)
               MOL_RXN_IS_FINAL_ION(K,1:4)=MOL_RXN_IS_FINAL_ION_RD(J,1:4)
               MOL_RXN_LEV_NAME(K,1:4)=MOL_RXN_LEV_NAME_RD(J,1:4)
               MOL_RXN_COEF(K,1:N_COEF_MAX_MOL)=MOL_RXN_COEF_RD(J,1:N_COEF_MAX_MOL)
               IF (MOL_RXN_TYPE_RD(J) .EQ. 3) THEN
                  MOL_RXN_RATE_SCALE(K) = MOL_RXN_G(K,3)*MOL_RXN_G(K,4)/G_SUM(J,3)/G_SUM(J,4)
               ELSE
                  MOL_RXN_COEF(K,1)=MOL_RXN_COEF(K,1)*MOL_RXN_G(K,3)*MOL_RXN_G(K,4)/G_SUM(J,3)/G_SUM(J,4)
                  MOL_RXN_RATE_SCALE(K)=1.0_LDP
               END IF
!
! Now correct for the fact that we have split some reactions based on super levels in the products
! Since my rates do not depend on the products (I'm not calculating reverse rates), we need to avoid
! overcounting the contribution from reactions which have multiple product states. I will assume that
! all such reactions go to all available product states with probability proportional to the statistical weight
!               
!               MOL_RXN_RATE_SCALE(K)=1.0/(LEV_CNT(J,3)*LEV_CNT(J,4))
            END DO
            LST=L
!
         END IF
      END DO
!
! Now check that the split levels included in the reactions correspond to a single LS state
!
      WRITE(LUOUT,'(A)')
      WRITE(LUOUT,'(A)')
      WRITE(LUOUT,'(A)')'Summary of molecular reaction info'
      DO J=1,N_MOL_RXN
         WRITE(LUOUT,'(/,4(A),7X,A,3X,A,8X,A,3X,A)')'   Type','     ID','  Level','   Lpop','z','MOL_ID','g','Species'
         WRITE(LUOUT,'(4(3X,I4),3X,F5.2,5X,I4,2X,F7.2,3X,A)') &
              (MOL_RXN_TYPE(J),ID_ION_MOL(J,K),LEV_IN_ION_MOL(J,K),LEV_IN_POPS_MOL(J,K),MOL_Z_CHG(J,K),MOL_ID(J,K),&
               MOL_RXN_G(J,K),TRIM(MOL_RXN_SPECIESID(J,K)),K=1,4)
      END DO
      CLOSE(LUOUT)
!
      CALL GEN_ASCI_OPEN(LUOUT,'MOL_RXN_CHK','UNKNOWN',' ','WRITE',IZERO,IOS)
      IF (IOS .NE. 0) THEN
         WRITE(LUER1,*) 'Unable to open LUOUT for the reaction check'
         WRITE(LUER1,*) 'IOS=',IOS
      END IF
!
      DO J=1,N_MOL_RXN
         STRING1=' '
         DO K=1,4
            ID=ID_ION_MOL(J,K)
            IF ((ID .GT. 0) .AND. (MOL_RXN_IS_FINAL_ION(J,K) .EQV. .FALSE.)) THEN
               IF(ATM(ID)%XzV_PRES) THEN
                  I_S=LEV_IN_ION_MOL(J,K)
                  NAME1=' '
                  NAME2=' '
                  DO I_F=1,ATM(ID)%NXzV_F
                     IF(ATM(ID)%F_TO_S_XzV(I_F) .EQ. I_S) THEN
                        IF (NAME1 .EQ. ' ') THEN
                           NAME1=ATM(ID)%XzVLEVNAME_F(I_F)
                           L=LEN_TRIM(STRING1)
                           IF(K .NE. 1)L=MAX(L,(K-1)*25+5*MOD(K+1,2))
                           STRING1(L+1:)=TRIM(MOL_RXN_SPECIESID(J,K))//'{'//TRIM(NAME1)//'}'
                           IPOS=INDEX(NAME1,'[')
                           IF(IPOS .NE. 0) THEN 
                              NAME1=NAME1(1:IPOS-1)
                           ELSE
                              IPOS=INDEX(NAME1,'|')
                              IF (IPOS .NE. 0) NAME1=NAME1(1:IPOS-1)
                           END IF
!
                           IF (ATM(ID)%NXzV .EQ. 1) THEN !To address the LTE scenario
                              IPOS=INDEX(NAME1,'|')
                              IF (IPOS .NE. 0) NAME1=NAME1(1:IPOS-1)
                           END IF
                              
!*
!                           IPOS=INDEX(NAME1,'|')
!*
!                           IF(IPOS .NE. 0)NAME1=NAME1(1:IPOS-1)
                        ELSE
                           NAME2=ATM(ID)%XzVLEVNAME_F(I_F)
                           IPOS=INDEX(NAME2,'[')
                           IF(IPOS .NE. 0) THEN 
                              NAME2=NAME2(1:IPOS-1)
                           ELSE
                              IPOS=INDEX(NAME2,'|')
                              IF (IPOS .NE. 0) NAME2=NAME2(1:IPOS-1)
                           END IF
!
                           IF (ATM(ID)%NXzV .EQ. 1) THEN !To address the LTE scenario
                              IPOS=INDEX(NAME2,'|')
                              IF (IPOS .NE. 0) NAME2=NAME2(1:IPOS-1)
                           END IF
!*
!                           IPOS=INDEX(NAME2,'|')
!*
!                           IF(IPOS .NE. 0)NAME2=NAME2(1:IPOS-1)
                           IF(NAME2 .NE. NAME1) THEN
                              WRITE(LUER1,*)'Error in SET_MOL_RXN_LEV'
                              WRITE(LUER1,*)'Level-superlevel mismatch in reaction'
                              WRITE(LUER1,*)J,K,ID,I_S
                              WRITE(LUER1,*)'NAME1=',NAME1
                              WRITE(LUER1,*)'NAME2=',NAME2
                              STOP
                           END IF
                           
                           STRING1(L+1:)=TRIM(MOL_RXN_SPECIESID(J,K))//'{'//TRIM(NAME2)//'}'
                        END IF
                     END IF
                  END DO
                  IF (NAME1 .EQ. ' ') THEN
                     WRITE(LUER1,*)'Error in SET_MOL_RXN_LEV: No match for level name'
                     WRITE(LUER1,*)J,K,ID,I_S
                     STOP
                  END IF
               ELSE
                  L=LEN_TRIM(STRING1)
                  IF(K .NE. 1) L=MAX(L,(K-1)*25+5*MOD(K+1,2))
                  STRING1(L+1:)=TRIM(MOL_RXN_SPECIESID(J,K))//'{ion}'
               END IF
            END IF
         END DO
         WRITE(LUOUT,'(1X,I4,3X,A,3X,2ES12.4)')J,TRIM(STRING1),MOL_RXN_COEF(J,1),MOL_RXN_COEF_RD(MOL_ID(J,1),1)
      END DO

!
! Set the number of species involved, as well as the list of species
!
      ALLOCATE(TMP_SP_ID_LIST(N_MOL_RXN))
      ALLOCATE(TMP_SP_NAME_LIST(N_MOL_RXN))
      TMP_SP_ID_LIST(:)='NONE'
      TMP_SP_NAME_LIST(:)='NONE'
      I=0
      DO J=1,N_MOL_RXN
         DO K=1,4
            DO L=1,N_MOL_RXN
               IF (MOL_RXN_SPECIESID(J,K) .EQ. TMP_SP_ID_LIST(L)) EXIT
               IF (TMP_SP_ID_LIST(L) .EQ. 'NONE') EXIT
            END DO
            IF (TMP_SP_ID_LIST(L) .EQ. 'NONE') THEN
               IF (MOL_RXN_SPECIESID(J,K) .EQ. 'NOATOM') THEN
                  CONTINUE
               ELSE IF (MOL_RXN_SPECIESID(J,K) .EQ. 'ELEC') THEN
                  CONTINUE
               ELSE
                  I=I+1
                  TMP_SP_ID_LIST(I) = MOL_RXN_SPECIESID(J,K)
                  DO ID=1,NUM_IONS-1
                     IF ( MOL_RXN_SPECIESID(J,K) .EQ. ION_ID(ID) ) THEN
                        IF (.NOT. ANY(TMP_SP_NAME_LIST .EQ. SPECIES(SPECIES_LNK(ID)))) THEN
                           TMP_SP_NAME_LIST(I) = SPECIES(SPECIES_LNK(ID))
                        END IF
                        EXIT
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END DO
      N_SPECIES=I
      ALLOCATE(SPECIESID_LIST(N_SPECIES))
      ALLOCATE(SPECIESNAME_LIST(N_SPECIES))
      SPECIESID_LIST(1:N_SPECIES)=TMP_SP_ID_LIST(1:N_SPECIES)
      SPECIESNAME_LIST(1:N_SPECIES)=TMP_SP_NAME_LIST(1:N_SPECIES)
      ALLOCATE(N_LVLS_SPECIES(N_SPECIES))
      DEALLOCATE(TMP_SP_ID_LIST)
      DEALLOCATE(TMP_SP_NAME_LIST)
      WRITE(*,*) 'IN SET_MOL_RXN, HERE ARE THE SPECIES NAME LISTS:'
      WRITE(*,*) 'ID LIST:'
      DO I=1,N_SPECIES
         WRITE(*,*) SPECIESID_LIST(I)
      END DO
      WRITE(*,*) 'NAME LIST:'
      DO I=1,N_SPECIES
         WRITE(*,*) SPECIESNAME_LIST(I)
      END DO
      CALL FLUSH(6)
!
! Now set the number of levels relevant for molecules
!
      N_MOL_LEVS=0
      DO I=1,N_SPECIES
         RXNLOOP: DO J=1,N_MOL_RXN_RD
            IF (PRODUCT(LEV_CNT(J,:)) .NE. 0) THEN
               DO K=1,4
                  IF (MOL_RXN_SPECIESID_RD(J,K) .EQ. SPECIESID_LIST(I)) THEN
                     N_MOL_LEVS = N_MOL_LEVS + LEV_CNT(J,K)
!                     WRITE(LUER1,*) 'n_mol_levs in set_mol_rxn_lev.f90: '
!                     WRITE(LUER1,*) 'J,K: ',J,K
!                     WRITE(LUER1,*) 'LEV_CNT: ',LEV_CNT(J,1:4)
!                     WRITE(LUER1,*) 'N_MOL_LEVS: ',N_MOL_LEVS
!                     WRITE(*,*) TRIM(ADJUSTL(MOL_RXN_SPECIESID_RD(J,1))),' ',TRIM(ADJUSTL(MOL_RXN_SPECIESID_RD(J,2))),' ',&
!                          TRIM(ADJUSTL(MOL_RXN_SPECIESID_RD(J,3))),' ',TRIM(ADJUSTL(MOL_RXN_SPECIESID_RD(J,4)))
                     N_LVLS_SPECIES(I)=LEV_CNT(J,K)
                     EXIT RXNLOOP
                  END IF
               END DO
            END IF
         END DO RXNLOOP
      END DO
!
! Write out some results of the level identifications
!
      CALL GEN_ASCI_OPEN(LUOUT,'MOL_RXN_RD_CHK','UNKNOWN',' ','WRITE',IZERO,IOS)
      CALL SET_LINE_BUFFERING(LUOUT)
      WRITE(LUOUT,*) 'Number of involved species=',N_SPECIES
      WRITE(LUOUT,*) 'List of involved species:'
      DO I=1,N_SPECIES
         WRITE(LUOUT,*) SPECIESID_LIST(I)
      END DO
      WRITE(LUOUT,*) 'Number of individual levels involved in molecular reactions: ',N_MOL_LEVS
!
      WRITE(LUOUT,'(/,A)')'LHS are included molecular reactions'
      WRITE(LUOUT,  '(A)')'RHS are excluded molecular reactions'
      WRITE(LUOUT,  '(A)')'A final ionization stage has same ID as the previous ionization stage'
      WRITE(LUOUT,  '(A)')'An ID of zero is used for electrons, and -1 for photons'
      WRITE(LUOUT,'(/,2(3X,A),3X,A,T20,4X,A,2X,A,3X,A)')'J','K','Species','ID','NLev','nGsum'
      DO J=1,N_MOL_RXN_RD
         L=LEV_CNT(J,1)*LEV_CNT(J,2)*LEV_CNT(J,3)*LEV_CNT(J,4)
         IF(L .NE. 0) THEN
            WRITE(LUOUT,'(A)')' '
            DO K=1,4
               WRITE(LUOUT,'(2I4,3X,A,T20,2I6,F10.1)')J,K,TRIM(MOL_RXN_SPECIESID_RD(J,K)),ID_POINTER(J,K)&
,LEV_CNT(J,K),G_SUM(J,K) !Included reactions
            END DO
         ELSE
            WRITE(LUOUT,'(A)')' '
            DO K=1,4
               WRITE(LUOUT,'(40X,2I4,3X,A,T60,2I6,F10.1)')J,K,TRIM(MOL_RXN_SPECIESID_RD(J,K)),ID_POINTER(J,K)&
,LEV_CNT(J,K),G_SUM(J,K) !Excluded reactions
            END DO
         END IF
      END DO

      IF (N_RXN_OMITTED .NE. 0) THEN
         WRITE(LUOUT,'(A)')
         WRITE(LUOUT,'(A)') 'The following reactions were not included: '
         DO J=1,N_MOL_RXN_RD
            IF (.NOT. MOL_RXN_INCL_RD(J)) THEN
               STRING1 = ' '
               DO K=1,4
                  STRING1(1+(K-1)*30:)=TRIM(MOL_RXN_SPECIESID_RD(J,K))//'{'//TRIM(MOL_RXN_LEV_NAME_RD(J,K))//'}'
               END DO
               WRITE(LUOUT,'(1X,A)') TRIM(STRING1)
            END IF
         END DO
      END IF
      
      CLOSE(LUOUT)
      
      !Deallocate local variables for memory

      DEALLOCATE(LEV_CNT)
      DEALLOCATE(MOL_ID)
      DEALLOCATE(ID_POINTER)
      DEALLOCATE(TMP_ID_ION,TMP_LEV_IN_POPS,TMP_LEV_IN_ION,TMP_G,TMP_Z)

      RETURN
      END SUBROUTINE SET_MOL_RXN_LEV
