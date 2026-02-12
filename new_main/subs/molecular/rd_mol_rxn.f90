!  A file to read in molecular reaction data from a file (originally MOL_RXNS_DATA)
!  Created August 24, 2021 by Collin McLeod
!
!  File to read in should contain the line 
!  NN    !Number of molecular reactions
!
!  Unlike charge exchange reactions, there are no strict rules on the order in which
!  reactants are listed in the data file
! 
!  Hopefully I can write decent code which will take care of the reaction quirks
!
!
!  Altered 25-Aug-2021: added final steps to read in data according to molecular reaction rules

      SUBROUTINE RD_MOL_RXN(LUIN,INCL_MOL_RXN,SCALE_FAC,CO_SCALE)
      USE MOL_RXN_MOD
      IMPLICIT NONE

      INTEGER LUIN
      LOGICAL INCL_MOL_RXN
      REAL(8) SCALE_FAC
      REAL(8) CO_SCALE

      INTEGER, PARAMETER :: IZERO=0
      INTEGER IOS
      INTEGER N_COEF
      INTEGER I,J,K,L,LB,RB
      INTEGER ERROR_LU
      EXTERNAL ERROR_LU
      CHARACTER*132 STRING1,STRING2
      CHARACTER*11  FORMAT_DATE

! Check whether we are including molecular reactions

      IF(INCL_MOL_RXN) THEN
         DO_MOL_RXNS=.TRUE.
         WRITE(*,*) 'In rd_mol_rxn: DO_MOL_RXNS == TRUE'
         WRITE(*,*) 'INCL_MOL_RXN: ',INCL_MOL_RXN
      ELSE
         DO_MOL_RXNS=.FALSE.
         N_MOL_RXN=0
         RETURN
      END IF
!
! Set the scale factor for reaction rates
!
      RXN_SCALE_FAC=SCALE_FAC
      CO_SCALE_FAC = CO_SCALE
!
! Making sure the file can be opened

      LUER1=ERROR_LU()
      WRITE(LUER1,*) "Beginning to read molecular reaction data"
      WRITE(LUER1,*) 'LUER1:',LUER1,'LUIN:',LUIN
      CALL GEN_ASCI_OPEN(LUIN,'MOL_RXNS_DATA','OLD',' ','READ',IZERO,IOS)
      IF (IOS .NE. 0) THEN
         WRITE(LUER1,*)"Unable to open MOL_RXNS_DATA in routine RD_MOL_RXN"
         WRITE(LUER1,*)'IOS=',IOS
         STOP
      END IF

! Now read in the number of reactions to include (skipping comments at the top of 
! MOL_RXNS_DATA

! Check for the string which indicates the number of reactions
      IOS=0
      L=0
      DO WHILE(L .EQ. 0 .AND. IOS .EQ. 0)
         READ(LUIN,'(A)',IOSTAT=IOS)STRING1
         L=INDEX(STRING1,"!Number of molecular reactions")
      END DO
      WRITE(LUER1,*) STRING1
! Write out any errors that occurred in that step
      IF(IOS .NE. 0) THEN
         WRITE(LUER1,*)'Error in RD_MOL_RXN'
         WRITE(LUER1,*)'Number of molecular reactions not found'
         STOP
      ELSE
         READ(STRING1(1:L),*) N_MOL_RXN_RD
      END IF

! Check the file date
      FORMAT_DATE=' '
      READ(LUIN,'(A)')STRING1
      IF(INDEX(STRING1,'!Format date') .NE. 0) THEN
         FORMAT_DATE=STRING1(1:11)
         READ(LUIN,'(A)')STRING2
         IF(INDEX(STRING2,'!Modification date') .EQ. 0) THEN
            WRITE(LUER1,*)'Error in RD_MOL_RXN'
            WRITE(LUER1,*)'Unable to read Modification date in MOL_RXNS_DATA'
            STOP
         END IF
      ELSE
         WRITE(LUER1,*)"Error in RD_MOL_RXN"
         WRITE(LUER1,*)"Unable to read Format date in MOL_RXNS_DATA"
         BACKSPACE(LUIN)
      END IF
! Backspace command for consistency with old format (following the design of rd_chg_exch

! Now allocate arrays which will store all the reaction data

      IF (.NOT. ALLOCATED(MOL_RXN_TYPE_RD)) THEN
         ALLOCATE(MOL_RXN_TYPE_RD(N_MOL_RXN_RD),STAT=IOS)
         IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_TLO_RD(N_MOL_RXN_RD),STAT=IOS)
         IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_THI_RD(N_MOL_RXN_RD),STAT=IOS)
         IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_INCL_RD(N_MOL_RXN_RD),STAT=IOS)
         IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_SPECIESID_RD(N_MOL_RXN_RD,4),STAT=IOS)
         IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_LEV_NAME_RD(N_MOL_RXN_RD,4),STAT=IOS)
         IF(IOS .EQ. 0) ALLOCATE(MOL_RXN_COEF_RD(N_MOL_RXN_RD,N_COEF_MAX_MOL),STAT=IOS)
      END IF

      WRITE(*,*) 'Allocated arrays in rd_mol_rxn'

! Check for errors in allocation
      IF(IOS .NE. 0) THEN
         WRITE(LUER1,*)"Error in RD_MOL_RXN"
         WRITE(LUER1,*)"Unable to allocate arrays properly"
         STOP
      END IF

! Now read in the reaction information in a big loop

      DO I=1,N_MOL_RXN_RD
         STRING1 = ' '
! Make sure the information can be read in properly
         DO WHILE(STRING1 .EQ. ' ' .OR. STRING1(1:1) .EQ. '!')
            READ(LUIN,'(A)',IOSTAT=IOS)STRING1
            IF (IOS .NE. 0) THEN
               WRITE(LUER1,'(A,I3)')'Error reading molecular reaction number ',I
               WRITE(LUER1,*)'Last reaction read was'
               WRITE(LUER1,*)TRIM(STRING2)
               STOP
            ENDIF
         END DO
         STRING2=STRING1

         DO K=1,4
            STRING1=ADJUSTL(STRING1)
            IF(STRING1 .EQ. ' ')THEN
               WRITE(LUER1,*)'Insufficient information in reaction string (1)'
               WRITE(LUER1,*)TRIM(STRING2)
            ENDIF
            
            L = INDEX(STRING1,' ')
            IF(L .NE. INDEX(STRING1,'  '))THEN
               WRITE(LUER1,*)'Reaction information must be separated by two spaces'
               WRITE(LUER1,*)TRIM(STRING2)
            ENDIF
            
! Start putting the reaction information into arrays, starting with the reaction string

            MOL_RXN_SPECIESID_RD(I,K)=STRING1(1:L-1)
            STRING1(1:)=STRING1(L+1:)
            IF(STRING1 .EQ. ' ')THEN
               WRITE(LUER1,*)'Insufficient information in reaction string (2)'
               WRITE(LUER1,*)TRIM(STRING2)
               STOP
            ENDIF
            STRING1=ADJUSTL(STRING1)

            L = INDEX(STRING1,' ')
            IF (L .NE. INDEX(STRING1,'  '))THEN
               WRITE(LUER1,*)'Error in RD_MOL_RXN'
               WRITE(LUER1,*)'Reaction information must be separated by two spaces'
               STOP
            ENDIF

! I only have one format (as of now) for reaction information, and no alternate level names, so 
! I have not included several steps that check for those issues in rd_chg_exch_v3.f

!           Read in the level name
            MOL_RXN_LEV_NAME_RD(I,K)=STRING1(1:L-1)
            
            STRING1(1:)=STRING1(L+1:)
         END DO

!        Read in the reaction type (I am currently using a single reaction formula for all molecular reactions
         READ(LUIN,*,IOSTAT=IOS)MOL_RXN_TYPE_RD(I)
         IF(IOS .NE. 0)THEN
            WRITE(LUER1,*)'Error in RD_MOL_RXN'
            WRITE(LUER1,*)'Unable to read reaction type for reaction number ',I
            WRITE(LUER1,*)TRIM(STRING2)
            STOP
         ENDIF
         
!        Read number of coefficients
         READ(LUIN,*,IOSTAT=IOS)N_COEF
         IF(IOS .NE. 0)THEN
            WRITE(LUER1,*)'Error in RD_MOL_RXN'
            WRITE(LUER1,*)'Unable to read number of coefficients for reaction number ',I
            WRITE(LUER1,*)TRIM(STRING2)
            STOP
         ENDIF

!        Check the number of coefficients
         IF(N_COEF .GT. N_COEF_MAX_MOL)THEN
            WRITE(LUER1,*)'Error in RD_MOL_RXN: N_COEF exceeds N_COEF_MAX in reaction ',I
            WRITE(LUER1,*)TRIM(STRING2)
            STOP
         ENDIF

!        Read in all the coefficients
         DO J=1,N_COEF
            READ(LUIN,*,IOSTAT=IOS)MOL_RXN_COEF_RD(I,J)
            IF(IOS .NE. 0)THEN
               WRITE(LUER1,*)'Error in RD_MOL_RXN'
               WRITE(LUER1,*)'Unable to read coefficients for reaction number ',I
               WRITE(LUER1,*)TRIM(STRING2)
               STOP
            ENDIF
         END DO

         !Read in the temperature limits
         READ(LUIN,*,IOSTAT=IOS)MOL_RXN_TLO_RD(I)
         IF(IOS .NE. 0)THEN
               WRITE(LUER1,*)'Error in RD_MOL_RXN'
               WRITE(LUER1,*)'Unable to read temperature limits for reaction number ',I
               WRITE(LUER1,*)TRIM(STRING2)
               STOP
         ENDIF

         READ(LUIN,*,IOSTAT=IOS)MOL_RXN_THI_RD(I)
         IF(IOS .NE. 0)THEN
               WRITE(LUER1,*)'Error in RD_MOL_RXN'
               WRITE(LUER1,*)'Unable to read temperature limits for reaction number ',I
               WRITE(LUER1,*)TRIM(STRING2)
               STOP
         ENDIF

      END DO

      !Set the code to include all reactions which have been read in
      MOL_RXN_INCL_RD(:) = .TRUE.
      WRITE(LUER1,*)'Molecular reaction data successfully read'
      CLOSE(LUIN)

      RETURN
      END SUBROUTINE RD_MOL_RXN
