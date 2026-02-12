PROGRAM HITRAN_CONVERT
  
      IMPLICIT NONE

!     Program to convert the HITRAN format data for Carbon Monoxide to
!     the CMFGEN format Will read data from Table_S5.par (from Li et
!     al., 2015) and write to a file called carb_monox_osc (name tbd)
!     Date created: Apr 6, 2021
!     Original author: Collin McLeod
!     Date edited: Apr 6, 2021
!     Basic outline of program format: declare
!     variables, create types for energy levels and transitions, read
!     chosen values for max v, max J.  Read data from Table_S5, select
!     out the chosen isotopologue, select out the transitions within max
!     v and J, determine the energy levels, allocate chosen transitions
!     and energy levels to arrays and sort them, print to file

!     A type which will contain all the data needed for the energy level list in the atomic files
      TYPE ENERGY_LEVEL
      CHARACTER(33) :: LEVEL_NAME
      REAL(KIND=8) :: LEVEL_G
      REAL(KIND=8) :: ENERGY_CM1
      REAL(KIND=8) :: FREQ_1015HZ
      REAL(KIND=8) :: ENERGY_BELOW_ION_EV
      REAL(KIND=8) :: LAMBDA_ANG
      INTEGER :: LEVEL_ID
      REAL(8) :: ARAD
      END TYPE


!     A type which will contain all the data needed for the transition list in the atomic files
      TYPE TRANSITION
      CHARACTER(28) :: LOWER_LEVEL_NAME
      CHARACTER(32) :: UPPER_LEVEL_NAME
      REAL(KIND=8) :: OSCILLATOR_STR
      REAL(KIND=8) :: EINSTEIN_A
      REAL(KIND=8) :: LAMBDA_ANG
      INTEGER :: LOWER_LEVEL_ID, UPPER_LEVEL_ID, TRANS_ID
      END TYPE
      
      INTEGER :: MAXV, MAXJ, NLEVELS, CHOSEN_ISO, IN_FILE_STAT, OUT_FILE_STAT, IN_FILE_NLINES
      INTEGER :: CHOOSE_ID
      CHARACTER(50) :: INFILENAME, OUTFILENAME, TEMPSTR1, TEMPSTR2, F_TO_S_FILE
      REAL(KIND=8), PARAMETER :: CM1TOEV=0.00012398470543514755D0
      REAL(KIND=8), PARAMETER :: PLANCKJOULE=6.62607015D-34
      REAL(KIND=8), PARAMETER :: PLANCKEV=4.135667696D-15
      REAL(KIND=8), PARAMETER :: LIGHTCM=29979245800.0D0 !cm/second
      REAL(KIND=8), PARAMETER :: LIGHTM=299792458.0D0 !Meters/second
      REAL(8), PARAMETER :: PI=3.14159265358979D0
      REAL(8), PARAMETER :: E_CHARGE=1.602176634D-19 !Coulombs
      REAL(8), PARAMETER :: EPS_0=8.8541878128D-12 !F/m, or N m^2 / C^2
      REAL(8), PARAMETER :: ELEC_MASS=9.1093837015D-31 !kilograms
      REAL(8), PARAMETER :: ION_ENERGY_EV=14.014D0

      REAL(8) :: MAX_LEVEL_ENERGY, MAX_LEVEL_G, TEMP1, TEMP2, TEMP3, TEMP4
      INTEGER :: TEMPINT1, TEMPINT2, TEMPINT3
      INTEGER, DIMENSION(1) :: TEMPINTARR

      TYPE(ENERGY_LEVEL), DIMENSION(:),  ALLOCATABLE :: ENERGY_LEVEL_LIST
      TYPE(TRANSITION), DIMENSION(:), ALLOCATABLE :: TRANSITION_LIST

      TYPE(ENERGY_LEVEL) :: TEMPLEV1, TEMPLEV2, TEMPLEV3

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: WAVENUMBER_LIST
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: EINSTEIN_A_LIST, AIR_BROADENING_LIST
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SELF_BROADENING_LIST, LOWER_STATE_ENERGY_LIST
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TEMP_DEPENDENCE_LIST, AIR_SHIFT_LIST, UPPER_G_LIST
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: LOWER_G_LIST, USED_ENERGIES
      INTEGER, DIMENSION(:), ALLOCATABLE :: MOLECULE_IDS, ISO_IDS, UPPER_V_LIST, LOWER_V_LIST,&
           LOWER_J_LIST
      CHARACTER(18), DIMENSION(:), ALLOCATABLE :: UNCERTAINTIES_AND_REFS_LIST
      CHARACTER(10), DIMENSION(:), ALLOCATABLE :: LINE_INTENSITY_LIST
      CHARACTER(21), DIMENSION(:), ALLOCATABLE :: BRANCH_LIST

      LOGICAL, DIMENSION(:), ALLOCATABLE :: ISO_MASK, VJMASK, LEVELMASK
      INTEGER :: NUM_CHOSEN_ISO, IOREAD1, COUNTER_I, VJCOUNT, COUNTER_J, COUNTER_K, COUNTER_L, COUNTER_M

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ISO_WAVENUMBERS, ISO_EINSTEIN_AS, ISO_LOWER_STATE_ENERGIES, &
           ISO_UPPER_GS, ISO_LOWER_GS
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISO_UPPER_VS, ISO_LOWER_VS, ISO_LOWER_JS, ISO_UPPER_JS
      CHARACTER(21), DIMENSION(:), ALLOCATABLE :: ISO_BRANCHES
      REAL(8) :: RED_FAC
      

      !Indicate the file which contains HITRAN data
      INFILENAME="Table_S5.par"
      CHOSEN_ISO = 1

      !Open the HITRAN file, read the number of lines, and allocate space for the arrays of input values
      OPEN(UNIT=1,FILE=INFILENAME,ACTION="READ",IOSTAT=IN_FILE_STAT)
      IF (IN_FILE_STAT>0) STOP "Problem opening HITRAN file"

      IN_FILE_NLINES=0
      DO
         READ(1,*,IOSTAT=IOREAD1)
         IF (IOREAD1/=0) EXIT
         IN_FILE_NLINES = IN_FILE_NLINES+1
      END DO
      REWIND(1)
      
      ALLOCATE(ISO_MASK(IN_FILE_NLINES))
      ALLOCATE(MOLECULE_IDS(IN_FILE_NLINES))
      ALLOCATE(ISO_IDS(IN_FILE_NLINES))
      ALLOCATE(WAVENUMBER_LIST(IN_FILE_NLINES))
      ALLOCATE(LINE_INTENSITY_LIST(IN_FILE_NLINES))
      ALLOCATE(EINSTEIN_A_LIST(IN_FILE_NLINES))
      ALLOCATE(AIR_BROADENING_LIST(IN_FILE_NLINES))
      ALLOCATE(SELF_BROADENING_LIST(IN_FILE_NLINES))
      ALLOCATE(LOWER_STATE_ENERGY_LIST(IN_FILE_NLINES))
      ALLOCATE(TEMP_DEPENDENCE_LIST(IN_FILE_NLINES))
      ALLOCATE(AIR_SHIFT_LIST(IN_FILE_NLINES))
      ALLOCATE(UPPER_V_LIST(IN_FILE_NLINES))
      ALLOCATE(LOWER_V_LIST(IN_FILE_NLINES))
      ALLOCATE(BRANCH_LIST(IN_FILE_NLINES))
      ALLOCATE(LOWER_J_LIST(IN_FILE_NLINES))
      ALLOCATE(UNCERTAINTIES_AND_REFS_LIST(IN_FILE_NLINES))
      ALLOCATE(UPPER_G_LIST(IN_FILE_NLINES))
      ALLOCATE(LOWER_G_LIST(IN_FILE_NLINES))

      !Read the lines in the HITRAN format, and sort them into the right arrays (allocated above)
      
      DO COUNTER_I = 1, IN_FILE_NLINES
         READ (1,'(I2,I1,F12.6,A10,F10.3,F5.4,F5.3,F10.4,F4.2,F8.6,I15,I15,A21,I9,A18,F7.1,F7.1)',IOSTAT=IOREAD1) &
              MOLECULE_IDS(COUNTER_I), ISO_IDS(COUNTER_I), WAVENUMBER_LIST(COUNTER_I), &
              LINE_INTENSITY_LIST(COUNTER_I), EINSTEIN_A_LIST(COUNTER_I), &
              AIR_BROADENING_LIST(COUNTER_I), SELF_BROADENING_LIST(COUNTER_I), &
              LOWER_STATE_ENERGY_LIST(COUNTER_I), TEMP_DEPENDENCE_LIST(COUNTER_I), &
              AIR_SHIFT_LIST(COUNTER_I), UPPER_V_LIST(COUNTER_I), LOWER_V_LIST(COUNTER_I), &
              BRANCH_LIST(COUNTER_I), LOWER_J_LIST(COUNTER_I), UNCERTAINTIES_AND_REFS_LIST(COUNTER_I),&
              UPPER_G_LIST(COUNTER_I), LOWER_G_LIST(COUNTER_I)
         IF (IOREAD1/=0) THEN
            WRITE(*,*) "Line Number: ", COUNTER_I
            STOP "Error reading variables from file"
         END IF
         ! ISO_MASK will record the line numbers which correspond to the correct isotopologue (12C16O)
         IF (ISO_IDS(COUNTER_I) .EQ. 1) THEN
            ISO_MASK(COUNTER_I)=.TRUE.
         ELSE
            ISO_MASK(COUNTER_I)=.FALSE.
         END IF
      END DO

      CLOSE(1)

      NUM_CHOSEN_ISO = COUNT(ISO_MASK)

      !Allocate space for arrays from the correct isotopologue, only accounting for the information I want to keep

      ALLOCATE(ISO_WAVENUMBERS(NUM_CHOSEN_ISO))
      ALLOCATE(ISO_EINSTEIN_AS(NUM_CHOSEN_ISO))
      ALLOCATE(ISO_LOWER_STATE_ENERGIES(NUM_CHOSEN_ISO))
      ALLOCATE(ISO_UPPER_GS(NUM_CHOSEN_ISO))
      ALLOCATE(ISO_LOWER_GS(NUM_CHOSEN_ISO))
      ALLOCATE(ISO_UPPER_VS(NUM_CHOSEN_ISO))
      ALLOCATE(ISO_LOWER_VS(NUM_CHOSEN_ISO))
      ALLOCATE(ISO_LOWER_JS(NUM_CHOSEN_ISO))
      ALLOCATE(ISO_BRANCHES(NUM_CHOSEN_ISO))

      !Now pick out the values I want from the full lists, using the mask I set up and the intrinsic function PACK

      ISO_WAVENUMBERS = PACK(WAVENUMBER_LIST,ISO_MASK)
      ISO_EINSTEIN_AS = PACK(EINSTEIN_A_LIST,ISO_MASK)
      ISO_LOWER_STATE_ENERGIES = PACK(LOWER_STATE_ENERGY_LIST,ISO_MASK)
      ISO_UPPER_GS = PACK(UPPER_G_LIST,ISO_MASK)
      ISO_LOWER_GS = PACK(LOWER_G_LIST,ISO_MASK)
      ISO_UPPER_VS = PACK(UPPER_V_LIST,ISO_MASK)
      ISO_LOWER_VS = PACK(LOWER_V_LIST,ISO_MASK)
      ISO_LOWER_JS = PACK(LOWER_J_LIST,ISO_MASK)
      ISO_BRANCHES = PACK(BRANCH_LIST,ISO_MASK)
      
!     Read in chosen values for the maximum v and J to include in the atomic files
      DO
         WRITE(*,*) "Enter values for maximum v (vibrational) and J (rotational) quantum numbers: "
         READ(*,*) MAXV, MAXJ
         IF (MAXV <= 41 .AND. MAXJ <=150) EXIT
         WRITE(*,*) "Max v must be less than 41, max J must be less than 150"
      END DO

      WRITE(*,*) "Pick factor to reduce emission rates: "
      READ(*,*) RED_FAC

      NLEVELS = (MAXV+1)*(MAXJ+1)

      !The values in Table_S5 are already ordered by wavenumber (smallest to largest). It'll be rather
      !simple to make an ordered list of transitions within the v and J limits


      !Select out the transitions which will be possible given the chosen maxima for v and J
      VJCOUNT=0
      
      ALLOCATE(VJMASK(NUM_CHOSEN_ISO))

      DO COUNTER_I=1, NUM_CHOSEN_ISO
         IF (ISO_UPPER_VS(COUNTER_I)<=MAXV .AND. ((ISO_LOWER_JS(COUNTER_I)<=MAXJ .AND. &
              TRIM(ADJUSTL(ISO_BRANCHES(COUNTER_I)))=="P") .OR. (ISO_LOWER_JS(COUNTER_I)<=(MAXJ-1) &
              .AND. TRIM(ADJUSTL(ISO_BRANCHES(COUNTER_I)))=="R"))) THEN
            VJMASK(COUNTER_I)=.TRUE.
            VJCOUNT = VJCOUNT + 1
         ELSE
            VJMASK(COUNTER_I)=.FALSE.
         END IF
      END DO

      VJCOUNT = COUNT(VJMASK)

      WRITE(*,*) "Max v: ", MAXV, "Max J: ", MAXJ, "Total number of possible transitions: ", VJCOUNT

      !Make lists of energy levels--first need to check for the highest energy level (v=MAXV, J=MAXJ)

      DO COUNTER_I=1, NUM_CHOSEN_ISO
         IF (ISO_LOWER_VS(COUNTER_I)==MAXV .AND. ISO_LOWER_JS(COUNTER_I)==MAXJ) THEN
            MAX_LEVEL_ENERGY=ISO_LOWER_STATE_ENERGIES(COUNTER_I)
            MAX_LEVEL_G = ISO_LOWER_GS(COUNTER_I)
            EXIT
         END IF
      END DO

      WRITE(*,*) "Number of energy levels: ", NLEVELS, "Max energy level: ", MAX_LEVEL_ENERGY, "cm^-1"

      ALLOCATE(ENERGY_LEVEL_LIST(NLEVELS))
      ALLOCATE(USED_ENERGIES(NLEVELS))
      
      !Loop over the included transitions, which includes information about lower state energies.
      !If that energy level has not yet been added to the list, add its information.
      !Exit the loop when the maximum number of levels is reached.
      !NOTE: THE ENERGY LEVEL LIST PRODUCED BY THIS LOOP IS NOT ORDERED
      COUNTER_J=1
      DO COUNTER_I=1, NUM_CHOSEN_ISO
         IF (COUNTER_J > NLEVELS) THEN
            EXIT
         ELSE IF ((COUNTER_J == NLEVELS) .AND. (NLEVELS<10)) THEN
            ENERGY_LEVEL_LIST(COUNTER_J)%LEVEL_NAME=LEVEL_NAME_VJ(0,0)
            ENERGY_LEVEL_LIST(COUNTER_J)%LEVEL_G=1.0
            ENERGY_LEVEL_LIST(COUNTER_J)%ENERGY_CM1=0.0D0
            TEMP1 = ION_ENERGY_EV
            ENERGY_LEVEL_LIST(COUNTER_J)%FREQ_1015HZ=(TEMP1/(PLANCKEV))*1D-15
            ENERGY_LEVEL_LIST(COUNTER_J)%ENERGY_BELOW_ION_EV=TEMP1
            ENERGY_LEVEL_LIST(COUNTER_J)%LAMBDA_ANG= ((PLANCKEV*LIGHTM)/(TEMP1))*1D10
            COUNTER_J=COUNTER_J+1

        !Above stuff was added to rectify an issue that was later rectified elsewhere, deprecated
            
         ELSE IF (.NOT. ANY(USED_ENERGIES==ISO_LOWER_STATE_ENERGIES(COUNTER_I)) .AND. &
              (ISO_LOWER_VS(COUNTER_I)<=MAXV) .AND. (ISO_LOWER_JS(COUNTER_I)<=MAXJ)) THEN
            ENERGY_LEVEL_LIST(COUNTER_J)%LEVEL_NAME=LEVEL_NAME_VJ(ISO_LOWER_VS(COUNTER_I),ISO_LOWER_JS(COUNTER_I))
            ENERGY_LEVEL_LIST(COUNTER_J)%LEVEL_G=ISO_LOWER_GS(COUNTER_I)
            ENERGY_LEVEL_LIST(COUNTER_J)%ENERGY_CM1=ISO_LOWER_STATE_ENERGIES(COUNTER_I)
            TEMP1 = ION_ENERGY_EV - (CM1TOEV*ISO_LOWER_STATE_ENERGIES(COUNTER_I))
            ENERGY_LEVEL_LIST(COUNTER_J)%FREQ_1015HZ=(TEMP1/(PLANCKEV))*1D-15
            ENERGY_LEVEL_LIST(COUNTER_J)%ENERGY_BELOW_ION_EV=TEMP1
            ENERGY_LEVEL_LIST(COUNTER_J)%LAMBDA_ANG= ((PLANCKEV*LIGHTM)/(TEMP1))*1D10
            USED_ENERGIES(COUNTER_J)=ISO_LOWER_STATE_ENERGIES(COUNTER_I)
            COUNTER_J=COUNTER_J+1
         END IF
      END DO
      
      !Sort the list by level energy, by iterative swapping
      DO COUNTER_I=1, NLEVELS-1
         TEMPINTARR=MINLOC(ENERGY_LEVEL_LIST(COUNTER_I:)%ENERGY_CM1)
         TEMPINT1=TEMPINTARR(1)
         TEMPLEV1=ENERGY_LEVEL_LIST(TEMPINT1+COUNTER_I-1)
         TEMPLEV2=ENERGY_LEVEL_LIST(COUNTER_I)
         ENERGY_LEVEL_LIST(COUNTER_I)=TEMPLEV1
         ENERGY_LEVEL_LIST(TEMPINT1+COUNTER_I-1)=TEMPLEV2
      END DO


      !Check that the energy level list is properly ordered
      TEMPINT2=0
      DO COUNTER_I=1, NLEVELS-1
         IF (ENERGY_LEVEL_LIST(COUNTER_I+1)%ENERGY_CM1<=ENERGY_LEVEL_LIST(COUNTER_I)%ENERGY_CM1) THEN
            TEMPINT2 = TEMPINT2+1
            TEMPINT3 = COUNTER_I
         END IF
      END DO

      WRITE(*,*) "Number of energy levels out of order: ", TEMPINT2, "Problem ID: ", TEMPINT3

      !Now index the sorted energy levels

      DO COUNTER_J=1,NLEVELS
         ENERGY_LEVEL_LIST(COUNTER_J)%LEVEL_ID=COUNTER_J
      END DO
!
      !Fix the bug involving the first level
      ENERGY_LEVEL_LIST(1)%LEVEL_NAME=LEVEL_NAME_VJ(0,0)
      ENERGY_LEVEL_LIST(1)%LEVEL_G=1.0
      ENERGY_LEVEL_LIST(1)%ENERGY_CM1=0.0D0
      TEMP1 = ION_ENERGY_EV
      ENERGY_LEVEL_LIST(1)%FREQ_1015HZ=(TEMP1/(PLANCKEV))*1D-15
      ENERGY_LEVEL_LIST(1)%ENERGY_BELOW_ION_EV=TEMP1
      ENERGY_LEVEL_LIST(1)%LAMBDA_ANG= ((PLANCKEV*LIGHTM)/(TEMP1))*1D10

      !Make the list of transitions
      !Will consist of nested loops: loop over the energy levels--loop over listed transitions. If the energy level in question is the lower level of the transition, then add it to the list.

      ALLOCATE(TRANSITION_LIST(VJCOUNT))

      COUNTER_L=1
      DO COUNTER_I=1, NLEVELS
         IF (COUNTER_L>VJCOUNT) EXIT
         TEMP4 = ENERGY_LEVEL_LIST(COUNTER_I)%ENERGY_CM1
         DO COUNTER_K=1,NUM_CHOSEN_ISO
            IF (COUNTER_L>VJCOUNT) EXIT
            IF ((ISO_LOWER_STATE_ENERGIES(COUNTER_K)==TEMP4) .AND. VJMASK(COUNTER_K)) THEN
               TRANSITION_LIST(COUNTER_L)%LOWER_LEVEL_NAME=ENERGY_LEVEL_LIST(COUNTER_I)%LEVEL_NAME
               TRANSITION_LIST(COUNTER_L)%LOWER_LEVEL_ID=COUNTER_I
               TRANSITION_LIST(COUNTER_L)%TRANS_ID=COUNTER_L
! Scale the A values to reduce the total CO emission
               TRANSITION_LIST(COUNTER_L)%EINSTEIN_A=RED_FAC*ISO_EINSTEIN_AS(COUNTER_K)
               TRANSITION_LIST(COUNTER_L)%OSCILLATOR_STR=RED_FAC*OSC_FROM_EINSTEIN(ISO_EINSTEIN_AS(COUNTER_K),&
                    ISO_LOWER_GS(COUNTER_K),ISO_UPPER_GS(COUNTER_K),ISO_WAVENUMBERS(COUNTER_K))
               TRANSITION_LIST(COUNTER_L)%LAMBDA_ANG=(1.0D0/ISO_WAVENUMBERS(COUNTER_K))*1D8
               SELECT CASE (TRIM(ADJUSTL(ISO_BRANCHES(COUNTER_K))))
                  CASE ("R")
               TRANSITION_LIST(COUNTER_L)%UPPER_LEVEL_NAME=LEVEL_NAME_VJ(ISO_UPPER_VS(COUNTER_K),&
                    ISO_LOWER_JS(COUNTER_K)+1)
                  CASE ("P")
               TRANSITION_LIST(COUNTER_L)%UPPER_LEVEL_NAME=LEVEL_NAME_VJ(ISO_UPPER_VS(COUNTER_K),&
                    ISO_LOWER_JS(COUNTER_K)-1)
               END SELECT
            
               DO COUNTER_M=1, NLEVELS
                  IF (ENERGY_LEVEL_LIST(COUNTER_M)%LEVEL_NAME==TRANSITION_LIST(COUNTER_L)%UPPER_LEVEL_NAME) THEN
                     TEMPINT3=COUNTER_M
                     EXIT
                  END IF
               END DO
               TRANSITION_LIST(COUNTER_L)%UPPER_LEVEL_ID=COUNTER_M
               COUNTER_L=COUNTER_L+1
            END IF
         END DO
      END DO

      !Determine ARAD for each energy level:
      
      DO COUNTER_I=1,NLEVELS
         TEMP1=0.0D0
         TEMPSTR1 = ENERGY_LEVEL_LIST(COUNTER_I)%LEVEL_NAME
         DO COUNTER_J=1,VJCOUNT
            IF (TRANSITION_LIST(COUNTER_J)%UPPER_LEVEL_NAME==TEMPSTR1) THEN
               TEMP1 = TEMP1 + TRANSITION_LIST(COUNTER_J)%EINSTEIN_A
            END IF
         END DO
         ENERGY_LEVEL_LIST(COUNTER_I)%ARAD=TEMP1
      END DO   

      !For checking accuracy of the level/transition list:

      WRITE(*,*) "Check accuracy? Enter 0 for no, 1 for levels, 2 for transitions"
      READ(*,*) TEMPINT1

      SELECT CASE (TEMPINT1)

      CASE (1)
         DO
          WRITE(*,*) "Pick a level id, 0 to exit: "
          READ(*,*) CHOOSE_ID
          IF (CHOOSE_ID==0) EXIT
          WRITE(*,*) "LEVEL_NAME: ", ENERGY_LEVEL_LIST(CHOOSE_ID)%LEVEL_NAME
          WRITE(*,*) "lEVEL_G: ", ENERGY_LEVEL_LIST(CHOOSE_ID)%LEVEL_G
          WRITE(*,*) "ENERGY_CM1: ", ENERGY_LEVEL_LIST(CHOOSE_ID)%ENERGY_CM1
          WRITE(*,*) "FREQ_1015HZ: ", ENERGY_LEVEL_LIST(CHOOSE_ID)%FREQ_1015HZ
          WRITE(*,*) "ENERGY_BELOW_ION_EV: ", ENERGY_LEVEL_LIST(CHOOSE_ID)%ENERGY_BELOW_ION_EV
          WRITE(*,*) "EINSTEIN_A: ", TRANSITION_LIST(CHOOSE_ID)%EINSTEIN_A
          WRITE(*,*) "LAMBDA_ANG: ", ENERGY_LEVEL_LIST(CHOOSE_ID)%LAMBDA_ANG
          WRITE(*,*) "LEVEL_ID: ", ENERGY_LEVEL_LIST(CHOOSE_ID)%LEVEL_ID
          WRITE(*,*) "ARAD: ", ENERGY_LEVEL_LIST(CHOOSE_ID)%ARAD
         END DO

      CASE (2)
       DO
        WRITE(*,*) "Pick a transition ID, 0 to exit: "
        READ(*,*) CHOOSE_ID
        IF (CHOOSE_ID==0) EXIT
        WRITE(*,*) "LOWER_LEVEL_NAME: ", TRANSITION_LIST(CHOOSE_ID)%LOWER_LEVEL_NAME
        WRITE(*,*) "UPPER_LEVEL_NAME: ", TRANSITION_LIST(CHOOSE_ID)%UPPER_LEVEL_NAME
        WRITE(*,*) "LOWER_LEVEL_ID: ", TRANSITION_LIST(CHOOSE_ID)%LOWER_LEVEL_ID
        WRITE(*,*) "UPPER_LEVEL_ID: ", TRANSITION_LIST(CHOOSE_ID)%UPPER_LEVEL_ID
        WRITE(*,*) "OSCILLATOR_STR: ", TRANSITION_LIST(CHOOSE_ID)%OSCILLATOR_STR
        WRITE(*,*) "EINSTEIN_A: ", TRANSITION_LIST(CHOOSE_ID)%EINSTEIN_A
        WRITE(*,*) "LAMBDA_ANG: ", TRANSITION_LIST(CHOOSE_ID)%LAMBDA_ANG
        WRITE(*,*) "TRANS_ID: ", TRANSITION_LIST(CHOOSE_ID)%TRANS_ID
       END DO

      CASE DEFAULT
       CONTINUE
    END SELECT

   !Default string for the name of the output file
    WRITE(OUTFILENAME,'(A,I2)') "COMI_OSC_RED",INT(100*RED_FAC)
    
   !Read in a name for the output file, or use a default string
    WRITE(*,'(A,A,A)',ADVANCE="NO") "Choose a name for the output file [", TRIM(ADJUSTL(OUTFILENAME)),"]: "
    READ(*,'(A)') TEMPSTR1
    IF (.NOT. LEN_TRIM(ADJUSTL(TEMPSTR1))==0) THEN
       OUTFILENAME=TEMPSTR1
    END IF

    !Open the file to write
    OPEN(UNIT=12,FILE=OUTFILENAME,STATUS='REPLACE',ACTION='READWRITE',IOSTAT=IOREAD1)

    !Write header information to the file

    WRITE(12,*) "*************************************************************************************"
    WRITE(12,*)
    WRITE(12,'(a,i3,a,i3)') "Ro-vibrational Energy Levels and Transitions for Carbon Monoxide, up to v=",MAXV," and J=", MAXJ
    WRITE(12,*)
    WRITE(12,*) "Data taken from Li et al., 2015.  Accessed 2021. URL: &
https://ui.adsabs.harvard.edu/abs/2015ApJS..216...15L/abstract"
    WRITE(12,*)

    WRITE(12,'(A33,A5,A15,A12,A9,A13,A8,A12,A12,A12)') "Level Name ", "g", "E(cm^-1)", "10^15 Hz", &
         "eV", "Lam(A)", "ID", "ARAD", "GAM2", "GAM4"
    WRITE(12,*) "*************************************************************************************"
    WRITE(12,*)

    TEMPSTR2 = FDATE()
    TEMPSTR1 = "15-Apr-2021"

    WRITE(12,'(A,T50,A)') TRIM(ADJUSTL('17-Oct-2000'))," !Format date"

    WRITE(12,'(A50,A)') ADJUSTL(TEMPSTR1), " !Date Created"

    WRITE(TEMPSTR1,*) NLEVELS
    WRITE(12,'(A50,A)') ADJUSTL(TEMPSTR1), " !Number of energy levels"

    WRITE(TEMPSTR1,'(F15.4)') ION_ENERGY_EV*8065.73D0
    WRITE(12,'(A50,A)') ADJUSTL(TEMPSTR1), " !Ionization energy"

!    WRITE(TEMPSTR1,*) "N/A"
    WRITE(12,'(F3.1,T50,A)') 1.0," !Screened nuclear charge"

    WRITE(TEMPSTR1,*) VJCOUNT
    WRITE(12,'(A50,A)') ADJUSTL(TEMPSTR1), " !Number of transitions"

    WRITE(12,*)

    !Loop over energy levels and print information in the desired format

    DO COUNTER_I=1, NLEVELS
       WRITE(12,'(A33,F5.1,F15.4,F12.5,F9.3,ES13.3,I8,ES12.3,ES12.3,ES12.3)') ADJUSTL(ENERGY_LEVEL_LIST(COUNTER_I)%LEVEL_NAME), &
            ENERGY_LEVEL_LIST(COUNTER_I)%LEVEL_G, ENERGY_LEVEL_LIST(COUNTER_I)%ENERGY_CM1, &
            ENERGY_LEVEL_LIST(COUNTER_I)%FREQ_1015HZ, ENERGY_LEVEL_LIST(COUNTER_I)%ENERGY_BELOW_ION_EV,&
            ENERGY_LEVEL_LIST(COUNTER_I)%LAMBDA_ANG, ENERGY_LEVEL_LIST(COUNTER_I)%LEVEL_ID, &
            ENERGY_LEVEL_LIST(COUNTER_I)%ARAD, 0.0, 0.0
    END DO

    !Write additional header-like information for the transitions

    WRITE(12,*)
    WRITE(12,*)
    WRITE(12,*) "****************************************************************************************"
    WRITE(12,*)
    WRITE(12,*) "oscillator strengths"
    WRITE(12,*)
    WRITE(12,*) "****************************************************************************************"
    WRITE(12,*)
    WRITE(12,'(A12,A55,A12,A16,A9,A12)') "Transition", "f", "A", "Lam(A)", "i-j", "Trans. #"
    WRITE(12,*)   

    !Loop over transition list and print information in the desired format

    DO COUNTER_I=1,VJCOUNT
       WRITE(12,'(A28,A1,A32,ES10.4,ES12.4,ES12.3,I8,A1,I4,I9)') ADJUSTL(TRANSITION_LIST(COUNTER_I)%LOWER_LEVEL_NAME), "-", &
            ADJUSTL(TRANSITION_LIST(COUNTER_I)%UPPER_LEVEL_NAME), TRANSITION_LIST(COUNTER_I)%OSCILLATOR_STR, &
            TRANSITION_LIST(COUNTER_I)%EINSTEIN_A, TRANSITION_LIST(COUNTER_I)%LAMBDA_ANG, &
            TRANSITION_LIST(COUNTER_I)%LOWER_LEVEL_ID, "-", TRANSITION_LIST(COUNTER_I)%UPPER_LEVEL_ID,&
            TRANSITION_LIST(COUNTER_I)%TRANS_ID
    END DO

    !Close the file and we're done
       
    CLOSE(12)

    CALL FLUSH

    !Now create the file with superlevel assignments, very similar format 
    !I'm assigning all CO levels to the same superlevel for purposes related to the reaction rates
!
    WRITE(F_TO_S_FILE,'(A,I2)')'COMI_F_TO_S_RED',INT(100*RED_FAC)
!
    OPEN(UNIT=13,FILE=F_TO_S_FILE,STATUS='REPLACE',ACTION='READWRITE',IOSTAT=IOREAD1)
    IF(IOREAD1 .NE. 0) THEN
       WRITE(*,*) 'Error opening superlevel assignment file'
       STOP
    END IF
!
    WRITE(13,*) "*************************************************************************************"
    WRITE(13,*)
    WRITE(13,'(a,i3,a,i3)') "Simple superlevel assignments for Carbon Monoxide, up to v=",MAXV," and J=", MAXJ
    WRITE(13,*)
    WRITE(13,*) "Data taken from Li et al., 2015.  Accessed 2021. URL: &
https://ui.adsabs.harvard.edu/abs/2015ApJS..216...15L/abstract"
    WRITE(13,*)

    WRITE(13,'(A33,A5,A15,A12,A9,A10)') "Level Name ", "g", "E(cm^-1)", "10^15 Hz", &
         "Lam(A)","SL"
    WRITE(13,*) "*************************************************************************************"
    WRITE(13,*)
!
    TEMPSTR1 = "25-Sep-2021"

    WRITE(13,'(A50,A)') ADJUSTL(TEMPSTR1), " !Date Created"

    WRITE(TEMPSTR1,*) NLEVELS
    WRITE(13,'(A50,A)') ADJUSTL(TEMPSTR1), " !Number of energy levels"

    WRITE(13,'(A,T50,A)') ADJUSTL('6'),' !Entry number of link to super level'
    WRITE(13,'(A)') ' '
!
    DO COUNTER_I=1, NLEVELS
       WRITE(13,'(A33,F5.1,F15.4,F12.5,ES13.3,I8,I8,I8)') ADJUSTL(ENERGY_LEVEL_LIST(COUNTER_I)%LEVEL_NAME), &
            ENERGY_LEVEL_LIST(COUNTER_I)%LEVEL_G, ENERGY_LEVEL_LIST(COUNTER_I)%ENERGY_CM1, &
            ENERGY_LEVEL_LIST(COUNTER_I)%FREQ_1015HZ, ENERGY_LEVEL_LIST(COUNTER_I)%LAMBDA_ANG, &
            1,0,ENERGY_LEVEL_LIST(COUNTER_I)%LEVEL_ID
    END DO
!
    CLOSE(13)
    CALL FLUSH
      
          CONTAINS
            
            !A function to convert Einstein A values to oscillator strengths (input A must be in s^-1, wavenum in cm^-1). G_2 is upper state statistical weight, G_1 lower
            
            REAL(8) FUNCTION OSC_FROM_EINSTEIN(EINS_A,G_1,G_2,WAVENUM) 
              REAL(8) :: STAT_RATIO, WAVENUM_METER, NUMERATOR, DENOMINATOR
              REAL(8), INTENT(IN) :: EINS_A, G_1, G_2, WAVENUM

              STAT_RATIO=G_2/G_1
              WAVENUM_METER=WAVENUM*100
              NUMERATOR=STAT_RATIO*EPS_0*ELEC_MASS*LIGHTM
              DENOMINATOR=2*PI*(WAVENUM_METER**2)*(E_CHARGE**2)

              OSC_FROM_EINSTEIN=(NUMERATOR/DENOMINATOR)*EINS_A
              RETURN
            END FUNCTION OSC_FROM_EINSTEIN
            
            !Formats the level name for a given v, J level
            CHARACTER(33) FUNCTION LEVEL_NAME_VJ(V,J)
              INTEGER, INTENT(IN) :: V,J
              IF (V .GT. 9 .AND. J .GT. 9) THEN
                 WRITE(LEVEL_NAME_VJ,"(A15,I2,A1,I2,A1)") "Rovibration_vJ{",v,"|",J,"}"
              ELSE IF (V .GT. 9 .AND. J .LE. 9) THEN
                 WRITE(LEVEL_NAME_VJ,"(A15,I2,A1,I1,A1)") "Rovibration_vJ{",v,"|",J,"}"
              ELSE IF (V .LE. 9 .AND. J .GT. 9) THEN
                 WRITE(LEVEL_NAME_VJ,"(A15,I1,A1,I2,A1)") "Rovibration_vJ{",v,"|",J,"}"                 
              ELSE IF (V .LE. 9 .AND. J .LE. 9) THEN
                 WRITE(LEVEL_NAME_VJ,"(A15,I1,A1,I1,A1)") "Rovibration_vJ{",v,"|",J,"}"
              END IF
              RETURN
            END FUNCTION LEVEL_NAME_VJ
            
    END PROGRAM HITRAN_CONVERT
      
