!Program to write out files with very crude photoionzation cross-sections for Carbon Monoxide

PROGRAM MK_PHOT_FILE
  IMPLICIT NONE

  
  TYPE ENERGY_LEVEL
    CHARACTER(33) :: LEVEL_NAME
    REAL(8) :: LEVEL_G
    REAL(8) :: ENERGY_CMI
    REAL(8) :: FREQ_1015HZ
    REAL(8) :: ENERGY_BELOW_ION_EV
    REAL(8) :: LAMBDA_ANG
    INTEGER :: LEVEL_ID
    REAL(8) :: ARAD
  END TYPE ENERGY_LEVEL

  TYPE(ENERGY_LEVEL), DIMENSION(:), ALLOCATABLE :: CO_LEVEL_LIST
  TYPE(ENERGY_LEVEL) :: TLVL1, TLVL2

  REAL(8) :: TREAL1, TREAL2, TREAL3, GION1
  CHARACTER(200) :: TSTR1, TSTR2, CMDATAFILE
  CHARACTER(200) :: HEADER1, HEADER2, HEADER3, HEADER4, HEADER5, HEADER6, ELVLHEAD

  INTEGER :: NLVLSA, NLVLSB, COUNTI, COUNTJ, COUNTK, FILEIO, MAXV, MAXJ

  character(20) :: dateandtime

  call date_time(dateandtime)

  DO
    WRITE(*,*) "Enter values for maximum v and J quantum numbers: "
    READ(*,*) MAXV, MAXJ
    IF (MAXV <= 41 .AND. MAXJ <= 150) EXIT
    WRITE(*,*) "Max v must be less than 41, max J must be less than 150"
  END DO

  WRITE(CMDATAFILE,'(A16,I3.3,A6,I3.3,A4)') "CARB_MONOX_MAXV_", MAXV, "_MAXJ_",MAXJ,"_OSC"!Uses chosen max v and max J

  OPEN(UNIT=11,FILE=CMDATAFILE,ACTION='READ',IOSTAT=FILEIO)
  
  RDLOOP: DO !Read through header info for the number of energy levels
    READ(11,'(A)') TSTR1
    IF (INDEX(TSTR1,"!Number of energy levels") .NE. 0) THEN 
      READ(TSTR1,'(I5,A)') NLVLSA, TSTR2
      EXIT RDLOOP
    END IF
  END DO RDLOOP

  DO COUNTI=1,3
    READ(11,'(A)') TSTR1 !Skip rest of header info
  END DO

  NLVLSB = NLVLSA + 1

  ALLOCATE(CO_LEVEL_LIST(NLVLSB+1))

  GION1 = 0.0D0
  DO COUNTI=1,NLVLSB
    READ(11,'(A33,F5.1,F15.4,F12.5,F9.3,E13.3,I8,E12.3,E12.3,E12.3)') CO_LEVEL_LIST(COUNTI)%LEVEL_NAME,&
    CO_LEVEL_LIST(COUNTI)%LEVEL_G, CO_LEVEL_LIST(COUNTI)%ENERGY_CMI, CO_LEVEL_LIST(COUNTI)%FREQ_1015HZ,&
    CO_LEVEL_LIST(COUNTI)%ENERGY_BELOW_ION_EV, CO_LEVEL_LIST(COUNTI)%LAMBDA_ANG, &
    CO_LEVEL_LIST(COUNTI)%LEVEL_ID, CO_LEVEL_LIST(COUNTI)%ARAD, TREAL1, TREAL2
    GION1 = GION1+CO_LEVEL_LIST(COUNTI)%LEVEL_G
  END DO

  CO_LEVEL_LIST(:NLVLSA) = CO_LEVEL_LIST(2:NLVLSA+1)

  CLOSE(11)
  CALL FLUSH(11)

  CO_LEVEL_LIST(NLVLSB)%LEVEL_NAME="Ionized_CO"
  CO_LEVEL_LIST(NLVLSB)%LEVEL_G=100
  CO_LEVEL_LIST(NLVLSB)%ENERGY_CMI=113030.53768
  CO_LEVEL_LIST(NLVLSB)%FREQ_1015HZ=0.0
  CO_LEVEL_LIST(NLVLSB)%LEVEL_ID = NLVLSB
  CO_LEVEL_LIST(NLVLSB)%ARAD = 0

  OPEN(UNIT=12,FILE='PHOTCMI_A',ACTION='READWRITE',IOSTAT=FILEIO)

  IF (FILEIO>0) STOP "Error opening PHOTCMI_A"

  !Write the basic header information

  WRITE(12,*) '****************************************************************************************'
  WRITE(12,*) '****************************************************************************************'

  WRITE(12,*)
  WRITE(12,*) 'Crude photoionization cross-sections for Carbon Monoxide, based on CMI_PHOT_DATA'

  WRITE(12,*)
  WRITE(12,*) '****************************************************************************************'
  WRITE(12,*) '****************************************************************************************'

  WRITE(12,*)

  WRITE(12,'(A,t40,A)') dateandtime, "!Date"

  write(tstr1,*)nlvlsa
  WRITE(12,'(a,t40,A)') trim(adjustl(tstr1)), '!Number of energy levels'

  WRITE(12,'(a,t40,a)') '1','!Number of photoionization routes'
  WRITE(12,'(A,t40,A)') '1.0','!Screened nuclear charge'
  WRITE(12,'(A,t40,A)') 'UNKNOWN','!Final state in ion'
  WRITE(12,'(A,t40,A)') '0.0','!Excitation energy of final state'
  WRITE(12,'(F7.1,t40,A)') REAL(NINT(GION1*0.9)),'!Statistical weight of ion'
  WRITE(12,'(A,t40,A)') 'Megabarns','!Cross-section unit'
  WRITE(12,'(A,t40,A)') 'False','!Split J levels'
  write(tstr1,*)3*nlvlsa
  WRITE(12,'(a,t40,A)')  trim(adjustl(tstr1)), '!Total number of data pairs'

  !Write the basic info for low each level

  DO COUNTI=1,NLVLSA !Write out the photoionization information -- this is the simplest possible data, basically
    WRITE(12,*)
    WRITE(12,'(A33,t40,A)') CO_LEVEL_LIST(COUNTI)%LEVEL_NAME, '!Configuration name'
    WRITE(12,'(A,t40,A)') '1','!Type of cross-section'
    WRITE(12,'(A,t40,A)') '3','!Number of cross-section points'
    IF (COUNTI <= 3) THEN
      WRITE(12,'(a)') '1.0D0'
    ELSE
      WRITE(12,'(a)') '0.0000'
    END IF
    WRITE(12,'(a)') '1.0D0'
    WRITE(12,'(a)') '2.0D0'
  END DO

  close(12)
  call flush (12)  

  END PROGRAM MK_PHOT_FILE

  
