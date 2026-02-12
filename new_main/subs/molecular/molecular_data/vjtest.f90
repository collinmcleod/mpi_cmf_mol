PROGRAM TESTVJ
!
  IMPLICIT NONE
!
  INTEGER :: I,J,K,L,V
!
  CHARACTER(50) :: TEST_STR
!
  TEST_STR='Rovibration_vJ|1|3|'
!
  WRITE(*,*)
  CALL VJ_FROM_NAME(TEST_STR,V,J)
  WRITE(*,*) 'Name: ',TEST_STR,', v: ',V,', J: ',J
!
  TEST_STR='Rovibration_vJ|10|1|'
!
  WRITE(*,*)
  CALL VJ_FROM_NAME(TEST_STR,V,J)
  WRITE(*,*) 'Name: ',TEST_STR,', v: ',V,', J: ',J
!
  TEST_STR='Rovibration_vJ|1|15|'
!
  WRITE(*,*)
  CALL VJ_FROM_NAME(TEST_STR,V,J)
  WRITE(*,*) 'Name: ',TEST_STR,', v: ',V,', J: ',J
!
  TEST_STR='Rovibration_vJ|210|30|'
!
  WRITE(*,*)
  CALL VJ_FROM_NAME(TEST_STR,V,J)
  WRITE(*,*) 'Name: ',TEST_STR,', v: ',V,', J: ',J
!
  TEST_STR='Rovibration_vJ|110|135|'
!
  WRITE(*,*)
  CALL VJ_FROM_NAME(TEST_STR,V,J)
  WRITE(*,*) 'Name: ',TEST_STR,', v: ',V,', J: ',J
!
  TEST_STR='Rovibration_vJ|10|13|'
!
  WRITE(*,*)
  CALL VJ_FROM_NAME(TEST_STR,V,J)
  WRITE(*,*) 'Name: ',TEST_STR,', v: ',V,', J: ',J
!
END PROGRAM TESTVJ
