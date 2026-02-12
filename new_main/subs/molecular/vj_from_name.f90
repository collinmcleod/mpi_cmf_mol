!
! Simple subroutine to extract v and J quantum numbers from a level name string
!
! Created 10-May-2022
!
SUBROUTINE VJ_FROM_NAME(NAME,V,J)
  IMPLICIT NONE
!
  CHARACTER(*) :: NAME
  INTEGER :: V,J
!
  CHARACTER(300) :: STR1,STR2,STR3
  INTEGER :: I,K,L
!
  K = INDEX(NAME,'|')
  STR1 = TRIM(ADJUSTL(NAME(K+1:)))
!
  IF ( K .EQ. 0 ) THEN
     WRITE(*,*) 'Error in VJ_FROM_NAME: name does not match template'
     WRITE(*,*) 'Name: ',NAME
     RETURN
  END IF
!
  L = INDEX(NAME(K+1:),'[')
!
  IF ( L .EQ. 2 ) THEN
     READ(STR1(1:L-1),'(I1)') V
  ELSE IF ( L .EQ. 3 ) THEN
     READ(STR1(1:L-1),'(I2)') V
  ELSE IF ( L .EQ. 4 ) THEN
     READ(STR1(1:L-1),'(I3)') V
  END IF
!
  STR2 = TRIM(ADJUSTL(STR1(L+1:)))
  L = INDEX(TRIM(ADJUSTL(STR2)),']')
!
  IF ( L .EQ. 2 ) THEN
     READ(STR2(1:L-1),'(I1)') J
  ELSE IF ( L .EQ. 3 ) THEN
     READ(STR2(1:L-1),'(I2)') J
  ELSE IF ( L .EQ. 4 ) THEN
     READ(STR2(1:L-1),'(I3)') J
  END IF
!
  RETURN
END SUBROUTINE VJ_FROM_NAME
