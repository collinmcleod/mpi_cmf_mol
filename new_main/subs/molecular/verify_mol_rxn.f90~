!
! Simple routine to check that the molecular reaction data read in obeys certain rules
! Most importantly, they must not violate charge balance
!
  SUBROUTINE VERIFY_MOL_RXN()
    USE MOL_RXN_MOD
    IMPLICIT NONE
!
! Created 22-Sep-2021
!
    INTEGER J
    INTEGER COUNT
!
    COUNT=0
    DO J=1,N_MOL_RXN
!       
       MOL_REACTION_AVAILABLE(J) = .TRUE.
       IF ((LEV_IN_POPS_MOL(J,1) .EQ. 0) .OR. (LEV_IN_POPS_MOL(J,2) .EQ. 0) .OR. &
       (LEV_IN_POPS_MOL(J,3) .EQ. 0) .OR. (LEV_IN_POPS_MOL(J,4) .EQ. 0)) THEN
          MOL_REACTION_AVAILABLE(J) = .FALSE.
          IF (COUNT .EQ. 0) WRITE(LUER1,*) ' '
          IF (COUNT .EQ. 0) WRITE(LUER1,*) 'Warning from VERIFY_MOL_RXN'
          WRITE(LUER1,*) 'Molecular reaction unavailable, Number is ',J
          COUNT = COUNT+1
       ELSE
          IF ((MOL_Z_CHG(J,1) + MOL_Z_CHG(J,2)) .NE. (MOL_Z_CHG(J,3) + MOL_Z_CHG(J,4)) .AND. (MOL_RXN_TYPE(J) .NE. 3)) THEN
             WRITE(LUER1,*) 'mol_rxn_type:', MOL_RXN_TYPE(J)
             WRITE(LUER1,*) 'MOL_RXN_TYPE(J) .NE. 3:', MOL_RXN_TYPE(J) .NE. 3
             WRITE(LUER1,*)'Error in VERIFY_MOL_RXN:'
             WRITE(LUER1,*)'Charge balance violated in molecular reaction number ',J
             WRITE(LUER1,*) MOL_RXN_SPECIESID(J,1),'+',MOL_RXN_SPECIESID(J,2),'->',MOL_RXN_SPECIESID(J,3),'+',MOL_RXN_SPECIESID(J,4)
             WRITE(LUER1,*) 'chgs:'
             WRITE(LUER1,*) MOL_Z_CHG(J,1), '---',MOL_Z_CHG(J,2),'---',MOL_Z_CHG(J,3),'---',MOL_Z_CHG(J,4)
!             IF (J .EQ. N_MOL_RXN) THEN
!                STOP
          END IF
       END IF
       
    END DO
    RETURN
    END SUBROUTINE VERIFY_MOL_RXN
!
