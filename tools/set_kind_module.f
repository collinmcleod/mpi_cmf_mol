	MODULE SET_KIND_MODULE
	IMPLICIT NONE
!
! This statment should not be changed.
!
	REAL(8) REAL_STAR_EIGHT
!
! Define the diuble precion variable to be used by the CMFGEN
! caculations.
!
	REAL(8), PRIVATE :: KVAR
        INTEGER, PARAMETER :: LDP=KIND(KVAR)
	INTEGER MYPE
	INTEGER NTHREAD
	INTEGER MY_MPI_DP
!
!        INTEGER, PARAMETER :: LDP=SELECTED_REAL_KIND(15,3000)
!
	SAVE
	END MODULE SET_KIND_MODULE
	
