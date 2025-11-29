	MODULE ALLOCATE_SHARED_MEM
	PUBLIC ALLOCATE_MPI_MEM
!
	CONTAINS
	SUBROUTINE ALLOCATE_MPI_MEM(OWNER,ARRAY_SIZE,VAR_SIZE,BASEPTR,SHARE_COMM,WIN_ID)
	USE SET_KIND_MODULE
	USE MPI
	USE ISO_C_BINDING, ONLY : C_PTR
	IMPLICIT NONE
!
! Arguments
!
! This line, and the declaration of X, are the two type dependent lines
! in this subroutine.
!
	INTEGER, INTENT(inout) :: WIN_ID            !Identifies the shared momory
	INTEGER, INTENT(in) :: SHARE_COMM	    !Communicator
	TYPE(C_PTR),    INTENT(out) :: BASEPTR
	INTEGER(KIND=MPI_ADDRESS_KIND) :: WINDOW_SIZE
!
	INTEGER, INTENT(in) :: ARRAY_SIZE	    !Size of array to be created
	INTEGER, INTENT(in) :: VAR_SIZE	            !Size (in bytes) of each element in array
	INTEGER, INTENT(in) :: OWNER                !Thread on which memory will be located
!
	INTEGER :: IER
!
! Memory will be allocated to the thread OWNER. For all other threads,
! we request shared memoery of length zero.
!
! The size of the shared memorey must be in bytes.
!
	WINDOW_SIZE = ARRAY_SIZE*VAR_SIZE
	IF(OWNER .NE. MYPE)WINDOW_SIZE=0
!
	CALL MPI_WIN_ALLOCATE_SHARED(WINDOW_SIZE, VAR_SIZE, MPI_INFO_NULL, SHARE_COMM, BASEPTR, WIN_ID ,IER)
!	CALL MPI_WIN_ALLOCATE_SHARED(WINDOW_SIZE, VAR_SIZE, MPI_INFO_NULL, MPI_COMM_WORLD, BASEPTR, WIN_ID ,IER)
	IF (MYPE /= OWNER) CALL MPI_WIN_SHARED_QUERY(WIN_ID, 0, WINDOW_SIZE, var_size, BASEPTR, IER)

	RETURN
	END SUBROUTINE ALLOCATE_MPI_MEM
	END MODULE ALLOCATE_SHARED_MEM
