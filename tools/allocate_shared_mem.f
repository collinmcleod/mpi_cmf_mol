!
! Routine to allocate shared memory. The type of variable is not needed for this routine --
! only the array size and the variable size is needed.
!
! After calling this routine you need to call C_F_POINTER:
!
!       CALL C_F_POINTER(BASEPTR,VEC,[NX,NY] where VEC[NX,NY...] is a pointer ot the appropiate type.
!
! In the calling routine you aso need: 
!
!       USE ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
!       TYPE(C_PTR) BASEPTR
!
! Altered: 03-Feb-2025 : Intent deleted for VAR_SIZE. Now compiles on MAC.
!
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
	INTEGER, INTENT(inout) :: WIN_ID            !Identifies the shared momory
	INTEGER, INTENT(in) :: SHARE_COMM	    !Communicator
	TYPE(C_PTR) :: BASEPTR
	INTEGER(KIND=MPI_ADDRESS_KIND) :: WINDOW_SIZE
!
!NB: Giving this variable INTENT(in) cause the mpif90 on my mac to fail (no valid MPI_WIN_SHARED_QUERY)
!
	INTEGER :: VAR_SIZE		            !Size (in bytes) of each element in array
!
	INTEGER, INTENT(in) :: ARRAY_SIZE	    !Size of array to be created
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
	IF (MYPE /= OWNER) CALL MPI_WIN_SHARED_QUERY(WIN_ID, 0, WINDOW_SIZE, var_size, BASEPTR, IER)

	RETURN
	END SUBROUTINE ALLOCATE_MPI_MEM
	END MODULE ALLOCATE_SHARED_MEM
