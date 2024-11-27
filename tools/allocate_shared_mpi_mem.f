	MODULE ALLOCATE_SHARED_MPI_MEM
	IMPLICIT NONE
	PUBLIC SUB_ALLOCATE_DP_MPI_MEM
!
	CONTAINS
	SUBROUTINE SUB_ALLOCATE_DP_MPI_MEM(array,n,subid,mymaster,share_comm,win)
	USE SET_KIND_MODULE
	USE MPI
	USE ISO_C_BINDING

	IMPLICIT NONE
!
! Arguments
!
	DOUBLE PRECISION, POINTER, INTENT(inout) :: ARRAY(:)
	INTEGER, INTENT(in) :: N                                  !Lenth of array
	INTEGER, INTENT(in) :: SUBID                              !Thread ID
	INTEGER, INTENT(in) :: MYMASTER                           !Master thread ID
	INTEGER, INTENT(in) :: SHARE_COMM                         !Communicator (can be world)
	INTEGER, INTENT(inout) :: WIN                              !Window/memory identifier
!
! Local variables
!
	REAL(KIND=LDP) X
	INTEGER :: VAR_SIZE, IER
	INTEGER :: ARRAY_SHAPE(1)
	INTEGER(KIND=MPI_ADDRESS_KIND) :: WINDOW_SIZE
	TYPE(C_PTR) :: BASEPTR
!
! Set window size in bytes.
!
	VAR_SIZE = SIZEOF(X)
	WINDOW_SIZE = N*VAR_SIZE
	ARRAY_SHAPE=N
	CALL MPI_BARRIER(SHARE_COMM, IER)
!
! Allocate memory.
!
	CALL MPI_WIN_ALLOCATE_SHARED(window_size, VAR_SIZE, MPI_INFO_NULL, SHARE_COMM, BASEPTR, WIN ,IER)
!
! Get C pointer for other threads, and the convert from C pointer to Fortran pointer.
!
	IF (SUBID /= MYMASTER) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, VAR_SIZE, baseptr, ier)
	CALL C_F_POINTER(BASEPTR, ARRAY, ARRAY_SHAPE)
	CALL MPI_WIN_FENCE(0, WIN, IER)

	RETURN
	END SUBROUTINE SUB_ALLOCATE_DP_MPI_MEM
	END MODULE ALLOCATE_SHARED_MPI_MEM
