! 
! This routine calls a C routine to get the maximum amount of memory
! that has been allocated to a FORTRAN process.
!
! If luout is non zero, the peak amount of memory is written to luout.
!
subroutine get_max_memory(peak_mem_bytes,peak_mem,mem_unit,desc,luout)
    use iso_c_binding
    use set_kind_module
    implicit none

    ! Interface to the C function
    interface
        subroutine get_max_rss(max_kb) bind(c, name="get_max_rss")
            import :: c_long
            integer(c_long), intent(out) :: max_kb
        end subroutine get_max_rss
    end interface

    real(kind=ldp) :: peak_mem_bytes		!Memory in bytes (retruned)
    real(kind=ldp) :: peak_mem                  !Memory is useful units.
    integer(c_long) :: peak_bytes               !Memory in bytes retruned by C function
!
    integer, save :: error_count=0
    integer :: luout
    character(len=*) mem_unit
    character(len=*) desc      			!Used as identified if output.

! Call the C function
    call get_max_rss(peak_bytes)
    peak_mem_bytes=peak_bytes
    if (peak_bytes >= 1024.0_LDP**4) then
      peak_mem=peak_bytes/(1024.0_LDP**4)
      mem_unit=' TB'
    else if (peak_bytes >= 1024**3) then
      peak_mem=peak_bytes/1.073741824E+09
      mem_unit=' GB'
    else if (peak_bytes >= 1024*1024) then
      peak_mem=peak_bytes/(1048576.0_LDP)
      mem_unit=' MB'
    else if (peak_bytes >= 1024) then
      peak_mem=peak_bytes/1024.0_LDP
      mem_unit=' kB'
    else if (peak_bytes >= 0) then
      peak_mem=peak_bytes
      mem_unit='  B'
    else
      error_count=error_count+1
      if(mype .eq. 0 .and. error_count .lt. 2)then
        write(6,*)'Error retrieving memory usage -- peak_mem=',peak_bytes
      end if
      peak_mem=0.0_ldp; mem_unit=' '
    end if
    if(luout .ne. 0)then
      write(luout,'(1X,A,F9.3,A)')'Memeory usage for '//TRIM(desc)//' is',peak_mem,trim(mem_unit)
    end if
    return

end subroutine get_max_memory

!
! Simple function to convert from bytes to a more convenient unit.
! Mainly for external use -- it was written after code used above.
!

subroutine cnvt_mem_unit(peak_mem_bytes,peak_mem,mem_unit)
    use set_kind_module
    implicit none
    real(kind=ldp) peak_mem_bytes
    real(kind=ldp) peak_mem
    character(len=3) mem_unit
!
    integer i,k
    integer, save :: cnt=0
    logical no_match
    character(len=3) unit_list(5)
    data unit_list/'  B',' kB',' MB',' GB',' TB'/
!
    peak_mem=peak_mem_bytes; k=1
    mem_unit=unit_list(1)
    do while(peak_mem .gt. 1000.0_ldp)
      k=k+1
      peak_mem=peak_mem/1000.0_ldp
      mem_unit=unit_list(k)
    end do
!
    return
    end
