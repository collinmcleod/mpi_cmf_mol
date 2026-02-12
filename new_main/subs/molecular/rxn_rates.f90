!Module for variables and functions important to the carbon monoxide chemical creation/destruction processes

module rxn_rates

  type reaction
    character(8) :: q1, q2, p1, p2, p3
    real(8) :: alpha, beta, gamma
    integer :: id
  end type

  real(8), allocatable :: rate_ks(:)
  type(reaction), allocatable :: reaction_list(:)
  type(reaction) :: rxn1, rxn2, rxn3

  character(20) :: rxn_file, tempstr
  integer :: nrxns

  contains

    subroutine make_reaction_list() !subroutine that reads the reaction information from a file (usu. REACTIONS_DATA) and packages it into the array reaction_list
      implicit none
      integer :: fileio, counti
      character(100) :: tstr1, tstr2

      rxn_file = "REACTIONS_DATA"

      open(unit=11,file=rxn_file,action='READ',iostat=fileio)
      if (fileio .ne. 0) stop "Error opening reaction_data file"

      do !Skip the header information in the file
        read(11,'(a)') tstr1
        if (index(tstr1,'^**') .ne. 0) exit
      end do

      read(11,*) nrxns, tstr1 !Read the number of reactions

      allocate(reaction_list(nrxns))

      read(11,'(a)') tstr1 !Skip a blank line

      do counti=1,nrxns !Loop to read reaction information from the file
        read(11,'(a)') tstr2
        read(tstr2,'(i8,a8,a8,a8,a8,a8,f12.2,f12.2,f12.5)') reaction_list(counti)%id,&
reaction_list(counti)%q1,reaction_list(counti)%q2,reaction_list(counti)%p1,reaction_list(counti)%p2,&
reaction_list(counti)%p3,reaction_list(counti)%alpha,reaction_list(counti)%beta,&
reaction_list(counti)%gamma
      end do

      call flush(11) !close the reactions_data file
      close(11)
    
      return

    end subroutine make_reaction_list

    subroutine make_k_list(T)

      real(8) :: T
      integer :: counti

      allocate(rate_ks(nrxns))

      do counti=1,nrxns
        call rate_k(reaction_list(counti)%alpha,reaction_list(counti)%beta,reaction_list(counti)%gamma,T,rate_ks(counti))
      end do

      return

    end subroutine make_k_list

    subroutine rate_k(alpha,beta,gamma,T,k_rate) !calculates the reaction rate constant k from reaction data specific to each reaction, as well as the temperature dependence
      implicit none
      real(8), intent(in) :: alpha, beta, gamma, T !T in units of 10^4 K
      real(8) :: tr1, tr2, k_rate
      
      if (T .le. 0.0) stop "Temperature must be positive"
     
      tr1 = (T/(0.03D0))**beta
      tr2 = exp(-gamma/(T*(10**4)))

      k_rate = alpha*tr1*tr2
      return

    end subroutine rate_k

end module rxn_rates
