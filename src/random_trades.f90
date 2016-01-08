module random_trades
  use constants,only:dp
  implicit none

  contains

  subroutine init_random_seed_input(nx,input_seed)
    integer,intent(in)::nx,input_seed
    integer :: i, n!, clock
    integer, dimension(:), allocatable :: seed

    n = nx
    call random_seed(size = n)
    allocate(seed(n))
    seed = input_seed + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(PUT = seed)
    deallocate(seed)

    return
  end subroutine init_random_seed_input

  
    ! it initializes the seed/s, it uses a clock based generator
  subroutine init_random_seed_clock(nx)
    integer,intent(in)::nx
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    n = nx
    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
!     write(*,'(10000(1x,i10))')seed
    call random_seed(PUT = seed)
    deallocate(seed)
    write(*,*)

    return
  end subroutine init_random_seed_clock

  ! random number from subroutine to function
  function random_scalar() result(xr)
    real(dp)::xr
    call random_number(xr)
  end function random_scalar
  
end module random_trades
