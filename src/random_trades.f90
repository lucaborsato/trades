module random_trades
    use constants, only: dp
    implicit none

contains

! ==============================================================================

    subroutine init_random_seed_input(nx, input_seed)
        integer, intent(in)::nx, input_seed
        integer :: i, n!, clock
        integer, dimension(:), allocatable :: seed

        n = nx
        call random_seed(size=n)
        allocate (seed(n))
        seed = input_seed+37*(/(i-1, i=1, n)/)
        call random_seed(PUT=seed)
        deallocate (seed)

        return
    end subroutine init_random_seed_input

! ==============================================================================

    ! it initializes the seed/s, it uses a clock based generator
    subroutine init_random_seed_clock(nx)
        integer, intent(in)::nx
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        n = nx
        call random_seed(size=n)
        allocate (seed(n))
        call system_clock(COUNT=clock)
        seed = clock+37*(/(i-1, i=1, n)/)
!     write(*,'(10000(1x,i10))')seed
        call random_seed(PUT=seed)
        deallocate (seed)
        write (*, *)

        return
    end subroutine init_random_seed_clock

! ==============================================================================

    ! random number from subroutine to function
    function random_scalar() result(xr)
        real(dp)::xr
        call random_number(xr)
    end function random_scalar

    ! Knuth Shuffle from http://rosettacode.org/wiki/Knuth_shuffle#Fortran
    subroutine shuffle_int(a)
        integer, dimension(:), intent(inout)::a
        integer::iter, randpos, temp
        real(dp)::r

        do iter = size(a), 2, -1
            call random_number(r)
            randpos = int(r*iter)+1
            temp = a(randpos)
            a(randpos) = a(iter)
            a(iter) = temp
        end do

    end subroutine shuffle_int

! ==============================================================================

    ! Numerical recepis
    subroutine random_int(mini, maxi, ri)
        integer, intent(in)::mini, maxi
        integer, intent(out)::ri

        real(dp)::r

        call random_number(r)
        ri = mini+int(maxi*r)

        return
    end subroutine random_int

    subroutine random_choice_int_with_repeatition(xvals, nchoice, xchoice)
        integer, dimension(:), intent(in)::xvals
        integer, intent(in)::nchoice
        integer, dimension(:), allocatable, intent(out)::xchoice

        integer::ichoice, nvals, idx

        nvals = size(xvals)
        if (.not. allocated(xchoice)) allocate (xchoice(nchoice))

        do ichoice = 1, nchoice
            call random_int(1, nvals, idx)
            xchoice(ichoice) = xvals(idx)
        end do

        return
    end subroutine random_choice_int_with_repeatition

    subroutine random_choice_int_without_repeatition(xvals, nchoice, xchoice)
        integer, dimension(:), intent(in)::xvals
        integer, intent(in)::nchoice
        integer, dimension(:), allocatable, intent(out)::xchoice

        logical, dimension(:), allocatable::used
        integer::ichoice, nvals, idx

        nvals = size(xvals)
        if (.not. allocated(xchoice)) allocate (xchoice(nchoice))
        allocate (used(nvals))
        used = .false.

        do ichoice = 1, nchoice
            inner: do
                call random_int(1, nvals, idx)
                if (.not. used(idx)) then
                    xchoice(ichoice) = xvals(idx)
                    used(idx) = .true.
                    exit inner
                end if
            end do inner
        end do

        return
    end subroutine random_choice_int_without_repeatition

! ==============================================================================
end module random_trades
