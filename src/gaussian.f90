module gaussian
    use constants, only: dp, zero, one, two, dpi
    implicit none

    interface gaussdev
        ! module procedure gaussdev_s, gaussdev_v
        module procedure random_normal_scalar, random_normal_1Darray, random_normal_2Darray
    end interface gaussdev

contains

    subroutine random_normal_scalar(x)
        real(dp), intent(out) :: x

        real(dp) :: r1, r2, u1, u2

        call random_number(r1)
        call random_number(r2)
        u1 = one-r1
        u2 = one-r2
        x = sqrt(-two*log(u1))*cos(dpi*u2)

        return
    end subroutine random_normal_scalar

    subroutine random_normal_1Darray(x)
        real(dp), dimension(:), intent(out) :: x

        integer::n
        real(dp), dimension(:), allocatable :: u1, u2

        n = size(x)
        allocate (u1(n), u2(n))
        call random_number(u1)
        call random_number(u2)
        u1 = one-u1
        u2 = one-u2
        x = sqrt(-two*log(u1))*cos(dpi*u2)
        deallocate (u1, u2)

        return
    end subroutine random_normal_1Darray

    subroutine random_normal_2Darray(x)
        real(dp), dimension(:, :), intent(out) :: x

        integer, dimension(2)::n
        real(dp), dimension(:, :), allocatable :: u1, u2

        n = shape(x)
        allocate (u1(n(1), n(2)), u2(n(1), n(2)))
        call random_number(u1)
        call random_number(u2)
        u1 = one-u1
        u2 = one-u2
        x = sqrt(-two*log(u1))*cos(dpi*u2)
        deallocate (u1, u2)

        return
    end subroutine random_normal_2Darray

end module gaussian
