module gaussian
    use constants, only: dp, zero, one, two, dpi
    implicit none

    interface gaussdev
        ! module procedure gaussdev_s, gaussdev_v
        module procedure random_normal_scalar, random_normal_1Darray, random_normal_2Darray
    end interface gaussdev

contains

    ! ! THESES SUBROUTINES COULD FAIL WITH N >= 1000000!!!

    ! ! Based on Numerical Recipes
    ! ! it calculates a gaussian number of kind N(0,1)
    ! subroutine gaussdev_s(gauss)
    !     real(dp), intent(out)::gauss
    !     real(dp), dimension(2)::v1
    !     real(dp)::rsq, lrsq, lrsq2, slrsq2

    !     do
    !         call random_number(v1)
    !         v1 = two*v1-one
    !         rsq = v1(1)**2+v1(2)**2
    !         if ((rsq .gt. zero) .and. (rsq .lt. one)) exit
    !     end do
    !     lrsq = log(rsq)
    !     lrsq2 = lrsq/rsq
    !     slrsq2 = sqrt(-two*lrsq2)
    !     gauss = v1(1)*slrsq2

    !     return
    ! end subroutine gaussdev_s

    ! ! Based on Numerical Recipes
    ! ! it calculates a gaussian vector of kind N(0,1)
    ! subroutine gaussdev_v(gauss)
    !     real(dp), dimension(:), intent(out)::gauss
    !     real(dp), dimension(size(gauss))::v1, v2, rsq, lrsq, lrsq2, slrsq2
    !     integer::n, ng, nn, m
    !     logical, dimension(size(gauss))::mask

    !     n = size(gauss)
    !     ng = 1
    !     do
    !         if (ng .gt. n) exit
    !         call random_number(v1(ng:n))
    !         call random_number(v2(ng:n))
    !         v1(ng:n) = two*v1(ng:n)-one
    !         v2(ng:n) = two*v2(ng:n)-one
    !         rsq(ng:n) = v1(ng:n)**2+v2(ng:n)**2
    !         mask(ng:n) = ((rsq(ng:n) .gt. zero) .and. (rsq(ng:n) .lt. one))
    !         call array_copy(pack(v1(ng:n), mask(ng:n)), v1(ng:), nn, m)
    !         v2(ng:ng+nn-1) = pack(v2(ng:n), mask(ng:n))
    !         rsq(ng:ng+nn-1) = pack(rsq(ng:n), mask(ng:n))
    !         ng = ng+nn
    !     end do
    !     lrsq = log(rsq)
    !     lrsq2 = lrsq/rsq
    !     slrsq2 = sqrt(-two*lrsq2)
    !     gauss = v1*slrsq2

    !     return
    ! end subroutine gaussdev_v

    ! ! Numerical Recipes
    ! ! it copies a masked array into another, needed by gaussdev_v
    ! subroutine array_copy(src, dest, n_copied, n_not_copied)
    !     real(dp), dimension(:), intent(in)::src
    !     real(dp), dimension(:), intent(out)::dest
    !     integer, intent(out)::n_copied, n_not_copied

    !     n_copied = min(size(src), size(dest))
    !     n_not_copied = size(src)-n_copied
    !     dest(1:n_copied) = src(1:n_copied)

    !     return
    ! end subroutine array_copy

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
