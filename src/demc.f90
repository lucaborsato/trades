module demc
    use constants
    use gaussian
    implicit none
    public

contains

    subroutine de_proposal(x, xp)
        real(dp), dimension(:, :), intent(in)::x
        real(dp), dimension(:, :), intent(out)::xp

        real(dp), parameter::b = 1.0e-5_dp
        real(dp)::g
        real(dp), dimension(:), allocatable::e, dx
        integer, dimension(:), allocatable::idx_walkers, idx_r

        integer, dimension(2)::nshape, r12
        integer::i_walker

        nshape = shape(x) ! ndim, nw

        g = 2.38_dp/sqrt(two*nshape(1))

        idx_walkers = (/(i_walker, i_walker=1, nshape(2), 1)/)

        allocate (e(nshape(1)))
        allocate (dx(nshape(1)))

        dowalk: do i_walker = 1, nshape(2)
            call gaussdev(e)
            e = e*b ! G(0, b)
            idx_r = pack(idx_walkers, idx_walkers /= i_walker)
            call random_choice_int_without_repeatition(idx_r, 2, r12)
            dx = x(:, r12(1))-x(:, r12(2))
            xp(:, i_walker) = x(:, i_walker)+g*dx+e
        end do dowalk

        deallocate (e, dx)

        return
    end subroutine de_proposal

end module demc
