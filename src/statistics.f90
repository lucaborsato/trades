module statistics
    use constants, only: sp, dp, zero, half
    use sorting, only: sort
    implicit none

contains

    ! calculates the mean of a vector
    function mean(x) result(xm)
        real(dp)::xm
        real(dp), dimension(:), intent(in)::x
        xm = sum(x)/size(x)

        return
    end function mean

    ! calculates the weighted mean
    function wmean(x, wx) result(xm)
        real(dp)::xm
        real(dp), dimension(:), intent(in)::x, wx
        xm = sum(x*wx)/sum(wx)

        return
    end function wmean

    ! =========================================================================
    ! computes median
    function median_dp(x) result(median)
        real(dp)::median
        real(dp), dimension(:), intent(in)::x

        integer::n_x, h_n
        ! integer,dimension(:),allocatable::idx
        real(dp), dimension(:), allocatable::xsorted

        n_x = size(x)
        h_n = n_x/2
        median = zero
        if (n_x == 1) then
            median = x(1)
        else if (n_x > 1) then
            ! allocate(idx(n_x))
            ! call indexx(x, idx)
            allocate (xsorted(n_x))
            xsorted = x
            call sort(xsorted)
            if (mod(n_x, 2) .eq. 0) then
                ! median = half*(x( idx( h_n+1 ) ) + x( idx( h_n ) ))
                median = half*(xsorted(h_n+1)+xsorted(h_n))
            else
                ! median = x( idx( h_n ) )
                median = xsorted(h_n)
            end if
            ! deallocate(idx)
            deallocate (xsorted)
        end if

        return
    end function median_dp

    ! =========================================================================
    ! computes residuals at 68.27-th percentile w.r.t. input value val
    function std_equivalent_dp(x, val) result(std_eq)
        real(dp)::std_eq
        real(dp), dimension(:), intent(in)::x
        real(dp), intent(in)::val

        real(dp), dimension(:), allocatable::res
        integer::n_x, p

        n_x = size(x)
        p = int(0.6828_sp*real(n_x, sp))
        if (p .eq. 0) then
            p = 1
        end if

        allocate (res(n_x))
        res = abs(x-val)
        call sort(res)
        std_eq = res(p)
        deallocate (res)

        return
    end function std_equivalent_dp

end module statistics
