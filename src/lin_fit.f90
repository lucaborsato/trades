module lin_fit
    use constants, only: dp, one, two
    implicit none

    interface linfit
        module procedure linfit_noerrors_0, linfit_noerrors_1, linfit_0, linfit_1
    end interface linfit

contains

    ! ------------------------------------------------------------------ !
    ! linear fit with vector x and y => y = m*x + q
    subroutine linfit_noerrors_0(x, y, m, q)
        integer, dimension(:), intent(in)::x
        real(dp), dimension(:), intent(in)::y
        real(dp), intent(out)::m, q

        integer::n
        real(dp), dimension(:), allocatable::xr
        real(dp)::Sxy, Sx, Sx2, Sy, nr, num, den

        n = size(x)
        allocate (xr(n))
        nr = real(n, dp)
        xr = real(x, dp)
        Sxy = sum(xr*y)
        Sx = sum(xr)
        Sx2 = sum(xr*xr)
        Sy = sum(y)
        den = nr*Sx2-(Sx*Sx)

        num = (nr*Sxy)-(Sx*Sy)
        m = num/den

        num = (Sx2*Sy)-(Sx*Sxy)
        q = num/den

        deallocate (xr)

        return
    end subroutine linfit_noerrors_0

    ! ------------------------------------------------------------------ !
    ! linear fit with vector x and y => y = m*x + q
    subroutine linfit_noerrors_1(x, y, m, q)
        real(dp), dimension(:), intent(in)::x
        real(dp), dimension(:), intent(in)::y
        real(dp), intent(out)::m, q

        integer::n
        real(dp)::Sxy, Sx, Sx2, Sy, nr, num, den

        n = size(x)
        nr = real(n, dp)

        Sxy = sum(x*y)
        Sx = sum(x)
        Sx2 = sum(x*x)
        Sy = sum(y)
        den = nr*Sx2-(Sx*Sx)

        num = (nr*Sxy)-(Sx*Sy)
        m = num/den

        num = (Sx2*Sy)-(Sx*Sxy)
        q = num/den

        return
    end subroutine linfit_noerrors_1

    ! ------------------------------------------------------------------ !
    ! linear fit with vector x and y and error on y (ey) => y = m*x + q
    subroutine linfit_0(x, y, ey, m, em, q, eq)
        integer, dimension(:), intent(in)::x
        real(dp), dimension(:), intent(in)::y, ey
        real(dp), intent(out)::m, em, q, eq

        integer::n
        real(dp), dimension(:), allocatable::xf, w, wyx, wy, wx, wx2
        real(dp)::Sw, Swyx, Swy, Swx, Swx2, delta, num, eq2, em2, nf

        n = size(x)
        nf = real(n, dp)

        allocate (xf(n), w(n), wyx(n), wy(n), wx(n), wx2(n))

        xf = real(x, dp)
        w = one/(ey*ey)
        Sw = sum(w)
        wyx = w*y*xf
        Swyx = sum(wyx)
        wy = w*y
        Swy = sum(wy)
        wx = w*xf
        Swx = sum(wx)
        wx2 = w*(xf*xf)
        Swx2 = sum(wx2)
        delta = Sw*Swx2-(Swx*Swx)
        num = Swy*Swx2-Swx*Swyx
        q = num/delta
        !write(*,*)"q => delta",delta,"num",num,"num/delta",q
        num = Sw*Swyx-Swy*Swx
        m = num/delta
        !write(*,*)"m => delta",delta,"num",num,"num/delta",m
        eq2 = Swx2/delta
        eq = sqrt(eq2/(nf-two))
        !write(*,*)"eq => eq2=Swx2/delta",eq2,"eq2/(n-2)",eq2/(n-2),"eq",eq
        em2 = Sw/delta
        em = sqrt(em2/(nf-two))
        !write(*,*)"em => em2=Sw/delta",em2,"em2/(n-2)",em2/(n-2),"em",em
        deallocate (xf, w, wyx, wy, wx, wx2)

        return
    end subroutine linfit_0

    ! linear fit with vector x and y and error on y (ey) => y = m*x + q
    subroutine linfit_1(x, y, ey, m, em, q, eq)
        real(dp), dimension(:), intent(in)::x
        real(dp), dimension(:), intent(in)::y, ey
        real(dp), intent(out)::m, em, q, eq
        integer::n
        real(dp), dimension(:), allocatable::w, wyx, wy, wx, wx2
        real(dp)::Sw, Swyx, Swy, Swx, Swx2, delta, num, eq2, em2, nf

        n = size(x)
        nf = real(n, dp)

        allocate (w(n), wyx(n), wy(n), wx(n), wx2(n))
        w = one/(ey*ey)
        Sw = sum(w)
        wyx = w*y*x
        Swyx = sum(wyx)
        wy = w*y
        Swy = sum(wy)
        wx = w*x
        Swx = sum(wx)
        wx2 = w*(x*x)
        Swx2 = sum(wx2)
        delta = Sw*Swx2-(Swx*Swx)
        num = Swy*Swx2-Swx*Swyx
        q = num/delta
        num = Sw*Swyx-Swy*Swx
        m = num/delta
        eq2 = Swx2/delta
        eq = sqrt(eq2/(nf-two))
        em2 = Sw/delta
        em = sqrt(em2/(nf-two))
        deallocate (w, wyx, wy, wx, wx2)

        return
    end subroutine linfit_1
    ! ------------------------------------------------------------------ !

end module lin_fit
