module numerical_integrator
    use constants, only: dp, zero, one, two, three, TOLERANCE, TOL_dp
    use parameters
    use eq_motion
    use celestial_mechanics, only: separation_mutual_Hill_check, statevector_amd_hill_stability
    implicit none

contains

    ! ------------------------------------------------------------------ !
    ! RUNGE-KUTTA-CASH-KARP by Numerical Recepis
    ! rkck integrator that calls only eq. of motion in astrocentric coordinates
    !
    ! subroutine rkck_a(m, rin, drdt, h, rout, error)
    !     real(dp), dimension(:), intent(in)::m, rin, drdt
    !     real(dp), intent(in)::h
    !     real(dp), dimension(:), intent(out)::rout, error
    !     real(dp), dimension(:), allocatable::rtemp

    !     !coefficients k(6*NB)[s], s=1,6
    !     real(dp), dimension(:), allocatable::k2, k3, k4, k5, k6

    !     !coefficients A[s,s-1]
    !     !used in the intermediate state vector:
    !     !ks:=[t+Cs*h,y[t]+sum(A[s,s-1]*k[s-1])]
    !     real(dp), parameter::A21 = one/5.0_dp
    !     real(dp), parameter::A31 = three/40.0_dp, A32 = 9.0_dp/40.0_dp
    !     real(dp), parameter::A41 = three/10.0_dp, A42 = -9.0_dp/10.0_dp,&
    !         &A43 = 6.0_dp/5.0_dp
    !     real(dp), parameter::A51 = -11.0_dp/54.0_dp, A52 = 5.0_dp/2.0_dp,&
    !         &A53 = -70.0_dp/27.0_dp, A54 = 35.0_dp/27.0_dp
    !     real(dp), parameter::A61 = 1631.0_dp/55296.0_dp, A62 = 175.0_dp/512.0_dp,&
    !         &A63 = 575.0_dp/13824.0_dp, A64 = 44275.0_dp/110592.0_dp, A65 = 253.0_dp/4096.0_dp

    !     !coefficients B[s]
    !     !used in the final state vector (4th order):
    !     !y[n+1]=y[n]+sum(B[s]*k[s])
    !     real(dp), parameter::B1 = 37.0_dp/378.0_dp, B3 = 250.0_dp/621.0_dp,&
    !         &B4 = 125.0_dp/594.0_dp, B6 = 512.0_dp/1771.0_dp

    !     !coefficients C[s]
    !     !used in the intermediate state vector:
    !     !ks:=[t+Cs*h,y[t]+sum(A[s,s-1]*k[s-1])]
    !     real(dp), parameter::C2 = one/5.0_dp, C3 = three/10.0_dp, C4 = three/5.0_dp,&
    !         &C6 = 7.0_dp/8.0_dp

    !     !coefficients E[s], s=1,6 for the error
    !     real(dp), parameter::E1 = B1-2825.0_dp/27648.0_dp,&
    !         &E3 = B3-18575.0_dp/48384.0_dp, E4 = B4-13525.0_dp/55296.0_dp,&
    !         &E5 = -277.0_dp/14336.0_dp, E6 = B6-0.25_dp

    !     allocate (rtemp(NBDIM), k2(NBDIM), k3(NBDIM), k4(NBDIM), k5(NBDIM), k6(NBDIM))
    !     !Numerical Recipes algorithm
    !     rtemp = rin+A21*h*drdt
    !     call eqmastro(m, rtemp, k2)
    !     rtemp = rin+h*(A31*drdt+A32*k2)
    !     call eqmastro(m, rtemp, k3)
    !     rtemp = rin+h*(A41*drdt+A42*k2+A43*k3)
    !     call eqmastro(m, rtemp, k4)
    !     rtemp = rin+h*(A51*drdt+A52*k2+A53*k3+A54*k4)
    !     call eqmastro(m, rtemp, k5)
    !     rtemp = rin+h*(A61*drdt+A62*k2+A63*k3+A64*k4+A65*k5)
    !     call eqmastro(m, rtemp, k6)
    !     rout = rin+h*(B1*drdt+B3*k3+B4*k4+B6*k6)
    !     error = h*(E1*drdt+E3*k3+E4*k4+E5*k5+E6*k6)
    !     deallocate (rtemp, k2, k3, k4, k5, k6)

    !     return
    ! end subroutine rkck_a

    subroutine rkck_a(m, rin, drdt, h, rout, error)
        real(dp), dimension(:), intent(in)::m, rin, drdt
        real(dp), intent(in)::h
        real(dp), dimension(:), intent(out)::rout, error
        real(dp), dimension(:), allocatable::rtemp

        !coefficients k(6*NB)[s], s=1,6
        real(dp), dimension(:), allocatable::k2, k3, k4, k5, k6

        ! !coefficients A[2:6]
        ! real(dp), parameter::A2 = one/5.0_dp
        ! real(dp), parameter::A3 = three/10.0_dp
        ! real(dp), parameter::A4 = three/5.0_dp
        ! real(dp), parameter::A5 = one
        ! real(dp), parameter::A6 = 7.0_dp/8.0_dp

        ! coefficients BXY, X=2:6, Y=1:5
        real(dp), parameter::B21 = one/5.0_dp
        real(dp), parameter::B31 = three/40.0_dp
        real(dp), parameter::B32 = 9.0_dp/40.0_dp
        real(dp), parameter::B41 = three/10.0_dp
        real(dp), parameter::B42 = -9.0_dp/10.0_dp
        real(dp), parameter::B43 = 6.0_dp/5.0_dp
        real(dp), parameter::B51 = -11.0_dp/54.0_dp
        real(dp), parameter::B52 = 5.0_dp/two
        real(dp), parameter::B53 = -70.0_dp/27.0_dp
        real(dp), parameter::B54 = 35.0_dp/27.0_dp
        real(dp), parameter::B61 = 1631.0_dp/55296.0_dp
        real(dp), parameter::B62 = 175.0_dp/512.0_dp
        real(dp), parameter::B63 = 575.0_dp/13824.0_dp
        real(dp), parameter::B64 = 44275.0_dp/110592.0_dp
        real(dp), parameter::B65 = 253.0_dp/4096.0_dp

        !coefficients C[1:6]
        real(dp), parameter::C1 = 37.0_dp/378.0_dp
        ! real(dp), parameter::C2 = zero
        real(dp), parameter::C3 = 250.0_dp/621.0_dp
        real(dp), parameter::C4 = 125.0_dp/594.0_dp
        ! real(dp), parameter::C5 = zero
        real(dp), parameter::C6 = 512.0_dp/1771.0_dp

        !coefficients E[s], s=1,6 for the error ==> DC in NR
        real(dp), parameter::E1 = C1-(2825.0_dp/27648.0_dp)
        ! real(dp), parameter::E2 = zero
        real(dp), parameter::E3 = C3-(18575.0_dp/48384.0_dp)
        real(dp), parameter::E4 = C4-(13525.0_dp/55296.0_dp)
        real(dp), parameter::E5 = -277.0_dp/14336.0_dp
        real(dp), parameter::E6 = C6-0.25_dp

        allocate (rtemp(NBDIM), k2(NBDIM), k3(NBDIM), k4(NBDIM), k5(NBDIM), k6(NBDIM))
        !Numerical Recipes algorithm
        rtemp = rin+B21*h*drdt
        call eqmastro(m, rtemp, k2)
        rtemp = rin+h*(B31*drdt+B32*k2)
        call eqmastro(m, rtemp, k3)
        rtemp = rin+h*(B41*drdt+B42*k2+B43*k3)
        call eqmastro(m, rtemp, k4)
        rtemp = rin+h*(B51*drdt+B52*k2+B53*k3+B54*k4)
        call eqmastro(m, rtemp, k5)
        rtemp = rin+h*(B61*drdt+B62*k2+B63*k3+B64*k4+B65*k5)
        call eqmastro(m, rtemp, k6)
        rout = rin+h*(C1*drdt+C3*k3+C4*k4+C6*k6)
        error = h*(E1*drdt+E3*k3+E4*k4+E5*k5+E6*k6)
        deallocate (rtemp, k2, k3, k4, k5, k6)

        return
    end subroutine rkck_a

!     function get_emax(emold, error, rscal) result(emax)
!         real(dp)::emax
!         real(dp), intent(in)::emold
!         real(dp), dimension(:), intent(in)::error, rscal
!         real(dp), dimension(:), allocatable::ertemp
!         real(dp)::maxtemp
!         integer::i

!         allocate (ertemp(size(error)))
!         ! ertemp = zero
!         do i = 7, NBDIM
!             if (abs(rscal(i)) .le. TOL_dp) then
!                 ertemp(i) = TOL_dp
!             else
!                 ertemp(i) = abs(error(i))/rscal(i)
!             end if
!         end do
!         maxtemp = maxval(ertemp(7:NBDIM))
!         emax = max(emold, maxtemp)
!         deallocate (ertemp)

!         return
!     end function get_emax

!     ! it calls the integrator and select the right step for the integration
!     subroutine int_rk_a(m, rin, drdt, initial_step, ok_step, next_step, rout, err)
!         real(dp), dimension(:), intent(in)::m, rin, drdt
!         real(dp), intent(in)::initial_step
!         real(dp), intent(out)::ok_step, next_step
!         real(dp), dimension(:), intent(out)::rout, err
!         real(dp), dimension(:), allocatable::rscal
!         real(dp)::emax, working_step, scale_factor
!         ! safety factor = sfac ;
!         real(dp), parameter::sfac = 0.9_dp
!         real(dp), parameter::exp1 = one/5.0_dp

!         working_step = initial_step !uses a temporary variable
!         scale_factor = one
!         allocate (rscal(NBDIM))
!         rscal = abs(rin)+abs(working_step*drdt)
!         sel: do
! !       emax=0._dp
!             emax = TOL_dp
!             call rkck_a(m, rin, drdt, working_step, rout, err)
!             emax = get_emax(emax, err, rscal)
!             scale_factor = sfac*((tol_int/emax)**(exp1))
!             ok_step = working_step
!             working_step = working_step*scale_factor
!             if (tol_int .ge. emax) exit sel
!         end do sel
!         next_step = working_step
!         deallocate (rscal)

!         return
!     end subroutine int_rk_a

    subroutine int_rk_a(mass, rin, drdt, initial_step, ok_step, next_step, rout, rerr)
        ! Input
        real(dp), dimension(:), intent(in) :: mass
        real(dp), dimension(:), intent(in) :: rin, drdt
        real(dp), intent(in) :: initial_step
        ! Output
        real(dp), intent(out) :: ok_step, next_step
        real(dp), dimension(:), intent(out) :: rout, rerr
        ! Local variables
        real(dp) :: err_old
        real(dp), parameter :: safe_factor = 0.9_dp
        real(dp), parameter :: minscale = 0.2_dp, maxscale = 5.0_dp
        real(dp), parameter :: atol = 1.0e-15_dp, rtol = 1.0e-13_dp
        real(dp), parameter :: err_default = 1.0e-10_dp
        ! | Metodo    | Soluzione accettata | Errore stimato | k corretto |
        ! | --------- | ------------------- | -------------- | ---------- |
        ! | RKCK5(4)  | ordine 5            | differenza 5−4 | **4**      |
        ! | DOPRI5(4) | ordine 5            | differenza 5−4 | **4**      |
        ! Valori moderni (Hairer / Söderlind / DOPRI)
        ! Per RK embedded di ordine 5(4), i valori più usati sono:
        ! alpha ≈ 0.7 / (k+1)
        ! beta  ≈ 0.4 / (k+1)
        ! con 
        ! 𝑘 +  1 = 5
        ! k+1=5:
        ! alpha = 0.14
        ! beta  = 0.08
        ! Oppure (molto comune):
        ! alpha = 0.2
        ! beta  = 0.1
        ! | Fonte             | α           | β     | Carattere          |
        ! | ----------------- | ----------- | ----- | ------------------ |
        ! | Numerical Recipes | 0.10        | 0.175 | Molto conservativo |
        ! | Hairer–Wanner     | 0.14        | 0.08  | Bilanciato         |
        ! | Söderlind         | 0.2         | 0.1   | Moderno, stabile   |
        ! | Solo P controller | 1/(k+1)=0.2 | 0     | Semplice           |
        ! real(dp), parameter :: k = 5.0_dp, beta = 0.4_dp/k, alpha = (one/k) - 0.75_dp*beta
        real(dp), parameter :: alpha = 0.2_dp, beta = 0.1_dp
        real(dp) :: working_step
        real(dp) :: sk, scale_factor, errval
        logical :: accept_step
        integer :: i, ndim

        err_old = 1.0e-12_dp
        ndim = size(rin)
        working_step = initial_step
        do 
            call rkck_a(mass, rin, drdt, working_step, rout, rerr)
            ! scaled norm of the rerr
            errval = zero
            doerr: do i = 7,ndim
                sk = atol + rtol * max(abs(rin(i)), abs(rout(i)))
                errval = errval + (rerr(i)/sk)**2
            end do doerr
            errval = sqrt(errval/real(ndim-6, dp))

            ! accept step if errval <= 1
            if (errval <= one)then
                accept_step = .true.
            else
                accept_step = .false.
            end if

            ! propose new step 
            if (errval > zero) then
                scale_factor = safe_factor * errval**(-alpha) * err_old**(beta)
            else
                scale_factor = maxscale
            end if

            ! new step
            scale_factor = min(maxscale, max(minscale, scale_factor))
            next_step = working_step * scale_factor

            ! ==== Aggiorna errore precedente ====
            if (accept_step) then
                ok_step = working_step
                err_old = max(errval, err_default)
                exit
            else
                working_step = next_step
            end if

        end do

        return
    end subroutine int_rk_a

    ! ------------------------------------------------------------------ !

    subroutine integrates_rk(m, rin, dt, rout)
        real(dp), dimension(:), intent(in)::m, rin
        real(dp), intent(in)::dt
        real(dp), dimension(:), intent(out)::rout

        real(dp)::itime, working_step, ok_step, next_step
        real(dp), dimension(:), allocatable::r1, r2, drdt, error

        rout = zero
        allocate (r1(NBDIM), r2(NBDIM), drdt(NBDIM), error(NBDIM))
        r1 = rin
        r2 = rin

        working_step = half*dt
        itime = zero

        loopint: do
            call eqmastro(m, r1, drdt)
            if (abs(itime+working_step) .gt. abs(dt)) then

                working_step = dt-itime
                call rkck_a(m, r1, drdt, working_step, r2, error)
                itime = dt

            else

                call int_rk_a(m, r1, drdt, working_step, ok_step, next_step, r2, error)
                itime = itime+ok_step
                working_step = next_step

            end if
            if (abs(itime-dt) .le. TOL_dp) exit loopint
            r1 = r2
        end do loopint
        rout = r2
        deallocate (r1, r2, drdt, error)

        return
    end subroutine integrates_rk

    ! ================================================================================
    ! performs the forward or backward integration in one direction
    ! from one time to another, computing only the final orbit
    ! none RVs or T0s are computed with this subroutine
    !
    ! **Input**
    ! mass     == masses of all bodies (star == 1)
    ! radius   == radius of all bodies
    ! rin      == stave vector of all bodies in astrocentric coordinates XYZVxVyVz
    ! time     == integration time
    ! stepsize == initial or proposed integration stepsize
    !
    ! **Output**
    ! rout   == final state vector of all bodies in astrocentric coordinates at time
    ! Hc     == Hill chek, .true. if stable, .false. if unstable for Hill's condition
    subroutine one_forward_step(mass, radius, rin, time, stepsize, Hc, rout)
        ! **Input**
        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::time
        ! **Input/Output**
        real(dp), intent(inout)::stepsize
        logical, intent(inout)::Hc
        ! **Output**
        real(dp), dimension(:), intent(out)::rout
        ! **local**
        integer::n_body, nb_dim
        real(dp)::working_step, ok_step, next_step, itime
        real(dp), dimension(:), allocatable::dr, r1, r2, err
        ! if you see a variable non listed here, probably it is a global variable
        ! defined in constants.f90 or parameters.f90

        n_body = size(mass)
        nb_dim = n_body*6
        working_step = sign(stepsize, time)

        itime = zero
        allocate (dr(nb_dim), r1(nb_dim), r2(nb_dim), err(nb_dim))
        r1 = rin

        if (abs(time) .gt. TOL_dp) then
            integration: do

                if (abs(itime+working_step) .gt. abs(time)) working_step = time-itime
                call eqmastro(mass, r1, dr)
                call int_rk_a(mass, r1, dr, working_step, ok_step, next_step, r2, err)
                
                if (close_encounter_check) then
                    Hc = separation_mutual_Hill_check(mass, radius, r2, do_hill_check)
                end if
                if (amd_hill_check) then
                    call statevector_amd_hill_stability(mass, r1, Hc)
                end if
                itime = itime+ok_step

                if (.not. Hc) then
                    exit integration
                end if
                if (abs(itime) .ge. abs(time)) then
                    exit integration
                end if

                working_step = next_step
                r1 = r2

            end do integration
            rout = r2

        else

            rout = rin
            next_step = working_step

        end if

        stepsize = next_step

        deallocate (dr, r1, r2, err)

        return
    end subroutine one_forward_step

    subroutine one_forward_step_no_check(mass, rin, time, stepsize, rout)
        ! **Input**
        real(dp), dimension(:), intent(in)::mass, rin
        real(dp), intent(in)::time
        ! **Input/Output**
        real(dp), intent(inout)::stepsize
        ! **Output**
        real(dp), dimension(:), intent(out)::rout
        ! **local**
        real(dp)::working_step, ok_step, next_step, itime
        real(dp), dimension(:), allocatable::dr, r1, r2, err
        ! if you see a variable non listed here, probably it is a global variable
        ! defined in constants.f90 or parameters.f90

        working_step = sign(stepsize, time)

        itime = zero
        allocate (dr(NBDIM), r1(NBDIM), r2(NBDIM), err(NBDIM))
        r1 = rin

        if (abs(time) .gt. TOL_dp) then
            integration: do

                if (abs(itime+working_step) .gt. abs(time)) working_step = time-itime
                call eqmastro(mass, r1, dr)
                call int_rk_a(mass, r1, dr, working_step, ok_step, next_step, r2, err)
                itime = itime+ok_step

                if (abs(itime) .ge. abs(time)) then
                    exit integration
                end if

                working_step = next_step
                r1 = r2

            end do integration
            rout = r2

        else

            rout = rin
            next_step = working_step

        end if

        stepsize = next_step

        deallocate (dr, r1, r2, err)

        return
    end subroutine one_forward_step_no_check

end module numerical_integrator
