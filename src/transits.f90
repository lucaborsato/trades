module transits
    use constants
    use custom_type
    use parameters
    use celestial_mechanics, only: rsky, barycenter, fgfunctions, elements_one_body, eleMD => elem_mer
    use linear_ephem
    use convert_type, only: string
    use numerical_integrator, only: one_forward_step_no_check, one_forward_step
    implicit none

    interface check_T0
        module procedure check_T0_1, check_T0_2
    end interface check_T0

contains

    ! ------------------------------------------------------------------ !

    !
    ! Impact parameter of the transit
    !
    function impact_parameter(Rs, sma_p, inc_p, ecc_p, arg_p, R_p) result(b)
        real(dp)::b
        real(dp), intent(in)::Rs, sma_p, inc_p
        real(dp), optional, intent(in)::ecc_p, arg_p, R_p

        real(dp)::Rsum, aRs, rhoc

        Rsum = Rs*RsunAU
        if (present(R_p)) Rsum = (Rs+R_p)*RsunAU
        aRs = sma_p/Rsum

        if (present(ecc_p) .and. present(arg_p)) then

            rhoc = (one-(ecc_p*ecc_p))/(one+ecc_p*sin(arg_p*deg2rad))

        else

            rhoc = one

        end if

        b = aRs*rhoc*cos(inc_p*deg2rad)

        return
    end function impact_parameter

    ! ------------------------------------------------------------------ !
    ! move the whole state vectors (for all the bodies) of time = dt
    ! using fgfunctions subroutine
    subroutine advancefg(m, rw, dt, Hc)
        real(dp), dimension(:), intent(in)::m
        real(dp), dimension(:), intent(inout)::rw
        real(dp), intent(in)::dt
        logical, intent(inout)::Hc

        integer::j1, j2, i1, i6
        real(dp)::mu

        integer::n_body
        n_body = size(m)

        do j1 = 2, n_body
            j2 = (j1-1)*6
            i1 = 1+j2
            i6 = 6+j2
            mu = Giau*(m(1)+m(j1))
            call fgfunctions(mu, rw(i1:i6), dt, Hc)
            if (.not. Hc) exit
        end do

        return
    end subroutine advancefg
    ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
    ! a bisection step, it updates a boundary and the step for the next iteration
    ! - transit -
    subroutine onetra_bis(itra, m, A, B, rw, dt, Hc)
        integer::itra
        real(dp), dimension(:), intent(in)::m
        real(dp), intent(inout)::A, B, dt
        real(dp), dimension(:), intent(inout)::rw
        logical, intent(inout)::Hc

        integer::i1, i2, i4, i5
        real(dp)::C

        ! let's use the F&G functions to get closer to the transit
        call advancefg(m, rw, dt, Hc)

        i1 = 1+(itra-1)*6
        i2 = 2+(itra-1)*6
        i4 = 4+(itra-1)*6
        i5 = 5+(itra-1)*6
        C = rw(i1)*rw(i4)+rw(i2)*rw(i5)
        if ((A*C) .lt. zero) then
            B = C
            dt = -dt
        else
            A = C
        end if

        return
    end subroutine onetra_bis

    ! given a state vector and the itra of the body to check it computes the next step
    ! for the Newton-Raphson method
    ! - transit -
    subroutine onetra_nr(itra, mass, rw, dt)
        use eq_motion, only: eqmastro
        integer, intent(in)::itra
        real(dp), dimension(:), intent(in)::mass, rw
        real(dp), intent(out)::dt
        integer::ix, iy, ivx, ivy
        real(dp)::x, y, vx, vy, ax, ay, vx2, xax, vy2, yay, fk, dfk
        real(dp), dimension(:), allocatable::drw

        dt = zero
        
        ix  = 1+(itra-1)*6
        iy  = 2+(itra-1)*6
        ivx = 4+(itra-1)*6
        ivy = 5+(itra-1)*6
        
        x  = rw(ix)
        y  = rw(iy)
        vx = rw(ivx)
        vy = rw(ivy)
        
        ! Newton-Raphson
        ! f(k) = x*vx + y*vy
        fk = x*vx+y*vy
        allocate (drw(NBDIM))
        call eqmastro(mass, rw, drw)
        ! df(k)/dt = vx^2 + x*ax + vy^2 + y*ay
        ax = drw(ivx)
        ay = drw(ivy)
        vx2 = vx*vx
        xax = x*ax
        vy2 = vy*vy
        yay = y*ay
        dfk = vx2+xax+vy2+yay
        dt = -(fk/dfk)
        deallocate (drw)

        return
    end subroutine onetra_nr

    ! given a state vector and the itra of the body to check,
    ! it computes the next step
    ! for the Newton-Raphson method
    ! - contact -
    subroutine onecont_nr(itra, Rcheck, rw, dt)
        integer, intent(in)::itra
        real(dp), intent(in)::Rcheck
        real(dp), dimension(:), intent(in)::rw
        real(dp), intent(out)::dt
        integer::ix, iy, ivx, ivy, jtra
        real(dp)::rs2, hk, dhk, xvx, yvy

        dt = zero

        jtra = (itra-1)*6
        ix = 1+jtra
        iy = 2+jtra
        ivx = 4+jtra
        ivy = 5+jtra

        ! Newton-Raphson
        rs2 = rw(ix)*rw(ix)+rw(iy)*rw(iy)
        hk = rs2-Rcheck

        xvx = rw(ix)*rw(ivx)
        yvy = rw(iy)*rw(ivy)
        dhk = two*(xvx+yvy)
        if (abs(dhk) .le. TOLERANCE) dhk = TOLERANCE
        dt = -(hk/dhk)

        return
    end subroutine onecont_nr

    ! a bisection step, it updates a boundary and the step for the next iteration
    ! - contact -
    subroutine onecont_bis(itra, Rcheck, m, A, B, rw, dt, Hc)
        integer::itra
        real(dp), intent(in)::Rcheck
        real(dp), dimension(:), intent(in)::m
        real(dp), intent(inout)::A, B, dt
        real(dp), dimension(:), intent(inout)::rw

        logical, intent(inout)::Hc
        integer::i1, i2, i4, i5
        real(dp)::C

        ! let's use the F&G functions to get closer to the transit
        call advancefg(m, rw, dt, Hc)

        i1 = 1+(itra-1)*6
        i2 = 2+(itra-1)*6
        i4 = 4+(itra-1)*6
        i5 = 5+(itra-1)*6
        C = (rsky(rw(i1:i2))**2)-Rcheck
        if ((A*C) .lt. zero) then
            B = C
            dt = -half*dt
        else
            A = C
            dt = half*dt
        end if

        return
    end subroutine onecont_bis

    ! ------------------------------------------------------------------ !

    ! IT DEFINES THE RIGHT DIMENSION FOR THE RADIUS CHECK IN TRANSIT AND CONTACT TIMES
    subroutine Rbounds(itra, radii, Rs, Rp, Rmin, Rmax)
        integer, intent(in)::itra
        real(dp), dimension(:), intent(in)::radii
        real(dp), intent(out)::Rs, Rp, Rmin, Rmax

        Rs = radii(1)*RsunAU
        Rp = radii(itra)*RsunAU
        Rmax = Rs+Rp
        Rmin = Rs-Rp

        return
    end subroutine Rbounds

    ! ------------------------------------------------------------------ !
    ! seeks the root of rsky = 0 using an hybrid seeker: Bisection + Newton-Raphson
    ! subroutine find_transit(itra, mass, radius, r1, r2, time_r1, int_step, tmidtra, lte, ro, status)
    subroutine find_transit(t_epoch, itra, mass, r1, r2, time_r1, int_step, tmidtra, lte, ro, status)
        ! Input
        real(dp),intent(in)::t_epoch
        integer, intent(in)::itra
        real(dp), dimension(:), intent(in)::mass
        ! real(dp), dimension(:), intent(in)::radius
        real(dp), dimension(:), intent(in)::r1, r2
        real(dp), intent(in)::time_r1, int_step
        ! Output
        real(dp), intent(out)::tmidtra, lte
        ! Input/Output
        real(dp), dimension(:), intent(inout)::ro
        logical, intent(inout)::status
        ! Local
        real(dp)::stepsize
        real(dp), dimension(:), allocatable::rw, rwbar
        real(dp), dimension(6)::bar
        real(dp)::A, B, dt1, dt2
        integer::ix, iy, ivx, ivy
        integer::loop, many_iter
        integer::n_body, nb_dim

        status = .true.

        ix  = 1+(itra-1)*6
        iy  = 2+(itra-1)*6
        ivx = 4+(itra-1)*6
        ivy = 5+(itra-1)*6
        ! A
        A = r1(ix)*r1(ivx)+r1(iy)*r1(ivy)
        ! B
        B = r2(ix)*r2(ivx)+r2(iy)*r2(ivy)

        stepsize = sign(step_0, int_step)
        dt1 = int_step
        dt2 = zero

        
        tmidtra = dt1
        
        n_body = size(mass)
        nb_dim = n_body*6
        allocate (rw(nb_dim))
        rw = r2

        if (abs(dt1) .gt. TOLERANCE) then
            loop = 0
            many_iter = 1000

            traloop: do
                loop = loop+1
                if (abs(dt1) .le. TOLERANCE) then
                    exit traloop
                end if
                if (loop .ge. many_iter) then
                    exit traloop
                end if
                call onetra_nr(itra, mass, rw, dt2)

                if (abs(dt2) .lt. abs(dt1)) then
                    ! continue with N-R
                    dt1 = dt2
                else
                    ! -- BISECTION
                    if (A*B .lt. zero)then
                        dt1 = -dt1*half
                    else
                        dt1 = dt1*half
                    end if
                end if
                call one_forward_step_no_check(mass, rw, dt1, stepsize, ro) ! internally: sign(stepsize, dt1)
                rw = ro
                tmidtra = tmidtra+dt1
                ! update B and A
                A = B
                B = rw(ix)*rw(ivx)+rw(iy)*rw(ivy)
            end do traloop

        end if
        flush(6)

        allocate(rwbar(nb_dim))
        call barycenter(mass, rw, bar, rwbar)
        lte = -rwbar(3)/speedaud
        ro = rw
        tmidtra = tmidtra+t_epoch+time_r1+lte
        deallocate (rw, rwbar)

        return
    end subroutine find_transit

    ! IT CALCULATES A CONTACT TIME
    subroutine one_contact(icon, itra, mass, radii, rtra, ttra, tcont)
        integer, intent(in)::icon, itra
        real(dp), dimension(:), intent(in)::mass, radii, rtra
        real(dp), intent(in)::ttra
        real(dp), intent(out)::tcont

        real(dp)::stepsize
        real(dp), dimension(:), allocatable::ro, rw, rwbar
        real(dp), dimension(6)::bar
        real(dp)::vmid, dt1, dt2, A, B
        real(dp)::Rp, Rs, Rmin, Rmax, Rcheck, tmid, tt, lte
        integer::jtra, ix, iy, ivx, ivy
        ! logical::Hc
        integer::loop, many_iter

        jtra = (itra-1)*6
        ix = 1+jtra
        iy = 2+jtra
        ivx = 4+jtra
        ivy = 5+jtra

        ! icon = 0 ==> t_1.5
        ! icon = 1 ==> t_1
        ! icon = 2 ==> t_2
        ! icon = 3 ==> t_3
        ! icon = 4 ==> t_4
        ! icon = 5 ==> t_3.5

        call Rbounds(itra, radii, Rs, Rp, Rmin, Rmax)
        vmid = rsky(rtra(ivx:ivy))

        stepsize = step_0
        dt1 = -Rs/vmid ! t_1.5, t_1, t_2: - sign
        dt2 = zero
        if (icon .ge. 3) dt1 = -dt1 ! t_3, t_4, t_3.5: + sign

        if ((icon .eq. 2) .or. (icon .eq. 3)) then
            Rcheck = Rmin*Rmin
        else if ((icon .eq. 0) .or. (icon .eq. 5)) then
            Rcheck = Rs*Rs
        else
            Rcheck = Rmax*Rmax
        end if

        allocate (ro(NBDIM), rw(NBDIM), rwbar(NBDIM))
        call barycenter(mass, rtra, bar, rwbar)
        lte = -rwbar(3)/speedaud
        tmid = ttra-lte

        bar = zero
        rwbar = zero
        tcont = zero
        tt = dt1

        ! move the state vector of dt1 as suggested by Fabricky 2010
        call one_forward_step_no_check(mass, rtra, dt1, stepsize, rw)

        ! --

        A = (rtra(ix)*rtra(ix)+rtra(iy)*rtra(iy))-Rcheck
        B = (rw(ix)*rw(ix)+rw(iy)*rw(iy))-Rcheck

        many_iter = 1000
        loop = 0
        contloop: do
            loop = loop+1
            if (abs(dt1) .le. TOLERANCE) exit contloop
            if (loop .ge. many_iter) then
                exit contloop
            end if
            call onecont_nr(itra, Rcheck, rw, dt2)
            if (abs(dt2) .le. abs(dt1)) then
                dt1 = dt2
            else ! bisection
                if (A*B .lt. zero) then
                    dt1 = -dt1*half
                else
                    dt1 = dt1*half
                end if
            end if
            call one_forward_step_no_check(mass, rw, dt1, stepsize, ro)
            rw = ro
            tt = tt+dt1
            A = B
            B = (rw(ix)*rw(ix)+rw(iy)*rw(iy))-Rcheck
        end do contloop

        call barycenter(mass, rw, bar, rwbar)
        lte = -rwbar(3)/speedaud
        deallocate (rw, ro, rwbar)
        tcont = tmid+tt+lte

        return
    end subroutine one_contact

    ! IT DETERMINES ALL CONTACT TIMES (IF THEY EXIST) OF TRANSIT
    subroutine find_contacts(itra, mass, radii, rtra, ttra, tcont)
        integer, intent(in)::itra
        real(dp), dimension(:), intent(in)::mass, radii, rtra
        real(dp), intent(in)::ttra
        real(dp), dimension(4), intent(out)::tcont
        real(dp)::r_sky, Rs, Rp, Rmin, Rmax
        integer::jcont, jtra, ix, iy, step

        jtra = (itra-1)*6
        ix = 1+jtra
        iy = 2+jtra
        r_sky = rsky(rtra(ix:iy))
        call Rbounds(itra, radii, Rs, Rp, Rmin, Rmax)
        tcont = zero
        if (r_sky .le. Rmax) then
            step = 3
            if (r_sky .lt. Rmin) step = 1
            do jcont = 1, 4, step
                call one_contact(jcont, itra, mass, radii, rtra, ttra, tcont(jcont))
            end do
            if (r_sky .ge. Rmin) then
                tcont(2) = tcont(1)
                tcont(3) = tcont(4)
            end if
        end if

        return
    end subroutine find_contacts

    ! computes transit duration as:
    ! duration = t_3.5 - t_1.5
    ! where
    ! t_1.5 = time between contact time 1 and 2,
    ! t_3.5 = time between contact time 3 and 4,
    ! when project distance of the centre of the planet is on the edge of the star:
    ! rsky == Rstar (ingress <-> t_1.5, egress <-> t_3.5)
    subroutine compute_transit_duration_c2c(id_body, mass, radii, rtra, ttra, duration)
        integer, intent(in)::id_body
        real(dp), dimension(:), intent(in)::mass, radii, rtra
        real(dp), intent(in)::ttra
        real(dp), intent(out)::duration

        integer::sel_r
        real(dp)::r_sky, Rs, Rp, Rmin, Rmax
        real(dp)::t_hing, t_hegr

        t_hing = zero
        t_hegr = zero

        call Rbounds(id_body, radii, Rs, Rp, Rmin, Rmax)
        sel_r = (id_body-1)*6
        r_sky = rsky(rtra(1+sel_r:2+sel_r))

        if (r_sky .le. Rmax) then
            ! computes the t_1.5 == t_hing = planet on the edge of the star
            call one_contact(0, id_body, mass, radii, rtra, ttra, t_hing)
            call one_contact(5, id_body, mass, radii, rtra, ttra, t_hegr)

            duration = t_hegr-t_hing

        end if

        return
    end subroutine compute_transit_duration_c2c

    ! computes transit duration as in Kipping 2010, eq. (15):
    ! duration = T1 = (P/pi) * (rhoc^2 / sqrt(1-ecc^2)) * arcsin( sqrt(1 - (a/Rs)^2 * rhoc^2 * cos(inc)^2) / (a/Rs) * rhoc * sin(inc) )
    ! where rhoc = (1-ecc^2) / (1 +/- ecc*sin(argp)) with + <-> transit, - <-> occultation (W11_eq7-8, K10 rhoc)
    subroutine compute_transit_duration_K10_15(id_body, mass, radii, rtra, duration)
        integer, intent(in)::id_body
        real(dp), dimension(:), intent(in)::mass, radii, rtra
        real(dp), intent(out)::duration

        integer::sel_r
        real(dp)::mu, P_p, sma_p, ecc_p, inc_p, mA_p, w_p, lN_p, f_p, dtau_p
        real(dp)::ome2, rhoc, aRs, num1, den1, asin_num, asin_den

        mu = Giau*(mass(1)+mass(id_body))
        sel_r = (id_body-1)*6
        call elements_one_body(id_body, mass, rtra, P_p, sma_p, ecc_p, inc_p, mA_p, w_p, lN_p, f_p, dtau_p)

        ! 1 - ecc^2
        ome2 = one-(ecc_p*ecc_p)
        rhoc = ome2/(one+ecc_p*sin(w_p*deg2rad)) ! in some cases the + would be change to - for occultation
        ! sma / Rstar
        aRs = sma_p/(radii(1)*RsunAU)
        ! T1 Kipping 2010 eq. 15: One-term expression
        num1 = P_p*rhoc*rhoc
        den1 = dpi*sqrt(ome2)
        asin_num = sqrt(one-aRs*aRs*rhoc*rhoc*cos(inc_p*deg2rad))
        asin_den = aRs*rhoc*sin(inc_p*deg2rad)
        duration = asin(asin_num/asin_den)*num1/den1

        return
    end subroutine compute_transit_duration_K10_15

    ! call find_transit to compute the transit time (TT) and it assigns the right place
    ! of the TT comparing with the observations
    ! Computes duration of the transit (does not assign it now) as T4 - T1
    subroutine check_T0_1(t_epoch, id_body, mass, radii, r1, r2, time_r1, integration_step, simT0, Hc)
        ! Input
        real(dp),intent(in)::t_epoch
        integer, intent(in)::id_body
        real(dp), dimension(:), intent(in)::mass, radii, r1, r2
        real(dp), intent(in)::time_r1, integration_step
        ! Input/Output
        type(dataT0), dimension(:), intent(inout)::simT0
        logical, intent(inout)::Hc

        call check_T0_2(t_epoch, id_body, mass, radii, r1, r2, time_r1, integration_step, do_transit, durcheck, obsData, simT0, Hc)

        return
    end subroutine check_T0_1

    ! same as check_T0_1, but the epochs of the transit have to be provided
    ! as the array T0_num
!   subroutine check_T0_2(id_body,mass,radii,r1,r2,time_r1,integration_step,transit_flag,dur_check,n_T0,T0_num,T0_stat,T0_sim,Hc)
    subroutine check_T0_2(t_epoch, id_body, mass, radii, r1, r2, time_r1, integration_step, transit_flag, dur_check, oDataIn, simT0, Hc)
        ! Input
        real(dp),intent(in)::t_epoch
        integer, intent(in)::id_body
        real(dp), dimension(:), intent(in)::mass, radii, r1, r2
        real(dp), intent(in)::time_r1, integration_step
        logical, dimension(:), intent(in)::transit_flag
        integer, intent(in)::dur_check
        type(dataObs), intent(in)::oDataIn
        ! Input/Output
        type(dataT0), dimension(:), intent(inout)::simT0
        logical, intent(inout)::Hc
        ! Local
        real(dp)::tmidtra, ttra_temp, lte, duration, r_sky, Rs, Rp, Rmin, Rmax
        real(dp), dimension(:), allocatable::rtra
        real(dp), dimension(4)::tcont
        integer::sel_r

        tmidtra = zero
        duration = zero
        tcont = zero
        allocate (rtra(NBDIM))
        rtra = r1


        call find_transit(t_epoch, id_body, mass, r1, r2, time_r1, integration_step, ttra_temp, lte, rtra, Hc)

        if (Hc) then

            call Rbounds(id_body, radii, Rs, Rp, Rmin, Rmax)
            sel_r = (id_body-1)*6
            r_sky = rsky(rtra(1+sel_r:2+sel_r))

            if (r_sky .le. Rmax) then ! planet transits the star (b <= 1)

                if (transit_flag(id_body)) then ! planet has to transit!

                    tmidtra = ttra_temp
                    if (dur_check .eq. 1) then
                        call find_contacts(id_body, mass, radii, rtra, ttra_temp, tcont) !tcont take into account LTE
                        duration = (tcont(4)-tcont(1))*1440.0_dp
                    end if
                    if (oDataIn%obsT0(id_body-1)%nT0 .gt. 0) then
                        call assign_T0_byNumber(id_body, mass, rtra, oDataIn, tmidtra, duration, simT0)
                    end if
                else ! ... but it should not!

                    Hc = .false.

                end if ! transit_flag

            end if ! r_sky

        end if ! Hc

        deallocate (rtra)

        return
    end subroutine check_T0_2

    ! IT DETERMINES WHICH IS THE RIGHT T_0,obs TO BE ASSOCIATED WITH
    ! THE SIMULATED T_0,sim = tmidtra
    ! v2
    subroutine assign_T0_byNumber(id_body, mass, rtra, oDataIn, tmidtra, duration, simT0)
        ! Input
        integer, intent(in)::id_body
        real(dp), dimension(:), intent(in)::mass, rtra
        type(dataObs), intent(in)::oDataIn
        real(dp), intent(in)::tmidtra, duration
        ! Input/Output
        type(dataT0), dimension(:), intent(inout)::simT0
        ! Local
        real(dp)::dT, dTP, Tref, Pref
        integer::ntmid, i_n, nTs, ibd, epo

        ibd = id_body-1
        nTs = oDataIn%obsT0(ibd)%nT0

        Tref = oDataIn%obsT0(ibd)%Tephem
        Pref = oDataIn%obsT0(ibd)%Pephem

        dT = tmidtra-Tref
        dTP = dT/Pref
        ntmid = nint(dTP)

        do i_n = 1, nTs
            epo = oDataIn%obsT0(ibd)%epo(i_n)
            if (ntmid .eq. epo) then
                simT0(ibd)%epo(i_n) = ntmid
                simT0(ibd)%T0(i_n) = tmidtra
                simT0(ibd)%T0_stat(i_n) = 1
                simT0(ibd)%nT0 = simT0(ibd)%nT0+1
                
                simT0(ibd)%dur(i_n) = duration
                simT0(ibd)%dur_stat(i_n) = 1
                simT0(ibd)%nDur = simT0(ibd)%nDur+1
                ! compute orbital elements of selected T0 and body
                call elements_one_body(id_body, mass, rtra,&
                    &simT0(ibd)%period(i_n), simT0(ibd)%sma(i_n),&
                    &simT0(ibd)%ecc(i_n), simT0(ibd)%inc(i_n), simT0(ibd)%meanA(i_n),&
                    &simT0(ibd)%argp(i_n), simT0(ibd)%trueA(i_n),&
                    &simT0(ibd)%longN(i_n), simT0(ibd)%dttau(i_n))

            end if
        end do

        return
    end subroutine assign_T0_byNumber

    ! it finds all transits of the selected planet (id_body) and store them
    ! in storetra variable, ready to be written into file
    subroutine all_transits(t_epoch, pos, id_body, mass, radii, r1, r2, time_r1, integration_step, stat_tra, storetra)
        ! Input
        real(dp), intent(in)::t_epoch
        integer, intent(in)::pos, id_body
        real(dp), dimension(:), intent(in)::mass, radii, r1, r2
        real(dp), intent(in)::time_r1, integration_step
        integer, dimension(:, :), intent(out)::stat_tra
        real(dp), dimension(:, :), intent(out)::storetra

        real(dp)::ttra, lte, Rs, Rp, Rmin, Rmax, r_sky
        real(dp), dimension(:), allocatable::rtra
        real(dp), dimension(4)::tcont
        logical::status
        integer::sel_r

        ttra = zero
        lte = zero
        tcont = zero
        status = .true.
        allocate (rtra(NBDIM))
        rtra = r1
        
        call find_transit(t_epoch, id_body, mass, r1, r2, time_r1, integration_step, ttra, lte, rtra, status)

        if (status) then

            call Rbounds(id_body, radii, Rs, Rp, Rmin, Rmax)
            sel_r = (id_body-1)*6
            r_sky = rsky(rtra(1+sel_r:2+sel_r))

            if (r_sky .le. Rmax) then ! planet transits the star (b <= 1)

                if (do_transit(id_body)) then ! planet has to transit!

                    call find_contacts(id_body, mass, radii, rtra, ttra, tcont)
                    stat_tra(id_body, pos) = 1
                    storetra(1, pos)       = ttra
                    storetra(2, pos)       = lte
                    storetra(3:6, pos)     = tcont
                    storetra(7:, pos)      = rtra

                else

                    status = .false.

                end if ! do_transit/transit_flag

            end if ! r_sky

        end if ! status

        deallocate (rtra)

        return
    end subroutine all_transits
    ! ------------------------------------------------------------------ !

! ==============================================================================
! compute the transit time and proper duration from K10 eq. 15
    subroutine transit_time(t_epoch, id_body, mass, radii, r1, r2, time_r1, integration_step, ttra, dur_tra, check_ttra)
        ! Input
        real(dp), intent(in)::t_epoch
        integer, intent(in)::id_body ! value: 2 to NB
        real(dp), dimension(:), intent(in)::mass, radii, r1, r2
        real(dp), intent(in)::time_r1, integration_step
        real(dp), intent(out)::ttra, dur_tra
        logical, intent(out)::check_ttra

        real(dp)::ttra_temp, lte
        real(dp), dimension(4)::tcont
        real(dp), dimension(:), allocatable::rtra

        real(dp)::Rs, Rp, Rmin, Rmax, r_sky
        integer::sel_r

        check_ttra = .true.
        ttra = -9.e10_dp
        dur_tra = -9.e10_dp

        allocate (rtra(NBDIM))
        rtra = one

        call find_transit(t_epoch, id_body, mass, r1, r2, time_r1, integration_step, ttra_temp, lte, rtra, check_ttra)

        if (check_ttra) then

            call Rbounds(id_body, radii, Rs, Rp, Rmin, Rmax)
            sel_r = (id_body-1)*6
            r_sky = rsky(rtra(1+sel_r:2+sel_r))

            if (r_sky .le. Rmax) then ! planet transits the star (b <= 1)

                if (do_transit(id_body)) then ! planet has to transit!

                    ttra = ttra_temp !+ lte
                    call find_contacts(id_body, mass, radii, rtra, ttra_temp, tcont) !tcont take into account LTE
                    dur_tra = (tcont(4)-tcont(1))*1440.0_dp

                else ! ... but it should not!

                    check_ttra = .false.

                end if ! do_transit
            else
                ! does not transit: Rsky > Rs+Rp ==> b > 1
                check_ttra = .false.

            end if ! r_sky

        end if ! check_ttra

        deallocate (rtra)

        return
    end subroutine transit_time

! ==============================================================================

end module transits
