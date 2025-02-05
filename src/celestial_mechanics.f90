module celestial_mechanics
    use constants
    use parameters
    use rotations
    use sorting, only: indexx
    implicit none

    interface dist
        module procedure dist_1, dist_2
    end interface dist

contains

! ------------------------------------------------------------------------------
!   function to compute the semi-major axis of an orbit from Period
    elemental subroutine period_to_sma(mass_star, mass_planet, period, sma)
        ! Input
        real(dp), intent(in)::mass_star, mass_planet, period
        ! Output
        real(dp), intent(out)::sma

        real(dp), parameter::twopi_square = dpi*dpi
        real(dp)::mu, P2

        mu = Giau*(mass_star+mass_planet)
        P2 = period*period
        sma = ((mu*P2)/twopi_square)**(onethird)

        return
    end subroutine period_to_sma
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   function to compute the period given the semi-major axis
    elemental subroutine sma_to_period(mass_star, mass_planet, sma, period)
        ! Input
        real(dp), intent(in)::mass_star, mass_planet, sma
        ! Output
        real(dp), intent(out)::period
        real(dp)::mu, sma3

        mu = Giau*(mass_star+mass_planet)
        sma3 = sma**3
        period = dpi*sqrt(sma3/mu)

        return
    end subroutine sma_to_period
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    ! time pericentre to mean anomaly in rad
    elemental subroutine pericenter_time_to_mean_anomaly_rad(tau, t_ref, period, mean_anomaly_rad)
        ! Input
        real(dp), intent(in)::tau, t_ref, period
        ! Output
        real(dp), intent(out)::mean_anomaly_rad

        mean_anomaly_rad = mod(((dpi/period)*(t_ref-tau)), dpi)

        return
    end subroutine pericenter_time_to_mean_anomaly_rad

    ! time pericentre to mean anomaly in deg
    elemental subroutine pericenter_time_to_mean_anomaly_deg(tau, t_ref, period, mean_anomaly)
        ! Input
        real(dp), intent(in)::tau, t_ref, period
        ! Output
        real(dp), intent(out)::mean_anomaly
        ! Local
        real(dp)::meana_rad

        call pericenter_time_to_mean_anomaly_rad(tau, t_ref, period, meana_rad)
        mean_anomaly = meana_rad*rad2deg

        return
    end subroutine pericenter_time_to_mean_anomaly_deg
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    ! mean anomaly in rad to time pericentre
    elemental subroutine mean_anomaly_rad_to_pericenter_time(mean_anomaly_rad, t_ref, period, tau)
        ! Input
        real(dp), intent(in)::mean_anomaly_rad, t_ref, period
        ! Output
        real(dp), intent(out)::tau

        tau = t_ref-(mean_anomaly_rad*period/dpi)

        return
    end subroutine mean_anomaly_rad_to_pericenter_time

    ! ------------------------------------------------------------------------------
    ! mean anomaly in deg to time pericentre in rad
    elemental subroutine mean_anomaly_deg_to_pericenter_time(mean_anomaly_deg, t_ref, period, tau)
        ! Input
        real(dp), intent(in)::mean_anomaly_deg, t_ref, period
        ! Output
        real(dp), intent(out)::tau
        ! Local
        real(dp)::mean_anomaly_rad

        mean_anomaly_rad = mean_anomaly_deg*deg2rad
        call mean_anomaly_rad_to_pericenter_time(mean_anomaly_rad, t_ref, period, tau)

        return
    end subroutine mean_anomaly_deg_to_pericenter_time
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   calculates Eccentric anomaly from meanAnomaly [deg] and eccentricity
    function EAnom(mA, ecc) result(EA)
        real(dp)::EA
        real(dp), intent(IN)::mA, ecc
        real(dp)::mArad, E, fE, dfE

        real(dp), parameter::TOL_check = TOL_dp ! TOLERANCE
        integer, parameter::maxcount = 100
        integer::icount

        EA = zero
        mArad = mod(mA*deg2rad+dpi, dpi)
        ! if (mArad .lt. zero) mArad = mArad+dpi
        icount = 0
        E = mArad
        if (ecc .gt. TOL_check) then
            ! if (ecc .gt. 0.6_dp) E = pi
            if (ecc .ge. 0.8_dp) E = pi
            loopo: do
                icount = icount+1
                fE = E-ecc*sin(E)-mArad
                dfE = one-ecc*cos(E)
                EA = E-(fE/dfE)
                if ((abs(E-EA) .le. TOL_check) .or. (icount .gt. maxcount)) then
                    exit loopo
                end if
                E = EA
            end do loopo
            EA = EA*rad2deg
        else
            ! ecc == 0
            EA = mA
        end if

        return
    end function EAnom
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    function calculate_true_anomaly(mean_anomaly, eccentricity) result(true_anomaly)
        real(dp)::true_anomaly

        real(dp), intent(in)::mean_anomaly, eccentricity

        real(dp)::EA, tan_EA, ecc_coeff

        if (eccentricity .le. TOLERANCE) then
            true_anomaly = mean_anomaly*deg2rad
        else
            EA = EAnom(mod(mean_anomaly+circ, circ), eccentricity)*deg2rad
            tan_EA = tan(EA*half)
            ecc_coeff = sqrt((one+eccentricity)/(one-eccentricity))
            true_anomaly = two*atan(ecc_coeff*tan_EA) ! output in rad!!!
        end if

        return
    end function calculate_true_anomaly
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    function trueAnom_ecc_to_eccAnom(trueAnom, ecc) result(eccAnom)
        real(dp)::eccAnom

        real(dp), intent(in)::trueAnom, ecc

        real(dp)::tan_htA, ecoeff

        if (ecc .le. TOLERANCE) then
            eccAnom = trueAnom
        else
            tan_htA = tan(half*trueAnom*deg2rad)
            ecoeff = sqrt((one-ecc)/(one+ecc))
            eccAnom = mod((two*atan(ecoeff*tan_htA))+dpi, dpi) ! output in rad!!!
        end if

        return
    end function trueAnom_ecc_to_eccAnom
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   calculate the module of a 3-D vector
    function dist_1(r) result(out)
        real(dp)::out
        real(dp), dimension(3), intent(in)::r
        real(dp)::x, y, z

        x = r(1)
        y = r(2)
        z = r(3)

        out = sqrt(x*x+y*y+z*z)

        return
    end function dist_1
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   function to compute the distance between to vectors
    function dist_2(r1, r2) result(out)
        real(dp)::out
        real(dp), dimension(3), intent(in)::r1, r2
        real(dp)::dx, dy, dz

        dx = r1(1)-r2(1)
        dy = r1(2)-r2(2)
        dz = r1(3)-r2(3)
        out = sqrt(dx*dx+dy*dy+dz*dz)

        return
    end function dist_2
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   module of a 2-D vector
    function rsky(r) result(out)
        real(dp)::out
        real(dp), dimension(2), intent(in)::r
        real(dp)::x, y

        x = r(1)
        y = r(2)
        out = sqrt(x*x+y*y)

        return
    end function rsky
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   calculate module of Rvec x Vvec
    function rvprod(r, v) result(out)
        real(dp)::out
        real(dp), dimension(3), intent(in)::r, v
        real(dp)::x, y, z, vx, vy, vz

        x = r(1)
        y = r(2)
        z = r(3)
        vx = v(1)
        vy = v(2)
        vz = v(3)

        out = x*vx+y*vy+z*vz

        return
    end function rvprod
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   function to compute the Hill radius w/o eccentricity
    function rHill_circ(ms, mp, sma) result(rH)
        real(dp)::rH
        real(dp), intent(in)::ms, mp, sma

        rH = sma*((onethird*(mp/ms))**onethird)

        return
    end function rHill_circ
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   function to compute the Hill radius w/ eccentricity
    function rHill_ecc(ms, mp, sma, ecc) result(rH)
        real(dp)::rH
        real(dp), intent(in)::ms, mp, sma, ecc
        real(dp)::e1

        e1 = one-ecc
        rH = sma*e1*((onethird*(mp/ms))**onethird)

        return
    end function rHill_ecc
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! mutual Hill radius
    function mutual_Hill_radius(ms, mp1, sma1, mp2, sma2) result(rH)
        real(dp)::rH
        real(dp), intent(in)::ms, mp1, sma1, mp2, sma2
        real(dp)::sma_mean, mass_ratio

        real(dp)::min_ratio

        rH = zero
        sma_mean = half*(sma1+sma2)
        mass_ratio = (mp1+mp2)/ms
        min_ratio = TOL_dp**onethird

        if (mass_ratio .le. TOL_dp) then

            rH = sma_mean*min_ratio
            write (*, '(2(a,es23.16),a)') ' mass_ratio = ', mass_ratio, ' <= ', TOL_dp, ' = TOL_dp ==> USING TOL_dp'
            write (*, '(3(a,es23.16))') ' Mstar [Msun] = ', ms, ' Mi [Msun]= ', mp1, ' Mj [Msun]= ', mp2
            flush (6)

        else if (ms .le. TOL_dp) then

            write (*, '(2(a,es23.16),a)') ' Mstar = ', ms, ' <= ', TOL_dp, ' = TOL_dp'
            write (*, '(3(a,es23.16))') ' Mstar [Msun] = ', ms, ' Mi [Msun]= ', mp1, ' Mj [Msun]= ', mp2
            flush (6)
            rH = 1.e10_dp

        else if (sma_mean .le. TOL_dp) then

            write (*, '(2(a,es23.16),a)') ' sma1 = ', sma1, ' au & sma2 = ', sma2, ' au'
            write (*, '(2(a,es23.16),a)') ' sma_mean = ', sma_mean, ' <= ', TOL_dp, ' = TOL_dp'
            write (*, '(3(a,es23.16))') ' Mstar [Msun] = ', ms, ' Mi [Msun]= ', mp1, ' Mj [Msun]= ', mp2
            flush (6)
            rH = 1.0e10_dp

        else

            rH = sma_mean*((onethird*mass_ratio)**onethird)

        end if

        return
    end function mutual_Hill_radius
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    subroutine close_encounters_collision(mass, radius, svin, hill_check, not_colliding)
        ! Input
        real(dp), dimension(:), intent(in)::mass, radius
        real(dp), dimension(:), intent(in)::svin
        logical, intent(in)::hill_check
        ! Output
        logical, intent(out)::not_colliding

        real(dp), dimension(3)::rvec1, rvec2
        ! real(dp),dimension(3)::vvec1,vvec2
        real(dp)::rsum, dr
        ! real(dp)::dv
        integer::i_body, j_body
        integer::nc1, nc2

        ! real(dp), dimension(:), allocatable::period, sma, ecc, inc, meanA, argp, trueA, longN, dttau
        real(dp)::sma1, ecc1, sma2, ecc2
        real(dp)::Hill_radius_ij, delta_ij, stability_criterion

        not_colliding = .true.

        body1: do i_body = 2, NB-1
            nc1 = (i_body-1)*6
            rvec1 = svin(nc1+1:nc1+3)
            rsum = (radius(1)+radius(i_body))*RsunAU
            ! |r1 - (Rstar+Rplanet1)| <= 1 cm?
            if ((abs(dist(rvec1))-rsum) .le. cm_au) then
                not_colliding = .false.
                exit body1
            end if

            if (hill_check) then
                call statevector_to_sma_ecc(mass(1), mass(i_body), svin(nc1+1:nc1+6), sma1, ecc1)
            end if

            body2: do j_body = i_body+1, NB
                nc2 = (j_body-1)*6
                rvec2 = svin(nc2+1:nc2+3)
                dr = dist(rvec1, rvec2)
                rsum = (radius(i_body)+radius(j_body))*RsunAU
                ! |(|r1-r2|) - (Rplanet1+Rplanet2)| <= 1 cm?
                if ((abs(dr)-rsum) .le. cm_au) then
                    not_colliding = .false.
                    exit body1
                end if
                if (hill_check) then
                    call statevector_to_sma_ecc(mass(1), mass(j_body), svin(nc2+1:nc2+6), sma2, ecc2)
                    Hill_radius_ij = mutual_Hill_radius(mass(1), mass(i_body), sma1, mass(j_body), sma2)
                    delta_ij = abs(sma2-sma1)
                    stability_criterion = (delta_ij/Hill_radius_ij)-sqrt_12
                    if (stability_criterion .le. TOL_dp) then
                        not_colliding = .false. ! if false is bad/unstable
                        exit body1
                    end if
                end if
            end do body2
        end do body1

        return
    end subroutine close_encounters_collision

    function separation_mutual_Hill_check(mass, radius, svin, do_hill_check) result(hill_check)
        ! Output
        logical::hill_check
        ! Input
        real(dp), dimension(:), intent(in)::mass, radius
        real(dp), dimension(:), intent(in)::svin
        logical, intent(in)::do_hill_check

        call close_encounters_collision(mass, radius, svin, do_hill_check, hill_check)

        return
    end function separation_mutual_Hill_check

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    ! other subroutine to compute the initial state vector in the sky frame
    ! from the keplerian elements
    subroutine kepelements2statevector(mass, sma, ecc, meanA, argp, inc, longN, rout)
        real(dp), dimension(:), intent(in)::mass, sma, ecc, meanA, argp, inc, longN
        real(dp), dimension(:), intent(out)::rout

        integer::ibd, nci
        real(dp)::mu, trueA, cosf, sinf, oneme2, rval, vval
        real(dp), dimension(:), allocatable::rtemp
        ! for the rotation matrices
        ! real(dp), dimension(3, 3)::Rargp, Rinc, RlN

        rout = zero
        allocate (rtemp(size(rout)))
        rtemp = zero

        do ibd = 2, NB
            ! mu = G(Mstar+Mplane)
            mu = Giau*(mass(1)+mass(ibd))
            trueA = calculate_true_anomaly(meanA(ibd), ecc(ibd)) ! true anomaly -> f
            cosf = cos(trueA)
            sinf = sin(trueA)
            oneme2 = one-(ecc(ibd)*ecc(ibd))
            rval = sma(ibd)*oneme2/(one+ecc(ibd)*cosf) ! r = a (1-e^2)/ (1+e*cosf)
!       vval = sqrt((mu/sma(ibd))/oneme2) ! v = sqrt(G(ms+mp)/a/(1-e^2))
            vval = sqrt(mu/(sma(ibd)*oneme2)) ! v = sqrt(G(ms+mp)/a/(1-e^2))

            ! Tx,y,z (vx,vy,vz) --> rotations --> X,Y,Z (VX,VY,VZ)
            nci = (ibd-1)*6
            rtemp(1+nci) = rval*cosf ! x = rcosf
            rtemp(2+nci) = rval*sinf ! y = rsinf
            rtemp(4+nci) = -vval*sinf ! vx = v*(-sinf)
            rtemp(5+nci) = vval*(ecc(ibd)+cosf) ! vy = v*(e+cosf)

            ! with matrix multiplication
            ! Rargp = zero
            ! call rotmat3(argp(ibd), Rargp)
            ! Rinc = zero
            ! call rotmat1(inc(ibd), Rinc)
            ! RlN = zero
            ! call rotmat3(longN(ibd), RlN)
            ! rout(1+nci:3+nci) = matmul(RlN, matmul(Rinc, matmul(Rargp, rtemp(1+nci:3+nci))))
            ! rout(4+nci:6+nci) = matmul(RlN, matmul(Rinc, matmul(Rargp, rtemp(4+nci:6+nci))))
            ! ! with pre-defined rotations
            call rotate_vector(rtemp(1+nci:3+nci), argp(ibd), inc(ibd), longN(ibd), rout(1+nci:3+nci))
            call rotate_vector(rtemp(4+nci:6+nci), argp(ibd), inc(ibd), longN(ibd), rout(4+nci:6+nci))
        end do
!     call orb2obs(rtemp,longN,inc,argp,rout)
        deallocate (rtemp)

        return
    end subroutine kepelements2statevector
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    ! computes the barycenter of the system
    subroutine barycenter(m, ri, rbar, ro)
        ! Input
        real(dp), dimension(:), intent(in)::m, ri
        ! Output
        real(dp), dimension(:), intent(out)::rbar, ro
        ! Local
        real(dp)::mtot
        integer::nbodies, j, nj

        nbodies  = size(m)
        mtot = sum(m)
        rbar = zero
        ro = zero
        do j = 1, nbodies
            nj = (j-1)*6
            rbar = rbar+ri(1+nj:6+nj)*m(j)
        end do
        rbar = rbar/mtot
        !compute the position of the star and other bodies respect to the barycenter
        do j = 1, nbodies
            nj = (j-1)*6
            ro(1+nj:6+nj) = ri(1+nj:6+nj)-rbar
        end do

        return
    end subroutine barycenter
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    ! move a body coordinates (x,y,z,vx,vy,vz) of time = dt using the
    ! f and g functions (Murray and Dermott 1999)
    subroutine fgfunctions(mu, rout, dt, Hc)
        real(dp), intent(in)::mu, dt
        real(dp), dimension(:), intent(inout)::rout
        logical, intent(inout)::Hc

        real(dp), dimension(3)::rtemp, vtemp
        real(dp)::r0, v0
        real(dp)::pxx, sma, ecc, ixx, mA0, wxx, tA0, lnxx, dtxx
        real(dp)::n, EA0, mA1, EA1, dEA
        real(dp)::ar0, art
        real(dp)::Ffunc, dFfunc, Gfunc, dGfunc


        call elem_mer(mu, rout, pxx, sma, ecc, ixx, mA0, wxx, lnxx, tA0, dtxx) ! in rad!
        if (sma .le. TOL_dp .or. ecc .gt. one) then
            Hc = .false.
            return
        end if
        if (ecc .le. TOLERANCE) then
            EA0 = tA0
        else
            EA0 = trueAnom_ecc_to_eccAnom(tA0*rad2deg, ecc)
        end if
        n = dpi/pxx
        mA1 = mA0+dt*n ! not checking if <0 or greater than 2pi, I want to know the 'direction'
        EA1 = EAnom(mA1*rad2deg, ecc)*deg2rad
        dEA = EA1-EA0

        r0 = dist(rout(1:3))
        v0 = dist(rout(4:6))

        ar0 = sma/r0
        Ffunc = ar0*(cos(dEA)-one)+one
        Gfunc = dt+(sin(dEA)-dEA)/n
        rtemp = Ffunc*rout(1:3)+Gfunc*rout(4:6)
        art = sma/dist(rtemp)
        dFfunc = -ar0*art*n*sin(dEA)
        dGfunc = art*(cos(dEA)-one)+one
        vtemp = dFfunc*rout(1:3)+dGfunc*rout(4:6)

        rout(1:3) = rtemp
        rout(4:6) = vtemp

        return
    end subroutine fgfunctions
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    ! computing the singular specific angular momentum
    subroutine angmom(rvin, Lvec)
        real(dp), intent(in)::rvin(6)
        real(dp), intent(out)::Lvec(3)
        real(dp)::x, y, z, vx, vy, vz

        x = rvin(1)
        y = rvin(2)
        z = rvin(3)
        vx = rvin(4)
        vy = rvin(5)
        vz = rvin(6)

        Lvec(1) = y*vz-z*vy
        Lvec(2) = z*vx-x*vz
        Lvec(3) = x*vy-y*vx

        return
    end subroutine angmom
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    ! computes the constant of the motion: Energy and Angular Momentum
    subroutine const_motion(m, rin, Etot, htot)
        real(dp), dimension(:), intent(in)::m
        real(dp), dimension(:), intent(in)::rin
        real(dp), intent(out)::Etot, htot
        real(dp), dimension(3)::hvec, htvec, rj1, rj2, r0
        real(dp), dimension(6)::rvj1
        real(dp)::Ekin, Epot, temp, rj10, vj1, rj21
        integer j1, j2, nj1, nj2

        htot = zero
        Etot = zero
        Ekin = zero
        Epot = zero
        hvec = zero
        htvec = zero
        r0 = rin(1:3)

        ! Total Angular Momentum & kinetic Energy partial
        j1do: do j1 = 1, NB
            nj1 = (j1-1)*6
            temp = zero
            rvj1 = rin(1+nj1:6+nj1)
            rj1 = rvj1(1:3)
            vj1 = dist(rvj1(4:6))
            call angmom(rvj1, hvec)
            htvec = htvec+m(j1)*hvec
            Ekin = Ekin+m(j1)*(vj1*vj1)
            if (j1 .gt. 1) then
                rj10 = dist(rj1, r0)
                j2do: do j2 = j1+1, NB
                    nj2 = (j2-1)*6
                    rj2 = rin(1+nj2:3+nj2)
                    rj21 = dist(rj2, rj1)
                    temp = temp+m(j2)/rj21
                end do j2do
                Epot = Epot-Giau*m(j1)*(temp+m(1)/rj10)
            end if
        end do j1do
        htot = dist(htvec)
        Etot = half*Ekin+Epot

        return
    end subroutine const_motion
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    ! small subroutine to call sequentially the subroutines
    ! needed to update the E and h constant of motion
    subroutine compute_con(m, rin, Etot, htot)
        real(dp), dimension(:), intent(in)::m, rin
        real(dp), intent(inout)::Etot, htot
        real(dp), dimension(6)::bar
        real(dp), dimension(:), allocatable::rbar

        allocate (rbar(NBDIM))
        call barycenter(m, rin, bar, rbar)
        call const_motion(m, rbar, Etot, htot)
        deallocate (rbar)

        return
    end subroutine compute_con
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    subroutine statevector_to_rvh(mu, svin, r, v2, rv, hx, hy, hz, h2, h, s)
        real(dp), intent(in)::mu
        real(dp), dimension(6), intent(in)::svin
        real(dp), intent(out)::r, v2, rv, hx, hy, hz, h2, h, s

        real(dp)::x, y, z, u, v, w

        x = svin(1)
        y = svin(2)
        z = svin(3)
        u = svin(4)
        v = svin(5)
        w = svin(6)

        r = sqrt(x*x+y*y+z*z)
        v2 = u*u+v*v+w*w
        rv = x*u+y*v+z*w

        hx = y*w-z*v
        hy = z*u-x*w
        hz = x*v-y*u
        h2 = hx*hx+hy*hy+hz*hz
        h = sqrt(h2)
        s = h2/mu

        return
    end subroutine statevector_to_rvh

    subroutine statevector_to_sma_ecc(mass_star, mass_planet, svin, sma, ecc)
        real(dp), intent(in)::mass_star, mass_planet
        real(dp), dimension(6), intent(in)::svin
        real(dp), intent(out)::sma, ecc

        ! real(dp)::x,y,z,u,v,w
        real(dp)::r, v2, rv
        real(dp)::hx, hy, hz, h2, h
        real(dp)::s
        real(dp)::mu, v2_to_mu, two_to_r

        mu = Giau*(mass_star+mass_planet)

        call statevector_to_rvh(mu, svin, r, v2, rv, hx, hy, hz, h2, h, s)

        v2_to_mu = v2/mu
        two_to_r = two/r
        ecc = one+s*(v2_to_mu-two_to_r)

        if (ecc .le. TOLERANCE) then
            ecc = zero
            sma = s
        else
            ecc = sqrt(ecc)
            sma = s/(one-(ecc*ecc))
        end if

        if (ecc .gt. one) sma = -sma ! hyperbola added by Luca

        return
    end subroutine statevector_to_sma_ecc
! ------------------------------------------------------------------------------

! ADAPTED FROM MERCURY
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_X2EL.FOR    (ErikSoft  20 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates Keplerian orbital elements given relative coordinates and
! velocities, and GM = G times the sum of the masses.
!
! The elements are: q = perihelion distance
!                   ecc = eccentricity
!                   inc = inclination
!                   p = longitude of perihelion (NOT argument of perihelion!!)
!                   ln = longitude of ascending node
!                   mA = mean anomaly (or mean longitude if e < 1.e-8)
!
!-------------------------------------------------------------------------------
!
    subroutine mco_x2el(mu, svec, q, ecc, inc, p, ln, mA, f)
        implicit none
        ! include 'mercury.inc'

        ! Input/Output
        real(dp), intent(in)::mu
        real(dp), intent(in), dimension(:)::svec
!     real(dp),intent(inout)::q,ecc,inc,p,ln,mA,f
        real(dp), intent(out)::q, ecc, inc, p, ln, mA, f

        ! Local
        real(dp)::x, y !, z, u, v, w
        real(dp)::hx, hy, hz, h2, h, v2, r, rv, s, true
        real(dp)::ci, t_o, temp, tmp2, bige, cf, ce

        ! init all output variables
        q = zero
        ecc = zero
        inc = zero
        p = zero
        ln = zero
        mA = zero
        f = zero

        x = svec(1)
        y = svec(2)

        call statevector_to_rvh(mu, svec, r, v2, rv, hx, hy, hz, h2, h, s)

        ! Inclination and node
        ci = hz/h
        if (abs(ci) .lt. one) then
            inc = acos(ci)
            ln = atan2(hx, -hy)
            if (ln .lt. zero) ln = ln+dpi
        else
!       if (ci.gt.zero) inc = zero
!       if (ci.lt.zero) inc = pi
            if (ci .gt. zero) then
                inc = zero
            else
                inc = pi
            end if
            ln = pi
        end if

        ! Eccentricity and perihelion distance
        ecc = zero
        temp = one+s*((v2/mu)-(two/r))

        ! if (temp.le.zero) then
        if (temp .le. TOL_dp) then
            ecc = zero
            q = s
        else
            ecc = sqrt(temp)
            q = s/(one+ecc)
            if (ecc .gt. one) q = -q ! hyperbola added by Luca
        end if

        ! True longitude
        if (hy .ne. zero) then
            t_o = -hx/hy
            temp = (one-ci)*t_o
            tmp2 = t_o*t_o
            true = mod(atan2((y*(one+tmp2*ci)-x*temp), (x*(tmp2+ci)-y*temp))+dpi, dpi)
        else
            true = mod(atan2(y*ci, x)+dpi, dpi)
        end if
        if (ci .lt. zero) true = mod(true+pi, dpi)

        if (ecc .le. TOLERANCE) then
            ! p = zero
            p = -half*pi
            mA = true
            f = mA
        else
            ce = (v2*r-mu)/(ecc*mu)

            ! Mean anomaly for ellipse
            if (ecc .lt. one) then
                if (abs(ce) .gt. one) ce = sign(one, ce)
                bige = acos(ce)
                if (rv .lt. zero) bige = dpi-bige
                mA = bige-ecc*sin(bige)
            else

                ! Mean anomaly for hyperbola
                if (ce .lt. one) ce = one
                bige = log(ce+sqrt(ce*ce-one))
                if (rv .lt. zero) bige = -bige
                mA = ecc*sinh(bige)-bige
            end if

            ! Longitude of perihelion
            cf = (s-r)/(ecc*r)
            if (abs(cf) .gt. one) cf = sign(one, cf)
            f = acos(cf)
            if (rv .lt. zero) f = dpi-f
            p = true-f
            p = mod(p+dpi+dpi, dpi)
        end if

        if (mA .lt. zero .and. ecc .lt. one) mA = mA+dpi
        if (mA .gt. dpi .and. ecc .lt. one) mA = mod(mA, dpi)

        return
    end subroutine mco_x2el
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! call the mco_x2el and 'completes' the kepler elements
    subroutine elem_mer(mu, svec, P, sma, ecc, inc, mA, w, lN, f, dt)
        real(dp), intent(in)::mu
        real(dp), intent(in), dimension(:)::svec
        real(dp), intent(out)::P, sma, ecc, inc, mA, w, lN, f, dt

        real(dp)::q, lperi
!     real(dp)::cf

        q = zero
        lperi = -half*pi
!     cf=zero

        call mco_x2el(mu, svec, q, ecc, inc, lperi, lN, mA, f)
        ! semi-major axis sma
        sma = q/(one-ecc)
        ! arg. pericentre w
        w = mod(lperi-lN+dpi, dpi)
        if (ecc .le. TOLERANCE) w = half*pi

        P = dpi*sqrt((sma**3)/mu)

        dt = mA*P/dpi

!     if(ecc.ge.one)then
!       write(*,*)'mu = ',mu
!       write(*,*)'svec = ',svec
!       write(*,*)'P,sma,ecc,inc,mA,w,lN,f,dt'
!       write(*,*)P,sma,ecc,inc,mA,w,lN,f,dt
!     end if

    end subroutine elem_mer
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! from state vector to orbital elements of a specified body
    subroutine elements_one_body(i_body, mass, rin, period, sma, ecc, inc, meanA, argp, trueA, longN, dttau)
        integer, intent(in)::i_body
        real(dp), dimension(:), intent(in)::mass, rin
        real(dp), intent(out)::period, sma, ecc, inc, meanA, argp, trueA, longN, dttau

        real(dp)::mu
        integer::sel

        ! init
        period = zero
        sma = zero
        ecc = zero
        inc = 90.0_dp
        meanA = zero ! mean anomaly
        argp = 90.0_dp ! argument of pericentre
        trueA = zero ! true anomaly
        longN = 180.0_dp ! longitude of node
        dttau = zero

        sel = (i_body-1)*6
        mu = Giau*(mass(1)+mass(i_body))

        call elem_mer(mu, rin(1+sel:6+sel),&
                &period, sma, ecc, inc,&
                &meanA, argp, longN, trueA,&
                &dttau)
        inc = inc*rad2deg
        meanA = meanA*rad2deg
        longN = longN*rad2deg
        argp = argp*rad2deg
        trueA = trueA*rad2deg

        return
    end subroutine elements_one_body
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    ! from state vector (astrocentric cartesian coordinates)
    ! to keplerian orbital elements
    ! for all the bodies
    subroutine elements(mass, rin, period, sma, ecc, inc, meanA, argp, trueA, longN, dttau)
        real(dp), dimension(:), intent(in)::mass, rin
        real(dp), dimension(:), intent(out)::period, sma, ecc, inc, meanA, argp, trueA, longN, dttau
        ! real(dp), dimension(6)::svec
        ! real(dp)::mu

        integer::n_body, i_body
        ! integer::ncj

        ! init all orbital elements
        ! period = zero
        ! sma = zero
        ! ecc = zero
        ! inc = 90.0_dp
        ! meanA = zero ! mean anomaly
        ! argp = 90.0_dp ! argument of pericentre
        ! trueA = zero ! true anomaly
        ! longN = 180.0_dp ! longitude of node
        ! dttau = zero

        n_body = size(mass)

        period(1) = zero
        sma(1) = zero
        ecc(1) = zero
        inc(1) = zero
        meanA(1) = zero ! mean anomaly
        argp(1) = zero ! argument of pericentre
        trueA(1) = zero ! true anomaly
        longN(1) = zero ! longitude of node
        dttau(1) = zero

        cicle: do i_body = 2, n_body
            call elements_one_body(i_body, mass, rin,&
                &period(i_body), sma(i_body), ecc(i_body), inc(i_body),&
                &meanA(i_body), argp(i_body), trueA(i_body), longN(i_body),&
                &dttau(i_body))
        end do cicle

        return
    end subroutine elements
! ------------------------------------------------------------------------------

    ! stability criterion based on dis/eq. 26 of Petit et al. 2018
    subroutine amd_hill(alpha, gamma, epsilon, c_hill)
        ! Input
        ! alpha = sma1 / sma2
        ! gamma = mass1 / mass2
        ! epsilon = (mass1 + mass2) / mass_star
        ! planet 1 inner, 2 outer
        real(dp), intent(in)::alpha, gamma, epsilon
        ! Output
        real(dp), intent(out)::c_hill
        ! Local
        real(dp)::term1,term2,term3,term4,term5


        term1 = gamma * sqrt(alpha) + one
        term2 = (one+gamma)**(three/two)
        term3 = alpha / (gamma + alpha)
        term4 = (three**(4.0_dp/three)) * (epsilon**(two/three)) * gamma
        term5 = (one + gamma)**two
        c_hill = term1 - term2*sqrt(term3*(one + (term4/term5)))

        return
    end subroutine amd_hill

    ! Relative AMD from eq. 23 of Petit et al. 2018
    subroutine relative_amd(alpha, gamma, ecc_in, inc_in, ecc_out, inc_out, amd_rel)
        ! Input
        ! alpha = sma_in / sma_out
        ! gamma = mass_in / mass_out
        real(dp), intent(in)::alpha, gamma, ecc_in, inc_in, ecc_out, inc_out
        ! Output
        real(dp), intent(out)::amd_rel
        ! Local
        real(dp)::cosin,cosout,term1,term2,term3

        cosin  = cos((inc_in -90.0_dp)*deg2rad)
        cosout = cos((inc_out-90.0_dp)*deg2rad)
        term1  = gamma*sqrt(alpha)
        term2  = sqrt(one-(ecc_in*ecc_in))*cosin
        term3  = sqrt(one-(ecc_out*ecc_out))*cosout

        amd_rel = term1*(one-term2) + one - term3

        return
    end subroutine relative_amd

    subroutine amd_hill_pair(mstar, mpl_inn, sma_in, ecc_in, inc_in, mpl_out, sma_out, ecc_out, inc_out, amd_h, amd_r)
        ! Input
        real(dp), intent(in)::mstar, mpl_inn, sma_in, ecc_in, inc_in, mpl_out, sma_out, ecc_out, inc_out
        ! Output
        real(dp), intent(out)::amd_h, amd_r
        ! Local
        real(dp)::gamma,alpha,epsilon

        gamma = mpl_inn / mpl_out
        alpha = sma_in / sma_out
        epsilon = (mpl_inn + mpl_out) / mstar
        call relative_amd(alpha, gamma, ecc_in, inc_in, ecc_out, inc_out, amd_r)
        call amd_hill(alpha, gamma, epsilon, amd_h)

        return
    end subroutine amd_hill_pair

    subroutine amd_hill_stability(mass, ecc, inc, periods, amd_r_pairs, amd_h_pairs, stable)
        ! Input
        real(dp), dimension(:), intent(in)::mass, ecc, inc, periods
        ! Output
        real(dp), dimension(:), intent(out)::amd_r_pairs, amd_h_pairs
        logical, intent(out)::stable
        ! Local
        real(dp)::amd_h, amd_r
        integer::n_body, i_sort, i_in, i_out, i_pairs
        integer, dimension(:), allocatable::sort_periods
        logical, dimension(:), allocatable::check

        n_body = size(mass)
        allocate(sort_periods(n_body), check(size(amd_r_pairs)))
        call indexx(periods, sort_periods)

        i_pairs = 1
        do i_sort = 2,n_body-1
            i_in = sort_periods(i_sort)
            i_out = sort_periods(i_sort+1)
            call amd_hill_pair(mass(1),&
                &mass(i_in), periods(i_in),ecc(i_in), inc(i_in),&
                &mass(i_out), periods(i_out),ecc(i_out), inc(i_out),&
                &amd_h, amd_r)
            amd_r_pairs(i_pairs) = amd_r
            amd_h_pairs(i_pairs) = amd_h
            check(i_pairs) = amd_r < amd_h
            i_pairs = i_pairs + 1
        end do

        stable = all(check)

        deallocate(sort_periods, check)

        return
    end subroutine amd_hill_stability

    subroutine statevector_amd_hill_stability(mass, sv, stable)
        ! Input
        real(dp), dimension(:), intent(in)::mass
        real(dp), dimension(:), intent(in)::sv
        ! Output
        logical, intent(out)::stable
        ! Local
        integer::n_body, n_pairs
        real(dp), dimension(:), allocatable::period, sma, ecc, inc, meanA, argp, trueA, longN, dttau
        real(dp), dimension(:), allocatable::ar, ah

        n_body = size(mass)
        n_pairs = n_body - 2
        allocate(period(n_body), sma(n_body), ecc(n_body), inc(n_body), meanA(n_body), argp(n_body), trueA(n_body), longN(n_body), dttau(n_body))
        call elements(mass, sv, period, sma, ecc, inc, meanA, argp, trueA, longN, dttau)
        allocate(ar(n_pairs), ah(n_pairs))
        call amd_hill_stability(mass, ecc, inc, period, ar, ah, stable)

        return
    end subroutine statevector_amd_hill_stability

    subroutine amd_one_planet(mstar, mplanet, ecc, inc, sma, lambda, amd_planet)
        ! Input
        real(dp), intent(in)::mstar, mplanet, ecc, inc, sma
        ! Output
        real(dp), intent(out)::lambda, amd_planet
        ! Local
        real(dp):: mratio, msum, cosi, oneecc2
        real(dp),parameter::d2yr = 365.25_dp
        real(dp),parameter::factor = one !10.0e10_dp ! Laskar 1997 Caption Fig.1

        msum = mstar + mplanet
        mratio = (mplanet*mstar)/msum

        oneecc2 = one - (ecc*ecc)
        cosi = cos((inc-90.0_dp)*deg2rad) ! Definition of inc on the reference plane

        lambda = mratio * sqrt(Giau * msum * sma) * d2yr * factor

        amd_planet = lambda * (one - (sqrt(oneecc2) * cosi))

        return 
    end subroutine amd_one_planet

    subroutine amd_system(masses, eccs, incs, periods, lambdas, amd_bodies, amd)
        ! Input
        real(dp), dimension(:), intent(in)::masses, eccs, incs, periods
        ! Output
        real(dp), dimension(:), intent(out):: lambdas, amd_bodies
        real(dp), intent(out)::amd
        ! Local
        integer::n_body, i_planet, i_sort
        real(dp)::amd_i, sma
        integer, dimension(:), allocatable::sort_periods

        !! WARNING !!
        ! to change with:
        ! Laskar & Petit 2017 eq. 29 with/without cos(inc)
        ! Compute Amd Hill stability
        ! Petit 2018 eq. 26 or 30
        ! AMD - AMD Hill stable if AMD < AMD Hill for each pairs!
        ! DO IT FOR PAIRS AND SORT IN PERIOD OR SMA

        n_body = size(masses)

        allocate(sort_periods(n_body))
        call indexx(periods, sort_periods)
        amd = zero
        lambdas(1) = zero
        amd_bodies(1) = zero
        do i_sort = 2,n_body
            i_planet = sort_periods(i_sort)
            call period_to_sma(masses(1), masses(i_planet), periods(i_planet), sma)
            call amd_one_planet(masses(1), masses(i_planet), eccs(i_planet),&
                &incs(i_planet), sma,&
                &lambdas(i_planet), amd_i)
            amd = amd + amd_i
            amd_bodies(i_planet) = amd_i
        end do
        deallocate(sort_periods)

        return
    end subroutine amd_system
! ------------------------------------------------------------------------------

end module celestial_mechanics
