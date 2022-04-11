module ode_run
    use constants
    use custom_type
    use parameters
    use parameters_conversion
    use linear_ephem
    use celestial_mechanics
    !   use rotations,only:orb2obs
    use eq_motion, only: eqmastro
    use numerical_integrator, only: int_rk_a, rkck_a
    use transits
    use radial_velocities
    use utils
    use gls_module, only: check_periodogram, check_periodogram_scale, check_and_write_periodogram
    use output_files
    use sorting, only: indexx
    use statistics, only: mean
    implicit none

contains

    ! subroutine to initialise the variables needed to check coordinates for the transit events.
    subroutine set_checking_coordinates(n_body,radius,idxX,idxY,idxZ,idxVX,idxVY,checkX,checkY,checkR,check_coord,XYmean)
        ! **Input**
        integer,intent(in)::n_body
        real(dp), dimension(:),intent(in)::radius
        ! **Output**
        integer, dimension(:), allocatable, intent(out)::idxX, idxY, idxZ
        integer, dimension(:), allocatable, intent(out)::idxVX, idxVY
        real(dp), dimension(:), allocatable, intent(out)::checkX, checkY, checkR ! variables to check if there is a transit or not
        real(dp), dimension(:), allocatable, intent(out)::check_coord, XYmean ! variables to check if there is a transit or not
        
        integer::i_body

        allocate (idxX(n_body), idxY(n_body), idxZ(n_body), idxVX(n_body), idxVY(n_body))
        allocate (checkX(n_body), checkY(n_body), checkR(n_body)) !, rmean(n_body))
        allocate (check_coord(n_body), XYmean(n_body))

        idxX = 0
        idxY = 0
        idxZ = 0
        idxVX = 0
        idxVY = 0
        do i_body = 2, n_body
            idxX(i_body) = 1 + (i_body - 1)*6
            idxY(i_body) = 2 + (i_body - 1)*6
            idxZ(i_body) = 3 + (i_body - 1)*6
            idxVX(i_body) = 4 + (i_body - 1)*6
            idxVY(i_body) = 5 + (i_body - 1)*6
        end do
        checkX = one
        checkY = one
        checkR = zero
        checkR(2:n_body) = 1.5_dp*(radius(1) + radius(2:n_body))*RsunAU
        check_coord = zero
        XYmean = one

        return
    end subroutine set_checking_coordinates

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
        real(dp),intent(inout)::stepsize
        logical, intent(inout)::Hc
        ! **Output**
        real(dp), dimension(:), intent(out)::rout
        ! **local**
        real(dp)::working_step, ok_step, next_step, itime
        real(dp), dimension(:), allocatable::dr, r1, r2, err
        ! if you see a variable non listed here, probably it is a global variable
        ! defined in constants.f90 or parameters.f90

        Hc = separation_mutual_Hill_check(mass, radius, rin, do_hill_check)
        if (.not. Hc) return


        working_step=stepsize
        if (time .lt. zero) working_step = -working_step
        itime = zero
        allocate (dr(NBDIM), r1(NBDIM), r2(NBDIM), err(NBDIM))
        r1 = rin
        r2 = zero
        dr = zero
        err = zero

        integration: do

            if (abs(itime + working_step) .gt. abs(time)) working_step = time - itime
            call eqmastro(mass, r1, dr)
            call int_rk_a(mass, r1, dr, working_step, ok_step, next_step, r2, err)
            Hc = separation_mutual_Hill_check(mass, radius, r2, do_hill_check)
            
            itime = itime + ok_step
            if (.not. Hc) exit integration
            if (abs(itime) .ge. abs(time)) exit integration
            working_step = next_step
            r1 = r2
        end do integration

        rout = r2
        stepsize = next_step
        deallocate(dr, r1, r2, err)

        return
    end subroutine one_forward_step


    ! ================================================================================
    ! performs the forward integration in one direction.
    ! computes the orbits and RVs and T0s (T14s) as in
    ! obsRV.dat and NB#_observations.dat files (stored in obsData).
    !
    ! **Input**
    ! mass   == masses of all bodies (star == 1)
    ! radius == radius of all bodies
    ! rin    == stave vector of all bodies in astrocentric coordinates XYZVxVyVz
    ! time   == integration time
    ! check_longN == for each body it specifies if check T0 along X o Y axis
    !
    ! **Input/Output**
    ! simRV  == dataRV type with simulated RV
    ! simT0  == array of dataT0 type with simulated T0 (T14), for each planet
    !
    ! **Output**
    ! Hc     == Hill chek, .true. if stable, .false. if unstable for Hill's condition
    subroutine ode_forward_data(mass, radius, rin, time, check_longN,simRV, simT0, Hc)
        ! **Input**
        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::time
        integer, dimension(:), intent(in)::check_longN
        ! **Input/Output**
        type(dataRV), intent(inout)::simRV
        type(dataT0), dimension(:), intent(inout)::simT0
        logical, intent(out)::Hc
        ! **local**
        real(dp), dimension(:), allocatable::dr, r1, r2, err
        integer, dimension(:), allocatable::X, Y, Z
        ! real(dp), dimension(:), allocatable::cX, cY, cR, rmean ! variables to check if there is a transit or not
        real(dp), dimension(:), allocatable::cX, cY, cR ! variables to check if there is a transit or not
        integer, dimension(:), allocatable::VX, VY
        real(dp), dimension(:), allocatable::check_coord, XYmean ! variables to check if there is a transit or not
        real(dp)::hw, hok, hnext, itime
        integer::i_body, iteration, body_transiting_start,body_transiting_end
        integer::nRV, nTTs, nDurs
        logical::do_transit_check

        real(dp)::Tr,Pr,Tx
        integer::epox
        
        ! if you see a variable non listed here, probably it is a global variable
        ! defined in constants.f90 or parameters.f90

        ! obsData global variable
        nRV   = obsData%obsRV%nRV
        nTTs  = obsData%nTTs
        nDurs = obsData%nDurs

        Hc = separation_mutual_Hill_check(mass, radius, rin, do_hill_check)
        if (.not. Hc) return

        ! == NEW ==
        call set_checking_coordinates(NB,radius,X,Y,Z,VX,VY,cX,cY,cR,check_coord,XYmean)
        ! ========

        allocate(dr(NBDIM), r1(NBDIM), r2(NBDIM), err(NBDIM))
        hw = step_0 ! working stepsize -> initial value based on global variable step_0
        if (time .lt. zero) hw = -hw
        itime = zero
        r1 = rin
        r2 = zero
        err = zero

        ! == NEW ==
        ! by default it will check the transit of all the planets (2->NB)
        do_transit_check=.true.
        body_transiting_start=2
        body_transiting_end=NB
        if((idtra.lt.1).or.(idtra.gt.NB))then
            do_transit_check=.false.
        else
            if(idtra.gt.1)then
                body_transiting_start=idtra
                body_transiting_end=idtra
            end if
        end if
        ! ========

        iteration = 0
        integration: do
            iteration = iteration + 1

            if (abs(itime + hw) .gt. abs(time)) hw = time - itime
            call eqmastro(mass, r1, dr)
            call int_rk_a(mass, r1, dr, hw, hok, hnext, r2, err)

            Hc = separation_mutual_Hill_check(mass, radius, r2, do_hill_check)
            if (.not. Hc) then
                return
            end if

            ! RV check
            if (nRV .gt. 0) then
                if (simRV%nRV .lt. nRV) then
                    call check_RV(mass, r1, dr, itime, hok, simRV)
                end if
            end if

            ! == NEW ==
            ! ! T0 check
            Tx = tepoch + itime + half*hok
            do i_body = body_transiting_start, body_transiting_end
                Tr = obsData%obsT0(i_body-1)%Tephem
                Pr = obsData%obsT0(i_body-1)%Pephem
                if (Pr .gt. zero) then ! only if we have Pephem > 0

                    cX(i_body) = r1(X(i_body))*r2(X(i_body))
                    cY(i_body) = r1(Y(i_body))*r2(Y(i_body))

                    if (check_longN(i_body) .eq. 0) then
                        check_coord(i_body) = cX(i_body)
                        XYmean(i_body) = half*(r1(Y(i_body))+r2(Y(i_body)))
                    else
                        check_coord(i_body) = cY(i_body)
                        XYmean(i_body) = half*(r1(X(i_body))+r2(X(i_body)))
                    end if
                    if(check_coord(i_body) .le. zero)then ! change of the sign of X or Y coordinate of two consecutive steps
                        if(XYmean(i_body).le.cR(i_body))then ! check if X or Y mean coordinates of two consecutive steps is less or equal the 1.5 x (Rplanet + Rstar)
                            if((r1(Z(i_body)).gt.zero) .or. (r2(Z(i_body)).gt.zero))then ! check if Z coordinates of one of two consecutives steps is positive, that is planet in front of the star
                                
                                epox = nint(((Tx-Tr)/Pr))
                                if(any(epox .eq. obsData%obsT0(i_body-1)%epo))then ! check if the epoch of the mean time of two consecutive steps is within the observed epochs
                                    call check_T0(i_body, mass, radius, r1, r2,&
                                        &itime, hok, simT0, Hc)
                                end if
                            end if
                        end if
                    end if

                end if
            end do ! i_body
            ! ========

            itime = itime + hok
            if (abs(itime) .ge. abs(time)) exit integration
            hw = hnext
            r1 = r2
        end do integration

        ! deallocate(X, Y, Z, cX, cY, cR, rmean)
        deallocate(X, Y, Z, VX, VY, cX, cY, cR)
        deallocate(dr, r1, r2, err)

        return
    end subroutine ode_forward_data
    ! ================================================================================

    ! ================================================================================
    ! subroutine that integrates orbit, computes, stores and
    ! writes orbits, orbital elements, constants of motion, transit time and RV into files.
    ! it computes all RV and T0 (T14) and those from obsRV.dat and NB#_observations.dat files (stored in obsData)
    ! it is called by ode_output, that check in which direction (in time) integrates
    ! called by:
    ! ode_output
    !
    ! **Input**
    ! uorb,ucon,uele,utra == units of the file to write: orbits, constants, elements, transits
    ! fmorb,fmcon,fmele   == format to write the files: orbits, constants, elements
    ! mass                == mass of the bodies
    ! radius              == radius of the bodies
    ! rin                 == state vector of all the bodies in astrocentric coordinates
    ! time                == integration time
    ! check_longN         == for each body it specifies if check T0 along X o Y axis
    !
    ! **Input/Output**
    ! simRV               == dataRV type with simulated RV
    ! simT0               == array of dataT0 type with simulated T0 (T14), for each planet
    ! Hc                  == Hill chek, .true. if stable, .false. if unstable for Hill's condition
    subroutine ode_forward_output(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
      &mass, radius, rin, time, check_longN, simRV, simT0, Hc)
        ! **Input**
        integer, intent(in)::uorb, ucon
        integer, dimension(:), intent(in)::uele, utra
        character(*), intent(in)::fmorb, fmcon, fmele
        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::time
        integer, dimension(:), intent(in)::check_longN
        ! **Input/Output**
        type(dataRV), intent(inout)::simRV
        type(dataT0), dimension(:), intent(inout)::simT0
        logical, intent(inout)::Hc
        ! local variables
        real(dp), dimension(:), allocatable::dr, r1, r2, err
        real(dp), dimension(:), allocatable::drdt_save, r_save
        integer, dimension(:), allocatable::X, Y, Z
        ! real(dp), dimension(:), allocatable::cX, cY, cR, rmean
        real(dp), dimension(:), allocatable::cX, cY, cR ! variables to check if there is a transit or not
        integer, dimension(:), allocatable::VX, VY
        real(dp), dimension(:), allocatable::check_coord, XYmean ! variables to check if there is a transit or not

        integer::i_body, body_transiting_start, body_transiting_end

        real(dp)::hw, hok, hnext, itime
        integer::iteration, save_iteration, jtra
        real(dp), dimension(:, :), allocatable::storeorb, storecon, storetra
        integer, dimension(:, :), allocatable::stat_tra
        integer::Norb
        real(dp)::itwrt, Etot, Eold, htot, hold

        real(dp)::step_save, step_write
        logical::do_transit_check

        real(dp)::Tr,Pr,Tx
        integer::epox

        integer::nRV
        ! if you see a variable non listed here, probably it is a global variable
        ! defined in constants.f90 or parameters.f90

        nRV = obsData%obsRV%nRV

        Hc = separation_mutual_Hill_check(mass, radius, rin, do_hill_check)
        if (.not. Hc) then
            return
        end if

        call set_checking_coordinates(NB,radius,X,Y,Z,VX,VY,cX,cY,cR,check_coord,XYmean)

        Norb = NBDIM + 3

        allocate (dr(NBDIM), r1(NBDIM), r2(NBDIM), err(NBDIM))
        allocate (drdt_save(NBDIM), r_save(NBDIM))
        hw = step_0
        if (time .lt. zero) hw = -hw
        itime = zero
        r1 = rin
        r2 = zero
        err = zero

        iteration = 0
        save_iteration = 1

        step_write = wrttime
        if (time .lt. zero) step_write = -step_write
        itwrt = step_write

        if ((wrtorb .eq. 1) .or. (wrtel .eq. 1)) then
            allocate (storeorb(Norb, DIMMAX))
            storeorb = zero
            call store_orb(save_iteration, itime, mass, r1, storeorb)
        end if

        if (wrtconst .eq. 1) then
            Etot = zero
            Eold = zero
            htot = zero
            hold = zero
            call compute_con(mass, r1, Eold, hold)
            allocate (storecon(5, DIMMAX))
            storecon = zero
            call store_con(save_iteration, itime, Eold, Eold, hold, hold, storecon)
        end if

        if ((idtra .ge. 1) .and. (idtra .le. NB)) then
            allocate (storetra(NBDIM + 6, DIMMAX))
            storetra = zero
            allocate (stat_tra(NB, DIMMAX))
            stat_tra = 0
            jtra = 0
        end if

        ! == NEW ==
        ! by default it will check the transit of all the planets (2->NB)
        do_transit_check=.true.
        body_transiting_start=2
        body_transiting_end=NB
        if((idtra.lt.1).or.(idtra.gt.NB))then
            do_transit_check=.false.
        else
            if(idtra.gt.1)then
                body_transiting_start=idtra
                body_transiting_end=idtra
            end if
        end if
        ! ========

        integration: do
            iteration = iteration + 1

            if (abs(itime + hw) .gt. abs(time)) hw = time - itime
            call eqmastro(mass, r1, dr) ! computes the eq. of motion
            call int_rk_a(mass, r1, dr, hw, hok, hnext, r2, err) ! computes the next orbit step

            ! check if it has to compute or not 'something'
            if ((wrtorb .eq. 1) .or. (wrtel .eq. 1) .or. (wrtconst .eq. 1)) then

                if (abs(itime + hok) .ge. (abs(itwrt))) then

                    saveloop: do

                        save_iteration = save_iteration + 1
                        ! computes the proper step -> itwrt will be updated in the loop
                        step_save = itwrt - itime
                        ! computes state vector from r1 to the new time step (smaller then hok)
                        drdt_save = dr
                        r_save = r1
                        call rkck_a(mass, r1, drdt_save, step_save, r_save, err)

                        ! check and store orbit or kep. elements
                        if ((wrtorb .eq. 1) .or. (wrtel .eq. 1)) then
                            call store_orb(save_iteration, itime + step_save, mass, r_save, storeorb) ! save r_save!!
                            if (save_iteration .eq. DIMMAX) then ! write into file if reach max dimension
                                if (wrtorb .eq. 1) call write_file(DIMMAX, uorb, fmorb, storeorb)
                                if (wrtel .eq. 1) call write_elem(DIMMAX, uele, fmele, mass, storeorb)
                                storeorb = zero ! reset the storeorb variable
                            end if
                        end if

                        ! check and store constants of motion
                        if (wrtconst .eq. 1) then
                            call compute_con(mass, r_save, Etot, htot) ! if selected it computes the Energy and Angular momentum
                            call store_con(save_iteration, itime + step_save, Etot, Eold, htot, hold, storecon)
                            if (save_iteration .eq. DIMMAX) then ! write into file if reach max dimension
                                call write_file(DIMMAX, ucon, fmcon, storecon)
                                storecon = zero
                            end if
                        end if

                        if (save_iteration .eq. DIMMAX) save_iteration = 0

                        itwrt = itwrt + step_write
                        if (abs(itwrt) .gt. abs(itime + hok)) exit saveloop

                    end do saveloop

                end if

            end if
            ! ===============================

            Hc = separation_mutual_Hill_check(mass, radius, r2, do_hill_check)

            if (.not. Hc) then
                return
            end if

            ! RV check
            if (nRV .gt. 0) then
                if (simRV%nRV .lt. nRV) then
                    call check_RV(mass, r1, dr, itime, hok, simRV)
                end if
            end if

            ! T0 check (to compare and all)
            Tx = tepoch + itime + half*hok
            do i_body = body_transiting_start, body_transiting_end

                Tr = obsData%obsT0(i_body-1)%Tephem
                Pr = obsData%obsT0(i_body-1)%Pephem
                if (Pr .gt. zero )then
                    cX(i_body) = r1(X(i_body))*r2(X(i_body))
                    cY(i_body) = r1(Y(i_body))*r2(Y(i_body))

                    if (check_longN(i_body) .eq. 0) then
                        check_coord(i_body) = cX(i_body)
                        XYmean(i_body) = half*(r1(Y(i_body))+r2(Y(i_body)))
                    else
                        check_coord(i_body) = cY(i_body)
                        XYmean(i_body) = half*(r1(X(i_body))+r2(X(i_body)))
                    end if
                    if(check_coord(i_body) .le. zero)then ! change of the sign of X or Y coordinate of two consecutive steps
                        if(XYmean(i_body).le.cR(i_body))then ! check if X or Y mean coordinates of two consecutive steps is less or equal the 1.5 x (Rplanet + Rstar)
                            if((r1(Z(i_body)).gt.zero) .or. (r2(Z(i_body)).gt.zero))then ! check if Z coordinates of one of two consecutives steps is positive, that is planet in front of the star
                                
                                epox = nint(((Tx-Tr)/Pr))
                                if(any(epox .eq. obsData%obsT0(i_body-1)%epo))then ! check if the epoch of the mean time of two consecutive steps is within the observed epochs
                                    call check_T0(i_body, mass, radius, r1, r2,&
                                        &itime, hok, simT0, Hc)
                                end if
                                jtra = jtra + 1
                                call all_transits(jtra, i_body, mass, radius, r1, r2,&
                                    &itime, hok, stat_tra, storetra)

                            end if
                        end if

                        if (jtra .eq. DIMMAX) then
                            call write_tra(jtra, utra, stat_tra, storetra)
                            jtra = 0
                            stat_tra = 0
                            storetra = zero
                        end if

                    end if
                end if

            end do ! i_body
            ! ========

            itime = itime + hok

            if (abs(itime) .ge. abs(time)) exit integration
            hw = hnext
            r1 = r2

        end do integration

        if ((idtra .ge. 1) .and. (idtra .le. NB)) call write_tra(jtra, utra, stat_tra, storetra)

        if ((wrtorb .eq. 1) .or. (wrtel .eq. 1)) then
            if (save_iteration .gt. 0) then
                if (wrtorb .eq. 1) call write_file(save_iteration, uorb, fmorb, storeorb)
                if (wrtel .eq. 1) call write_elem(save_iteration, uele, fmele, mass, storeorb)
            end if
            deallocate (storeorb)
        end if
        if (wrtconst .eq. 1) then
            if (save_iteration .gt. 0) then
                call write_file(save_iteration, ucon, fmcon, storecon)
            end if
            deallocate (storecon)
        end if

        ! deallocate(X, Y, Z, cX, cY, cR, rmean)
        deallocate(X, Y, Z, VX, VY, cX, cY, cR)
        deallocate (dr, r1, r2, err, drdt_save, r_save)

        return
    end subroutine ode_forward_output
    ! ================================================================================

    ! ================================================================================
    ! subroutine that integrates orbit, computes, stores and
    ! writes orbits, orbital elements, constants of motion, transit time and RV into files
    ! it doesn't need obsRV.dat or NB#_observations.dat files,
    ! it computes only all the possible transits
    ! and the RVs are stored in the #_#_rotorbit.dat file.
    !
    ! it is called by:
    ! ode_integrates, that check in which direction (in time) integrates
    !
    ! **Input**
    ! uorb,ucon,uele,utra == units of the file to write: orbits, constants, elements, transits
    ! fmorb,fmcon,fmele   == format to write the files: orbits, constants, elements
    ! mass                == mass of the bodies
    ! radius              == radius of the bodies
    ! rin                 == state vector of all the bodies in astrocentric coordinates
    ! time                == integration time
    ! check_longN         == for each body it specifies if check T0 along X o Y axis
    !
    ! **Input/Output**
    ! Hc                  == Hill chek, .true. if stable, .false. if unstable for Hill's condition
    subroutine ode_forward_output_nodata(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
        &mass, radius, rin, time, check_longN, Hc)
        integer, intent(in)::uorb, ucon
        integer, dimension(:), intent(in)::uele, utra
        character(*), intent(in)::fmorb, fmcon, fmele
        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::time
        integer, dimension(:), intent(in)::check_longN

        logical, intent(inout)::Hc

        real(dp), dimension(:), allocatable::dr, r1, r2, err
        real(dp), dimension(:), allocatable::drdt_save, r_save
        integer, dimension(:), allocatable::X, Y, Z
        ! real(dp), dimension(:), allocatable::cX, cY, cR, rmean
        real(dp), dimension(:), allocatable::cX, cY, cR ! variables to check if there is a transit or not
        integer, dimension(:), allocatable::VX, VY
        real(dp), dimension(:), allocatable::check_coord, XYmean ! variables to check if there is a transit or not

        integer::i_body, body_transiting_start, body_transiting_end

        real(dp)::hw, hok, hnext, itime
        integer::iteration, save_iteration, jtra
        real(dp), dimension(:, :), allocatable::storeorb, storecon, storetra
        integer, dimension(:, :), allocatable::stat_tra
        integer::Norb
        real(dp)::itwrt, Etot, Eold, htot, hold

        real(dp)::step_save, step_write
        logical::do_transit_check

        Hc = separation_mutual_Hill_check(mass, radius, rin, do_hill_check)
        if (.not. Hc) then
            return
        end if

        call set_checking_coordinates(NB,radius,X,Y,Z,VX,VY,cX,cY,cR,check_coord,XYmean)

        Norb = NBDIM + 3

        allocate (dr(NBDIM), r1(NBDIM), r2(NBDIM), err(NBDIM))
        allocate (drdt_save(NBDIM), r_save(NBDIM))
        hw = step_0
        if (time .lt. zero) hw = -hw
        itime = zero
        r1 = rin
        r2 = zero
        err = zero

        iteration = 0
        save_iteration = 1

        step_write = wrttime
        if (time .lt. zero) step_write = -step_write
        itwrt = step_write

        if ((wrtorb .eq. 1) .or. (wrtel .eq. 1)) then
            allocate (storeorb(Norb, DIMMAX))
            storeorb = zero
            call store_orb(save_iteration, itime, mass, r1, storeorb)
        end if

        if (wrtconst .eq. 1) then
            Etot = zero
            Eold = zero
            htot = zero
            hold = zero
            call compute_con(mass, r1, Eold, hold)
            allocate (storecon(5, DIMMAX))
            storecon = zero
            call store_con(save_iteration, itime, Eold, Eold, hold, hold, storecon)
        end if

        if ((idtra .ge. 1) .and. (idtra .le. NB)) then
            allocate (storetra(NBDIM + 6, DIMMAX))
            storetra = zero
            allocate (stat_tra(NB, DIMMAX))
            stat_tra = 0
            jtra = 0
        end if

        ! == NEW ==
        ! by default it will check the transit of all the planets (2->NB)
        do_transit_check=.true.
        body_transiting_start=2
        body_transiting_end=NB
        if((idtra.lt.1).or.(idtra.gt.NB))then
            do_transit_check=.false.
        else
            if(idtra.gt.1)then
                body_transiting_start=idtra
                body_transiting_end=idtra
            end if
        end if
        ! ========

        integration: do
            iteration = iteration + 1

            if (abs(itime + hw) .gt. abs(time)) hw = time - itime
            call eqmastro(mass, r1, dr) ! computes the eq. of motion
            call int_rk_a(mass, r1, dr, hw, hok, hnext, r2, err) ! computes the next orbit step

            ! ===============================
            ! NEW VERSION
            ! check if it has to compute or not 'something'
            if ((wrtorb .eq. 1) .or. (wrtel .eq. 1) .or. (wrtconst .eq. 1)) then

                if (abs(itime + hok) .ge. (abs(itwrt))) then

                    saveloop: do

                        save_iteration = save_iteration + 1
                        ! computes the proper step -> itwrt will be updated in the loop
                        step_save = itwrt - itime
                        ! computes state vector from r1 to the new time step (smaller then hok)
                        drdt_save = dr
                        r_save = r1
                        call rkck_a(mass, r1, drdt_save, step_save, r_save, err)

                        ! check and store orbit or kep. elements
                        if ((wrtorb .eq. 1) .or. (wrtel .eq. 1)) then
                            call store_orb(save_iteration, itime + step_save, mass, r_save, storeorb) ! save r_save!!
                            if (save_iteration .eq. DIMMAX) then ! write into file if reach max dimension
                                if (wrtorb .eq. 1) call write_file(DIMMAX, uorb, fmorb, storeorb)
                                if (wrtel .eq. 1) call write_elem(DIMMAX, uele, fmele, mass, storeorb)
                                storeorb = zero ! reset the storeorb variable
                            end if
                        end if

                        ! check and store constants of motion
                        if (wrtconst .eq. 1) then
                            call compute_con(mass, r_save, Etot, htot) ! if selected it computes the Energy and Angular momentum
                            call store_con(save_iteration, itime + step_save, Etot, Eold, htot, hold, storecon)
                            if (save_iteration .eq. DIMMAX) then ! write into file if reach max dimension
                                call write_file(DIMMAX, ucon, fmcon, storecon)
                                storecon = zero
                            end if
                        end if

                        if (save_iteration .eq. DIMMAX) save_iteration = 0

                        itwrt = itwrt + step_write
                        if (abs(itwrt) .gt. abs(itime + hok)) exit saveloop

                    end do saveloop

                end if

            end if
            ! ===============================

            Hc = separation_mutual_Hill_check(mass, radius, r2, do_hill_check)

            if (.not. Hc) then
                return
            end if

            ! T0 check (to compare and all)
            do i_body = body_transiting_start, body_transiting_end
                cX(i_body) = r1(X(i_body))*r2(X(i_body))
                cY(i_body) = r1(Y(i_body))*r2(Y(i_body))

                if (check_longN(i_body) .eq. 0) then
                    check_coord(i_body) = cX(i_body)
                    XYmean(i_body) = half*(r1(Y(i_body))+r2(Y(i_body)))
                else
                    check_coord(i_body) = cY(i_body)
                    XYmean(i_body) = half*(r1(X(i_body))+r2(X(i_body)))
                end if
                if(check_coord(i_body) .le. zero)then ! change of the sign of X or Y coordinate of two consecutive steps
                    if(XYmean(i_body).le.cR(i_body))then ! check if X or Y mean coordinates of two consecutive steps is less or equal the 1.5 x (Rplanet + Rstar)
                        if((r1(Z(i_body)).gt.zero) .or. (r2(Z(i_body)).gt.zero))then ! check if Z coordinates of one of two consecutives steps is positive, that is planet in front of the star
                            
                            jtra = jtra + 1
                            call all_transits(jtra, i_body, mass, radius, r1, r2,&
                                &itime, hok, stat_tra, storetra)

                        end if
                    end if

                    if (jtra .eq. DIMMAX) then
                        call write_tra(jtra, utra, stat_tra, storetra)
                        jtra = 0
                        stat_tra = 0
                        storetra = zero
                    end if

                end if
            end do ! i_body
            ! ========

            itime = itime + hok

            if (abs(itime) .ge. abs(time)) exit integration
            hw = hnext
            r1 = r2

        end do integration

        if ((idtra .ge. 1) .and. (idtra .le. NB)) call write_tra(jtra, utra, stat_tra, storetra)

        if ((wrtorb .eq. 1) .or. (wrtel .eq. 1)) then
            if (save_iteration .gt. 0) then
                if (wrtorb .eq. 1) call write_file(save_iteration, uorb, fmorb, storeorb)
                if (wrtel .eq. 1) call write_elem(save_iteration, uele, fmele, mass, storeorb)
            end if
            deallocate (storeorb)
        end if
        if (wrtconst .eq. 1) then
            if (save_iteration .gt. 0) then
                call write_file(save_iteration, ucon, fmcon, storecon)
            end if
            deallocate (storecon)
        end if

        ! deallocate(X, Y, Z, cX, cY, cR, rmean)
        deallocate(X, Y, Z, VX, VY, cX, cY, cR)
        deallocate (dr, r1, r2, err, drdt_save, r_save)

        return
    end subroutine ode_forward_output_nodata
    ! ================================================================================

! ================================================================================
    ! subroutine that given keplerian elements integrates orbits and computes RV and T0 (T14)
    ! based on global information, such tepoch, tstart, time_int and obsData (etc.)
    ! called by ode_lm
    !
    ! **Input**
    ! mass   == masses of the bodies in M_sun
    ! radius == radii of the bodies in R_sun
    ! period == period of the bodies in days
    ! ecc    == eccentricities of the bodies
    ! argp   == argument of pericenter (or periastro) of the bodies in deg
    ! meanA  == mean anomaly of the bodies in deg
    ! inc    == orbital inclination of the bodies in deg, a planet that transits exactly at mid
    !           of the stellar disk has an inclination of 90 deg
    ! longN  == longitude of the ascending Node of the bodies in deg.
    !           Suggested to 180 deg for a ref. planet.
    !
    ! **Output**
    ! simRV
    ! simT0
    ! Hc
    subroutine ode_keplerian_elements_to_data(mass, radius, period, ecc, argp, meanA, inc, longN,&
      &simRV, simT0, Hc)
        ! Input
        real(dp), dimension(:), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN
        ! Output
        type(dataRV), intent(out)::simRV
        type(dataT0), dimension(:), allocatable, intent(out)::simT0
        logical, intent(out)::Hc

        ! local variables
        integer::nRV, nTTs, nT0, ibd
        integer, dimension(:), allocatable::check_longN
        real(dp), dimension(:), allocatable::sma
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp)::dt1, dt2
        ! integer::ns,ne

        ! if not defined here, the variables are global and defined in
        ! constants.f90 or parameters.f90

        Hc = .true.

        ! set RV and T0s data type for simulated ones
        nRV = obsData%obsRV%nRV
        nTTs = obsData%nTTs

        ! RV
        if (nRV .gt. 0) then
            call init_dataRV(nRV, simRV)
            simRV%nRV = 0
            if (.not. allocated(simRV%jitter)) allocate (simRV%jitter(nRVset))
            simRV%jitter = zero
            if (.not. allocated(simRV%gamma)) allocate (simRV%gamma(nRVset))
            simRV%gamma = zero
            simRV%RVsetID = obsData%obsRV%RVsetID
        end if
        ! T0
        if (nTTs .gt. 0) then
            allocate (simT0(NB - 1))
            do ibd = 1, NB - 1
                nT0 = obsData%obsT0(ibd)%nT0
                call init_dataT0(nT0, simT0(ibd), durcheck)
                simT0(ibd)%nT0 = 0
                simT0(ibd)%nDur = 0
                !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
                simT0(ibd)%eT0 = one/s24h ! 1s in day
                simT0(ibd)%edur = one/60.0_dp ! 1s in min
            end do
        end if

        ! it is needed to define which is the alarm coordinate for the transit detection
        call lNset(longN, check_longN)

        NBDIM = 6*NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        allocate (sma(NB))
        sma = zero
        call get_semax_vec(mass(1), mass(2:NB), period(2:NB), sma(2:NB))
        call kepelements2statevector(mass, sma, ecc, meanA, argp, inc, longN, ra0)
        ra1 = ra0
        deallocate (sma)

        dt1 = tstart - tepoch
        dt2 = dt1 + tint
        if (dt1 .lt. zero) then

            call ode_forward_data(mass, radius, ra1, dt1, check_longN, simRV, simT0, Hc)

            if (Hc) then
                if (abs(dt1) .le. tint) then
                    call ode_forward_data(mass, radius, ra1, dt2, check_longN, simRV, simT0, Hc)
                end if
            end if

        else

            call ode_forward_data(mass, radius, ra1, dt2, check_longN, simRV, simT0, Hc)

        end if

        return
    end subroutine ode_keplerian_elements_to_data

    ! ================================================================================

    ! ================================================================================
    ! subroutine called by the L-M to fit the parameters
    !
    ! **Input**
    ! allpar  == all the parameters of the system as one vector of length nm (see also system_parameters, size npar)
    !            it is used set all the parameters fixed and not fitted
    ! nm, nn  == dimension of allpar and par, used previously by examples for lmfit function
    !            but not used here...it could be removed but I am not doing it to avoid issues with lmfit
    ! par     == fitted parameters (not only re-parameterization of keplerian elements, but also RV offset, jitter, and trend)
    !
    ! **Output**
    ! resw    == vector of weighted residuals with size ndata,
    !            it is sum(resw) = chi_square !! changed 2021-11-18
    ! iflag   == integer flag needed by lmfit
    subroutine ode_lm(allpar, nm, nn, par, resw, iflag)
        integer, intent(in)::nm, nn
        real(dp), dimension(:), intent(in)::allpar, par
        real(dp), dimension(:), intent(out)::resw
        integer, intent(inout)::iflag

        real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN
        integer, dimension(:), allocatable::check_longN
        real(dp), dimension(:), allocatable::ra0, ra1
        ! real(dp)::dt1,dt2
        logical::Hc

        type(dataRV)::simRV
        type(dataT0), dimension(:), allocatable::simT0

        logical::checkpar, gls_check

        integer::ndata
        integer::nRV
        integer::nTTs !,nT0
        integer::ibd
        integer::ns, ne

        ! write (*, *) " DEBUG: in ode_lm "

        ndata = obsData%ndata
        nRV = obsData%obsRV%nRV
        ! nRVset=obsData%obsRV%nRVset
        nTTs = obsData%nTTs

        Hc = .true.
        resw = zero
        allocate (mass(NB), radius(NB), period(NB), sma(NB), ecc(NB))
        allocate (argp(NB), meanA(NB), inc(NB), longN(NB), check_longN(NB))

        ! write (*, *) " DEBUG: convert_parameters "
        call convert_parameters(allpar, par, mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
        ! write (*, *) " DEBUG: checkpar ", checkpar
        ! write (*, *) " DEBUG: mass   ", mass
        ! write (*, *) " DEBUG: radius ", radius
        ! write (*, *) " DEBUG: period ", period
        ! write (*, *) " DEBUG: sma    ", sma
        ! write (*, *) " DEBUG: ecc    ", ecc
        ! write (*, *) " DEBUG: argp   ", argp
        ! write (*, *) " DEBUG: meanA  ", meanA
        ! write (*, *) " DEBUG: inc    ", inc
        ! write (*, *) " DEBUG: longN  ", longN

        if (.not. checkpar) then
            resw = set_max_residuals(ndata)
            return
        end if

        ! +++ FROM HERE CALL NEW (2021-11-18) SUBROUTINE ode_elements_to_data
        ! write (*, *) " DEBUG: ode_keplerian_elements_to_data"
        call ode_keplerian_elements_to_data(mass, radius, period, ecc, argp, meanA, inc, longN,&
          &simRV, simT0, Hc)

        ! write (*, *) " DEBUG: Hc ", Hc
        gls_check = .true.
        if (.not. Hc) then
            resw = set_max_residuals(ndata)
        else

            ! write (*, *) " DEBUG: Hc true"
            ! write (*, *) " DEBUG: nRVset ", nRVset
            if (nRVset .gt. 0) then
                ! set jitter into simRV
                ns = nkel + 1
                ne = nkel + nRVset
                simRV%jitter = two**par(ns:ne)
                ! set gamma into simRV
                ns = ne + 1
                ne = ne + nRVset
                simRV%gamma = par(ns:ne)
                call set_gamma_rv(simRV%gamma, simRV%RVsetID, simRV%gamma_rv)
                ! compute RV trend
                if (rv_trend_order .gt. 0) then
                    ! write (*, *) " DEBUG: addRVtrend"
                    ns = ne + 1
                    ne = ne + rv_trend_order
                    call addRVtrend(simRV%jd, tepoch, par(ns:ne), simRV%trend)
                end if
            end if

            ! write (*, *) " DEBUG: sum(resw*resw) = ", sum(resw*resw)
            ! write (*, *) " DEBUG: set_weighted_residuals"
            call set_weighted_residuals(obsData, simRV, simT0, resw, oc_fit)
            ! write (*, *) " DEBUG: sum(resw*resw) = ", sum(resw*resw)

            ! ! write (*, *) " DEBUG: check_periodogram"
            ! if (simRV%nRV .gt. 0) call check_periodogram(obsData%obsRV%jd,&
            !   &obsData%obsRV%RV - simRV%RV - simRV%gamma_rv - simRV%trend,&
            !   &obsData%obsRV%eRV,&
            !   &period, gls_check)

            ! ! write (*, *) " DEBUG: gls_check ", gls_check
            ! if (.not. gls_check) resw = set_max_residuals(ndata)

            ! ! write (*, *) " DEBUG: check_max_residuals"
            ! call check_max_residuals(resw, ndata)

        end if
        ! +++

        call deallocate_dataRV(simRV)
        if (nTTs .gt. 0) then
            do ibd = 1, NB - 1
                call deallocate_dataT0(simT0(ibd))
            end do
        end if

        if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN, check_longN)
        if (allocated(ra0)) deallocate (ra0, ra1)

        ! write (*, *) " DEBUG: LAST ode_lm: sum(resw*resw) = ", sum(resw*resw)

        return
    end subroutine ode_lm
    ! ================================================================================

    ! ================================================================================
    ! subroutines called by the bootstrap to calculate the new set of T_0 and RV from
    ! the fitted parameters
    subroutine ode_parok(allpar, par, simRV, simT0, resw,&
      &chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)
        real(dp), dimension(:), intent(in)::allpar, par

        type(dataRV), intent(inout)::simRV
        type(dataT0), dimension(:), intent(inout)::simT0

        real(dp), dimension(:), intent(out)::resw
        real(dp), intent(out)::chi_square, reduced_chi_square, lnLikelihood, ln_const, bic

        real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN
        integer, dimension(:), allocatable::check_longN
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp)::dt1, dt2
        logical::Hc

        logical::checkpar
        ! logical::gls_check

        integer::ndata, nRV, nTTs, nT0
        integer::ns, ne

        integer::ii, ibd

        ndata = obsData%ndata
        nRV = obsData%obsRV%nRV
        ! nRVset=obsData%obsRV%nRVset
        nTTs = obsData%nTTs

        Hc = .true.
        resw = zero
        chi_square = -one

        allocate (mass(NB), radius(NB), period(NB), sma(NB), ecc(NB),&
            &argp(NB), meanA(NB), inc(NB), longN(NB), check_longN(NB))

        call convert_parameters(allpar, par, mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)

        if (nRV .gt. 0) then
            call init_dataRV(nRV, simRV)
            simRV%nRV = 0
            if (.not. allocated(simRV%jitter)) allocate (simRV%jitter(nRVset))
            simRV%jitter = zero
            if (.not. allocated(simRV%gamma)) allocate (simRV%gamma(nRVset))
            simRV%gamma = zero
            simRV%RVsetID = obsData%obsRV%RVsetID
        end if
!
        if (nTTs .gt. 0) then
            do ibd = 1, NB - 1
                nT0 = obsData%obsT0(ibd)%nT0
                call init_dataT0(nT0, simT0(ibd), durcheck)
                simT0(ibd)%nT0 = 0
                simT0(ibd)%nDur = 0
                !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
                simT0(ibd)%eT0 = one/s24h ! 1s in day
                simT0(ibd)%edur = one/60.0_dp ! 1s in min
            end do
        end if

        if (.not. checkpar) then
            if (nRV .gt. 0) simRV%RV = zero

            if (nTTs .gt. 0) then
                do ibd = 1, NB - 1
                    simT0(ibd)%T0 = zero
                end do
            end if

            resw = set_max_residuals(ndata)
            chi_square = resmax
            call set_fitness_values(par, chi_square,&
              &reduced_chi_square, lnLikelihood, ln_const, bic)
            return
        end if

        do ii = 1, NB
            if (mass(ii) .lt. zero) then
                write (*, '(a,i3,a,10000(1x,f22.16))') ' neg mass(ii=', ii, ' ) : masses = ', mass
                stop
            end if
        end do

        ! it is needed to define the which is the alarm coordinate for the transit detection
        call lNset(longN, check_longN)

        NBDIM = 6*NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(mass, sma, ecc, argp, meanA, inc, longN, ra0)
        ra1 = ra0

        dt1 = tstart - tepoch
        dt2 = dt1 + tint
        if (dt1 .lt. zero) then

            call ode_forward_data(mass, radius, ra1, dt1, check_longN, simRV, simT0, Hc)

            if (Hc) then
                if (abs(dt1) .le. tint) then
                    dt2 = dt1 + tint
                    call ode_forward_data(mass, radius, ra1, dt2, check_longN, simRV, simT0, Hc)
                end if
            end if

        else

            dt2 = dt1 + tint
            call ode_forward_data(mass, radius, ra1, dt2, check_longN, simRV, simT0, Hc)

        end if
        if (.not. Hc) then
            if (nRV .gt. 0) simRV%RV = zero

            if (nTTs .gt. 0) then
                do ibd = 1, NB - 1
                    simT0(ibd)%T0 = zero
                end do
            end if

            resw = set_max_residuals(ndata)
            chi_square = resmax

        else

            if (nRVset .gt. 0) then
                ! set jitter into simRV
                ns = nkel + 1
                ne = nkel + nRVset
                simRV%jitter = two**par(ns:ne)
                ! set gamma into simRV
                ns = ne + 1
                ne = ne + nRVset
                simRV%gamma = par(ns:ne)
                call set_gamma_rv(simRV%gamma, simRV%RVsetID, simRV%gamma_rv)
                ! compute RV trend
                if (rv_trend_order .gt. 0) then
                    ns = ne + 1
                    ne = ne + rv_trend_order
                    call addRVtrend(simRV%jd, tepoch, par(ns:ne), simRV%trend)
                end if
            end if

            call set_weighted_residuals(obsData, simRV, simT0, resw, oc_fit)

            ! if (simRV%nRV .gt. 0) call check_periodogram(obsData%obsRV%jd,&
            !   &obsData%obsRV%RV - simRV%RV - simRV%gamma_rv - simRV%trend,&
            !   &obsData%obsRV%eRV,&
            !   &period, gls_check)

            ! call check_max_residuals(resw, ndata)

            chi_square = sum(resw*resw)
        end if
        call set_fitness_values(par, chi_square,&
            &reduced_chi_square, lnLikelihood, ln_const, bic)

        if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN, check_longN)
        if (allocated(ra0)) deallocate (ra0, ra1)

        return
    end subroutine ode_parok

    ! subroutine called by the L-M to fit the parameters with observations in input (RV and T0)
    subroutine ode_boot(allpar, nm, nn, par, oDataIn, resw, iflag)
        integer, intent(in)::nm, nn
        real(dp), dimension(:), intent(in)::allpar, par
        type(dataObs), intent(in)::oDataIn

        real(dp), dimension(:), intent(out)::resw
        integer, intent(inout)::iflag

        real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN
        integer, dimension(:), allocatable::check_longN
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp)::dt1, dt2
        logical::Hc

        type(dataRV)::simRV
        type(dataT0), dimension(:), allocatable::simT0

        integer::ndata, nRV, nTTs, nT0, ibd
        integer::ns, ne

        logical::checkpar
        ! logical::gls_check

        ndata = oDataIn%ndata
        nRV = oDataIn%obsRV%nRV
        ! nRVset=oDataIn%obsRV%nRVset
        nTTs = oDataIn%nTTs

        Hc = .true.
        resw = zero
        allocate (mass(NB), radius(NB), period(NB), sma(NB), ecc(NB),&
            &argp(NB), meanA(NB), inc(NB), longN(NB), check_longN(NB))

        call convert_parameters(allpar, par, mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
        if (.not. checkpar) then
            resw = set_max_residuals(ndata)
            return
        end if

        ! it is needed to define the which is the alarm coordinate for the transit detection
        call lNset(longN, check_longN)

        NBDIM = 6*NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(mass, sma, ecc, argp, meanA, inc, longN, ra0)
        ra1 = ra0

        if (nRV .gt. 0) then
            call init_dataRV(nRV, simRV)
            simRV%nRV = 0
            if (.not. allocated(simRV%jitter)) allocate (simRV%jitter(nRVset))
            simRV%jitter = zero
            if (.not. allocated(simRV%gamma)) allocate (simRV%gamma(nRVset))
            simRV%gamma = zero
            simRV%RVsetID = obsData%obsRV%RVsetID
        end if

        if (nTTs .gt. 0) then
            allocate (simT0(NB - 1))
            do ibd = 1, NB - 1
                nT0 = obsData%obsT0(ibd)%nT0
                call init_dataT0(nT0, simT0(ibd), durcheck)
                simT0(ibd)%nT0 = 0
                simT0(ibd)%nDur = 0
                !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
                simT0(ibd)%eT0 = one/s24h ! 1s in day
                simT0(ibd)%edur = one/60.0_dp ! 1s in min
            end do
        end if

        dt1 = tstart - tepoch
        if (dt1 .lt. zero) then
            call ode_forward_data(mass, radius, ra1, dt1, check_longN, simRV, simT0, Hc)
            if (Hc) then
                if (abs(dt1) .le. tint) then
                    dt2 = dt1 + tint
                    call ode_forward_data(mass, radius, ra1, dt2, check_longN, simRV, simT0, Hc)
                end if
            end if
        else
            dt2 = dt1 + tint
            call ode_forward_data(mass, radius, ra1, dt2, check_longN, simRV, simT0, Hc)
        end if
        if (.not. Hc) then
            resw = set_max_residuals(ndata)
        else

            if (nRVset .gt. 0) then
                ! set jitter into simRV
                ns = nkel + 1
                ne = nkel + nRVset
                simRV%jitter = two**par(ns:ne)
                ! set gamma into simRV
                ns = ne + 1
                ne = ne + nRVset
                simRV%gamma = par(ns:ne)
                call set_gamma_rv(simRV%gamma, simRV%RVsetID, simRV%gamma_rv)
                ! compute RV trend
                if (rv_trend_order .gt. 0) then
                    ns = ne + 1
                    ne = ne + rv_trend_order
                    call addRVtrend(simRV%jd, tepoch, par(ns:ne), simRV%trend)
                end if
            end if

            call set_weighted_residuals(oDataIn, simRV, simT0, resw, oc_fit)

            ! if (simRV%nRV .gt. 0) call check_periodogram(oDataIn%obsRV%jd,&
            !   &oDataIn%obsRV%RV - simRV%RV - simRV%gamma_rv - simRV%trend,&
            !   &obsData%obsRV%eRV,&
            !   &period, gls_check)

            ! call check_max_residuals(resw, ndata)

        end if

        call deallocate_dataRV(simRV)
        if (nTTs .gt. 0) then
            do ibd = 1, NB - 1
                call deallocate_dataT0(simT0(ibd))
            end do
        end if

        if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN, check_longN)
        if (allocated(ra0)) deallocate (ra0, ra1)

        return
    end subroutine ode_boot
    ! ================================================================================

    subroutine ode_full_args(mass, radius, rin, t_epoch, time, step_in, check_longN, simRV,&
      &id_transit_body, transit_flag, dur_check, simT0, Hc)

        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::t_epoch, time, step_in
        integer, dimension(:), intent(in)::check_longN

        type(dataRV), intent(inout)::simRV

        integer, intent(in)::id_transit_body
        logical, dimension(:), intent(in)::transit_flag
        integer, intent(in)::dur_check

        type(dataT0), dimension(:), intent(inout)::simT0

        logical, intent(inout)::Hc

        integer::n_body, nb_dim
        real(dp), dimension(:), allocatable::dr, r1, r2, err
        integer, dimension(:), allocatable::X, Y, Z
        real(dp), dimension(:), allocatable::cX, cY, cR ! variables to check if there is a transit or not
        integer, dimension(:), allocatable::VX, VY
        real(dp), dimension(:), allocatable::check_coord, XYmean ! variables to check if there is a trans
        real(dp)::hw, hok, hnext, itime
        real(dp)::Tr, Pr, Tx, epox
        integer::i_body, iteration, body_transiting_start,body_transiting_end
        logical::do_transit_check

        integer::nRV, nTTs

        Hc = separation_mutual_Hill_check(mass, radius, rin, do_hill_check)
        if (.not. Hc) return

        nRV = obsData%obsRV%nRV
        nTTs = obsData%nTTs

        n_body = size(mass)
        call set_checking_coordinates(NB,radius,X,Y,Z,VX,VY,cX,cY,cR,check_coord,XYmean)

        nb_dim = 6*n_body
        allocate (dr(nb_dim), r1(nb_dim), r2(nb_dim), err(nb_dim))
        hw = step_in
        if (time .lt. zero) hw = -hw
        itime = zero
        r1 = rin
        r2 = zero
        err = zero

        ! == NEW ==
        ! by default it will check the transit of all the planets (2->NB)
        do_transit_check=.true.
        body_transiting_start=2
        body_transiting_end=n_body
        if((id_transit_body.lt.1).or.(id_transit_body.gt.n_body))then
            do_transit_check=.false.
        else
            if(id_transit_body.gt.1)then
                body_transiting_start=id_transit_body
                body_transiting_end=id_transit_body
            end if
        end if
        ! ========

        iteration = 0
        integration: do
            iteration = iteration + 1

            if (abs(itime + hw) .gt. abs(time)) hw = time - itime
            call eqmastro(mass, r1, dr)
            call int_rk_a(mass, r1, dr, hw, hok, hnext, r2, err)

            Hc = separation_mutual_Hill_check(mass, radius, r2, do_hill_check)
            if (.not. Hc) return

            ! RV check
            if (nRV .gt. 0) then
                if (simRV%nRV .lt. nRV) then
                    call check_RV(mass, r1, dr, itime, hok, simRV)
                end if
            end if

            Tx = t_epoch + itime + half*hok
            do i_body = body_transiting_start, body_transiting_end
                Tr = obsData%obsT0(i_body-1)%Tephem
                Pr = obsData%obsT0(i_body-1)%Pephem
                if (Pr .gt. zero) then
                    
                    cX(i_body) = r1(X(i_body))*r2(X(i_body))
                    cY(i_body) = r1(Y(i_body))*r2(Y(i_body))

                    if (check_longN(i_body) .eq. 0) then
                        check_coord(i_body) = cX(i_body)
                        XYmean(i_body) = half*(r1(Y(i_body))+r2(Y(i_body)))
                    else
                        check_coord(i_body) = cY(i_body)
                        XYmean(i_body) = half*(r1(X(i_body))+r2(X(i_body)))
                    end if
                    if(check_coord(i_body) .le. zero)then ! change of the sign of X or Y coordinate of two consecutive steps
                        if(XYmean(i_body).le.cR(i_body))then ! check if X or Y mean coordinates of two consecutive steps is less or equal the 1.5 x (Rplanet + Rstar)
                            if((r1(Z(i_body)).gt.zero) .or. (r2(Z(i_body)).gt.zero))then ! check if Z coordinates of one of two consecutives steps is positive, that is planet in front of the star
                                
                                epox = nint(((Tx-Tr)/Pr))
                                if(any(epox .eq. obsData%obsT0(i_body-1)%epo))then ! check if the epoch of the mean time of two consecutive steps is within the observed epochs
                                    call check_T0(i_body, mass, radius, r1, r2,&
                                        &itime, hok,&
                                        &transit_flag, dur_check, obsData, simT0, Hc)
                                end if
                            
                            end if
                        end if
                    end if
                end if
            end do ! i_body
            ! ========

            itime = itime + hok
            if (abs(itime) .ge. abs(time)) exit integration
            hw = hnext
            r1 = r2
        end do integration

        deallocate(X, Y, Z, VX, VY, cX, cY, cR)
        deallocate (dr, r1, r2, err)

        return
    end subroutine ode_full_args

! =============================================================================
! =============================================================================
! subroutine that integrates orbits of the planets and computes
! the observables (RV, T0, T14).

! **Input**
! integration information:
! start time (t_start), reference time (t_epoch), initial stepsize (step_in),
! integration time in days (t_int)
! the keplerian elements of each body (star=1) at given time (tepoch):
! mass, radius, period, ecc, argp, meanA, inc, longN
! observations information:
! obsData_in (dataObs type), id number of transiting body to check (id_transit_body),
! flag for the duration T14 check (dur_check)
! **Output**
! simulated observables:
! simulated RV (simRV as dataRV type), simulated transit times (simT0 as dataT0 type for each planet)
! =============================================================================
! =============================================================================
    subroutine integration_info_to_observables(t_start, t_epoch, step_in, t_int,&
      &mass, radius, period, ecc, argp, meanA, inc, longN,&
      &obsData_in,&
      &id_transit_body, transit_flag, dur_check,&
      &simRV, simT0)

        real(dp), intent(in)::t_start, t_epoch, step_in, t_int
        real(dp), dimension(:), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN

        type(dataObs), intent(in)::obsData_in

        integer, intent(in)::id_transit_body
        logical, dimension(:), intent(in)::transit_flag
        integer, intent(in)::dur_check

        type(dataRV), intent(out)::simRV
        type(dataT0), dimension(:), allocatable, intent(out)::simT0

        integer::n_body, nb_dim
        integer, dimension(:), allocatable::check_longN
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp), dimension(:), allocatable::sma
        real(dp)::dt1, dt2
        logical::Hc

        integer::nRV, nTTs, nT14s
        integer::ibd, nT0

        nRV = obsData_in%obsRV%nRV
        nTTs = obsData_in%nTTs
        nT14s = obsData_in%nDurs

        Hc = .true.

        ! it is needed to define which is the alarm coordinate for the transit detection
        call lNset(longN, check_longN)

        n_body = size(mass)
        allocate (sma(n_body))
        sma = zero
        call get_semax_vec(mass(1), mass(2:n_body), period(2:n_body), sma(2:n_body))

        nb_dim = 6*n_body
        if (.not. allocated(ra0)) allocate (ra0(nb_dim), ra1(nb_dim))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(mass, sma, ecc, meanA, argp, inc, longN, ra0)
        ra1 = ra0

        if (nRV .gt. 0) then
            call init_dataRV(nRV, simRV)
            simRV%nRV = 0
            simRV%jd = obsData_in%obsRV%jd
            simRV%RVsetID = obsData_in%obsRV%RVsetID
        end if

        if (nTTs .gt. 0) then
            allocate (simT0(n_body - 1))
            do ibd = 1, n_body - 1
                nT0 = obsData_in%obsT0(ibd)%nT0
                call init_dataT0(nT0, simT0(ibd), dur_check)
                simT0(ibd)%nT0 = 0
                simT0(ibd)%nDur = 0
                !  let's set error on TT and Dur to 1 sec! It is needed to avoid division by zero
                simT0(ibd)%eT0 = one/s24h ! 1s in day
                simT0(ibd)%edur = one/60.0_dp ! 1s in min
            end do
        end if

        dt1 = t_start - t_epoch
        dt2 = dt1 + t_int
        if (dt1 .lt. zero) then
            call ode_full_args(mass, radius, ra1, t_epoch, dt1, step_in, check_longN,&
              &simRV,&
              &id_transit_body, transit_flag, dur_check, simT0,&
              &Hc)
            if (Hc) then
                if (abs(dt1) .le. t_int) then
                    dt2 = dt1 + t_int
                    call ode_full_args(mass, radius, ra1, t_epoch, dt2, step_in, check_longN,&
                      &simRV,&
                      &id_transit_body, transit_flag, dur_check, simT0,&
                      &Hc)
                end if
            end if
        else
            dt2 = dt1 + t_int
            call ode_full_args(mass, radius, ra1, t_epoch, dt2, step_in, check_longN,&
              &simRV,&
              &id_transit_body, transit_flag, dur_check, simT0,&
              &Hc)
        end if

        ! unstable or some bad parameters case!
        if (.not. Hc) then
            if (nRV .gt. 0) then
                simRV%RV = zero
                simRV%RV_stat = 0
            end if
            if (nTTs .gt. 0) then
                do ibd = 1, n_body - 1
                    if (simT0(ibd)%nT0 .gt. 0) then
                        simT0(ibd)%T0 = zero
                        simT0(ibd)%T0_stat = 0
                        if (dur_check .eq. 1) then
                            simT0(ibd)%dur = zero
                            simT0(ibd)%dur_stat = 0
                        end if
                    end if
                end do
            end if
        end if

        if (allocated(ra0)) deallocate (ra0, ra1)


        return
    end subroutine integration_info_to_observables

! =============================================================================
! =============================================================================
! subroutine that integrates orbits of the planets and computes
! the observables (RV, T0, T14).

! **Input**
! integration information:
! start time (t_start), reference time (t_epoch), initial stepsize (step_in),
! integration time in days (t_int)
! the keplerian elements of each body (star=1) at given time (tepoch):
! mass, radius, period, ecc, argp, meanA, inc, longN
! RV information:
! time of RV observations (tRV)
! T0 information:
! number id of the transiting body (1 for all), transiting flag (should transit or not?),
! duration check flag (check the transit duration T14 or not?)

! **Output**
! simulated observables:
! simulated RV (RV_sim), simulate transit times (T0_sim, in array with rows
! given by max number of T0 among all planets, and columns the number of planets)
! =============================================================================
! =============================================================================
    subroutine integration_info_to_data(t_start, t_epoch, step_in, t_int,&
      &mass, radius, period, ecc, argp, meanA, inc, longN,&
      &time_RV, RV_sim,&
      &id_transit_body, transit_flag, dur_check,&
      &T0_sim, T14_sim, kep_elem_sim)

        real(dp), intent(in)::t_start, t_epoch, step_in, t_int
        real(dp), dimension(:), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN

        real(dp), dimension(:), intent(in)::time_RV
        real(dp), dimension(:), allocatable, intent(out)::RV_sim

        integer, intent(in)::id_transit_body
        logical, dimension(:), intent(in)::transit_flag
        integer, intent(in)::dur_check

        real(dp), dimension(:, :), allocatable, intent(out)::T0_sim, T14_sim
        real(dp), dimension(:, :, :), allocatable, intent(out)::kep_elem_sim

        type(dataRV)::simRV
        type(dataT0), dimension(:), allocatable::simT0

        integer::n_body
        integer::nRV, nTTs
        integer::ibd, nT0
        integer::max_nT0
        integer,parameter::n_kepelem=8

        n_body = size(mass)

        nRV = size(time_RV)
        if (nRV .gt. 0) then
            allocate (RV_sim(nRV))
            RV_sim = zero
            ! needed to allocate obsData%obsRV here?
            if (.not. allocated(obsData%obsRV%jd)) then
                call init_dataRV(nRV, obsData%obsRV)
                obsData%obsRV%jd = time_rv
            end if
        end if

        call integration_info_to_observables(t_start, t_epoch, step_in, t_int,&
          &mass, radius, period, ecc, argp, meanA, inc, longN,&
          &obsData, id_transit_body, transit_flag, dur_check,&
          &simRV, simT0)

        if (nRV .gt. 0) then
            ! allocate (RV_sim(nRV))
            RV_sim = simRV%RV
        end if

        nTTs = obsData%nTTs
        if (nTTs .gt. 0) then
            max_nT0 = maxval(obsData%obsT0(:)%nT0)
            allocate (T0_sim(max_nT0, n_body), T14_sim(max_nT0, n_body))
            allocate (kep_elem_sim(max_nT0, n_body, n_kepelem))
            T0_sim = zero
            T14_sim = zero
            kep_elem_sim = zero
            do ibd = 1, n_body - 1
                nT0 = simT0(ibd)%nT0
                if (nT0 .gt. 0) then
                    T0_sim(1:nT0, ibd + 1) = simT0(ibd)%T0
                    T14_sim(1:nT0, ibd + 1) = simT0(ibd)%dur
                    kep_elem_sim(1:nT0, ibd + 1, 1) = simT0(ibd)%period
                    kep_elem_sim(1:nT0, ibd + 1, 2) = simT0(ibd)%sma
                    kep_elem_sim(1:nT0, ibd + 1, 3) = simT0(ibd)%ecc
                    kep_elem_sim(1:nT0, ibd + 1, 4) = simT0(ibd)%inc
                    kep_elem_sim(1:nT0, ibd + 1, 5) = simT0(ibd)%meana
                    kep_elem_sim(1:nT0, ibd + 1, 6) = simT0(ibd)%argp
                    kep_elem_sim(1:nT0, ibd + 1, 7) = simT0(ibd)%truea
                    kep_elem_sim(1:nT0, ibd + 1, 8) = simT0(ibd)%longn
                end if
            end do
        end if

        ! if(allocated(ra0)) deallocate(ra0,ra1)

        call deallocate_dataRV(simRV)
        if (nTTs .gt. 0) then
            do ibd = 1, n_body - 1
                call deallocate_dataT0(simT0(ibd))
            end do
        end if

        return
    end subroutine integration_info_to_data

    ! ================================================================================
    ! subroutine to write file and to screen the data ... what it writes
    ! depends on the option isim and wrtid/lmon in arg.in file
    subroutine ode_output(cpuid, isim, wrtid, allpar, par, resw, fit_scale, gls_scale)
        integer, intent(in)::cpuid, isim, wrtid
        real(dp), dimension(:), intent(in)::allpar, par
        real(dp), dimension(:), intent(out)::resw
        real(dp), optional, intent(out)::fit_scale, gls_scale

        real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN
        integer, dimension(:), allocatable::check_longN
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp)::dt1, dt2
        logical::Hc

        type(dataRV)::simRV
        type(dataT0), dimension(:), allocatable::simT0

        logical::checkpar, gls_check

        real(dp)::chi_square, reduced_chi_square, lnLikelihood, ln_const, bic

        integer::ndata, nRV, nTTs, nDurs, ibd, nT0
        integer::ns, ne

        ! units and file names to store and write to files
        integer::uorb, ucon
        character(512)::florb, fmorb, flcon, fmcon, fmele
        integer, dimension(:), allocatable::uele, utra
        character(512), dimension(:), allocatable::flele, fltra

        integer::ifit

        write (*, '(a)') ''
        write (*, '(a)') " EXECUTING SIMPLE INTEGRATION AND WRITING FINAL FILES"
        write (*, '(a)') ''

        ndata = obsData%ndata
        nRV = obsData%obsRV%nRV
        ! nRVset=obsDAta%obsRV%nRVset
        nTTs = obsData%nTTs
        nDurs = obsData%nDurs

        resw = zero
        Hc = .true.
        checkpar = .true.

        write (*, '(a)') ' fitting parameters'
        do ifit = 1, nfit
            write (*, '(a12,1x,es23.16)') trim(parid(ifit)), par(ifit)
        end do

        allocate (mass(NB), radius(NB), period(NB), sma(NB), ecc(NB))
        allocate (argp(NB), meanA(NB), inc(NB), longN(NB), check_longN(NB))

        ! write (*, *) " DEBUG: convert_parameters "
        if (present(fit_scale)) then
            fit_scale = one
            gls_scale = one
            call convert_parameters_scale(allpar, par,&
              &mass, radius, period, sma, ecc, argp, meanA, inc, longN, fit_scale)
        else
            ! call convert_parameters(allpar, par,&
            !   &mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
            call only_convert_parameters(allpar, par,&
                &mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
            if (.not. checkpar) then
                write (*, '(a)') ' checkpar after convert_parameters is FALSE ... BAD'
                if (ndata > 0) resw = set_max_residuals(ndata)
            end if
        end if
        ! write (*, *) " DEBUG: checkpar ", checkpar
        ! write (*, '(1X,a,100(1x,F12.6))') " DEBUG: mass   ", mass
        ! write (*, '(1X,a,100(1x,F12.6))') " DEBUG: radius ", radius
        ! write (*, '(1X,a,100(1x,F12.6))') " DEBUG: period ", period
        ! write (*, '(1X,a,100(1x,F12.6))') " DEBUG: sma    ", sma
        ! write (*, '(1X,a,100(1x,F12.6))') " DEBUG: ecc    ", ecc
        ! write (*, '(1X,a,100(1x,F12.6))') " DEBUG: argp   ", argp
        ! write (*, '(1X,a,100(1x,F12.6))') " DEBUG: meanA  ", meanA
        ! write (*, '(1X,a,100(1x,F12.6))') " DEBUG: inc    ", inc
        ! write (*, '(1X,a,100(1x,F12.6))') " DEBUG: longN  ", longN

        if (present(fit_scale)) write (*, '(a,es23.16)') ' fit_scale = ', fit_scale
        flush (6)

        if (nRV .gt. 0) then
            call init_dataRV(nRV, simRV)
            simRV%nRV = 0
            if (.not. allocated(simRV%jitter)) allocate (simRV%jitter(nRVset))
            simRV%jitter = zero
            if (.not. allocated(simRV%gamma)) allocate (simRV%gamma(nRVset))
            simRV%gamma = zero
            simRV%RVsetID = obsData%obsRV%RVsetID
        end if

        if ((idtra .ge. 1) .and. (idtra .le. NB)) call set_file_tra(cpuid, isim, wrtid, utra, fltra)

        allocate (simT0(NB - 1))
        if (nTTs .gt. 0) then
            do ibd = 1, NB - 1
                nT0 = obsData%obsT0(ibd)%nT0
                call init_dataT0(nT0, simT0(ibd), durcheck)
                simT0(ibd)%nT0 = 0
                simT0(ibd)%nDur = 0
                !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
                simT0(ibd)%eT0 = one/s24h ! 1s in day
                simT0(ibd)%edur = one/60.0_dp ! 1s in min
            end do
        end if

        fmorb = trim(adjustl(fmtorbit()))
        if (wrtorb .eq. 1) call set_file_orb(cpuid, isim, wrtid, uorb, florb)
        fmcon = trim(adjustl(fmtconst()))
        if (wrtconst .eq. 1) call set_file_con(cpuid, isim, wrtid, ucon, flcon)
        fmele = fmtele()
        if (wrtel .eq. 1) call set_file_elem(cpuid, isim, wrtid, uele, flele)

        ! write orbital elements into a file
        call outElements(isim, wrtid, mass, radius, period, sma, ecc, argp, meanA, inc, longN)

        ! it is needed to define the which is the alarm coordinate for the transit detection
        call lNset(longN, check_longN)

        NBDIM = 6*NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(mass, sma, ecc, meanA, argp, inc, longN, ra0)
        ra1 = ra0

        dt1 = tstart - tepoch
        dt2 = dt1 + tint
        if (dt1 .lt. zero) then

            call ode_forward_output(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
              &mass, radius, ra1, dt1, check_longN, simRV, simT0, Hc)

            if (Hc) then
                if (abs(dt1) .le. tint) then
                    call ode_forward_output(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
                      &mass, radius, ra1, dt2, check_longN, simRV, simT0, Hc)
                end if
            end if

        else

            call ode_forward_output(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
              &mass, radius, ra1, dt2, check_longN, simRV, simT0, Hc)

        end if

        write (*, *)
        if (.not. Hc) then
            write (*, '(a,l2,a)') ' flag Hc = ', Hc, ' : BAD ==> PROBLEM DURING INTEGRATION'
        else
            write (*, '(a,l2,a)') ' flag Hc = ', Hc, ' : OK  ==> NO PROBLEM DURING INTEGRATION'
        end if
        write (*, *)
        flush (6)

        ! write (*, *) " DEBUG: close_tra"
        if ((idtra .ge. 1) .and. (idtra .le. NB)) call close_tra(utra, fltra)
        ! write (*, *) " DEBUG: close(uorb)"
        if (wrtorb .eq. 1) close (uorb)
        ! write (*, *) " DEBUG: close(ucon)"
        if (wrtconst .eq. 1) close (ucon)
        ! write (*, *) " DEBUG: close_elem"
        if (wrtel .eq. 1) call close_elem(uele, flele)

        if (ndata .gt. 0) then

            ! write (*, *) " DEBUG: ndata > 0"
            chi_square = -one
            reduced_chi_square = -one
            lnLikelihood = half
            ln_const = zero
            bic = zero

            ! write (*, *) " DEBUG: Hc ", Hc
            if (.not. Hc) then
                resw = set_max_residuals(ndata)

            else

                ! write (*, *) " DEBUG: Hc true"
                ! write (*, *) " DEBUG: nRVset ", nRVset
                if (nRVset .gt. 0) then
                    ! set jitter into simRV
                    ns = nkel + 1
                    ne = nkel + nRVset
                    simRV%jitter = two**par(ns:ne)
                    ! set gamma into simRV
                    ns = ne + 1
                    ne = ne + nRVset
                    simRV%gamma = par(ns:ne)
                    ! write (*, *) " DEBUG: set_gamma_rv"
                    call set_gamma_rv(simRV%gamma, simRV%RVsetID, simRV%gamma_rv)
                    ! compute RV trend
                    if (rv_trend_order .gt. 0) then
                        ! write (*, *) " DEBUG: addRVtrend"
                        ns = ne + 1
                        ne = ne + rv_trend_order
                        call addRVtrend(simRV%jd, tepoch, par(ns:ne), simRV%trend)
                    end if
                end if

                ! write (*, *) " DEBUG: sum(resw*resw) = ", sum(resw*resw)
                ! write (*, *) " DEBUG: RV set_weighted_residuals"
                call set_weighted_residuals(obsData, simRV, simT0, resw, oc_fit)
                ! write (*, *) " DEBUG: sum(resw*resw) = ", sum(resw*resw)

                if (present(fit_scale)) then
                    if (simRV%nRV .gt. 0) call check_periodogram_scale(obsData%obsRV%jd,&
                      &obsData%obsRV%RV - simRV%RV - simRV%gamma_rv - simRV%trend,&
                      &obsData%obsRV%eRV,&
                      &period, gls_check, gls_scale)
                    write (*, *)
                    write (*, '(a)') ' CHECK GLS RESIDUALS-RV'
                    write (*, '(a,es23.16)') ' gls_scale = ', gls_scale
                    write (*, *)
                    flush (6)
                    resw = resw*fit_scale*gls_scale
                else
                    ! write (*, *) " DEBUG: check_periodogram"
                    if (simRV%nRV .gt. 0) call check_periodogram(obsData%obsRV%jd,&
                      &obsData%obsRV%RV - simRV%RV - simRV%gamma_rv - simRV%trend,&
                      &obsData%obsRV%eRV,&
                      &period, gls_check)
                end if
                ! write (*, *) " DEBUG: gls_check ", gls_check

                ! write (*, *) " DEBUG: check_max_residuals"
                call check_max_residuals(resw, ndata)

            end if

            if (simRV%nRV .gt. 0) then

                write (*, *) ""
                write (*, '(a,i4)') " RADIAL VELOCITIES found: ", simRV%nRV
                write (*, *) ""

                ! if (present(to_screen)) call write_RV(simRV)
                call write_RV(cpuid, isim, wrtid, simRV)

                ! write (*,*) ' DEBUG: check_and_write_periodogram'
                ! 2016-04-08: added gls check
                if (present(gls_scale)) then
                    call check_and_write_periodogram(cpuid, isim, wrtid,&
                      &obsData%obsRV%jd,&
                      &obsData%obsRV%RV - simRV%RV - simRV%gamma_rv - simRV%trend,&
                      &obsData%obsRV%eRV, period, gls_check, gls_scale)
                else
                    call check_and_write_periodogram(cpuid, isim, wrtid,&
                      &obsData%obsRV%jd,&
                      &obsData%obsRV%RV - simRV%RV - simRV%gamma_rv - simRV%trend,&
                      &obsData%obsRV%eRV, period, gls_check)
                end if
                ! if (.not. gls_check) resw = set_max_residuals(ndata)

            else

                write (*, *)
                write (*, '(a)') ' RADIAL VELOCITIES NOT FOUND'
                write (*, *)

            end if

            if (allocated(simT0) .and. sum(simT0(:)%nT0) .gt. 0) then
                write (*, *)
                write (*, '(a,i5)') " T0 SIM found ", sum(simT0(:)%nT0)
                write (*, *)
                ! if (present(to_screen)) call write_T0(simT0)
                call write_T0(cpuid, isim, wrtid, simT0)

            else
                write (*, *)
                write (*, '(a)') ' TRANSIT TIMES NOT FOUND'
                write (*, *)
            end if

            ! Reduced Chi Squares etc for report/summary
            ! write (*, *) " DEBUG: set_fitness_values"
            chi_square = sum(resw*resw)
            ! write(*,*)" DEBUG: chi_square         = ",chi_square
            call set_fitness_values(par, chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)
            ! write(*,*)" DEBUG: chi_square         = ",chi_square

        else ! ndata <= 0

            write (*, '(a)') ' ndata == 0: NO FITNESS SUMMARY (SCREEN AND FILES)'
            flush (6)

        end if  ! ndata

        ! write (*, *) ' DEBUG: deallocating simRV'
        call deallocate_dataRV(simRV)
        ! write (*, *) ' DEBUG: deallocating simT0'
        if (nTTs .gt. 0) then
            do ibd = 1, NB - 1
                call deallocate_dataT0(simT0(ibd))
            end do
        end if

        ! write (*, *) ' DEBUG: ode_out deallocated simRV and simT0...now deallocating kep/ra'

        if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN, check_longN)
        if (allocated(ra0)) deallocate (ra0, ra1)

        ! write (*, *) " DEBUG: LAST ode_output: sum(resw*resw) = ", sum(resw*resw)

        return
    end subroutine ode_output

    ! ================================================================================

    subroutine ode_integrates(cpuid, isim, wrtid, m, R, P, sma, ecc, w, mA, inc, lN)
        integer, intent(in)::cpuid, isim, wrtid
        real(dp), dimension(:), intent(in)::m, R, P, sma, ecc, w, mA, inc, lN

        integer, dimension(:), allocatable::check_longN
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp)::dt1, dt2
        logical::Hc

        ! units and file names to store and write to files
        integer::uorb, ucon
        character(512)::florb, fmorb, flcon, fmcon, fmele
        integer, dimension(:), allocatable::uele, utra
        character(512), dimension(:), allocatable::flele, fltra

        write (*, '(a)') ''
        write (*, '(a)') " EXECUTING SIMPLE INTEGRATION AND WRITING FINAL FILES"
        write (*, '(a)') ''

!     resw=zero
        Hc = .true.

        ! write orbital elements into a file
        call outElements(isim, wrtid, m, R, P, sma, ecc, w, mA, inc, lN)

        ! it is needed to define the which is the alarm coordinate for the transit detection
        call lNset(lN, check_longN)

        NBDIM = 6*NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(m, sma, ecc, mA, w, inc, lN, ra0)

        ra1 = ra0

        if ((idtra .ge. 1) .and. (idtra .le. NB)) call set_file_tra(cpuid, isim, wrtid, utra, fltra)

        fmorb = trim(adjustl(fmtorbit()))
        if (wrtorb .eq. 1) call set_file_orb(cpuid, isim, wrtid, uorb, florb)
        fmcon = trim(adjustl(fmtconst()))
        if (wrtconst .eq. 1) call set_file_con(cpuid, isim, wrtid, ucon, flcon)
        fmele = fmtele()
        if (wrtel .eq. 1) call set_file_elem(cpuid, isim, wrtid, uele, flele)

        write (*, '(a)') ' INITIALISED ORBIT AND OUTPUT FILES'
        write (*, '(a)') ' RUNNING INTEGRATION ...'
        flush (6)

        dt1 = tstart - tepoch
        dt2 = dt1 + tint
        if (dt1 .lt. zero) then
            call ode_forward_output_nodata(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
                &m, R, ra1, dt1, check_longN, Hc)
            if (abs(dt1) .le. tint) then
                call ode_forward_output_nodata(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
                    &m, R, ra1, dt2, check_longN, Hc)
            end if
        else
            call ode_forward_output_nodata(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
                &m, R, ra1, dt2, check_longN, Hc)
        end if

        write (*, '(a)') ' COMPLETED'
        if (.not. Hc) then
            write (*, '(a)') 'WARNING'
            write (*, '(a,l2,a)') ' flag Hc = ', Hc, ' : BAD ==> PROBLEM DURING INTEGRATION'
        else
            write (*, '(a,l2,a)') ' flag Hc = ', Hc, ' : OK  ==> NO PROBLEM DURING INTEGRATION'
        end if
        write (*, *)
        flush (6)

        if ((idtra .ge. 1) .and. (idtra .le. NB)) call close_tra(utra, fltra)
        if (wrtorb .eq. 1) close (uorb)
        if (wrtconst .eq. 1) close (ucon)
        if (wrtel .eq. 1) call close_elem(uele, flele)

        if (allocated(ra0)) deallocate (ra0, ra1)

        return
    end subroutine ode_integrates
    ! ================================================================================

! ==============================================================================
! ==============================================================================
! ORBITAL PARAMETERS TO ALL TRANSIT TIMES OF ALL PLANETS AND RV MODEL
! NO T0 FIT/DATA, NO RV FIT/DATA, NO OUTPUT FILES
! ==============================================================================
! ==============================================================================

    ! ================================================================================
    subroutine ode_b_orbit(m, R, rin, time_int, wrt_time, check_longN,&
      &last_tra, ttra_full, dur_full, id_ttra_full, stats_ttra,&
      &last_rv, time_rv_nmax, rv_nmax, stats_rv,&
      &Hc)
        real(dp), dimension(:), intent(in)::m, R, rin
        real(dp), intent(in)::time_int, wrt_time
        integer, dimension(:), intent(in)::check_longN

        ! transit times variables to be updated!!
        integer, intent(inout)::last_tra ! if first call it is zero, otherwise it is the last transit position
        real(dp), dimension(:), intent(inout)::ttra_full, dur_full
        integer, dimension(:), intent(inout)::id_ttra_full
        logical, dimension(:), intent(inout)::stats_ttra

        ! radial velocities variables to be updated!!
        integer, intent(inout)::last_rv
        real(dp), dimension(:), intent(inout)::time_rv_nmax, rv_nmax
        logical, dimension(:), intent(inout)::stats_rv

        logical, intent(out)::Hc

        integer::ntra_full

        real(dp)::step_rv, rv_temp, ttra_temp, dur_tra_temp
        logical::check_ttra
        integer::nrv_max

        real(dp), dimension(:), allocatable::dr, r1, r2, err
        integer, dimension(:), allocatable::X, Y, Z
        real(dp), dimension(:), allocatable::cX, cY, cR, rmean
        real(dp)::hw, hok, hnext, iter_time, iter_write, step_write
        integer::j, step_num

        ntra_full = size(ttra_full)
        nrv_max = size(rv_nmax)

        Hc = .true.
!     if(do_hill_check) Hc=mutual_Hill_check(m,rin)
        Hc = separation_mutual_Hill_check(m, R, rin, do_hill_check)
        if (.not. Hc) return

        ! init the state vector
        allocate (X(NB), Y(NB), Z(NB), cX(NB), cY(NB), cR(NB), rmean(NB))
        X = 0
        Y = 0
        Z = 0
        do j = 2, NB
            X(j) = 1 + (j - 1)*6
            Y(j) = 2 + (j - 1)*6
            Z(j) = 3 + (j - 1)*6
        end do
        cX = one
        cY = one
        cR = zero
        cR(2:NB) = 1.5_dp*(R(1) + R(2:NB))*RsunAU
        rmean = 9.0e9_dp

        allocate (dr(NBDIM), r1(NBDIM), r2(NBDIM), err(NBDIM))

        ! set initial stepsize to the value in the arg.in file
        hw = step_0
        if (time_int .lt. zero) hw = -hw ! reverse if backward integration

        ! set the iteration time: it will be updated with each step
        iter_time = zero

        ! set the initial state vector
        r1 = rin
        r2 = zero
        err = zero

        step_num = 0 ! step counter

        ! set the iteration write time: updated each wrt_time passed
        step_write = wrt_time
        if (time_int .lt. zero) step_write = -step_write
        iter_write = step_write

        integration: do
            step_num = step_num + 1

            if (abs(iter_time + hw) .gt. abs(time_int)) hw = time_int - iter_time ! if last step exceeds integration time create new stepsize
            call eqmastro(m, r1, dr) ! computes the eq. of motion
            call int_rk_a(m, r1, dr, hw, hok, hnext, r2, err) ! computes the next orbit step

            Hc = separation_mutual_Hill_check(m, R, r2, do_hill_check)
            if (.not. Hc) return

            ! check if it passes the iter_write and compute the rv at proper time and update last_rv
            if (abs(iter_time + hok) .ge. (abs(iter_write))) then

                rvloop: do

                    last_rv = last_rv + 1
                    if (last_rv .gt. nrv_max) then
                        Hc = .false.
                        return
                    end if
!           computes the proper step
                    step_rv = iter_write - iter_time
                    call calcRV(m, r1, dr, step_rv, rv_temp) ! computes the rv
!           save time rv status
                    time_rv_nmax(last_rv) = iter_write
                    rv_nmax(last_rv) = rv_temp
                    stats_rv(last_rv) = .true.
!           update next writing time
                    iter_write = iter_write + step_write
                    if (abs(iter_write) .gt. abs(iter_time + hok)) exit rvloop

                end do rvloop

            end if

            ! transit times!!
            ! loop on planets -> checks for all planets
            do j = 2, NB
                ! transit time check variables
                cX(j) = r1(X(j))*r2(X(j))
                cY(j) = r1(Y(j))*r2(Y(j))

                ! check the transit criteria for each body

                check_ttra = .false. ! initialise to .false.
                ttra_temp = -9.0e10_dp     !               zero
                dur_tra_temp = -9.0e10_dp  !               zero

                if (check_longN(j) .eq. 0) then ! condition to check X

                    if ((cX(j) .le. zero) .and. (r1(Z(j)) .gt. zero)) then
                        call transit_time(j, m, R, r1, r2, iter_time, hok, ttra_temp, dur_tra_temp, check_ttra)
                        if (check_ttra) then
                            last_tra = last_tra + 1
                            if (last_tra .gt. ntra_full) then
                                Hc = .false.
                                return
                            end if
                            ttra_full(last_tra) = ttra_temp
                            dur_full(last_tra) = dur_tra_temp
                            id_ttra_full(last_tra) = j
                            stats_ttra(last_tra) = .true.
                        end if
                    end if

                else ! condition to check Y

                    if ((cY(j) .le. zero) .and. (r1(Z(j)) .gt. zero)) then
                        call transit_time(j, m, R, r1, r2, iter_time, hok, ttra_temp, dur_tra_temp, check_ttra)
                        if (check_ttra) then
                            last_tra = last_tra + 1
                            if (last_tra .gt. ntra_full) then
                                Hc = .false.
                                return
                            end if
                            ttra_full(last_tra) = ttra_temp
                            dur_full(last_tra) = dur_tra_temp
                            id_ttra_full(last_tra) = j
                            stats_ttra(last_tra) = .true.
                        end if
                    end if

                end if ! end condition X,Y

!         end if ! end transit criteria

            end do

            ! update iteration time with the stepsize used (hok)
            iter_time = iter_time + hok

            ! check end of integration: time_int reached
            if (abs(iter_time) .ge. abs(time_int)) exit integration
            ! update step and state vector
            hw = hnext
            r1 = r2

        end do integration

        deallocate (X, Y, Z, cX, cY, cR, rmean)
        deallocate (dr, r1, r2, err)

        return
    end subroutine ode_b_orbit
    ! ================================================================================

    ! ================================================================================
    subroutine ode_all_ttra_rv(wrt_time, m, R, P, a, e, w, mA, inc, lN,&
      &ttra_full, dur_full, id_ttra_full, stats_ttra,&
      &time_rv_nmax, rv_nmax, stats_rv)
        real(dp), intent(in)::wrt_time
        real(dp), dimension(:), intent(in)::m, R, P, a, e, w, mA, inc, lN

        real(dp), dimension(:), intent(out)::ttra_full, dur_full
        integer, dimension(:), intent(out)::id_ttra_full
        logical, dimension(:), intent(out)::stats_ttra
        real(dp), dimension(:), intent(out)::time_rv_nmax, rv_nmax
        logical, dimension(:), intent(out)::stats_rv

        integer, dimension(:), allocatable::check_longN
        real(dp), dimension(:), allocatable::ra0, ra1

        integer::last_tra, last_rv

        real(dp)::dt1, dt2
        logical::Hc

        Hc = .true.

        ! it is needed to define the which is the alarm coordinate for the transit detection
        call lNset(lN, check_longN)

        NBDIM = 6*NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(m, a, e, mA, w, inc, lN, ra0)
        ra1 = ra0

        ! prepare rv arrays
        time_rv_nmax = zero
        rv_nmax = zero
        stats_rv = .false.
        last_rv = 0

        ! prepare ttra arrays
        ttra_full = -9.e10_dp
        dur_full = -9.e10_dp
        id_ttra_full = 0
        stats_ttra = .false.
        last_tra = 0

        dt1 = tstart - tepoch
        dt2 = dt1 + tint

        if (dt1 .lt. zero) then

            ! backward integration
            call ode_b_orbit(m, R, ra1, dt1, wrt_time, check_longN,&
              &last_tra, ttra_full, dur_full, id_ttra_full, stats_ttra,&
              &last_rv, time_rv_nmax, rv_nmax, stats_rv, Hc)
            if (.not. Hc) then
                if (allocated(ra0)) deallocate (ra0, ra1)
                return
            end if

            if (abs(dt1) .le. tint) then
                ! forward integration
                call ode_b_orbit(m, R, ra1, dt2, wrt_time, check_longN,&
                  &last_tra, ttra_full, dur_full, id_ttra_full, stats_ttra,&
                  &last_rv, time_rv_nmax, rv_nmax, stats_rv, Hc)
                if (.not. Hc) then
                    if (allocated(ra0)) deallocate (ra0, ra1)
                    return
                end if
            end if

        else

            ! only forward integration
            call ode_b_orbit(m, R, ra1, dt2, wrt_time, check_longN,&
              &last_tra, ttra_full, dur_full, id_ttra_full, stats_ttra,&
              &last_rv, time_rv_nmax, rv_nmax, stats_rv, Hc)

        end if

        if (allocated(ra0)) deallocate (ra0, ra1)

        return
    end subroutine ode_all_ttra_rv
    ! ================================================================================

end module ode_run

