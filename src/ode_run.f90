module ode_run
    use constants
    use custom_type
    use parameters
    use parameters_conversion
    use linear_ephem
    use celestial_mechanics
    use eq_motion, only: eqmastro
    use numerical_integrator, only: int_rk_a, rkck_a, one_forward_step
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
    subroutine set_checking_coordinates(n_body, idxX, idxY, idxZ, idxVX, idxVY, idxVZ)
        ! **Input**
        integer, intent(in)::n_body
        ! **Output**
        integer, dimension(:), allocatable, intent(out)::idxX, idxY, idxZ
        integer, dimension(:), allocatable, intent(out)::idxVX, idxVY, idxVZ
        ! Local
        real(dp), parameter::Rfactor = one
        integer::i_body

        allocate (idxX(n_body), idxY(n_body), idxZ(n_body))
        allocate (idxVX(n_body), idxVY(n_body), idxVZ(n_body))

        idxX = 0
        idxY = 0
        idxZ = 0
        idxVX = 0
        idxVY = 0
        do i_body = 2, n_body
            idxX(i_body) = 1 + (i_body - 1) * 6
            idxY(i_body) = 2 + (i_body - 1) * 6
            idxZ(i_body) = 3 + (i_body - 1) * 6
            idxVX(i_body) = 4 + (i_body - 1) * 6
            idxVY(i_body) = 5 + (i_body - 1) * 6
            idxVZ(i_body) = 6 + (i_body - 1) * 6
        end do

        return
    end subroutine set_checking_coordinates

    subroutine transit_conditions(rv1, rv2, A, B, AB, ABflag, Zflag)
        ! Input
        real(dp), dimension(6), intent(in)::rv1, rv2
        ! Output
        real(dp), intent(out)::A, B, AB
        logical, intent(out)::ABflag, Zflag

        ! C = X*VX + Y*VY
        A = rv1(1) * rv1(4) + rv1(2) * rv1(5)
        B = rv2(1) * rv2(4) + rv2(2) * rv2(5)
        AB = A * B
        ABflag = AB .le. zero
        Zflag = (rv1(3) .gt. zero) .or. (rv2(3) .gt. zero)

        return
    end subroutine transit_conditions

    ! subroutine check_if_excluded_timerange(i_body, mass, radii, t_epoch, itime1, itime2, r1, r2, integration_step,&
    !     &excl_body, excl_tranges, excl_status, transit_in_excluded)
    !     ! Input
    !     integer, intent(in)::i_body
    !     real(dp), dimension(:), intent(in)::mass, radii
    !     real(dp), intent(in)::t_epoch
    !     real(dp), intent(in)::itime1, itime2
    !     real(dp), dimension(:), intent(in)::r1, r2
    !     real(dp), intent(in)::integration_step

    !     integer, dimension(:), intent(in)::excl_body
    !     real(dp), dimension(:, :), intent(in)::excl_tranges
    !     ! Input/Output
    !     integer, dimension(:), intent(inout)::excl_status
    !     ! Output
    !     logical, intent(inout)::transit_in_excluded
        
    !     ! Local
    !     integer::n_excl, i_excl
    !     logical, dimension(:), allocatable::sel_body
    !     integer::n_sel
    !     ! real(dp), dimension(:, :), allocatable::excl_tranges_body
    !     ! integer, dimension(:), allocatable::excl_status_body
    !     real(dp)::tmin, tmax, tr1, tr2
    !     logical::overlap
    !     real(dp)::ttra, dur_tra, half_dur, alpha
    !     real(dp), dimension(4)::tcont
    !     logical::check_ttra, check_t0, check_t1, check_t4

    !     tr1 = t_epoch + itime1
    !     tr2 = t_epoch + itime2
    !     n_excl = size(excl_body)

    !     allocate (sel_body(size(excl_body)))
    !     sel_body = excl_body .eq. i_body
    !     n_sel = count(sel_body)

    !     if (n_sel .gt. 0) then

    !         do i_excl = 1, n_excl
    !             if (excl_status(i_excl) .eq. 0) then
    !                 tmin = excl_tranges(i_excl, 1)
    !                 tmax = excl_tranges(i_excl, 2)
    !                 call compare_inttime_timerange(tr1, tr2, tmin, tmax, overlap)
    !                 if (overlap) then
    !                     call transit_and_contact_time(t_epoch, i_body, mass, radii, r1, r2, tr1, integration_step,&
    !                         &ttra, dur_tra, tcont, alpha, check_ttra)
    !                     ! update CHECKS t0 t1 t4 ==> transit_in_excluded
    !                     ! transit_in_excluded = &
    !                     !     & ((tcont(1) .ge. excl_tranges_body(i_excl, 1)) .and. (tcont(1) .lt. excl_tranges_body(i_excl, 2)))&
    !                     !     & .or. &
    !                     !     & ((tcont(4) .gt. excl_tranges_body(i_excl, 1)) .and. (tcont(4) .le. excl_tranges_body(i_excl, 2)))
    !                     ! if (transit_in_excluded) write(*,*)"i_body ",i_body," i_excl ", i_excl, tcont(1), tcont(4), excl_tranges_body(i_excl, :)
    !                 end if
    !             end if
    !         end do
    !         ! deallocate (excl_status_body, excl_tranges_body)
    !     end if
    !     deallocate (sel_body)

    !     return
    ! end subroutine check_if_excluded_timerange

    subroutine check_if_excluded_timerange(i_body, mass, radii, tepoch, itime1, itime2, r1, r2, integration_step,&
        &T0,transit_in_excluded)
        ! Input
        integer, intent(in)::i_body
        real(dp), dimension(:), intent(in)::mass, radii
        real(dp), intent(in)::tepoch
        real(dp), intent(in)::itime1, itime2
        real(dp), dimension(:), intent(in)::r1, r2
        real(dp), intent(in)::integration_step
        ! Input/Output
        type(dataT0), intent(inout)::T0
        logical, intent(inout)::transit_in_excluded
        
        ! Local
        integer::nexcl, n_status, i_excl
        real(dp)::tmin, tmax, tr1, tr2
        logical::overlap
        real(dp)::ttra, dur_tra, half_dur, alpha
        real(dp), dimension(4)::tcont
        logical::check_ttra, check_t0, check_t1, check_t4, check_this

        tr1 = tepoch + itime1
        tr2 = tepoch + itime2

        nexcl = T0%n_excl
        n_status = sum(T0%excluded_status)

        if ((nexcl .gt. 0).and.(n_status .lt. nexcl)) then

            do i_excl = 1, nexcl
                if (T0%excluded_status(i_excl) .eq. 0) then
                    tmin = T0%excluded_time_ranges(i_excl, 1)
                    tmax = T0%excluded_time_ranges(i_excl, 2)
                    call compare_inttime_timerange(tr1, tr2, tmin, tmax, overlap)
                    if (overlap) then
                        call transit_and_contact_time(tepoch, i_body, mass, radii, r1, r2, itime1, integration_step,&
                        &ttra, dur_tra, tcont, alpha, check_ttra)
                        check_t0 = (ttra .ge. tmin) .and. (ttra .le. tmax)
                        check_t1 = (tcont(1) .ge. tmin) .and. (tcont(1) .lt. tmax)
                        check_t4 = (tcont(4) .gt. tmin) .and. (tcont(4) .le. tmax)
                        check_this = (check_t0 .or. check_t1 .or. check_t4)
                        transit_in_excluded = (transit_in_excluded .or. check_this)
                        T0%excluded_t0 = ttra
                        T0%excluded_t1 = tcont(1)
                        T0%excluded_t4 = tcont(4)
                        write(*,*)"i_body ",i_body," i_excl ", i_excl
                        write(*,*)"t1 = ",tr1, " t2 = ",tr2," tmin = ",tmin," tmax = ",tmax," overlap = ",overlap
                        write(*,*)"T0 = ",ttra," T1 = ",tcont(1)," T4 = ",tcont(4)
                        write(*,*)"check T0 - T1 - T4 & this: ",check_t0, check_t1, check_t4, check_this
                        write(*,*)"transit_in_excluded: ",transit_in_excluded
                        flush(6)
                        if (check_this) then
                            T0%excluded_status(i_excl) = 1
                        end if
                        write(*,*)"status = ",T0%excluded_status(i_excl)
                        flush(6)
                    end if
                end if
            end do
            ! deallocate (excl_status_body, excl_tranges_body)
        end if
        ! deallocate (sel_body)

        return
    end subroutine check_if_excluded_timerange

    ! ================================================================================
    ! performs the forward integration in one direction with explicit flags and args.
    ! computes the orbits and RVs and T0s (T14s) as in
    ! obsRV.dat and NB#_observations.dat files (stored in obsDataIn).
    !
    ! **Input**
    ! mass          == masses of all bodies (star == 1)
    ! radius        == radius of all bodies
    ! rin           == stave vector of all bodies in astrocentric coordinates XYZVxVyVz
    ! t_epoch       == reference time epoch of orbital parameters
    ! time_to_int   == integration time
    ! step_in       == integration step
    ! obsDataIn     == input data as dataObs type
    !
    ! id_transit_body == id of the transiting body/ies
    ! transit_flag == transit or not needed to check_T0 subroutine
    ! dur_check == flag to check or not the transit duration
    !
    ! **Input/Output**
    ! simRV  == dataRV type with simulated RV
    ! simT0  == array of dataT0 type with simulated T0 (T14), for each planet
    !
    ! **Output**
    ! Hc     == Hill chek, .true. if stable, .false. if unstable for Hill's condition
    subroutine ode_full_args(mass, radius, rin, t_epoch, time_to_int, step_in, obsDataIn, simRV,&
        &id_transit_body, transit_flag, dur_check, simT0, Hc)
        ! Input
        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::t_epoch, time_to_int, step_in
        type(dataObs), intent(in)::obsDataIn
        ! Input/Output
        type(dataRV), intent(inout)::simRV
        ! Input
        integer, intent(in)::id_transit_body
        logical, dimension(:), intent(in)::transit_flag
        integer, intent(in)::dur_check
        ! Input/Output
        type(dataT0), dimension(:), intent(inout)::simT0
        ! integer, dimension(:), intent(inout)::excl_status
        logical, intent(inout)::Hc
        ! Locals
        integer::n_body, nb_dim
        real(dp), dimension(:), allocatable::r1, r2
        integer, dimension(:), allocatable::X, Y, Z
        integer, dimension(:), allocatable::VX, VY, VZ

        real(dp)::A, B, AB
        logical::ABflag, Zflag

        real(dp)::integration_step, working_step

        real(dp), dimension(:), allocatable::t_all_steps
        integer, dimension(:), allocatable::all_idx
        integer::n_all
        real(dp)::trun1, trun2

        real(dp), dimension(:), allocatable::period, sma, ecc, inc, meanA, argp, trueA, longN, dttau
        real(dp), dimension(:), allocatable::Tref, Pref
        integer, dimension(:), allocatable::epo_obs
        real(dp)::Tr, Pr, Tx

        integer::i_body, iteration, body_transiting_start, body_transiting_end
        logical::do_transit_check, search_obs_transit, transit_in_excluded

        integer::nRV, nTTs, nexcl_all

        nRV = obsDataIn % obsRV % nRV
        nTTs = obsDataIn % nTTs
        transit_in_excluded = .false.
        nexcl_all = obsDataIn%n_excluded

        Hc = .true.
        if (close_encounter_check) then
            Hc = separation_mutual_Hill_check(mass, radius, rin, do_hill_check)
        end if
        if (amd_hill_check) call statevector_amd_hill_stability(mass, rin, Hc)
        if (.not. Hc) then
            return
        end if

        n_body = size(mass)
        call set_checking_coordinates(NB, X, Y, Z, VX, VY, VZ)

        integration_step = sign(step_in, time_to_int)
        call set_check_steps(t_epoch, time_to_int, integration_step, obsDataIn, n_all, t_all_steps, all_idx)

        nb_dim = 6 * n_body
        allocate (r1(nb_dim), r2(nb_dim))
        r1 = rin
        working_step = step_0

        call set_transiting_bodies(id_transit_body, do_transit_check, body_transiting_start, body_transiting_end)
        ! Define the reference transit time for each body
        call set_transit_references_all_bodies(t_epoch, obsDataIn, Tref)
        allocate (period(n_body))
        allocate (sma(n_body), ecc(n_body), inc(n_body), meanA(n_body), argp(n_body), trueA(n_body), longN(n_body), dttau(n_body))
        call elements(mass, rin, period, sma, ecc, inc, meanA, argp, trueA, longN, dttau)
        allocate (Pref(n_body - 1))
        Pref = period(2:n_body)

        trun1 = t_all_steps(1)
        integration: do iteration = 1, n_all

            trun2 = t_all_steps(iteration)
            integration_step = trun2 - trun1
            call one_forward_step(mass, radius, r1, integration_step, working_step, Hc, r2)
            if (.not. Hc) then
                exit integration
            end if

            ! ! ===============================
            ! RV computation
            if ((nRV .gt. 0) .and. (all_idx(iteration) .gt. 0)) then
                call get_and_set_RV_data(mass, r2, all_idx(iteration), obsDataIn, simRV)
            end if

            if (do_transit_check) then
                Tx = t_epoch + trun1 + half * integration_step

                do i_body = body_transiting_start, body_transiting_end
                    call transit_conditions(r1(X(i_body):VZ(i_body)), r2(X(i_body):VZ(i_body)),&
                    &A, B, AB, ABflag, Zflag)

                    if (ABflag .and. Zflag) then

                        call check_if_observed_transit(i_body, Tx, Tref, Pref, obsDataIn, search_obs_transit)
                        if (search_obs_transit) then
                            Tr = Tref(i_body - 1)
                            Pr = Pref(i_body - 1)
                            call check_T0(t_epoch, Tr, Pr, i_body, mass, radius, r1, r2,&
                                &trun1, integration_step,&
                                &transit_flag, dur_check, obsDataIn, simT0, Hc)
                        end if

                        if (nexcl_all .gt. 0) then
                            call check_if_excluded_timerange(i_body, mass, radius,&
                                &t_epoch, trun1, trun2, r1, r2, integration_step,&
                                &simT0(i_body-1),transit_in_excluded)
                            if (transit_in_excluded) then
                                ! Hc = .false.
                                ! exit integration
                                write(*,*)" in ode_full_args: transit_in_excluded is ", transit_in_excluded, " and Hc is ", Hc
                                flush(6)
                            end if
                        end if

                    end if ! ABflag, Zflag
                end do ! i_body
            end if ! do_transit_check
            ! ========

            r1 = r2
            trun1 = trun2

            if (.not. Hc) then
                exit integration
            end if
        end do integration

        deallocate (r1, r2)
        deallocate (X, Y, Z)
        deallocate (VX, VY, VZ)
        deallocate (period, sma, ecc, inc, meanA, argp, trueA, longN, dttau)
        deallocate (Tref, Pref)
        deallocate (t_all_steps)
        deallocate (all_idx)

        return
    end subroutine ode_full_args

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
    ! step   == step of check
    !
    ! **Input/Output**
    ! simRV  == dataRV type with simulated RV
    ! simT0  == array of dataT0 type with simulated T0 (T14), for each planet
    !
    ! **Output**
    ! Hc     == Hill chek, .true. if stable, .false. if unstable for Hill's condition
    subroutine ode_forward_data(mass, radius, rin, time_to_int, step, obsDataIn, simRV, simT0, Hc)
        ! **Input**
        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::time_to_int, step
        type(dataObs), intent(in)::obsDataIn
        ! **Input/Output**
        type(dataRV), intent(inout)::simRV
        type(dataT0), dimension(:), intent(inout)::simT0
        ! Output
        logical, intent(out)::Hc
        !

        Hc = .true.
        call ode_full_args(mass, radius, rin, tepoch, time_to_int, step, obsDataIn, simRV,&
            &idtra, do_transit, durcheck, simT0, Hc)

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
    ! time_to_int         == integration time
    ! obsDataIn           == RV and T0s/Durs dataset of type dataObs
    !
    ! **Input/Output**
    ! simRV               == dataRV type with simulated RV
    ! simT0               == array of dataT0 type with simulated T0 (T14), for each planet
    ! Hc                  == Hill chek, .true. if stable, .false. if unstable for Hill's condition
    subroutine ode_forward_output_data(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
      &mass, radius, rin, time_to_int, obsDataIn, simRV, simT0, Hc)
        ! **Input**
        integer, intent(in)::uorb, ucon
        integer, dimension(:), intent(in)::uele, utra
        character(*), intent(in)::fmorb, fmcon, fmele
        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::time_to_int
        type(dataObs), intent(in)::obsDataIn
        ! **Input/Output**
        type(dataRV), intent(inout)::simRV
        type(dataT0), dimension(:), intent(inout)::simT0
        logical, intent(inout)::Hc
        ! Local variables
        real(dp), dimension(:), allocatable::r1, r2

        integer, dimension(:), allocatable::X, Y, Z
        integer, dimension(:), allocatable::VX, VY, VZ

        real(dp)::A, B, AB
        logical::ABflag, Zflag

        integer::n_body, i_body, body_transiting_start, body_transiting_end
        integer::nRV, nTTs, nDurs

        real(dp)::integration_step, working_step !, ok_step, next_step, itime

        integer::iteration, save_iteration, jtra
        real(dp), dimension(:, :), allocatable::storeorb, storecon, storetra
        integer, dimension(:, :), allocatable::stat_tra
        integer::Norb
        real(dp)::Etot, Eold, htot, hold

        real(dp)::step_write
        logical::do_transit_check

        real(dp), dimension(:), allocatable::t_all_steps
        integer, dimension(:), allocatable::all_idx
        integer::n_all
        real(dp)::trun1, trun2

        real(dp), dimension(:), allocatable::period, sma, ecc, inc, meanA, argp, trueA, longN, dttau
        real(dp), dimension(:), allocatable::Tref, Pref
        ! integer, dimension(:), allocatable::epo_obs
        real(dp)::Tr, Pr, Tx
        ! integer::epox, ntx
        logical::search_obs_transit, transit_in_excluded

        ! if you see a variable non listed here, probably it is a global variable
        ! defined in constants.f90 or parameters.f90

        flush (6)

        nRV = obsDataIn % obsRV % nRV
        nTTs = obsDataIn % nTTs
        nDurs = obsDataIn % nDurs

        Hc = .true.
        if (close_encounter_check) then
            Hc = separation_mutual_Hill_check(mass, radius, rin, do_hill_check)
        end if
        if (amd_hill_check) call statevector_amd_hill_stability(mass, rin, Hc)
        if (.not. Hc) then
            return
        end if

        call set_checking_coordinates(NB, X, Y, Z, VX, VY, VZ)

        step_write = sign(wrttime, time_to_int)
        call set_check_steps(tepoch, time_to_int, step_write, obsDataIn, n_all, t_all_steps, all_idx)

        Norb = NBDIM + 3

        allocate (r1(NBDIM), r2(NBDIM))
        r1 = rin
        working_step = step_0

        ! ==================
        ! PREPARE VARIABLE TO STORE ORBIT/ENERGY/MOMENTUM/TRANSITS TO WRITE
        iteration = 0
        save_iteration = 1

        if ((wrtorb .eq. 1) .or. (wrtel .eq. 1)) then
            allocate (storeorb(Norb, DIMMAX))
            storeorb = zero
            call store_orb(save_iteration, zero, mass, r1, storeorb)
        end if

        if (wrtconst .eq. 1) then
            Etot = zero
            Eold = zero
            htot = zero
            hold = zero
            call compute_con(mass, r1, Eold, hold)
            allocate (storecon(5, DIMMAX))
            storecon = zero
            call store_con(save_iteration, zero, Eold, Eold, hold, hold, storecon)
        end if

        if ((idtra .ge. 1) .and. (idtra .le. NB)) then
            allocate (storetra(NBDIM + 6, DIMMAX))
            storetra = zero
            allocate (stat_tra(NB, DIMMAX))
            stat_tra = 0
            jtra = 0
        end if
        ! ==================

        call set_transiting_bodies(idtra, do_transit_check, body_transiting_start, body_transiting_end)
        n_body = size(mass)
        call set_transit_references_all_bodies(tepoch, obsDataIn, Tref)
        allocate (period(n_body), sma(n_body), ecc(n_body), inc(n_body), meanA(n_body), argp(n_body), trueA(n_body), longN(n_body))
        allocate (dttau(n_body))
        call elements(mass, rin, period, sma, ecc, inc, meanA, argp, trueA, longN, dttau)
        allocate (Pref(n_body - 1))
        Pref = period(2:n_body)

        trun1 = t_all_steps(1)
        integration: do iteration = 1, n_all

            trun2 = t_all_steps(iteration)
            integration_step = trun2 - trun1
            call one_forward_step(mass, radius, r1, integration_step, working_step, Hc, r2)
            if (.not. Hc) then
                exit integration
            end if

            ! ! ===============================
            ! RV computation
            if ((nRV .gt. 0) .and. (all_idx(iteration) .gt. 0)) then
                call get_and_set_RV_data(mass, r2, all_idx(iteration), obsDataIn, simRV)
            end if

            if (do_transit_check) then
                ! T0 check (to compare and all)
                Tx = tepoch + trun1 + half * integration_step

                ! T0 check
                do i_body = body_transiting_start, body_transiting_end
                    ! change of the sign of X or Y coordinate of two consecutive steps
                    call transit_conditions(r1(X(i_body):VZ(i_body)), r2(X(i_body):VZ(i_body)),&
                        &A, B, AB, ABflag, Zflag)

                    if (ABflag .and. Zflag) then

                        jtra = jtra + 1
                        call all_transits(tepoch, jtra, i_body, mass, radius, r1, r2,&
                            &trun1, integration_step, stat_tra, storetra)

                        if (jtra .eq. DIMMAX) then
                            call write_tra(jtra, utra, stat_tra, storetra)
                            jtra = 0
                            stat_tra = 0
                            storetra = zero
                        end if

                        ! if (allocated(obsDataIn % obsT0)) then
                        !     ntx = obsDataIn % obsT0(i_body - 1) % nT0
                        !     if (ntx .gt. 0) then
                        !         ! Tr = obsDataIn%obsT0(i_body-1)%Tephem
                        !         ! Pr = obsDataIn%obsT0(i_body-1)%Pephem
                        !         Tr = Tref(i_body - 1)
                        !         Pr = Pref(i_body - 1)
                        !         if (Pr .gt. zero) then ! only if we have Pephem > 0

                        !             epox = nint(((Tx - Tr) / Pr))
                        !             allocate (epo_obs(ntx))
                        !             epo_obs = nint(((obsDataIn % obsT0(i_body - 1) % T0 - Tr) / Pr))
                        !             ! check if the epoch of the mean time of two consecutive steps is within the observed epochs
                        !             if (any(epox .eq. epo_obs)) then
                        !                 call check_T0(tepoch, Tr, Pr, i_body, mass, radius, r1, r2,&
                        !                     &trun1, integration_step, simT0, Hc)
                        !             end if
                        !             if (allocated(epo_obs)) deallocate (epo_obs)
                        !         end if ! Pr
                        !     end if ! ntx
                        ! end if ! allocated(obsDataIn%obsT0)
                        call check_if_observed_transit(i_body, Tx, Tref, Pref, obsDataIn, search_obs_transit)
                        if (search_obs_transit) then
                            Tr = Tref(i_body - 1)
                            Pr = Pref(i_body - 1)
                            call check_T0(tepoch, Tr, Pr, i_body, mass, radius, r1, r2,&
                                &trun1, integration_step, simT0, Hc)
                        end if

                        if (n_excluded .gt. 0) then
                            call check_if_excluded_timerange(i_body, mass, radius,&
                                &tepoch, trun1, trun2, r1, r2, integration_step,&
                                &simT0(i_body-1),transit_in_excluded)
                                if (transit_in_excluded) Hc = .false.
                        end if

                    end if ! ABflag, Zflag
                end do ! i_body
            end if ! end do_transit_check
            ! ========

            ! check if it has to compute or not 'something' like orbits/constants
            if ((wrtorb .eq. 1) .or. (wrtel .eq. 1) .or. (wrtconst .eq. 1)) then

                if (all_idx(iteration) .eq. 0) then
                    save_iteration = save_iteration + 1
                    if ((wrtorb .eq. 1) .or. (wrtel .eq. 1)) then
                        call store_orb(save_iteration, trun2, mass, r2, storeorb) ! save r_save!!
                        if (save_iteration .eq. DIMMAX) then ! write into file if reach max dimension
                            if (wrtorb .eq. 1) call write_file(DIMMAX, uorb, fmorb, storeorb)
                            if (wrtel .eq. 1) call write_elem(DIMMAX, uele, fmele, mass, storeorb)
                            storeorb = zero ! reset the storeorb variable
                        end if
                    end if

                    ! check and store constants of motion
                    if (wrtconst .eq. 1) then
                        call compute_con(mass, r2, Etot, htot) ! if selected it computes the Energy and Angular momentum
                        call store_con(save_iteration, trun2, Etot, Eold, htot, hold, storecon)
                        if (save_iteration .eq. DIMMAX) then ! write into file if reach max dimension
                            call write_file(DIMMAX, ucon, fmcon, storecon)
                            storecon = zero
                        end if
                    end if

                    if (save_iteration .eq. DIMMAX) save_iteration = 0

                end if

            end if

            r1 = r2
            trun1 = trun2

            if (.not. Hc) then
                exit integration
            end if

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

        deallocate (r1, r2)
        deallocate (X, Y, Z)
        deallocate (VX, VY, VZ)
        deallocate (period, sma, ecc, inc, meanA, argp, trueA, longN, dttau)
        deallocate (Tref, Pref)
        deallocate (t_all_steps)
        deallocate (all_idx)
        flush (6)

        return
    end subroutine ode_forward_output_data
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
    !
    ! **Input/Output**
    ! Hc                  == Hill chek, .true. if stable, .false. if unstable for Hill's condition
    subroutine ode_forward_output_nodata(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
        &mass, radius, rin, time_to_int, Hc)
        ! Input
        integer, intent(in)::uorb, ucon
        integer, dimension(:), intent(in)::uele, utra
        character(*), intent(in)::fmorb, fmcon, fmele
        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::time_to_int
        ! Input/Output
        logical, intent(inout)::Hc
        ! Local
        type(dataObs)::empty_data
        type(dataRV)::empty_simrv
        type(dataT0), dimension(:), allocatable::empty_simt0

        ! nRV, nTTs, and nDurs already set to 0 by default in empty data and sim

        call ode_forward_output_data(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
            &mass, radius, rin, time_to_int, empty_data, empty_simrv, empty_simt0, Hc)

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
      &obsDataIn, simRV, simT0, Hc)
        ! Input
        real(dp), dimension(:), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN
        type(dataObs), intent(in)::obsDataIn
        ! Output
        type(dataRV), intent(out)::simRV
        type(dataT0), dimension(:), allocatable, intent(out)::simT0
        logical, intent(out)::Hc

        ! local variables
        integer::nRV, nTTs, nT0, ibd, nexcl
        real(dp), dimension(:), allocatable::sma
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp)::dt1, dt2, step

        ! if not defined here, the variables are global and defined in
        ! constants.f90 or parameters.f90

        Hc = .true.

        ! set RV and T0s data type for simulated ones
        nRV = obsDataIn % obsRV % nRV
        nTTs = obsDataIn % nTTs

        ! RV
        if (nRV .gt. 0) then
            call init_dataRV(nRV, simRV)
            simRV % nRV = 0
            if (.not. allocated(simRV % jitter)) allocate (simRV % jitter(nRVset))
            simRV % jitter = zero
            if (.not. allocated(simRV % gamma)) allocate (simRV % gamma(nRVset))
            simRV % gamma = zero
            simRV % RVsetID = obsDataIn % obsRV % RVsetID
        end if
        ! T0
        if (nTTs .gt. 0) then
            allocate (simT0(NB - 1))
            do ibd = 1, NB - 1
                nT0 = obsDataIn % obsT0(ibd) % nT0
                call init_dataT0(nT0, simT0(ibd), durcheck)
                simT0(ibd) % nT0 = 0
                simT0(ibd) % nDur = 0
                !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
                simT0(ibd) % eT0 = one / s24h ! 1s in day
                simT0(ibd) % edur = one / 60.0_dp ! 1s in min
                if (n_excluded .gt. 0)then
                    nexcl = obsDataIn % obsT0(ibd) % n_excl
                    if (nexcl .gt. 0)then
                        call init_dataT0_excluded(nexcl, obsData%obsT0(ibd)%excluded_time_ranges, simT0(ibd))
                    end if
                end if
            end do
        end if

        NBDIM = 6 * NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        allocate (sma(NB))
        sma = zero
        ! call get_semax_vec(mass(1), mass(2:NB), period(2:NB), sma(2:NB))
        call period_to_sma(mass(1), mass(2:NB), period(2:NB), sma(2:NB))

        call kepelements2statevector(mass, sma, ecc, meanA, argp, inc, longN, ra0)
        ra1 = ra0
        deallocate (sma)

        dt1 = tstart - tepoch
        dt2 = dt1 + tint
        step = minval(period(2:NB)) / ten

        if (dt1 .lt. zero) then

            call ode_forward_data(mass, radius, ra1, dt1, -step, obsDataIn, simRV, simT0, Hc)

            if (Hc) then
                if (abs(dt1) .le. tint) then
                    call ode_forward_data(mass, radius, ra1, dt2, step, obsDataIn, simRV, simT0, Hc)
                end if
            end if
        else
            call ode_forward_data(mass, radius, ra1, dt2, step, obsDataIn, simRV, simT0, Hc)

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
        real(dp), dimension(:), allocatable::ra0, ra1
        ! real(dp)::dt1,dt2
        logical::stable

        type(dataRV)::simRV
        type(dataT0), dimension(:), allocatable::simT0

        logical::checkpar

        integer::ndata
        integer::nRV
        integer::nTTs, nDurs !,nT0
        integer::ibd
        integer::ns, ne

        ! ! DEBUG
        ! logical::debugging
        ! character(512)::debug

        ! debugging=.true.

        ndata = obsData % ndata
        nRV = obsData % obsRV % nRV
        nTTs = obsData % nTTs
        nDurs = obsData % nDurs

        stable = .true.
        ! resw = zero
        checkpar = .true.

        allocate (mass(NB), radius(NB), period(NB), sma(NB), ecc(NB))
        allocate (argp(NB), meanA(NB), inc(NB), longN(NB))

        call convert_parameters(allpar, par, mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)

        ! write(debug, '(a,l)')"----------- ode_lm: checkpar: ",checkpar

        if (.not. checkpar) then
            resw = set_max_residuals(ndata)
            return
        end if

        ! +++ FROM HERE CALL NEW (2021-11-18) SUBROUTINE ode_elements_to_data
        call ode_keplerian_elements_to_data(mass, radius, period, ecc, argp, meanA, inc, longN,&
            &obsData, simRV, simT0, stable)

        if (.not. stable) then
            resw = set_max_residuals(ndata)
        else

            if (nRVset .gt. 0) then
                ! set jitter into simRV
                ns = nkel + 1
                ne = nkel + nRVset
                simRV % jitter = two**par(ns:ne)
                ! set gamma into simRV
                ns = ne + 1
                ne = ne + nRVset
                simRV % gamma = par(ns:ne)
                call set_gamma_rv(simRV % gamma, simRV % RVsetID, simRV % gamma_rv)
                ! compute RV trend
                if (rv_trend_order .gt. 0) then
                    ns = ne + 1
                    ne = ne + rv_trend_order
                    call addRVtrend(simRV % jd, tepoch, par(ns:ne), simRV % trend)
                end if
            end if

            ! call set_weighted_residuals(obsData, simRV, simT0, resw, oc_fit)
            ! write(debug, '(a,a,4(i4,1x))')trim(debug)," obsData%obsT0(:)nT0 = ",obsData%obsT0(:)%nT0
            call set_weighted_residuals(obsData, simRV, simT0, resw)

        end if

        ! -- DEBUG
        ! call write_T0_s(simT0)

        call deallocate_dataRV(simRV)
        if (nTTs .gt. 0) then
            ! -- DEBUG
            do ibd = 1, NB - 1
                call deallocate_dataT0(simT0(ibd))
            end do
        end if

        if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)
        if (allocated(ra0)) deallocate (ra0, ra1)

        ! write(debug, '(a,a,es25.15E3)')trim(debug)," sum(resw*resw) = ",sum(resw*resw)
        ! if (debugging)then
        !     write(*, *)trim(debug)
        !     flush(6)
        ! end if

        return
    end subroutine ode_lm
    ! ================================================================================

    ! ================================================================================
    ! subroutines called by the bootstrap to calculate the new set of T_0 and RV from
    ! the fitted parameters
    subroutine ode_parok(allpar, par, simRV, simT0, resw,&
      &chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)
        ! Input
        real(dp), dimension(:), intent(in)::allpar, par
        ! Input/Output
        type(dataRV), intent(inout)::simRV
        type(dataT0), dimension(:), intent(inout)::simT0
        ! Output
        real(dp), dimension(:), intent(out)::resw
        real(dp), intent(out)::chi_square, reduced_chi_square, lnLikelihood, ln_const, bic
        ! Local
        real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp)::dt1, dt2, step
        logical::Hc

        logical::checkpar
        ! logical::gls_check

        integer::ndata, nRV, nTTs, nT0, nexcl
        integer::ns, ne

        integer::ibd

        ndata = obsData % ndata
        nRV = obsData % obsRV % nRV
        ! nRVset=obsData%obsRV%nRVset
        nTTs = obsData % nTTs

        Hc = .true.
        resw = zero
        chi_square = -one

        allocate (mass(NB), radius(NB), period(NB), sma(NB), ecc(NB),&
            &argp(NB), meanA(NB), inc(NB), longN(NB))

        call convert_parameters(allpar, par, mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)

        if (nRV .gt. 0) then
            call init_dataRV(nRV, simRV)
            simRV % nRV = 0
            if (.not. allocated(simRV % jitter)) allocate (simRV % jitter(nRVset))
            simRV % jitter = zero
            if (.not. allocated(simRV % gamma)) allocate (simRV % gamma(nRVset))
            simRV % gamma = zero
            simRV % RVsetID = obsData % obsRV % RVsetID
        end if
!
        if (nTTs .gt. 0) then
            do ibd = 1, NB - 1
                nT0 = obsData % obsT0(ibd) % nT0
                call init_dataT0(nT0, simT0(ibd), durcheck)
                simT0(ibd) % nT0 = 0
                simT0(ibd) % nDur = 0
                !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
                simT0(ibd) % eT0 = one / s24h ! 1s in day
                simT0(ibd) % edur = one / 60.0_dp ! 1s in min
                if (n_excluded .gt. 0)then
                    nexcl = obsData % obsT0(ibd) % n_excl
                    if (nexcl .gt. 0)then
                        call init_dataT0_excluded(nexcl, obsData%obsT0(ibd)%excluded_time_ranges, simT0(ibd))
                    end if
                end if
            end do
        end if

        if (.not. checkpar) then
            if (nRV .gt. 0) simRV % RV = zero

            if (nTTs .gt. 0) then
                do ibd = 1, NB - 1
                    simT0(ibd) % T0 = zero
                end do
            end if

            resw = set_max_residuals(ndata)
            chi_square = resmax
            call set_fitness_values(par, chi_square,&
              &reduced_chi_square, lnLikelihood, ln_const, bic)
            return
        end if

        ! it is needed to define the which is the alarm coordinate for the transit detection

        NBDIM = 6 * NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(mass, sma, ecc, argp, meanA, inc, longN, ra0)
        ra1 = ra0

        dt1 = tstart - tepoch
        dt2 = dt1 + tint
        step = minval(period(2:NB)) / ten

        if (dt1 .lt. zero) then

            call ode_forward_data(mass, radius, ra1, dt1, -step, obsData, simRV, simT0, Hc)

            if (Hc) then
                if (abs(dt1) .le. tint) then
                    dt2 = dt1 + tint
                    call ode_forward_data(mass, radius, ra1, dt2, step, obsData, simRV, simT0, Hc)
                end if
            end if

        else

            dt2 = dt1 + tint
            call ode_forward_data(mass, radius, ra1, dt2, step, obsData, simRV, simT0, Hc)

        end if

        if (.not. Hc) then
            if (nRV .gt. 0) simRV % RV = zero

            if (nTTs .gt. 0) then
                do ibd = 1, NB - 1
                    simT0(ibd) % T0 = zero
                end do
            end if

            resw = set_max_residuals(ndata)
            chi_square = resmax

        else

            if (nRVset .gt. 0) then
                ! set jitter into simRV
                ns = nkel + 1
                ne = nkel + nRVset
                simRV % jitter = two**par(ns:ne)
                ! set gamma into simRV
                ns = ne + 1
                ne = ne + nRVset
                simRV % gamma = par(ns:ne)
                call set_gamma_rv(simRV % gamma, simRV % RVsetID, simRV % gamma_rv)
                ! compute RV trend
                if (rv_trend_order .gt. 0) then
                    ns = ne + 1
                    ne = ne + rv_trend_order
                    call addRVtrend(simRV % jd, tepoch, par(ns:ne), simRV % trend)
                end if
            end if

            ! call set_weighted_residuals(obsData, simRV, simT0, resw, oc_fit)
            call set_weighted_residuals(obsData, simRV, simT0, resw)

            chi_square = sum(resw * resw)
        end if
        call set_fitness_values(par, chi_square,&
            &reduced_chi_square, lnLikelihood, ln_const, bic)

        if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)
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
        real(dp), dimension(:), allocatable::ra0, ra1
        logical::Hc

        type(dataRV)::simRV
        type(dataT0), dimension(:), allocatable::simT0

        integer::ndata, nRV, nTTs, ibd
        integer::ns, ne

        logical::checkpar
        ! logical::gls_check

        ndata = oDataIn % ndata
        nRV = oDataIn % obsRV % nRV
        ! nRVset=oDataIn%obsRV%nRVset
        nTTs = oDataIn % nTTs

        Hc = .true.
        ! resw = zero
        allocate (mass(NB), radius(NB), period(NB), sma(NB), ecc(NB))
        allocate (argp(NB), meanA(NB), inc(NB), longN(NB))

        call convert_parameters(allpar, par, mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)

        if (.not. checkpar) then
            resw = set_max_residuals(ndata)
            return
        end if

        call ode_keplerian_elements_to_data(mass, radius, period, ecc, argp, meanA, inc, longN,&
            &oDataIn, simRV, simT0, Hc)

        if (.not. Hc) then
            resw = set_max_residuals(ndata)
        else

            if (nRVset .gt. 0) then
                ! set jitter into simRV
                ns = nkel + 1
                ne = nkel + nRVset
                simRV % jitter = two**par(ns:ne)
                ! set gamma into simRV
                ns = ne + 1
                ne = ne + nRVset
                simRV % gamma = par(ns:ne)
                call set_gamma_rv(simRV % gamma, simRV % RVsetID, simRV % gamma_rv)
                ! compute RV trend
                if (rv_trend_order .gt. 0) then
                    ns = ne + 1
                    ne = ne + rv_trend_order
                    call addRVtrend(simRV % jd, tepoch, par(ns:ne), simRV % trend)
                end if
            end if

            ! call set_weighted_residuals(oDataIn, simRV, simT0, resw, oc_fit)
            call set_weighted_residuals(oDataIn, simRV, simT0, resw)

        end if

        if (rv_res_gls) then
            call check_periodogram(simRV % jd,&
                &oDataIn % obsRV % RV - (simRV % RV + simRV % gamma_rv + simRV % trend),&
                &simRV % eRV, period, Hc)
            if (.not. Hc) resw = set_max_residuals(ndata)
        end if

        call deallocate_dataRV(simRV)
        if (nTTs .gt. 0) then
            do ibd = 1, NB - 1
                call deallocate_dataT0(simT0(ibd))
            end do
        end if

        if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)
        if (allocated(ra0)) deallocate (ra0, ra1)

        return
    end subroutine ode_boot
    ! ================================================================================

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
    subroutine integration_info_to_observables(t_epoch, t_start, step_in, t_int,&
        &mass, radius, period, ecc, argp, meanA, inc, longN,&
        &obsData_in,&
        &id_transit_body, transit_flag, dur_check,&
        &simRV, simT0,&
        &stable)

        real(dp), intent(in)::t_epoch, t_start, step_in, t_int
        real(dp), dimension(:), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN

        type(dataObs), intent(in)::obsData_in

        integer, intent(in)::id_transit_body
        logical, dimension(:), intent(in)::transit_flag
        integer, intent(in)::dur_check

        type(dataRV), intent(out)::simRV
        type(dataT0), dimension(:), allocatable, intent(out)::simT0
        logical, intent(out)::stable

        integer::n_body, nb_dim
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp), dimension(:), allocatable::sma
        real(dp)::dt1, dt2

        integer::nRV, nTTs, nT14s, nexcl
        integer::ibd, nT0

        nRV = obsData_in % obsRV % nRV
        nTTs = obsData_in % nTTs
        nT14s = obsData_in % nDurs

        stable = .true.

        n_body = size(mass)
        allocate (sma(n_body))
        sma = zero
        ! call get_semax_vec(mass(1), mass(2:n_body), period(2:n_body), sma(2:n_body))
        call period_to_sma(mass(1), mass(2:n_body), period(2:n_body), sma(2:n_body))

        nb_dim = 6 * n_body
        if (.not. allocated(ra0)) allocate (ra0(nb_dim), ra1(nb_dim))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(mass, sma, ecc, meanA, argp, inc, longN, ra0)
        ra1 = ra0

        if (nRV .gt. 0) then
            call init_dataRV(nRV, simRV)
            simRV % nRV = 0
            simRV % jd = obsData_in % obsRV % jd
            simRV % RVsetID = obsData_in % obsRV % RVsetID
        end if

        if (nTTs .gt. 0) then
            allocate (simT0(n_body - 1))
            do ibd = 1, n_body - 1
                nT0 = obsData_in % obsT0(ibd) % nT0
                call init_dataT0(nT0, simT0(ibd), dur_check)
                simT0(ibd) % nT0 = 0
                simT0(ibd) % nDur = 0
                !  let's set error on TT and Dur to 1 sec! It is needed to avoid division by zero
                simT0(ibd) % eT0 = one / s24h ! 1s in day
                simT0(ibd) % edur = one / 60.0_dp ! 1s in min
                if (n_excluded .gt. 0)then
                    nexcl = obsData_in % obsT0(ibd) % n_excl
                    if (nexcl .gt. 0)then
                        call init_dataT0_excluded(nexcl, obsData_in%obsT0(ibd)%excluded_time_ranges, simT0(ibd))
                    end if
                end if
            end do
        end if

        dt1 = t_start - t_epoch
        dt2 = dt1 + t_int

        if (dt1 .lt. zero) then
            call ode_full_args(mass, radius, ra1, t_epoch, dt1, -step_in, obsData_in,&
              &simRV,&
              &id_transit_body, transit_flag, dur_check, simT0,&
              &stable)
            if (stable) then
                if (abs(dt1) .le. t_int) then
                    dt2 = dt1 + t_int
                    call ode_full_args(mass, radius, ra1, t_epoch, dt2, step_in, obsData_in,&
                      &simRV,&
                      &id_transit_body, transit_flag, dur_check, simT0,&
                      &stable)
                end if
            end if
        else
            dt2 = dt1 + t_int
            call ode_full_args(mass, radius, ra1, t_epoch, dt2, step_in, obsData_in,&
              &simRV,&
              &id_transit_body, transit_flag, dur_check, simT0,&
              &stable)
        end if

        ! unstable or some bad parameters case!
        if (.not. stable) then
            if (nRV .gt. 0) then
                simRV % RV = zero
                simRV % RV_stat = 0
            end if
            if (nTTs .gt. 0) then
                do ibd = 1, n_body - 1
                    if (simT0(ibd) % nT0 .gt. 0) then
                        simT0(ibd) % T0 = zero
                        simT0(ibd) % T0_stat = 0
                        if (dur_check .eq. 1) then
                            simT0(ibd) % dur = zero
                            simT0(ibd) % dur_stat = 0
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
    subroutine integration_info_to_data(t_epoch, t_start, step_in, t_int,&
        &mass, radius, period, ecc, argp, meanA, inc, longN,&
        &id_transit_body, transit_flag, dur_check,&
        &RV_sim,&
        &body_T0_sim, epo_sim,&
        &T0_sim, T14_sim, lambda_rm_sim, kep_elem_sim,&
        &stable)
        ! Input
        real(dp), intent(in)::t_epoch, t_start, step_in, t_int
        real(dp), dimension(:), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN
        integer, intent(in)::id_transit_body
        logical, dimension(:), intent(in)::transit_flag
        integer, intent(in)::dur_check
        ! Output
        real(dp), dimension(:), allocatable, intent(out)::RV_sim
        integer, dimension(:), allocatable, intent(out)::body_T0_sim, epo_sim
        real(dp), dimension(:), allocatable, intent(out)::T0_sim, T14_sim, lambda_rm_sim
        real(dp), dimension(:, :), allocatable, intent(out)::kep_elem_sim
        logical, intent(out)::stable
        ! Local
        type(dataRV)::simRV
        type(dataT0), dimension(:), allocatable::simT0

        integer::n_body
        integer::nRV, nTTs
        integer::ibd, nT0, nstart, nend
        integer, parameter::n_kepelem = 8

        n_body = size(mass)

        nRV = obsData % obsRV % nRV
        if (nRV .gt. 0) then
            allocate (RV_sim(nRV))
            RV_sim = zero
        end if

        nTTs = obsData % nTTs

        call integration_info_to_observables(t_epoch, t_start, step_in, t_int,&
            &mass, radius, period, ecc, argp, meanA, inc, longN,&
            &obsData, id_transit_body, transit_flag, dur_check,&
            &simRV, simT0, stable)

        if (nRV .gt. 0) then
            ! allocate (RV_sim(nRV))
            RV_sim = simRV % RV
        end if

        nstart = 1
        nend = 0
        if (nTTs .gt. 0) then
            allocate (body_T0_sim(nTTs), epo_sim(nTTs))
            allocate (T0_sim(nTTs), T14_sim(nTTs), lambda_rm_sim(nTTs))
            allocate (kep_elem_sim(nTTs, n_kepelem))
            do ibd = 1, n_body - 1
                nT0 = obsData % obsT0(ibd) % nT0
                if (nT0 .gt. 0) then
                    nend = nend + nT0
                    epo_sim(nstart:nend) = obsData % obsT0(ibd) % epo
                    T0_sim(nstart:nend) = simT0(ibd) % T0
                    T14_sim(nstart:nend) = simT0(ibd) % dur
                    lambda_rm_sim(nstart:nend) = simT0(ibd) % lambda_rm
                    body_T0_sim(nstart:nend) = ibd + 1
                    kep_elem_sim(nstart:nend, 1) = simT0(ibd) % period
                    kep_elem_sim(nstart:nend, 2) = simT0(ibd) % sma
                    kep_elem_sim(nstart:nend, 3) = simT0(ibd) % ecc
                    kep_elem_sim(nstart:nend, 4) = simT0(ibd) % inc
                    kep_elem_sim(nstart:nend, 5) = simT0(ibd) % meana
                    kep_elem_sim(nstart:nend, 6) = simT0(ibd) % argp
                    kep_elem_sim(nstart:nend, 7) = simT0(ibd) % truea
                    kep_elem_sim(nstart:nend, 8) = simT0(ibd) % longn
                    nstart = nstart + nT0
                end if
            end do
        end if

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
    subroutine ode_output(cpuid, isim, wrtid, allpar, par, resw)
        ! Input
        integer, intent(in)::cpuid, isim, wrtid
        real(dp), dimension(:), intent(in)::allpar, par
        ! Output
        real(dp), dimension(:), intent(out)::resw
        ! Local
        real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN
        real(dp), dimension(:), allocatable::ra0, ra1
        real(dp)::dt1, dt2
        logical::Hc

        type(dataRV)::simRV
        type(dataT0), dimension(:), allocatable::simT0

        logical::checkpar, gls_check

        real(dp)::chi_square, reduced_chi_square, lnLikelihood, ln_const, bic

        integer::ndata, nRV, nTTs, nDurs, ibd, nT0, nexcl
        integer::ns, ne

        ! units and file names to store and write to files
        integer::uorb, ucon
        character(512)::florb, fmorb, flcon, fmcon, fmele
        integer, dimension(:), allocatable::uele, utra
        character(512), dimension(:), allocatable::flele, fltra

        write (*, '(a)') ''
        write (*, '(a)') " EXECUTING SIMPLE INTEGRATION AND WRITING FINAL FILES"
        write (*, '(a)') ''

        ndata = obsData % ndata
        nRV = obsData % obsRV % nRV
        nTTs = obsData % nTTs
        nDurs = obsData % nDurs

        Hc = .true.
        resw = zero
        checkpar = .true.

        allocate (mass(NB), radius(NB), period(NB), sma(NB), ecc(NB))
        allocate (argp(NB), meanA(NB), inc(NB), longN(NB))

        call only_convert_parameters(allpar, par,&
            &mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
        if (.not. checkpar) then
            write (*, '(a)') ' checkpar after convert_parameters is FALSE ... BAD'
            if (ndata > 0) resw = set_max_residuals(ndata)
        end if

        flush (6)

        if (nRV .gt. 0) then
            call init_dataRV(nRV, simRV)
            simRV % nRV = 0
            if (.not. allocated(simRV % jitter)) allocate (simRV % jitter(nRVset))
            simRV % jitter = zero
            if (.not. allocated(simRV % gamma)) allocate (simRV % gamma(nRVset))
            simRV % gamma = zero
            simRV % RVsetID = obsData % obsRV % RVsetID
        end if

        if ((idtra .ge. 1) .and. (idtra .le. NB)) call set_file_tra(cpuid, isim, wrtid, utra, fltra)

        allocate (simT0(NB - 1))
        if (nTTs .gt. 0) then
            do ibd = 1, NB - 1
                nT0 = obsData % obsT0(ibd) % nT0
                call init_dataT0(nT0, simT0(ibd), durcheck)
                simT0(ibd) % nT0 = 0
                simT0(ibd) % nDur = 0
                !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
                simT0(ibd) % eT0 = one / s24h ! 1s in day
                simT0(ibd) % edur = one / 60.0_dp ! 1s in min
                if (n_excluded .gt. 0)then
                    nexcl = obsData % obsT0(ibd) % n_excl
                    if (nexcl .gt. 0)then
                        call init_dataT0_excluded(nexcl, obsData%obsT0(ibd)%excluded_time_ranges, simT0(ibd))
                    end if
                end if
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

        NBDIM = 6 * NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(mass, sma, ecc, meanA, argp, inc, longN, ra0)
        ra1 = ra0

        dt1 = tstart - tepoch
        dt2 = dt1 + tint
        if (dt1 .lt. zero) then

            call ode_forward_output_data(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
                &mass, radius, ra1, dt1, obsData, simRV, simT0, Hc)

            if (Hc) then
                if (abs(dt1) .le. tint) then
                    call ode_forward_output_data(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
                        &mass, radius, ra1, dt2, obsData, simRV, simT0, Hc)
                end if
            end if

        else

            call ode_forward_output_data(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
                &mass, radius, ra1, dt2, obsData, simRV, simT0, Hc)

        end if

        write (*, *)
        if (.not. Hc) then
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

        if (ndata .gt. 0) then

            chi_square = -one
            reduced_chi_square = -one
            lnLikelihood = half
            ln_const = zero
            bic = zero

            if (.not. Hc) then
                resw = set_max_residuals(ndata)

            else

                if (nRVset .gt. 0) then
                    ! set jitter into simRV
                    ns = nkel + 1
                    ne = nkel + nRVset
                    simRV % jitter = two**par(ns:ne)
                    ! set gamma into simRV
                    ns = ne + 1
                    ne = ne + nRVset
                    simRV % gamma = par(ns:ne)
                    call set_gamma_rv(simRV % gamma, simRV % RVsetID, simRV % gamma_rv)
                    ! compute RV trend
                    if (rv_trend_order .gt. 0) then
                        ns = ne + 1
                        ne = ne + rv_trend_order
                        call addRVtrend(simRV % jd, tepoch, par(ns:ne), simRV % trend)
                    end if
                end if

                ! call set_weighted_residuals(obsData, simRV, simT0, resw, oc_fit)
                call set_weighted_residuals(obsData, simRV, simT0, resw)

                if (simRV % nRV .gt. 0) then
                    call check_periodogram(obsData % obsRV % jd,&
                        &obsData % obsRV % RV - simRV % RV - simRV % gamma_rv - simRV % trend,&
                        &obsData % obsRV % eRV,&
                        &period, gls_check)
                end if

                call check_max_residuals(resw, ndata)

            end if

            if (simRV % nRV .gt. 0) then

                write (*, *) ""
                write (*, '(a,i4)') " RADIAL VELOCITIES found: ", simRV % nRV
                write (*, *) ""

                ! if (present(to_screen)) call write_RV(simRV)
                call write_RV(cpuid, isim, wrtid, simRV)

                call check_and_write_periodogram(cpuid, isim, wrtid,&
                    &obsData % obsRV % jd,&
                    &obsData % obsRV % RV - simRV % RV - simRV % gamma_rv - simRV % trend,&
                    &obsData % obsRV % eRV, period, gls_check)

            else

                write (*, *)
                write (*, '(a)') ' RADIAL VELOCITIES NOT FOUND'
                write (*, *)

            end if

            if (allocated(simT0) .and. sum(simT0(:) % nT0) .gt. 0) then
                write (*, *)
                write (*, '(a,i5)') " T0 SIM found ", sum(simT0(:) % nT0)
                if (durcheck .eq. 1) write (*, '(a,i5)') " T41 SIM found ", sum(simT0(:) % nDur)
                write (*, *)
                call write_T0(cpuid, isim, wrtid, simT0)
                flush (6)
                call write_T0(simT0)
                flush (6)

            else
                write (*, *)
                write (*, '(a)') ' TRANSIT TIMES NOT FOUND'
                write (*, *)
            end if

            ! Reduced Chi Squares etc for report/summary
            chi_square = sum(resw * resw)
            call set_fitness_values(par, chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)

        else ! ndata <= 0

            write (*, '(a)') ' ndata == 0: NO FITNESS SUMMARY (SCREEN AND FILES)'
            flush (6)

        end if  ! ndata

        call deallocate_dataRV(simRV)
        if (nTTs .gt. 0) then
            do ibd = 1, NB - 1
                call deallocate_dataT0(simT0(ibd))
            end do
        end if

        if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)
        if (allocated(ra0)) deallocate (ra0, ra1)

        return
    end subroutine ode_output

    ! ================================================================================

    subroutine ode_integrates(cpuid, isim, wrtid, mass, radius, period, sma, ecc, argp, meana, inc, longn)
        ! Input
        integer, intent(in)::cpuid, isim, wrtid
        real(dp), dimension(:), intent(in)::mass, radius, period, sma, ecc, argp, meana, inc, longn
        ! Local
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
        call outElements(isim, wrtid, mass, radius, period, sma, ecc, argp, meana, inc, longn)

        NBDIM = 6 * NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(mass, sma, ecc, meana, argp, inc, longn, ra0)

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
                &mass, radius, ra1, dt1, Hc)
            if (abs(dt1) .le. tint) then
                call ode_forward_output_nodata(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
                    &mass, radius, ra1, dt2, Hc)
            end if
        else
            call ode_forward_output_nodata(uorb, ucon, uele, utra, fmorb, fmcon, fmele,&
                &mass, radius, ra1, dt2, Hc)
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
    subroutine ode_orbit_to_full_observables(mass, radius, rin, time_to_int, wrt_time,&
        &last_tra, ttra_full, dur_full, lambda_rm_full, id_ttra_full, stats_ttra,&
        &last_rv, time_rv_nmax, rv_nmax, stats_rv,&
        &Hc)
        ! Input
        real(dp), dimension(:), intent(in)::mass, radius, rin
        real(dp), intent(in)::time_to_int, wrt_time
        ! Input/Output
        ! transit times variables to be updated!!
        integer, intent(inout)::last_tra ! if first call it is zero, otherwise it is the last transit position
        real(dp), dimension(:), intent(inout)::ttra_full, dur_full, lambda_rm_full
        integer, dimension(:), intent(inout)::id_ttra_full
        logical, dimension(:), intent(inout)::stats_ttra
        ! radial velocities variables to be updated!!
        integer, intent(inout)::last_rv
        real(dp), dimension(:), intent(inout)::time_rv_nmax, rv_nmax
        logical, dimension(:), intent(inout)::stats_rv
        ! Output
        logical, intent(out)::Hc
        ! Local
        type(dataObs)::empty_data
        integer::ntra_full

        real(dp)::rv_temp, ttra_temp, dur_tra_temp, alpha_temp
        logical::check_ttra
        integer::nrv_max

        real(dp), dimension(:), allocatable::r1, r2

        integer, dimension(:), allocatable::X, Y, Z
        integer, dimension(:), allocatable::VX, VY, VZ

        real(dp)::A, B, AB
        logical::ABflag, Zflag

        integer::i_body, body_transiting_start, body_transiting_end
        logical::dummy_transit

        real(dp)::integration_step, working_step
        real(dp), dimension(:), allocatable::t_all_steps
        integer, dimension(:), allocatable::all_idx
        integer::n_all
        real(dp)::trun1, trun2
        integer::iteration

        ntra_full = size(ttra_full)
        nrv_max = size(rv_nmax)

        Hc = .true.
        Hc = .true.
        if (close_encounter_check) then
            Hc = separation_mutual_Hill_check(mass, radius, rin, do_hill_check)
        end if
        if (amd_hill_check) call statevector_amd_hill_stability(mass, rin, Hc)
        if (.not. Hc) then
            return
        end if

        ! == NEW ==
        call set_checking_coordinates(NB, X, Y, Z, VX, VY, VZ)
        ! ========

        integration_step = sign(wrt_time, time_to_int)
        call set_check_steps(tepoch, time_to_int, integration_step, empty_data, n_all, t_all_steps, all_idx)

        allocate (r1(NBDIM), r2(NBDIM))
        r1 = rin
        working_step = step_0

        call set_transiting_bodies(1, dummy_transit, body_transiting_start, body_transiting_end)

        trun1 = t_all_steps(1)
        integration: do iteration = 1, n_all

            trun2 = t_all_steps(iteration)
            integration_step = trun2 - trun1
            call one_forward_step(mass, radius, r1, integration_step, working_step, Hc, r2)
            if (.not. Hc) then
                exit integration
            end if

            ! RV
            last_rv = last_rv + 1
            if (last_rv .gt. nrv_max) then
                Hc = .false.
                ! return
                exit integration
            end if
            call get_RV(mass, r2, rv_temp) ! computes the rv
            ! save time rv status
            time_rv_nmax(last_rv) = trun2
            rv_nmax(last_rv) = rv_temp
            stats_rv(last_rv) = .true.

            ! transit times!!
            ! loop on planets -> checks for all planets
            do i_body = body_transiting_start, body_transiting_end
                !
                ! change of the sign of X or Y coordinate of two consecutive steps
                call transit_conditions(r1(X(i_body):VZ(i_body)), r2(X(i_body):VZ(i_body)),&
                &A, B, AB, ABflag, Zflag)

                if (ABflag .and. Zflag) then
                    call transit_time(tepoch, i_body, mass, radius, r1, r2, trun1, integration_step,&
                        &ttra_temp, dur_tra_temp, alpha_temp, check_ttra)
                    if (check_ttra) then
                        last_tra = last_tra + 1
                        if (last_tra .gt. ntra_full) then
                            Hc = .false.
                            exit integration
                        end if
                        ttra_full(last_tra) = ttra_temp
                        dur_full(last_tra) = dur_tra_temp
                        lambda_rm_full(last_tra) = alpha_temp
                        id_ttra_full(last_tra) = i_body
                        stats_ttra(last_tra) = .true.
                    end if ! check_ttra
                end if ! ABflag, Zflag
            end do ! i_body

            r1 = r2
            trun1 = trun2

            if (.not. Hc) then
                exit integration
            end if
        end do integration

        deallocate (r1, r2)
        deallocate (X, Y, Z)
        deallocate (VX, VY, VZ)
        deallocate (t_all_steps)
        deallocate (all_idx)

        return
    end subroutine ode_orbit_to_full_observables
    ! ================================================================================

    ! ================================================================================
    subroutine ode_all_ttra_rv(wrt_time, mass, radius, period, sma, ecc, argp, meana, inc, longn,&
      &ttra_full, dur_full, lambda_rm_full, id_ttra_full, stats_ttra,&
      &time_rv_nmax, rv_nmax, stats_rv)
        ! Input
        real(dp), intent(in)::wrt_time
        real(dp), dimension(:), intent(in)::mass, radius, period, sma, ecc, argp, meana, inc, longn
        ! Output
        real(dp), dimension(:), intent(out)::ttra_full, dur_full, lambda_rm_full
        integer, dimension(:), intent(out)::id_ttra_full
        logical, dimension(:), intent(out)::stats_ttra
        real(dp), dimension(:), intent(out)::time_rv_nmax, rv_nmax
        logical, dimension(:), intent(out)::stats_rv
        ! Local
        real(dp), dimension(:), allocatable::ra0, ra1

        integer::last_tra, last_rv

        real(dp)::dt1, dt2
        logical::Hc

        Hc = .true.

        NBDIM = 6 * NB
        if (.not. allocated(ra0)) allocate (ra0(NBDIM), ra1(NBDIM))

        ! NEW VERSION 2017-11-21
        call kepelements2statevector(mass, sma, ecc, meana, argp, inc, longn, ra0)
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
            call ode_orbit_to_full_observables(mass, radius, ra1, dt1, -wrt_time,&
              &last_tra, ttra_full, dur_full, lambda_rm_full, id_ttra_full, stats_ttra,&
              &last_rv, time_rv_nmax, rv_nmax, stats_rv, Hc)
            if (.not. Hc) then
                if (allocated(ra0)) deallocate (ra0, ra1)
                return
            end if

            if (abs(dt1) .le. tint) then
                ! forward integration
                call ode_orbit_to_full_observables(mass, radius, ra1, dt2, wrt_time,&
                  &last_tra, ttra_full, dur_full, lambda_rm_full, id_ttra_full, stats_ttra,&
                  &last_rv, time_rv_nmax, rv_nmax, stats_rv, Hc)
                if (.not. Hc) then
                    if (allocated(ra0)) deallocate (ra0, ra1)
                    return
                end if
            end if

        else

            ! only forward integration
            call ode_orbit_to_full_observables(mass, radius, ra1, dt2, wrt_time,&
              &last_tra, ttra_full, dur_full, lambda_rm_full, id_ttra_full, stats_ttra,&
              &last_rv, time_rv_nmax, rv_nmax, stats_rv, Hc)

        end if

        if (allocated(ra0)) deallocate (ra0, ra1)

        return
    end subroutine ode_all_ttra_rv
    ! ================================================================================

    ! ================================================================================
    subroutine ode_keplerian_elements_to_orbits(time_steps, mass, radius, period, ecc, argp, meanA, inc, longN, orbits, check)
        ! Input
        real(dp), dimension(:), intent(in)::time_steps
        real(dp), dimension(:), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN
        ! Output
        ! real(dp), dimension(:, :), intent(out), allocatable::orbits
        real(dp), dimension(:, :), intent(out)::orbits
        logical, intent(out)::check
        ! Local
        integer::n_steps, iteration
        integer::n_body, nb_dim
        real(dp), dimension(:), allocatable::sma
        real(dp), dimension(:), allocatable::sv1, sv2
        real(dp)::trun1, trun2, integration_step
        real(dp)::working_step

        n_steps = size(time_steps)
        n_body = size(mass)
        nb_dim = n_body * 6

        if (.not. allocated(sma)) allocate (sma(n_body))
        sma = zero
        call period_to_sma(mass(1), mass(2:n_body), period(2:n_body), sma(2:n_body))

        if (.not. allocated(sv1)) allocate (sv1(nb_dim), sv2(nb_dim))
        call kepelements2statevector(mass, sma, ecc, meanA, argp, inc, longN, sv1)

        ! if (.not. allocated(orbits)) allocate (orbits(n_steps, nb_dim))
        orbits(1, :) = sv1

        check = .true.
        working_step = step_0 ! default value of the stepsize
        trun1 = time_steps(1)
        integration: do iteration = 2, n_steps

            ! update integration step
            trun2 = time_steps(iteration)
            integration_step = trun2 - trun1

            ! compute the orbit at from trun1 to trun2
            call one_forward_step(mass, radius, sv1, integration_step, working_step, check, sv2)
            ! save statevector
            orbits(iteration, :) = sv2
            ! update statevecto and time for the next iteration
            sv1 = sv2
            trun1 = trun2
            ! check if stable otherwise exit
            if (.not. check) then
                exit integration
            end if
        end do integration

        deallocate (sv1, sv2)
        deallocate (sma)

        return
    end subroutine ode_keplerian_elements_to_orbits

    ! ================================================================================

end module ode_run

