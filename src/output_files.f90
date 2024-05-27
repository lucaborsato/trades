module output_files
    use constants
    use custom_type
    use parameters
    use parameters_conversion
    use init_trades, only: get_unit, state2string, get_lnec_full, get_lnec
    use convert_type, only: string
    use celestial_mechanics, only: barycenter
    use utils, only: set_fitness_values, ln_priors
    implicit none

    interface write_RV
        module procedure write_RV_s, write_RV_f
    end interface write_RV

    interface write_T0
        module procedure write_T0_s, write_T0_f
    end interface write_T0

contains

    ! THIS IS VERY OLD....
    ! it writes to screen and into file the inital combination of the grid search and the output
    ! of the LM, in particular:
    ! simID = number that identify the simulation
    ! LM_info = flag for LM convergence
    ! Chi2, Chi2r = Chi square and reduced chi square of the simulation with those initial parameters
    subroutine write_simlst(cpuid, Ngrid, grid, infos)
        integer, intent(in)::cpuid, Ngrid
        real(dp), dimension(:, :), intent(in)::grid
        integer, dimension(:), intent(in)::infos
        integer::ulst, j
        character(512)::fllst
        character(80)::fmt
        real(dp)::rdof

        rdof = real(obsData%dof, dp)

        fllst = trim(path)//"sim.lst"
        fllst = trim(adjustl(fllst))
        ulst = get_unit(cpuid)
        open (ulst, file=trim(fllst))
        write (ulst, '(a,2x,a,8x,a,22x,a,22x,a,21x,a,25x,a,22x,a,20x,a)')&
            &"#    simID", "LM_info", "m[Msun]", "P[d]", "a[AU]", "e", "w[°]",&
            &"fitness", "fitness*dof"
        write (*, '(a,2x,a,8x,a,22x,a,22x,a,21x,a,25x,a,22x,a,20x,a)') &
            "#     simID", "LM_info", "m[Msun]", "P[d]", "a[AU]", "e", "w[°]",&
            &"fitness", "fitness*dof"
        fmt = adjustl("(i6,1x,i3,1x,1000("//trim(sprec)//",1x))")
        do j = 1, Ngrid
            write (ulst, trim(fmt)) j, infos(j), grid(:, j), (grid(6, j)*rdof)
            write (*, trim(fmt)) j, infos(j), grid(:, j), (grid(6, j)*rdof)
        end do
        close (ulst)

        return
    end subroutine write_simlst

    ! GOOD GRID WRITE SUMMARY
    subroutine write_grid_summary(cpuid, n_grid, lm_flag, grid_summary)
        integer, intent(in)::cpuid, n_grid, lm_flag
        real(dp), dimension(:, :), intent(in)::grid_summary

        character(512)::output_file, header
        character(80)::fmt_wrt
        integer::u_wrt, i_sim

        u_wrt = get_unit(cpuid)

        fmt_wrt = "(i6,1x,1000("//trim(adjustl(sprec))//",1x))"
        output_file = trim(path)//"summary_grid_sims_"//&
          &trim(adjustl(string(lm_flag)))//".dat"
        header = "# id_sim "//trim(adjustl(all_names_str))//" fitness_x_dof fitness"

        open (u_wrt, file=trim(output_file))
        write (u_wrt, '(a)') trim(header)
        do i_sim = 1, n_grid
            write (u_wrt, trim(fmt_wrt)) i_sim, grid_summary(i_sim, :)
        end do
        close (u_wrt)
        write (*, '(a,a)') " WRITTEN FILE:", trim(output_file)
        flush (6)

        return
    end subroutine write_grid_summary

    ! format string to write orbit file
    function fmtorbit() result(fmt)
        integer::ncol
        character(128)::fmt
        character(72)::col
        character(5)::nstr

        col = "(1x,"//trim(sprec)//")"
        ncol = 2+(NB*6)+1
        nstr = trim(adjustl(string(ncol)))
        fmt = "("//trim(nstr)//trim(col)//")"
        fmt = trim(adjustl(fmt))

        return
    end function fmtorbit

    ! format string to write const. of motion file
    function fmtconst() result(fmt)
        implicit none
        character(128)::fmt

        fmt = "(5(1x,"//trim(sprec)//"))"

        return
    end function fmtconst

    ! format string to write Keplerian elements file
    function fmtele() result(fmt)
        implicit none
        character(128)::fmt

        fmt = "(10(1x,"//trim(sprec)//"))"

        return
    end function fmtele

    ! ------------------------------------------------------------------ !
    ! write to screen input parameters
    subroutine write_lm_inputpar(cpuid, par)
        integer, intent(in)::cpuid
        real(dp), dimension(:), intent(in)::par
        integer::i
        character(80)::fmt

        write (*, '(a,i5,a)') " CPU ", cpuid, " PARAM: "
        fmt = adjustl("(1x,a6,"//trim(sprec)//")")
        do i = 1, nfit
            write (*, trim(fmt)) parid(i), par(i)
        end do

        return
    end subroutine write_lm_inputpar
    ! ------------------------------------------------------------------ !

    ! subroutine to write initial parameters of planets in a file
    ! during the call ode_out(...)
    subroutine outElements(isim, wrtid, mass, radius, period, sma, ecc, argp, meana, inc, longn)
!$      use omp_lib
        integer, intent(in)::isim, wrtid
        real(dp), dimension(:), intent(in)::mass, radius, period, sma, ecc, argp, meana, inc, longn
        integer::cpuid, uwrt, ii
        character(512)::fwrt, fmtw

        fmtw = adjustl('(1000('//trim(sprec)//'))')
        fwrt = trim(path)//trim(adjustl(string(isim)))//"_"//&
          &trim(adjustl(string(wrtid)))//"_initialElements.dat"
        fwrt = trim(adjustl(fwrt))
        cpuid = 1
!$      cpuid = omp_get_thread_num()+1
        uwrt = get_unit(cpuid)
        open (uwrt, file=trim(fwrt))
        write (uwrt, '(a)') "# M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg"
        do ii = 2, NB
            write (uwrt, trim(fmtw)) mass(ii), radius(ii), period(ii), sma(ii), ecc(ii),&
                &argp(ii), meana(ii), inc(ii), longn(ii)
        end do
        close (uwrt)

        return
    end subroutine outElements
    ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
    ! subroutine to store the orbital state vectors, ready to be written into file
    subroutine store_orb(pos, itime, mass, rw, storeorb)
        integer, intent(in)::pos
        real(dp), intent(in)::itime
        real(dp), dimension(:), intent(in)::mass, rw
        real(dp), dimension(:, :), intent(inout)::storeorb
        real(dp), dimension(6)::bar
        real(dp), dimension(:), allocatable::rbar

        allocate (rbar(NBDIM))
        storeorb(1, pos) = tepoch+itime
        call barycenter(mass, rw, bar, rbar)
        storeorb(2, pos) = -rbar(3)/speedaud
        storeorb(3:NBDIM+2, pos) = rw
        storeorb(NBDIM+3, pos) = (-rbar(6)*AU/s24h)
        deallocate (rbar)

        return
    end subroutine store_orb

    ! subroutine to store constants of motion, ready to be written into file
    subroutine store_con(pos, itime, Etot, Eold, htot, hold, storecon)
        integer, intent(in)::pos
        real(dp), intent(in)::itime, Etot, Eold, htot, hold
        real(dp), dimension(:, :), intent(inout)::storecon

        storecon(1, pos) = tepoch+itime
        storecon(2, pos) = htot
        storecon(3, pos) = htot-hold
        storecon(4, pos) = Etot
        storecon(5, pos) = Etot-Eold

        return
    end subroutine store_con
    ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
    ! write orbit state vector into file
    subroutine write_orb(maxpos, uorb, fmorb, storeorb)
        integer, intent(in)::uorb, maxpos
        character(*), intent(in)::fmorb
        real(dp), dimension(:, :), intent(inout)::storeorb
        integer::j

        do j = 1, maxpos
            write (uorb, trim(adjustl(fmorb))) storeorb(:, j)
        end do

        return
    end subroutine write_orb

    ! generic subroutine to write something into file
    subroutine write_file(maxpos, uwrt, fmwrt, store)
        integer, intent(in)::uwrt, maxpos
        character(*), intent(in)::fmwrt
        real(dp), dimension(:, :), intent(inout)::store
        integer::j

        do j = 1, maxpos
            write (uwrt, trim(adjustl(fmwrt))) store(:, j)
        end do

        return
    end subroutine write_file

    ! write keplerian orbital elements into file
    subroutine write_elem(pos, uele, fmele, mass, storeorb)
        use celestial_mechanics, only: elements
        integer, intent(in)::pos
        integer, dimension(:), intent(in)::uele
        character(*), intent(in)::fmele
        real(dp), dimension(:), intent(in)::mass
        real(dp), dimension(:, :), intent(in)::storeorb
        real(dp), dimension(:), allocatable::period, sma, ecc, inc, meana, argp, truea, longn, dttau
        integer::j1, j2

        allocate (period(NB), sma(NB), ecc(NB), inc(NB), meana(NB), argp(NB), truea(NB), longn(NB), dttau(NB))
        j1loop: do j1 = 1, pos
            call elements(mass, storeorb(3:NBDIM+2, j1), period, sma, ecc, inc, meana, argp, truea, longn, dttau)
            j2loop: do j2 = 2, NB
                write (uele(j2), trim(adjustl(fmele))) storeorb(1, j1),&
                    &period(j2), sma(j2), ecc(j2), inc(j2), meana(j2), argp(j2), longn(j2),&
                    &truea(j2), dttau(j2)
            end do j2loop
        end do j1loop
        deallocate (period, sma, ecc, inc, meana, argp, truea, longn, dttau)

        return
    end subroutine write_elem

    ! write transit time, light travel-time effect and state vector (at transit time)
    ! into file
    subroutine write_tra(pos, utra, stat_tra, storetra)
        integer, intent(in)::pos
        integer, dimension(:), intent(in)::utra
        integer, dimension(:, :), intent(in)::stat_tra
        real(dp), dimension(:, :), intent(in)::storetra
        integer::j1, j2
        character(80)::fmt

        fmt = adjustl("(1000("//trim(sprec)//",1x))") ! '(1000(es23.16,1x))'
        j1loop: do j1 = 1, pos
            j2loop: do j2 = 2, NB
                if (stat_tra(j2, j1) .eq. 1) then
                    write (utra(j2), trim(fmt)) storetra(:, j1)
                end if
            end do j2loop
        end do j1loop

        return
    end subroutine write_tra
    ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
    ! defines the unit and the name of orbit file
    subroutine set_file_orb(cpuid, isim, wrtid, uorb, florb)
        integer, intent(in)::cpuid, isim, wrtid
        integer, intent(inout)::uorb
        character(512), intent(inout)::florb
        character(512)::strState

        uorb = 0
        uorb = get_unit(cpuid)
        florb = trim(path)//trim(adjustl(string(isim)))//"_"//&
          &trim(adjustl(string(wrtid)))//"_rotorbit.dat"
        florb = trim(adjustl(florb))
        open (uorb, file=trim(florb))
        strState = trim(adjustl(state2string(NB)))
        write (uorb, '(a,a,a)') "# Time_JD LTE_d ", trim(strState), " RV_ms^-1"

        return
    end subroutine set_file_orb

    ! defines the unit and the name of constants file
    subroutine set_file_con(cpuid, isim, wrtid, ucon, flcon)
        integer, intent(in)::cpuid, isim, wrtid
        integer, intent(inout)::ucon
        character(512), intent(inout)::flcon

        ucon = 0
        ucon = get_unit(cpuid)
        flcon = trim(path)//trim(adjustl(string(isim)))//"_"//&
          &trim(adjustl(string(wrtid)))//"_constants.dat"
        flcon = trim(adjustl(flcon))
        open (ucon, file=trim(flcon))
        write (ucon, '(a)') "# Time htot dh Etot dE"

        return
    end subroutine set_file_con

    ! defines the unit and the name of keplerian orbital elements files
    subroutine set_file_elem(cpuid, isim, wrtid, uele, flele)
        integer, intent(in)::cpuid, isim, wrtid
        integer, dimension(:), allocatable, intent(inout)::uele
        character(512), dimension(:), allocatable, intent(inout)::flele
        character(512)::strElem
        integer::je

        allocate (uele(NB), flele(NB))
        uele = 0
        do je = 2, NB
            uele(je) = get_unit(cpuid)
            flele(je) = trim(path)//trim(adjustl(string(isim)))//"_"//&
              &trim(adjustl(string(wrtid)))//"_NB"//&
              &trim(adjustl(string(je)))//"_elements.dat"
            flele(je) = trim(adjustl(flele(je)))
            open (uele(je), file=trim(flele(je)))
            strElem = "# Time "//&
                &"P"//trim(adjustl(string(je)))//"_d "//&
                &"a"//trim(adjustl(string(je)))//"_AU "//&
                &"e"//trim(adjustl(string(je)))//" "//&
                &"i"//trim(adjustl(string(je)))//"_deg "//&
                &"mA"//trim(adjustl(string(je)))//"_deg "//&
                &"w"//trim(adjustl(string(je)))//"_deg "//&
                &"lN"//trim(adjustl(string(je)))//"_deg "//&
                &"f"//trim(adjustl(string(je)))//"_deg "//&
                &"dtau"//trim(adjustl(string(je)))//"_d "
            write (uele(je), '(a)') trim(adjustl(strElem))
        end do

        return
    end subroutine set_file_elem

    ! closes and deallocates units and files for orbital elements
    subroutine close_elem(uele, flele)
        integer, dimension(:), allocatable, intent(inout)::uele
        character(512), dimension(:), allocatable, intent(inout)::flele
        integer::je

        do je = 2, NB
            close (uele(je))
        end do
        deallocate (uele, flele)

        return
    end subroutine close_elem

    ! defines the unit and the name of transit files
    subroutine set_file_tra(cpuid, isim, wrtid, utra, fltra)
        integer, intent(in)::cpuid, isim, wrtid
        integer, dimension(:), allocatable, intent(inout)::utra
        character(512), dimension(:), allocatable, intent(inout)::fltra
        character(512)::strState, strTra
        integer::j

        strState = state2string(NB)
        allocate (utra(NB), fltra(NB))
        utra = 0
        fltra = ""
        do j = 2, NB
            utra(j) = get_unit(cpuid)
            fltra = trim(path)//trim(adjustl(string(isim)))//"_"//&
              &trim(adjustl(string(wrtid)))//"_NB"//&
              &trim(adjustl(string(j)))//"_tra.dat"
            fltra(j) = trim(adjustl(fltra(j)))
            open (utra(j), file=trim(fltra(j)))
            strTra = "# ttra_"//trim(adjustl(string(j)))//&
                &" LTE_"//trim(adjustl(string(j)))//&
                &" t1_"//trim(adjustl(string(j)))//&
                &" t2_"//trim(adjustl(string(j)))//&
                &" t3_"//trim(adjustl(string(j)))//&
                &" t4_"//trim(adjustl(string(j)))//&
                &" "//trim(adjustl(strState))
            write (utra(j), '(a)') trim(adjustl(strTra))
        end do

        return
    end subroutine set_file_tra

    ! closes and deallocates units and files for transit times
    subroutine close_tra(utra, fltra)
        integer, dimension(:), allocatable, intent(inout)::utra
        character(512), dimension(:), allocatable, intent(inout)::fltra
        integer::j

        do j = 2, NB
            close (utra(j))
        end do
        deallocate (utra, fltra)

        return
    end subroutine close_tra
    ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
    ! write to screen the RV simulated compared to the RV observated
!   subroutine write_RV_s(gamma,RV_sim,RV_stat)
    subroutine write_RV_s(simRV)
!     real(dp),dimension(:,:),intent(in)::gamma
!     real(dp),dimension(:),intent(in)::RV_sim
!     integer,dimension(:),intent(in)::RV_stat
        type(dataRV), intent(in)::simRV

        real(dp), dimension(:), allocatable::RV_simwrt
        ! real(dp),dimension(:,:),allocatable::gamma_wrt
        integer::j, nRV
        character(80)::fmt

        real(dp),dimension(:),allocatable::jitter
        integer::k, nj


        nRV = simRV%nRV
        ! allocate(RV_simwrt(nRV),gamma_wrt(nRV,2))
        ! call setWriteRV(simRV,RV_simwrt,gamma_wrt)
        allocate (RV_simwrt(nRV), jitter(nRV))
        RV_simwrt = simRV%RV+simRV%gamma_rv+simRV%trend

        jitter = zero
        nj = size(simRV%jitter)
        do k=1,nj
            do j=1,nRV
                if (obsData%obsRV%RVsetID(j) .eq. k) then
                    jitter(j) = simRV%jitter(k)
                end if
            end do
        end do

        write (*, '(a)') "# JD RVobs eRVobs rv_sim RV_sim gamma_rv trend jitter &
          &RVsetID RV_stat"
        fmt = adjustl("(7("//trim(sprec)//",1x),i3,1x,i3)")
        do j = 1, nRV
            write (*, trim(fmt)) obsData%obsRV%jd(j), obsData%obsRV%RV(j),&
              &obsData%obsRV%eRV(j), simRV%RV(j), RV_simwrt(j),&
              &simRV%gamma_rv(j), simRV%trend(j), jitter(j),&
              &obsData%obsRV%RVsetID(j),&
              &simRV%RV_stat(j)
        end do
        ! deallocate(RV_simwrt,gamma_wrt)
        deallocate (RV_simwrt,jitter)

        return
    end subroutine write_RV_s

    ! --------------
    ! write into file the RV simulated compared to the RV observated
!   subroutine write_RV_f(cpuid,isim,wrtid,gamma,RV_sim,RV_stat)
    subroutine write_RV_f(cpuid, isim, wrtid, simRV)
        integer, intent(in)::cpuid, isim, wrtid ! wrtid = if 0 original parameters before LM, 1 parameters after LM
        type(dataRV), intent(in)::simRV
        real(dp), dimension(:), allocatable::RV_simwrt
        ! real(dp),dimension(:,:),allocatable::gamma_wrt
        integer::uRV, j
        character(512)::flRV
        character(80)::fmt
        integer::nRV

        real(dp),dimension(:),allocatable::jitter
        integer::k, nj

        nRV = simRV%nRV
        ! allocate(RV_simwrt(nRV),gamma_wrt(nRV,2))
        ! call setWriteRV(simRV,RV_simwrt,gamma_wrt)
        allocate (RV_simwrt(nRV), jitter(nRV))
        RV_simwrt = simRV%RV+simRV%gamma_rv+simRV%trend

        jitter = zero
        nj = size(simRV%jitter)
        do k=1,nj
            do j=1,nRV
                if (obsData%obsRV%RVsetID(j) .eq. k) then
                    jitter(j) = simRV%jitter(k)
                end if
            end do
        end do

        uRV = get_unit(cpuid)
        flRV = ""
        flRV = trim(path)//trim(adjustl(string(isim)))//"_"//&
          &trim(adjustl(string(wrtid)))//"_simRV.dat"
        flRV = trim(adjustl(flRV))
        open (uRV, file=trim(flRV))
!     fmt=adjustl("(a,28x,a,19x,a,18x,a,10x,2(a,"//trim(sprec)//"))")
        write (uRV, '(a)') "# JD RVobs eRVobs rv_sim RV_sim gamma_rv trend jitter &
          &RVsetID RV_stat"
        fmt = ""
        ! fmt = adjustl("(7("//trim(sprec)//",1x),i3,1x,i3)")
        fmt = adjustl("(8("//trim(sprec)//",1x),i3,1x,i3)")
        do j = 1, nRV
            write (uRV, trim(fmt)) obsData%obsRV%jd(j), obsData%obsRV%RV(j),&
              &obsData%obsRV%eRV(j), simRV%RV(j), RV_simwrt(j),&
              &simRV%gamma_rv(j), simRV%trend(j), jitter(j),&
              &obsData%obsRV%RVsetID(j),&
              &simRV%RV_stat(j)
        end do
        close (uRV)
        ! deallocate(RV_simwrt,gamma_wrt)
        deallocate (RV_simwrt,jitter)
        write (*, '(a,a)') " WRITTEN RV INTO FILE: ", trim(flRV)

        return
    end subroutine write_RV_f
    ! ------------------------------------------------------------------ !

    !=++=
    !=++= HERE
    !=++=

    ! ------------------------------------------------------------------ !
    ! write to screen the T0 simulated compared to the T0 observated
    subroutine write_T0_s(simT0)
        type(dataT0), dimension(:), intent(in)::simT0
!     real(dp),dimension(:,:),intent(in)::T0_sim
!     integer,dimension(:,:),intent(in)::T0_stat
        integer::i_body, i_T0, nT0o, nT0s, nT0
        character(128)::fmt,hea

        ! fmt = adjustl("(i6,2(4("//trim(sprec)//",1x),i3))")
        fmt = adjustl("i6,1x,4("//trim(sprec)//",1x),i3,1x,")
        hea = adjustl("a,5(1x,a),")
        do i_body = 1, NB-1
            nT0s = simT0(i_body)%nT0
            nT0o = obsData%obsT0(i_body)%nT0
            nT0 = nT0o

            if (nT0 .gt. 0) then

                if (durcheck .eq. 0) then ! do not check duration
                    fmt = adjustl(trim(fmt))//"i3"
                    hea = adjustl(trim(hea))//"1x,a"
                    write (*, "("//trim(hea)//")")&
                        &"# epoT0obs",&
                        &"T0obs", "eT0obs", "T0_sim", "T0obs-T0_sim",&
                        &"T0_stat",&
                        &"source_id"
                    ! write (*, '(a,1x,a,19x,a,18x,a,18x,a,5x,a)') &
                    !   "# epoT0obs ", " T0obs ", " eT0obs ", " T0_sim ",&
                    !   &" T0obs-T0_sim ", " T0_stat ", " source_id "
                    do i_T0 = 1, nT0
                        write (*, "("//trim(fmt)//")")&
                            &obsData%obsT0(i_body)%epo(i_T0),&
                            &obsData%obsT0(i_body)%T0(i_T0),&
                            &obsData%obsT0(i_body)%eT0(i_T0),&
                            &simT0(i_body)%T0(i_T0),&
                            &(obsData%obsT0(i_body)%T0(i_T0)-simT0(i_body)%T0(i_T0)),&
                            &simT0(i_body)%T0_stat(i_T0),&
                            &obsData%obsT0(i_body)%source_id(i_T0)
                    end do

                else ! do check duration

                    fmt = adjustl(trim(fmt))//"4("//trim(sprec)//",1x),2(i3,1x)"
                    hea = adjustl(trim(hea))//"7(1x,a)"
                    write (*, "("//trim(hea)//")")&
                       &"# epoT0obs",&
                       &"T0obs", "eT0obs", "T0_sim", "T0obs-T0_sim",&
                       &"T0_stat",&
                       &"Dur_obs", "eDur_obs", "Dur_sim", "Dur_obs-Dur_sim",&
                       &"Dur_stat ",&
                       &"source_id"
                    ! write (*, '(a,2(1x,a,19x,a,18x,a,18x,a,5x,a))')&
                    !     &"# epoT0obs ", " T0obs ", " eT0obs ", " T0_sim ",&
                    !     &" T0obs-T0_sim ", " T0_stat ",&
                    !     &" Dur_obs ", " eDur_obs ", " Dur_sim ", " Dur_obs-Dur_sim",&
                    !     &" Dur_stat ", " source_id "
                    do i_T0 = 1, nT0
                        write (*, "("//trim(fmt)//")")&
                            &obsData%obsT0(i_body)%epo(i_T0),&
                            &obsData%obsT0(i_body)%T0(i_T0), obsData%obsT0(i_body)%eT0(i_T0),&
                            &simT0(i_body)%T0(i_T0),&
                            &(obsData%obsT0(i_body)%T0(i_T0)-simT0(i_body)%T0(i_T0)),&
                            &simT0(i_body)%T0_stat(i_T0),&
                            &obsData%obsT0(i_body)%dur(i_T0), obsData%obsT0(i_body)%edur(i_T0),&
                            &simT0(i_body)%dur(i_T0), obsData%obsT0(i_body)%dur(i_T0)-simT0(i_body)%dur(i_T0),&
                            &simT0(i_body)%dur_stat(i_T0),&
                            &obsData%obsT0(i_body)%source_id(i_T0)
                    end do

                end if
                write (*, *) ""

            end if
        end do

        return
    end subroutine write_T0_s

    ! write into file the T0 simulated compared to the T0 observated
    subroutine write_T0_f(cpuid, isim, wrtid, simT0)
        integer, intent(in)::cpuid, isim, wrtid ! wrtid = if 0 original parameters before LM, 1 parameters after LM
        type(dataT0), dimension(:), intent(in)::simT0
!     real(dp),dimension(:,:),intent(in)::T0_sim
!     integer,dimension(:,:),intent(in)::T0_stat
        character(512)::flT0
        integer::uT0, i_body, i_T0, nT0o, nT0s, nT0
        character(128)::fmt, hea

        
        ! fmt = adjustl("i6,1x,4("//trim(sprec)//",1x),i3,1x,8("//trim(sprec)//",1x),i3")
        ! hea = adjustl("a,5(1x,a),9(1x,a)")
        fmt = adjustl("i6,1x,4("//trim(sprec)//",1x),i3,1x,")
        hea = adjustl("a,5(1x,a),")
        do i_body = 1, NB-1
            nT0s = simT0(i_body)%nT0
            nT0o = obsData%obsT0(i_body)%nT0

            ! write(*,*)" DEBUG: ** body id = ",i_body+1
            ! write(*,*)" DEBUG: ** sim nT0 = ",nT0s
            ! write(*,*)" DEBUG: ** obs nT0 = ",nT0o

            nT0 = nT0o

            if (nT0 .gt. 0) then

                uT0 = get_unit(cpuid)
                flT0 = ""
                flT0 = trim(path)//trim(adjustl(string(isim)))//"_"//&
                  &trim(adjustl(string(wrtid)))//"_NB"//&
                  &trim(adjustl(string(i_body+1)))//"_simT0.dat"
                flT0 = trim(adjustl(flT0))
                open (uT0, file=trim(flT0))

                if (durcheck .eq. 0) then ! do not check duration
                    fmt = adjustl(trim(fmt))//"8("//trim(sprec)//",1x),i3"
                    hea = adjustl(trim(hea))//"9(1x,a)"
                    write (uT0, "("//trim(hea)//")")&
                        &"# epoT0obs",&
                        &"T0obs", "eT0obs", "T0_sim", "T0obs-T0_sim",&
                        &"T0_stat",&
                        &"period_d", "sma_au", "ecc", "inc_deg",&
                        &"meana_deg", "argp_deg", "truea_deg", "longn_deg",&
                        &"source_id"
                    do i_T0 = 1, nT0
                        write (uT0, "("//trim(fmt)//")")&
                            &obsData%obsT0(i_body)%epo(i_T0),&
                            &obsData%obsT0(i_body)%T0(i_T0), obsData%obsT0(i_body)%eT0(i_T0),&
                            &simT0(i_body)%T0(i_T0),&
                            &(obsData%obsT0(i_body)%T0(i_T0)-simT0(i_body)%T0(i_T0)),&
                            &simT0(i_body)%T0_stat(i_T0),&
                            &simT0(i_body)%period(i_T0), simT0(i_body)%sma(i_T0), simT0(i_body)%ecc(i_T0),&
                            &simT0(i_body)%inc(i_T0), simT0(i_body)%meana(i_T0), simT0(i_body)%argp(i_T0),&
                            &simT0(i_body)%truea(i_T0), simT0(i_body)%longn(i_T0),&
                            &obsData%obsT0(i_body)%source_id(i_T0)
                    end do

                else ! do check duration
                    fmt = adjustl(trim(fmt))//"12("//trim(sprec)//",1x),2(i3,1x)"
                    hea = adjustl(trim(hea))//"14(1x,a)"
                    write (uT0, "("//trim(hea)//")")&
                       &"# epoT0obs",&
                       &"T0obs", "eT0obs", "T0_sim", "T0obs-T0_sim",&
                       &"T0_stat",&
                       &"period_d", "sma_au", "ecc", "inc_deg",&
                       &"meana_deg", "argp_deg", "truea_deg", "longn_deg",&
                       &"Dur_obs", "eDur_obs", "Dur_sim", "Dur_obs-Dur_sim",&
                       &"Dur_stat ",&
                       &"source_id"
                    do i_T0 = 1, nT0
                        write (uT0, "("//trim(fmt)//",4("//trim(sprec)//",1x),i3)")&
                          &obsData%obsT0(i_body)%epo(i_T0),&
                          &obsData%obsT0(i_body)%T0(i_T0), obsData%obsT0(i_body)%eT0(i_T0),&
                          &simT0(i_body)%T0(i_T0),&
                          &(obsData%obsT0(i_body)%T0(i_T0)-simT0(i_body)%T0(i_T0)),&
                          &simT0(i_body)%T0_stat(i_T0),&
                          &simT0(i_body)%period(i_T0), simT0(i_body)%sma(i_T0), simT0(i_body)%ecc(i_T0),&
                          &simT0(i_body)%inc(i_T0), simT0(i_body)%meana(i_T0), simT0(i_body)%argp(i_T0),&
                          &simT0(i_body)%truea(i_T0), simT0(i_body)%longn(i_T0),&
                          &obsData%obsT0(i_body)%dur(i_T0), obsData%obsT0(i_body)%edur(i_T0),&
                          &simT0(i_body)%dur(i_T0), obsData%obsT0(i_body)%dur(i_T0)-simT0(i_body)%dur(i_T0),&
                          &simT0(i_body)%dur_stat(i_T0),&
                          &obsData%obsT0(i_body)%source_id(i_T0)
                    end do

                end if

                flush (uT0)
                close (uT0)
                write (*, '(a,a)') " WRITTEN T0 INTO FILE: ", trim(flT0)
!         flush(6)

            end if
        end do

        return
    end subroutine write_T0_f
    ! ------------------------------------------------------------------ !

    ! ---
    ! adapted subroutine to write only final parameters
    subroutine write_parameters(cpuid, isim, wrtid, par, resw, to_screen)
        ! Input
        integer, intent(in)::cpuid, isim, wrtid
        real(dp), dimension(:), intent(in)::par
        real(dp), dimension(:), intent(in)::resw
        logical, optional, intent(in)::to_screen
        ! Local parameters
        real(dp)::chi_square, reduced_chi_square, lnLikelihood, ln_const, lnprior, lnprob, bic
        integer::upar, j
        character(512)::flpar
        character(80)::fmt
        ! real(dp)::fitness,fitxdof,bic,chi2,chi2r,lnL,ln_const
        ! integer::ns,ne

        chi_square = sum(resw*resw)
        call set_fitness_values(par, chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)
        lnprior = zero
        if (n_priors .gt. 0) then
            call ln_priors(system_parameters, par,&
                &priors_names, priors_values, lnprior)
        end if
        lnprob = lnLikelihood+lnprior
        bic = bic-two*lnprior

        upar = get_unit(cpuid)
        flpar = ""
        flpar = trim(path)//trim(adjustl(string(isim)))//"_"//&
          &trim(adjustl(string(wrtid)))//"_final"//&
          &trim(adjustl(string(nfit)))//"par.dat"
        flpar = trim(adjustl(flpar))
        open (upar, file=trim(flpar))

        write (upar, '(a,i3,a)') "# CPU ", cpuid, " PARAMETERS "

        if (present(to_screen)) then

            fmt = adjustl("(a,i5)")
            write (*, trim(fmt)) "# ndata   = ", obsData%ndata
            write (*, trim(fmt)) "# nfit    = ", nfit
            write (*, trim(fmt)) "# nfree   = ", obsData%nfree
            write (*, trim(fmt)) "# nRVset  = ", nRVset
            write (*, trim(fmt)) "# nRV     = ", obsData%obsRV%nRV
            write (*, trim(fmt)) "# nT0     = ", obsData%nTTs
            write (*, trim(fmt)) "# nT14    = ", obsData%nDurs

            fmt = adjustl("(2(a,"//trim(sprec)//"),a,i6)")
            write (*, trim(fmt)) "# Chi Square = ",&
            &chi_square, " => Reduced Chi Square = ", reduced_chi_square,&
            &" dof = ", obsData%dof

            fmt = adjustl("(a,"//trim(sprec)//",a)")
            write (*, trim(fmt)) "# ln_const(err, jitter) = ", ln_const, " ( = -ndata/2 *ln(2pi) - sum( ln(err^2 + jitter^2)) )"
            write (*, trim(fmt)) "# lnLikelihood = ", lnLikelihood
            write (*, trim(fmt)) "# lnPrior      = ", lnprior
            write (*, trim(fmt)) "# lnProb = lnLikelihood + lnPrior = ", lnprob
            write (*, trim(fmt)) "# BIC = -2xlnProb + nfit x ln(ndata) = ", bic
            write (*, '(a,5x,a)') "# parameter ", " value "
        end if

        fmt = adjustl("(a,i5)")
        write (upar, trim(fmt)) "# ndata   = ", obsData%ndata
        write (upar, trim(fmt)) "# nfit    = ", nfit
        write (upar, trim(fmt)) "# nfree   = ", obsData%nfree
        write (upar, trim(fmt)) "# nRVset  = ", nRVset
        write (upar, trim(fmt)) "# nRV     = ", obsData%obsRV%nRV
        write (upar, trim(fmt)) "# nT0     = ", obsData%nTTs
        write (upar, trim(fmt)) "# nT14    = ", obsData%nDurs

        fmt = adjustl("(2(a,"//trim(sprec)//"),a,i6)")
        write (upar, trim(fmt)) "# Chi Square = ",&
        &chi_square, " => Reduced Chi Square = ", reduced_chi_square,&
        &" dof = ", obsData%dof

        fmt = adjustl("(a,"//trim(sprec)//",a)")
        write (upar, trim(fmt)) "# ln_const(err, jitter) = ", ln_const, " ( = -ndata/2 *ln(2pi) - sum( ln(err^2 + jitter^2)) )"
        write (upar, trim(fmt)) "# lnLikelihood = ", lnLikelihood
        write (upar, trim(fmt)) "# lnPrior      = ", lnprior
        write (upar, trim(fmt)) "# lnProb = lnLikelihood + lnPrior = ", lnprob
        write (upar, trim(fmt)) "# BIC = -2xlnProb + nfit x ln(ndata) = ", bic
        write (upar, '(a,5x,a)') "# parameter ", " value "

        fmt = adjustl("(a10,4x,"//trim(sprec)//")")
        do j = 1, nfit
            if (present(to_screen)) write (*, trim(fmt)) trim(parid(j)), par(j)
            write (upar, trim(fmt)) trim(parid(j)), par(j)
        end do
        flush (6)
        flush (upar)
        close (upar)

        return
    end subroutine write_parameters

end module output_files
