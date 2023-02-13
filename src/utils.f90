module utils
    use constants
    use parameters
    use custom_type
    use sorting, only: indexx

contains

    function get_ln_err_const(eRV, eT0, edur) result(ln_const)
        real(dp)::ln_const
        real(dp), dimension(:), intent(in)::eRV
        real(dp), dimension(:, :), intent(in)::eT0, edur

        real(dp)::ln_eRV, ln_eT0, ln_edur
        integer::nRV, nTTs

        nRV = size(eRV)
        nTTs = obsData%nTTs

        ln_eRV = zero
        ln_eT0 = zero
        ln_edur = zero

        if (nRV .gt. 0) then
            ln_eRV = sum(log(pack(eRV, eRV /= zero)*pack(eRV, eRV /= zero)))
        else
            ln_eRV = zero
        end if
        if (nTTs .gt. 0) then
            ln_eT0 = sum(log(pack(eT0, eT0 /= zero)*pack(eT0, eT0 /= zero)))
            if (durcheck .eq. 1) then
                ln_edur = sum(log(pack(edur, edur /= zero)*pack(edur, edur /= zero)))
            else
                ln_edur = zero
            end if
        else
            ln_eT0 = zero
        end if
        ln_const = -(half*real(obsData%dof, dp)*log(dpi))-(half*(ln_eRV+ln_eT0+ln_edur))

    end function get_ln_err_const

    ! ln errors constant with derived data type
    function get_lnec(obsData_in) result(ln_const)
        real(dp)::ln_const
        type(dataObs), intent(in)::obsData_in

        real(dp)::ln_eRV, ln_eT0, ln_edur
        integer::ipl

        ln_eRV = zero
        ln_eT0 = zero
        ln_edur = zero

        if (obsData_in%obsRV%nRV .gt. 0) then

            ln_eRV = sum(log(&
              &pack(obsData_in%obsRV%eRV, obsData_in%obsRV%eRV /= zero)*&
              &pack(obsData_in%obsRV%eRV, obsData_in%obsRV%eRV /= zero)&
              &))

        end if

        if (obsData_in%nTTs .gt. 0) then

            do ipl = 1, NB-1
                if (obsData_in%obsT0(ipl)%nT0 .gt. 0) then
                    ln_eT0 = ln_eT0+sum(log(pack(&
                      &obsData_in%obsT0(ipl)%eT0, obsData_in%obsT0(ipl)%eT0 /= zero)*&
                      &pack(obsData_in%obsT0(ipl)%eT0, obsData_in%obsT0(ipl)%eT0 /= zero)&
                      &))
                    if (durcheck .eq. 1) ln_edur = ln_edur+sum(log(&
                      &pack(obsData_in%obsT0(ipl)%edur, obsData_in%obsT0(ipl)%edur /= zero)*&
                      &pack(obsData_in%obsT0(ipl)%edur, obsData_in%obsT0(ipl)%edur /= zero)&
                      &))
                end if

            end do

        end if

        ln_const = -(half*real(obsData_in%ndata, dp)*log(dpi))&
          &-(half*(ln_eRV+ln_eT0+ln_edur))

    end function get_lnec

    ! ln errors constant with derived data type and taking into account RV jitter
    function get_lnec_full(obsData_in, jitter) result(ln_const)
        real(dp)::ln_const
        type(dataObs), intent(in)::obsData_in
        real(dp), dimension(:), intent(in)::jitter

        real(dp)::ln_eRV, ln_eT0, ln_edur
        integer::nset, iset, aRV, bRV, ipl

        ln_eRV = zero
        ln_eT0 = zero
        ln_edur = zero
        ln_const = zero

        if (obsData_in%obsRV%nRV .gt. 0) then

            nset = nRVset !obsData_in%obsRV%nRVset
            aRV = 0
            bRV = 0
            do iset = 1, nset
                aRV = bRV+1
                bRV = bRV+obsData_in%obsRV%nRVsingle(iset)
                ln_eRV = ln_eRV+sum(log(obsData_in%obsRV%eRV(aRV:bRV)**2+jitter(iset)**2))
            end do

        end if

        if (obsData_in%nTTs .gt. 0) then

            do ipl = 1, NB-1
                if (obsData_in%obsT0(ipl)%nT0 .gt. 0) then
                    ln_eT0 = ln_eT0+sum(log(pack(&
                      &obsData_in%obsT0(ipl)%eT0, obsData_in%obsT0(ipl)%eT0 /= zero)*&
                      &pack(obsData_in%obsT0(ipl)%eT0, obsData_in%obsT0(ipl)%eT0 /= zero)&
                      &))
                    if (durcheck .eq. 1) ln_edur = ln_edur+sum(log(&
                      &pack(obsData_in%obsT0(ipl)%edur, obsData_in%obsT0(ipl)%edur /= zero)*&
                      &pack(obsData_in%obsT0(ipl)%edur, obsData_in%obsT0(ipl)%edur /= zero)&
                      &))
                end if

            end do

        end if

        ln_const = -(half*real(obsData_in%ndata, dp)*log(dpi))&
          &-(half*(ln_eRV+ln_eT0+ln_edur))

    end function get_lnec_full

    ! set only RV to the weighted residuals
!   subroutine set_RV_resw(RV_obs,RV_sim,e_RVobs,gamma,resw)
    subroutine set_RV_resw(obsRV, simRV, resw)
        use statistics, only: wmean
!     real(dp),dimension(:),intent(in)::RV_obs,RV_sim,e_RVobs
!     real(dp),dimension(:,:),allocatable,intent(out)::gamma
        type(dataRV), intent(in)::obsRV
        type(dataRV), intent(inout)::simRV

        real(dp), dimension(:), intent(inout)::resw

        integer::iRV, nRV, RVset
        real(dp)::xRV

        nRV = obsRV%nRV
        if (nRV .gt. 0) then
            do iRV = 1, nRV
                xRV = obsRV%RV(iRV)-(simRV%RV(iRV)+simRV%gamma_rv(iRV)+simRV%trend(iRV))
                RVset = simRV%RVsetID(iRV)
                resw(iRV) = xRV/sqrt(obsRV%eRV(iRV)**2+simRV%jitter(RVset)**2)
                ! write(*,*)iRV, RVset, obsRV%RV(iRV),&
                !   &simRV%RV(iRV), simRV%gamma_rv(iRV), simRV%trend(iRV),&
                !   &simRV%jitter(RVset),resw(iRV)
            end do
        end if

        return
    end subroutine set_RV_resw

    ! set only T0 to the weighted residuals
    subroutine set_T0_resw(obsT0, simT0, resw)
!     real(dp),dimension(:,:),intent(in)::T0_obs,T0_sim,eT0_obs
        type(dataT0), dimension(:), intent(in)::obsT0
        type(dataT0), dimension(:), intent(inout)::simT0
        real(dp), dimension(:), intent(out)::resw

        integer::j, j1, a, b, nT0

        ! NB is a common variable
        resw = zero

        a = 0
        b = 0
        do j = 2, NB
            j1 = j-1
            nT0 = obsT0(j1)%nT0
            if (nT0 .gt. 0) then
                a = a+1
                b = b+nT0
                resw(a:b) = (obsT0(j1)%T0-simT0(j1)%T0)/obsT0(j1)%eT0
                a = b
            end if
        end do

        return
    end subroutine set_T0_resw

    ! set only dur to the weighted residuals
    subroutine set_dur_resw(obsT0, simT0, resw)
!     real(dp),dimension(:,:),intent(in)::T0_obs,T0_sim,eT0_obs
        type(dataT0), dimension(:), intent(in)::obsT0
        type(dataT0), dimension(:), intent(inout)::simT0
        real(dp), dimension(:), intent(out)::resw

        integer::j, j1, a, b, nd

        ! NB is a common variable
        resw = zero

        a = 0
        b = 0
        do j = 2, NB
            j1 = j-1
            nd = obsT0(j1)%nDur
            if (nd .gt. 0) then
                a = a+1
                b = b+nd
                resw(a:b) = (obsT0(j1)%dur-simT0(j1)%dur)/obsT0(j1)%edur
                a = b
            end if
        end do

        return
    end subroutine set_dur_resw

    ! compute linear ephemeris for both obs and sim
    ! use O-Cobs/sim to compute residuals
!   subroutine set_oc_resw(epoT0_obs,T0_obs,eT0_obs,T0_sim,resw)
    subroutine set_oc_resw(obsT0, simT0, resw)
        use linear_ephem
        type(dataT0), dimension(:), intent(in)::obsT0
        type(dataT0), dimension(:), intent(inout)::simT0
        real(dp), dimension(:), intent(out)::resw

        integer::i_body, nTx, a, b

        ! NB is a common variables
        resw = zero

        a = 0
        b = 0
        do i_body = 2, NB
            nTx = obsT0(i_body-1)%nT0

            if (nTx .gt. 0) then
                a = a+1
                b = b+nTx
                ! it computes the linear ephemeris from simulated data
                ! call set_ephem(simT0(i_body-1))
                call set_ephem_simT0(simT0(i_body-1))
                call compute_oc_one_planet(simT0(i_body-1))
                resw(a:b) = (obsT0(i_body-1)%oc-simT0(i_body-1)%oc)/obsT0(i_body-1)%eT0
                a = b
            end if
        end do

        return
    end subroutine set_oc_resw

! ===============================================================================================
    subroutine set_fitness_values(fit_parameters, chi_square, &
      &reduced_chi_square, lnLikelihood, ln_const, bic)
        ! Input
        real(dp), intent(in), dimension(:)::fit_parameters
        ! Input/Output
        real(dp), intent(inout)::chi_square
        ! Output
        real(dp), intent(out)::reduced_chi_square, lnLikelihood, ln_const, bic

        ! Local variables
        real(dp), dimension(:), allocatable::jitter
        real(dp)::inv_dof

        if (chi_square .ge. resmax) then
            chi_square = resmax
        end if

        inv_dof = obsData%inv_dof !one/real(obsData%dof, dp)
        reduced_chi_square = chi_square*inv_dof

        ln_const = zero
        if (nRVset .gt. 0) then
            ns = nkel+1
            ne = nkel+nRVset
            allocate (jitter(nRVset))
            jitter = two**fit_parameters(ns:ne)
            ln_const = get_lnec_full(obsData, jitter)
            deallocate (jitter)
        else
            ln_const = get_lnec(obsData)
        end if

        lnLikelihood = -half*chi_square+ln_const

        ! bic = -two*lnLikelihood + real(nfit,dp)*log(real(ndata,dp))
        bic = -two*lnLikelihood+bic_const ! bic_const global variable

        return
    end subroutine set_fitness_values

    ! setting properly the weighted residuals: RV and T0 or OC
    !   subroutine set_weighted_residuals(RV_obs,RV_sim,e_RVobs,gamma,epoT0_obs,T0_obs,eT0_obs,T0_sim,resw,ocfit)
    subroutine set_weighted_residuals(oDataIn, simRV, simT0, resw, ocfit)
        type(dataObs), intent(in)::oDataIn
!     type(dataRV),intent(in)::simRV
!     type(dataT0),dimension(:),intent(in)::simT0
        type(dataRV), intent(inout)::simRV
        type(dataT0), dimension(:), intent(inout)::simT0
        real(dp), dimension(:), intent(out)::resw
        integer, intent(in)::ocfit

        real(dp), dimension(:), allocatable::val, val_T0, val_dur, val_oc
        integer::nRV, nTTs, nDurs

        nRV = oDataIn%obsRV%nRV
        nTTs = oDataIn%nTTs
        nDurs = oDataIn%nDurs
        allocate (val(oDataIn%ndata))

        resw = zero
        val = zero

        ! write(*,*)" DEBUG: || in set_weighted_residuals"
        ! write(*,*)" DEBUG: || sum(resw*resw) = ",sum(resw*resw)
        ! write(*,*)" DEBUG: || nRV = ",nRV
        ! write(*,*)" DEBUG: || set_RV_resw"
        if (nRV .gt. 0) call set_RV_resw(oDataIn%obsRV, simRV, val(1:nRV))
        ! write(*,*)" DEBUG: || sum(val*val) = ",sum(val*val)

        if (ocfit .eq. 1) then
            ! fit only O-C

            allocate (val_oc(nTTs))
            val_oc = zero
            call set_oc_resw(oDataIn%obsT0, simT0, val_oc)
            val(nRV+1:nRV+nTTs) = val_oc
            deallocate (val_oc)

            if (durcheck .eq. 1) then
                allocate (val_dur(nDurs))
                val_dur = zero
                call set_dur_resw(oDataIn%obsT0, simT0, val_dur)
                val(nRV+nTTs+1:nRV+nTTs+nDurs) = val_dur
                deallocate (val_dur)
            end if

        else if (ocfit .eq. 2) then
            ! fit T0 and O-C and weight it half

            allocate (val_T0(nTTs), val_oc(nTTs))
            val_T0 = zero
            call set_T0_resw(oDataIn%obsT0, simT0, val_T0)
            val_oc = zero
            call set_oc_resw(oDataIn%obsT0, simT0, val_oc)
            val(nRV+1:nRV+nTTs) = sqrt_half*sqrt(val_T0*val_T0+val_oc*val_oc)
            deallocate (val_T0, val_oc)

            if (durcheck .eq. 1) then
                allocate (val_dur(nDurs))
                val_dur = zero
                call set_dur_resw(oDataIn%obsT0, simT0, val_dur)
                val(nRV+nTTs+1:nRV+nTTs+nDurs) = val_dur
                deallocate (val_dur)
            end if

        else
            ! fit only T0, default way

            allocate (val_T0(nTTs))
            val_T0 = zero
            ! write(*,*)" DEBUG: || set_T0_resw"
            call set_T0_resw(oDataIn%obsT0, simT0, val_T0)
            ! write(*,*)" DEBUG: || sum(val_T0*val_T0) = ",sum(val_T0*val_T0)
            val(nRV+1:nRV+nTTs) = val_T0
            ! write(*,*)" DEBUG: || sum(val*val) = ",sum(val*val)
            deallocate (val_T0)

            if (durcheck .eq. 1) then
                allocate (val_dur(nDurs))
                val_dur = zero
                call set_dur_resw(oDataIn%obsT0, simT0, val_dur)
                val(nRV+nTTs+1:nRV+nTTs+nDurs) = val_dur
                deallocate (val_dur)
            end if

        end if
        ! call set_fitness(obsData,val,resw)
        resw = val
        ! write(*,*)" DEBUG: || sum(resw*resw) = ",sum(resw*resw)
        deallocate (val)

        return
    end subroutine set_weighted_residuals

    ! ================================================================================
    function set_max_residuals(ndata) result(resw_max)
        integer, intent(in)::ndata
        real(dp)::resw_max

        resw_max = sqrt(resmax/real(ndata, dp))

        return
    end function set_max_residuals

    subroutine check_max_residuals(resw, ndata)
        real(dp), dimension(:), intent(inout)::resw
        integer, intent(in)::ndata
        real(dp)::resw_max, resw_max_possible

        if (ndata .gt. 0) then
            resw_max = maxval(resw)
            resw_max_possible = set_max_residuals(ndata)
            if (resw_max .gt. resw_max_possible) resw = resw_max_possible
        end if

        return
    end subroutine check_max_residuals

    subroutine linrange(start_val, end_val, step, vector, nvec)
        ! **Input**
        real(dp), intent(in)::start_val, end_val, step
        ! **Output**
        real(dp), dimension(:), allocatable, intent(out)::vector
        integer, intent(out)::nvec
        ! **local variables**
        real(dp)::delta
        integer::n0, n1, ivec

        n0 = int((end_val-start_val)/step)+1
        n1 = n0+1

        delta = end_val-step*real(n0-1, dp)
        if (abs(delta) .le. TOL_dp) then
            nvec = n0
        else
            nvec = n1
        end if
        allocate (vector(nvec))
        ! do ivec=1,nvec
        !     vector(ivec) = start_val + step*ivec
        !     if(ivec.eq.nvec)then
        !         vector(ivec) = end_val
        !     end if
        ! end do
        do ivec = 1, nvec-1
            vector(ivec) = start_val+step*(ivec-1)
        end do
        vector(nvec) = end_val

        return
    end subroutine linrange

        ! subroutine to prepare steps of the integration, for each of this step
    ! a check of RV and T0 will be performed.
    subroutine set_check_steps(t_epoch, time_to_int, step, observed_data, n_all, all_steps, all_idx)
        ! **Input**
        real(dp), intent(in)::t_epoch, time_to_int, step
        type(dataObs), intent(in)::observed_data
        ! **Output**
        integer, intent(out)::n_all
        real(dp), dimension(:), allocatable, intent(out)::all_steps
        integer, dimension(:), allocatable, intent(out)::all_idx
        ! **local variables**
        integer::nRV
        real(dp), dimension(:), allocatable::t_check_steps, t_rv, sel_t_rv
        integer, dimension(:), allocatable::check_idx, rv_idx, sel_rv_idx, sort_idx
        logical, dimension(:), allocatable::sel_rv
        integer::i_rv, n_checks, i_all

        nRV = observed_data%obsRV%nRV

        ! prepare vector with steps to check, based only on time and step
        call linrange(zero, time_to_int, step, t_check_steps, n_checks)
        allocate (check_idx(n_checks))
        check_idx = 0 ! all zero, not needed to check RV
        if (nRV .gt. 0) then
            ! prepare RV times and index, based on the sign of time
            allocate (t_rv(nRV), rv_idx(nRV), sel_rv(nRV))
            t_rv = obsData%obsRV%jd-t_epoch
            rv_idx = (/(i_rv, i_rv=1, nRV)/)
            if (time_to_int .lt. zero) then
                sel_rv = t_rv .lt. zero
            else
                sel_rv = t_rv .ge. zero
            end if
            sel_t_rv = pack(t_rv, sel_rv)
            sel_rv_idx = pack(rv_idx, sel_rv)
            deallocate (t_rv, rv_idx, sel_rv)
            ! concatenate check and rv steps
            all_steps = [t_check_steps, sel_t_rv]
            all_idx = [check_idx, sel_rv_idx]
            n_all = size(all_idx)
            deallocate (sel_rv_idx)
        else
            n_all = n_checks
            allocate (all_steps(n_all))
            all_steps = t_check_steps
            all_idx = check_idx
        end if

        deallocate (t_check_steps, check_idx)
        ! prepare sorting vector
        allocate (sort_idx(n_all))
        call indexx(all_steps, sort_idx) ! get the sort in ascending order!
        ! needed to reverse the order if time_to_int < 0
        if (time_to_int .lt. zero) then
            sort_idx = sort_idx((/(i_all, i_all=n_all, 1, -1)/))
        end if
        ! sort the vectors with all the steps and idx
        all_steps = all_steps(sort_idx)
        all_idx = all_idx(sort_idx) ! 0 if no RV, > 1 if RV to check, idx of RV position in obsData%obsRV

        return
    end subroutine set_check_steps

    subroutine set_transiting_bodies(transiting_body, do_transit_check, body_transiting_start, body_transiting_end)
        !Input
        integer,intent(in)::transiting_body
        !Output
        logical,intent(out)::do_transit_check
        integer,intent(out)::body_transiting_start, body_transiting_end

        ! == NEW ==
        ! by default it will check the transit of all the planets (2->NB)
        do_transit_check = .true.
        body_transiting_start = 2
        body_transiting_end = NB
        if ((transiting_body .lt. 1) .or. (transiting_body .gt. NB)) then
            do_transit_check = .false.
        else
            if (transiting_body .gt. 1) then
                body_transiting_start = transiting_body
                body_transiting_end = transiting_body
            end if
        end if
        ! ========

        return
    end subroutine set_transiting_bodies

end module utils
