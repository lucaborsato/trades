module utils
    use constants
    use parameters
    use custom_type

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
        ln_const = -(half*real(obsData%dof, dp)*log(dpi)) - (half*(ln_eRV + ln_eT0 + ln_edur))

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

            do ipl = 1, NB - 1
                if (obsData_in%obsT0(ipl)%nT0 .gt. 0) then
                    ln_eT0 = ln_eT0 + sum(log(pack(&
                      &obsData_in%obsT0(ipl)%eT0, obsData_in%obsT0(ipl)%eT0 /= zero)*&
                      &pack(obsData_in%obsT0(ipl)%eT0, obsData_in%obsT0(ipl)%eT0 /= zero)&
                      &))
                    if (durcheck .eq. 1) ln_edur = ln_edur + sum(log(&
                      &pack(obsData_in%obsT0(ipl)%edur, obsData_in%obsT0(ipl)%edur /= zero)*&
                      &pack(obsData_in%obsT0(ipl)%edur, obsData_in%obsT0(ipl)%edur /= zero)&
                      &))
                end if

            end do

        end if

        ln_const = -(half*real(obsData_in%ndata, dp)*log(dpi))&
          &- (half*(ln_eRV + ln_eT0 + ln_edur))

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
                aRV = bRV + 1
                bRV = bRV + obsData_in%obsRV%nRVsingle(iset)
                ln_eRV = ln_eRV + sum(log(obsData_in%obsRV%eRV(aRV:bRV)**2 + jitter(iset)**2))
            end do

        end if

        if (obsData_in%nTTs .gt. 0) then

            do ipl = 1, NB - 1
                if (obsData_in%obsT0(ipl)%nT0 .gt. 0) then
                    ln_eT0 = ln_eT0 + sum(log(pack(&
                      &obsData_in%obsT0(ipl)%eT0, obsData_in%obsT0(ipl)%eT0 /= zero)*&
                      &pack(obsData_in%obsT0(ipl)%eT0, obsData_in%obsT0(ipl)%eT0 /= zero)&
                      &))
                    if (durcheck .eq. 1) ln_edur = ln_edur + sum(log(&
                      &pack(obsData_in%obsT0(ipl)%edur, obsData_in%obsT0(ipl)%edur /= zero)*&
                      &pack(obsData_in%obsT0(ipl)%edur, obsData_in%obsT0(ipl)%edur /= zero)&
                      &))
                end if

            end do

        end if

        ln_const = -(half*real(obsData_in%ndata, dp)*log(dpi))&
          &- (half*(ln_eRV + ln_eT0 + ln_edur))

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
                xRV = obsRV%RV(iRV) - (simRV%RV(iRV) + simRV%gamma_rv(iRV) + simRV%trend(iRV))
                RVset = simRV%RVsetID(iRV)
                resw(iRV) = xRV/sqrt(obsRV%eRV(iRV)**2 + simRV%jitter(RVset)**2)
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
            j1 = j - 1
            nT0 = obsT0(j1)%nT0
            if (nT0 .gt. 0) then
                a = a + 1
                b = b + nT0
                resw(a:b) = (obsT0(j1)%T0 - simT0(j1)%T0)/obsT0(j1)%eT0
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
            j1 = j - 1
            nd = obsT0(j1)%nDur
            if (nd .gt. 0) then
                a = a + 1
                b = b + nd
                resw(a:b) = (obsT0(j1)%dur - simT0(j1)%dur)/obsT0(j1)%edur
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
            nTx = obsT0(i_body - 1)%nT0

            if (nTx .gt. 0) then
                a = a + 1
                b = b + nTx
                ! it computes the linear ephemeris from simulated data
                ! call set_ephem(simT0(i_body-1))
                call set_ephem_simT0(simT0(i_body - 1))
                call compute_oc_one_planet(simT0(i_body - 1))
                resw(a:b) = (obsT0(i_body - 1)%oc - simT0(i_body - 1)%oc)/obsT0(i_body - 1)%eT0
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
            ns = nkel + 1
            ne = nkel + nRVset
            allocate (jitter(nRVset))
            jitter = two**fit_parameters(ns:ne)
            ln_const = get_lnec_full(obsData, jitter)
            deallocate (jitter)
        else
            ln_const = get_lnec(obsData)
        end if

        lnLikelihood = -half*chi_square + ln_const

        ! bic = -two*lnLikelihood + real(nfit,dp)*log(real(ndata,dp))
        bic = -two*lnLikelihood + bic_const ! bic_const global variable

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
            val(nRV + 1:nRV + nTTs) = val_oc
            deallocate (val_oc)

            if (durcheck .eq. 1) then
                allocate (val_dur(nDurs))
                val_dur = zero
                call set_dur_resw(oDataIn%obsT0, simT0, val_dur)
                val(nRV + nTTs + 1:nRV + nTTs + nDurs) = val_dur
                deallocate (val_dur)
            end if

        else if (ocfit .eq. 2) then
            ! fit T0 and O-C and weight it half

            allocate (val_T0(nTTs), val_oc(nTTs))
            val_T0 = zero
            call set_T0_resw(oDataIn%obsT0, simT0, val_T0)
            val_oc = zero
            call set_oc_resw(oDataIn%obsT0, simT0, val_oc)
            val(nRV + 1:nRV + nTTs) = sqrt_half*sqrt(val_T0*val_T0 + val_oc*val_oc)
            deallocate (val_T0, val_oc)

            if (durcheck .eq. 1) then
                allocate (val_dur(nDurs))
                val_dur = zero
                call set_dur_resw(oDataIn%obsT0, simT0, val_dur)
                val(nRV + nTTs + 1:nRV + nTTs + nDurs) = val_dur
                deallocate (val_dur)
            end if

        else
            ! fit only T0, default way

            allocate (val_T0(nTTs))
            val_T0 = zero
            ! write(*,*)" DEBUG: || set_T0_resw"
            call set_T0_resw(oDataIn%obsT0, simT0, val_T0)
            ! write(*,*)" DEBUG: || sum(val_T0*val_T0) = ",sum(val_T0*val_T0)
            val(nRV + 1:nRV + nTTs) = val_T0
            ! write(*,*)" DEBUG: || sum(val*val) = ",sum(val*val)
            deallocate (val_T0)

            if (durcheck .eq. 1) then
                allocate (val_dur(nDurs))
                val_dur = zero
                call set_dur_resw(oDataIn%obsT0, simT0, val_dur)
                val(nRV + nTTs + 1:nRV + nTTs + nDurs) = val_dur
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

end module utils
