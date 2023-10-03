module radial_velocities
    use constants
    use custom_type
    use parameters
    use celestial_mechanics, only: barycenter
    implicit none

    interface check_RV
        module procedure check_RV_1, check_RV_2
    end interface check_RV

contains

    subroutine get_RV(mass, ri, rv_sim)
        real(dp), dimension(:), intent(in)::mass
        real(dp), dimension(:), intent(in)::ri
        real(dp), intent(out)::rv_sim

        real(dp), dimension(:), allocatable::rbarRV
        real(dp), dimension(:), allocatable::barRV

        allocate (barRV(6), rbarRV(NBDIM))
        call barycenter(mass, ri, barRV, rbarRV) ! astrocentric 2 barycentric
        rv_sim = -rbarRV(6)*AU/s24h ! rv as -Zstar,bar m/s
        deallocate (barRV, rbarRV)

        return
    end subroutine get_RV

    subroutine get_and_set_RV_data(mass, rin, idx, obsDataIn, simRV)
        ! **Input**
        real(dp), dimension(:), intent(in)::mass
        real(dp), dimension(:), intent(in)::rin
        integer, intent(in)::idx
        type(dataObs), intent(in)::obsDataIn
        ! **Input/Output**
        type(dataRV), intent(inout)::simRV
        ! **local variables**
        real(dp)::rv

        call get_RV(mass, rin, rv)
        simRV%jd(idx) = obsDataIn%obsRV%jd(idx)
        simRV%RV(idx) = rv
        simRV%nRV = simRV%nRV+1
        simRV%RVsetID(idx) = obsDataIn%obsRV%RVsetID(idx)
        simRV%RV_stat(idx) = 1

        return
    end subroutine get_and_set_RV_data

    subroutine get_and_set_RV(mass, rin, idx, simRV)
        ! **Input**
        real(dp), dimension(:), intent(in)::mass
        real(dp), dimension(:), intent(in)::rin
        integer, intent(in)::idx
        ! **Input/Output**
        type(dataRV), intent(inout)::simRV

        call get_and_set_RV_data(mass, rin, idx, obsData, simRV)

        return
    end subroutine get_and_set_RV

    ! ------------------------------------------------------------------ !
    ! uses the integrator to move to the obs RV jd and computed the RV (sim in m/s)
    subroutine calcRV(m, ri, drdt, stepRV, RVsim)
        use numerical_integrator, only: rkck_a
        ! use celestial_mechanics, only: barycenter
        real(dp), dimension(:), intent(in)::m
        real(dp), dimension(:), intent(in)::ri, drdt
        real(dp), intent(in)::stepRV
        real(dp), intent(out)::RVsim

        real(dp), dimension(:), allocatable::ro, err !, rbarRV
        ! real(dp), dimension(:), allocatable::barRV

        allocate (ro(NBDIM), err(NBDIM))
        call rkck_a(m, ri, drdt, stepRV, ro, err) ! call integrator
        call get_RV(m, ro, RVsim)

        return
    end subroutine calcRV

    ! calculation of the trend of RV, zeroth order == gamma, not taken into account
    subroutine addRVtrend(t, tref, coeff, trend)
        real(dp), dimension(:), intent(in)::t
        real(dp), intent(in)::tref
        real(dp), dimension(:), intent(in)::coeff
        real(dp), dimension(:), intent(out)::trend

        real(dp), dimension(size(t))::w
        integer::o, order

        w = t-tref ! tepoch is a global variable
        trend = zero
        order = size(coeff)
        do o = 1, order
            trend = trend+coeff(o)*w**o
        end do ! o -> order

        return
    end subroutine addRVtrend

    subroutine set_gamma_rv(gamma, rv_set_id, gamma_rv)
        real(dp), dimension(:), intent(in)::gamma
        integer, dimension(:), intent(in)::rv_set_id
        real(dp), dimension(:), intent(out)::gamma_rv
        integer::iset, iRV, nRV !,aRV,bRV

        ! write(*,*)" ====== set_gamma_rv ====="
        nRV = size(rv_set_id)
        ! write(*,*)"nRV = ",nRV
        do iRV = 1, nRV
            iset = rv_set_id(iRV)
            gamma_rv(iRV) = gamma(iset)
            ! write(*,*)iRV, iset, gamma_rv(iRV)
        end do
        ! write(*,*)" ========================="

        return
    end subroutine set_gamma_rv

    ! checks the RV, it uses the calcRV subroutine
    subroutine check_RV_1(m, ri, drdt, ttemp, hok, simRV)
        real(dp), dimension(:), intent(in)::m, ri, drdt
        real(dp), intent(in)::ttemp
        !     integer,intent(inout)::cntRV
        !     real(dp),dimension(:),intent(inout)::RV_sim
        !     integer,dimension(:),intent(inout)::RV_stat
        type(dataRV), intent(inout)::simRV
        real(dp), intent(inout)::hok

        real(dp)::stepRV, xRV
        integer::nRV, j

        nRV = obsData%obsRV%nRV
        xRV = zero
        ! RV check
        do j = 1, nRV
            if (simRV%RV_stat(j) .eq. 0) then
                stepRV = obsData%obsRV%jd(j)-tepoch-ttemp
                if (stepRV*(stepRV-hok) .le. zero) then
                    ! if ((stepRV .ge. zero) .and. (hok .ge. stepRV)) then
                    call calcRV(m, ri, drdt, stepRV, xRV)
                    simRV%RV(j) = xRV
                    simRV%jd(j) = obsData%obsRV%jd(j)
                    simRV%nRV = simRV%nRV+1
                    simRV%RVsetID(j) = obsData%obsRV%RVsetID(j)
                    simRV%RV_stat(j) = 1
                end if
            end if
        end do

        return
    end subroutine check_RV_1

!   subroutine check_RV_2(m,ri,drdt,ttemp,hok,cntRV,tRV,RV_stat,RV_sim)
    subroutine check_RV_2(m, ri, drdt, ttemp, hok, obsjd, simRV)
        real(dp), dimension(:), intent(in)::m, ri, drdt
        real(dp), intent(in)::ttemp
!     integer,intent(inout)::cntRV
!     real(dp),dimension(:),intent(in)::tRV
!     integer,dimension(:),intent(inout)::RV_stat
!     real(dp),dimension(:),intent(inout)::RV_sim
        real(dp), dimension(:), intent(in)::obsjd
        type(dataRV), intent(inout)::simRV
        real(dp), intent(inout)::hok
        real(dp)::stepRV
        integer::j, n_RV

        n_RV = size(obsjd)
        ! RV check
        do j = 1, n_RV
            if (simRV%RV_stat(j) .eq. 0) then
                stepRV = obsjd(j)-tepoch-ttemp
                if (stepRV*(stepRV-hok) .le. zero) then
                    ! if ((stepRV .ge. zero) .and. (hok .ge. stepRV)) then
                    call calcRV(m, ri, drdt, stepRV, simRV%RV(j))
                    simRV%jd(j) = obsjd(j)
                    simRV%nRV = simRV%nRV+1
                    simRV%RV_stat(j) = 1
                end if
            end if
        end do

        return
    end subroutine check_RV_2

end module radial_velocities
