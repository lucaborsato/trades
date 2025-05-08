module custom_type
    use constants, only: dp, zero, one
    implicit none

    ! ==============================================================================

!DEC$ OPTIONS /ALIGN=(RECORDS=NATURAL)
    ! OBS DATA AS STRUCTURE
    type dataRV
        ! sequence ! needed to store variables contiguous
        ! 
        integer::nRVset
        ! nRV will be used as the number of RV datapoints for the observed one,
        ! and as a counter for the simulated one
        integer::nRV = 0
        real(dp), dimension(:), allocatable::jd, RV, eRV, trend
        integer, dimension(:), allocatable::RV_stat
        real(dp), dimension(:), allocatable::gamma, gamma_rv
        real(dp), dimension(:), allocatable::jitter
        integer, dimension(:), allocatable::RVsetID, nRVsingle
    end type dataRV
!DEC$ END OPTIONS

! ==============================================================================

    ! T0 data for one body!
    type dataT0
        ! sequence ! needed to store variables contiguous

        ! nT0 will be used as the number of TTs for the observed one,
        ! and as a counter for the simulated one
        integer::nT0 = 0, nDur = 0
        logical::transiting = .true.
        integer, dimension(:), allocatable::epo
        real(dp), dimension(:), allocatable::T0, eT0, dur, edur, oc, lambda_rm
        integer, dimension(:), allocatable::source_id
        integer, dimension(:), allocatable::T0_stat, dur_stat
        real(dp)::Tephem = zero, eTephem = zero, Pephem = zero, ePephem = zero
        real(dp), dimension(:), allocatable::period, sma, ecc, inc, meana
        real(dp), dimension(:), allocatable::argp, truea, longn, dttau
        ! adding excluded time ranges for a specific body
        integer::n_excl = 0
        real(dp), dimension(:,:), allocatable::excluded_time_ranges
        real(dp), dimension(:), allocatable::excluded_t0, excluded_t1, excluded_t4
        integer, dimension(:), allocatable::excluded_status
    end type dataT0
! ==============================================================================

    ! full data type
    type dataObs
        ! sequence ! needed to store variables contiguous
        integer::ndata, nfree, dof = 1
        real(dp)::inv_dof
        ! RV data
        type(dataRV)::obsRV
        ! T0 and duration data
        ! and as a counter for the simulated ones
        integer::nTTs = 0, nDurs = 0, n_excluded = 0 ! default set to zero
        type(dataT0), dimension(:), allocatable::obsT0
    end type dataObs
! ==============================================================================

    ! define new data type for grid search
    type parameter_grid
        sequence
        character(15)::name ! name/id of the parameter
        real(dp), dimension(3)::input_values = zero ! min, max, step as read from input file
        character(2)::step_type = 'rn' ! input step type as read from the input file: 'ss/rn/sn'
        real(dp)::step_grid = one ! calculated step size of the parameter
        integer::n_steps = 1 ! number of steps calculated, it take into account the min and max values
        real(dp), dimension(:), allocatable::grid_values ! values of the grid parameter with dimension n_steps
    end type parameter_grid
! ==============================================================================

contains

! ==============================================================================

    subroutine init_dataRV(nRV, RV)
        integer, intent(in)::nRV
        type(dataRV), intent(inout)::RV

        if (allocated(RV%jd)) then
            call deallocate_dataRV(RV)
        end if
        RV%nRV = nRV
        allocate (RV%jd(nRV), RV%RV(nRV), RV%eRV(nRV), RV%gamma_rv(nRV), RV%trend(nRV))
        allocate (RV%RV_stat(nRV), RV%RVsetID(nRV))
        RV%jd = zero
        RV%RV = zero
        RV%eRV = zero
        RV%gamma_rv = zero
        RV%trend = zero
        RV%RV_stat = 0
        RV%RVsetID = 0

    !!
        ! nRVset & RVsetID & nRVsingle & gamma & jitter HAVE TO BE SET BY 'HAND'
    !!

        return
    end subroutine init_dataRV

! ==============================================================================

    subroutine deallocate_dataRV(RV)
        type(dataRV), intent(inout)::RV

        if (allocated(RV%jd)) deallocate (RV%jd)
        if (allocated(RV%RV)) deallocate (RV%RV)
        if (allocated(RV%eRV)) deallocate (RV%eRV)
        if (allocated(RV%trend)) deallocate (RV%trend)
        if (allocated(RV%RV_stat)) deallocate (RV%RV_stat)
        if (allocated(RV%gamma)) deallocate (RV%gamma)
        if (allocated(RV%gamma_rv)) deallocate (RV%gamma_rv)
        if (allocated(RV%jitter)) deallocate (RV%jitter)
        if (allocated(RV%RVsetID)) deallocate (RV%RVsetID)
        if (allocated(RV%nRVsingle)) deallocate (RV%nRVsingle)
        RV%nRV = 0

        return
    end subroutine deallocate_dataRV

! ==============================================================================

    subroutine init_dataT0(nT0, T0, dur_check)
        integer, intent(in)::nT0
        type(dataT0), intent(inout)::T0
        integer, intent(in)::dur_check

        if (allocated(T0%T0)) then
            call deallocate_dataT0(T0)
        end if
        T0%nT0 = nT0
        allocate (T0%epo(nT0), T0%T0(nT0), T0%eT0(nT0), T0%oc(nT0), T0%lambda_rm(nT0), T0%T0_stat(nT0), T0%source_id(nT0))
        T0%epo = 0
        T0%T0 = zero
        T0%eT0 = zero
        T0%oc = zero
        T0%lambda_rm = zero
        T0%T0_stat = 0
        T0%source_id = 1

        ! duration
        allocate (T0%dur(nT0), T0%edur(nT0), T0%dur_stat(nT0))
        T0%dur = zero
        T0%edur = zero
        T0%dur_stat = 0
        if (dur_check .eq. 1) T0%nDur = nT0

        allocate (T0%period(nT0), T0%sma(nT0), T0%ecc(nT0), T0%inc(nT0), T0%meana(nT0))
        allocate (T0%argp(nT0), T0%truea(nT0), T0%longn(nT0), T0%dttau(nT0))
        ! init keplerian elements to fake values:
        T0%period = -one
        T0%sma = -one
        T0%ecc = -one
        T0%inc = 90.0_dp
        T0%meana = zero
        T0%argp = 90.0_dp
        T0%truea = zero
        T0%longn = 180.0_dp
        T0%dttau = zero

        return
    end subroutine init_dataT0

    subroutine init_dataT0_excluded(n_excl, excluded_ranges, T0)
        ! Input
        integer, intent(in)::n_excl
        real(dp), dimension(:,:), intent(in):: excluded_ranges
        ! Input/Output
        type(dataT0), intent(inout)::T0

        ! write(*,*)"init excluded"
        T0%n_excl = n_excl
        ! write(*,*)"T0%n_excl",T0%n_excl
        allocate (T0%excluded_time_ranges(n_excl, 2), T0%excluded_status(n_excl))
        allocate (T0%excluded_t0(n_excl), T0%excluded_t1(n_excl), T0%excluded_t4(n_excl))
        T0%excluded_time_ranges = excluded_ranges
        ! write(*,*)T0%excluded_time_ranges(:,1)
        ! write(*,*)T0%excluded_time_ranges(:,2)
        T0%excluded_status = 0
        T0%excluded_t0 = zero
        T0%excluded_t1 = zero
        T0%excluded_t4 = zero
        ! flush(6)

        return
    end subroutine init_dataT0_excluded

! ==============================================================================

    subroutine deallocate_dataT0(T0)
        type(dataT0), intent(inout)::T0

        if (allocated(T0%epo)) deallocate (T0%epo)
        if (allocated(T0%T0)) deallocate (T0%T0)
        if (allocated(T0%eT0)) deallocate (T0%eT0)
        if (allocated(T0%oc)) deallocate (T0%oc)
        if (allocated(T0%lambda_rm)) deallocate (T0%lambda_rm)
        if (allocated(T0%T0_stat)) deallocate (T0%T0_stat)
        if (allocated(T0%source_id)) deallocate (T0%source_id)
        if (allocated(T0%dur)) deallocate (T0%dur)
        if (allocated(T0%edur)) deallocate (T0%edur)
        if (allocated(T0%dur_stat)) deallocate (T0%dur_stat)
        if (allocated(T0%period)) then
            deallocate (T0%period, T0%sma, T0%ecc, T0%inc, T0%meana)
            deallocate (T0%argp, T0%truea, T0%longn, T0%dttau)
        end if
        if (allocated(T0%excluded_time_ranges)) deallocate (T0%excluded_time_ranges)
        if (allocated(T0%excluded_status)) deallocate (T0%excluded_status)
        if (allocated(T0%excluded_t0)) deallocate(T0%excluded_t0)
        if (allocated(T0%excluded_t1)) deallocate(T0%excluded_t1)
        if (allocated(T0%excluded_t4)) deallocate(T0%excluded_t4)
        T0%nT0 = 0
        T0%nDur = 0
        T0%Tephem = zero
        T0%eTephem = zero
        T0%Pephem = zero
        T0%ePephem = zero
        T0%n_excl = 0

        return
    end subroutine deallocate_dataT0

! ==============================================================================

    subroutine deallocate_dataObs(oData)
        type(dataObs), intent(inout)::oData

        integer::n, i

        if (allocated(oData%obsRV%jd)) call deallocate_dataRV(oData%obsRV)
        if (allocated(oData%obsT0)) then
            n = size(oData%obsT0)
            if (oData%ntts .gt. 0) then
                do i = 1, n
                    call deallocate_dataT0(oData%obsT0(i))
                end do
            end if
            deallocate (oData%obsT0)
        end if
        oData%nTTs = 0
        oData%nDurs = 0

    end subroutine deallocate_dataObs

! ==============================================================================

end module custom_type
