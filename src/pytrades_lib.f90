! TRADES fortran module to be called fron python

module f90trades
!$  use omp_lib
    use constants
    use parameters
    use parameters_conversion
    use custom_type
    use convert_type, only: string
    use celestial_mechanics
    use init_trades
    use derived_parameters_mod
    use radial_velocities, only: addRVtrend, set_gamma_rv, get_RV
    use lin_fit, only: linfit
    use linear_ephem
    use fitness_module
    use utils
    use ode_run, only: integration_info_to_data, ode_all_ttra_rv,&
        &ode_keplerian_elements_to_orbits, set_checking_coordinates, transit_conditions
    use gls_module, only: check_periodogram
    use transits, only: compute_transit_time, spin_orbit_misalignment
    use grid_search

    implicit none

    ! exposing variables in parameters module to trades_lib
    !f2py integer,parameter::dp=selected_real_kind(8)
    !f2py character(512)::path

    !f2py integer::lmon
    !f2py integer::npar,nkel, nfit
    !f2py integer::nRVset
    !f2py integer::rv_trend_order
    !f2py real(dp),parameter::resmax
    !f2py real(dp)::tepoch,tstart,tint
    !f2py integer::idpert
    !f2py real(dp)::wrttime
    !f2py integer::durcheck
    !f2py logical::do_hill_check,amd_hill_check
    !f2py real(dp),dimension(2,2)::MR_star
    !f2py integer::n_priors
    !f2py integer::seed_pso,np_pso,nit_pso,wrt_pso
    !f2py integer::ncpu_in

    ! needed because TRADES now uses these variables in data type
    integer::ndata, nfree, dof
    integer::nRV=0 !,nRVset it is global now
    ! nT0, Pephem, Tephem returned by a subroutine now
    integer::nTTs=0, nDurs=0

!   real(dp)::ln_err_const,inv_dof
    real(dp)::inv_dof
    ! !f2py real(dp)::ln_err_const

    ! variables:  parameters to fit
    ! fitting_parameters and parameters_minmax returned by a subroutine now
    integer::str_len
    integer::n_global, n_bodies

contains

    ! ============================================================================

    ! SET ARGS AND PARAMETERS WITH SINGLE SUBROUTINES TO BE CALLED BY PYTHON IF NEEDED

    ! --- subroutine useful to modify the working path of TRADES from python
    subroutine path_change(new_path)
        character(512), intent(in)::new_path

        path = trim(adjustl(new_path))

        return
    end subroutine path_change

    subroutine get_path(path_out)
        character(512), intent(out)::path_out

        path_out = trim(adjustl(path))

        return
    end subroutine get_path
    ! ---

    ! --- set number of bodies
    subroutine set_n_bodies(n_body)
        integer, intent(in)::n_body

        NB = n_body
        NBDIM = NB*6

        return
    end subroutine set_n_bodies

    subroutine get_n_bodies(n_body)
        integer, intent(out)::n_body

        n_body = NB

        return
    end subroutine get_n_bodies
    ! ---

    ! --- set epcoh/reference time
    subroutine set_epoch_time(t_epoch)
        real(dp), intent(in)::t_epoch

        tepoch = t_epoch

        return
    end subroutine set_epoch_time

    ! --- set starting time
    subroutine set_starting_time(t_starting)
        real(dp), intent(in)::t_starting

        tstart = t_starting

        return
    end subroutine set_starting_time

    ! --- set integration time
    subroutine set_integration_time(t_integration)
        real(dp), intent(in)::t_integration

        tint = t_integration

        return
    end subroutine set_integration_time
    ! ============================================================================

    ! ============================================================================
    ! ============================================================================
    subroutine initialize_trades(path_in, sub_folder, n_threads_in)
        !f2py real(dp),dimension(:),allocatable::eRVobs
        !f2py real(dp),dimension(:,:),allocatable::eT0obs
        ! 
        ! character*(*), intent(in)::path_in, sub_folder
        character*(*), intent(in)::path_in, sub_folder
        integer, intent(in),optional::n_threads_in
        ! Local
        real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN
        integer, dimension(:), allocatable::nset

        integer:: ipar

        write (*, *) "INITIALISING TRADES ..."

        call initu(nfiles, 1)

        ! IT READS THE COMMAND ARGUMENT THAT DEFINE THE PATH OF THE FOLDER WITH THE FILES
        path_0 = trim(adjustl(path_in))
        path = trim(adjustl(path_in))

        ! IT DEFINES THE STRING TO WRITE THE REAL WITH RIGHT DECIMAL: PRECISION
        sprec = 'es23.16'

        ! IT READS THE ARGUMENTS OF INTEGRATION AND STORE IN COMMON MODULE PARAMETERS.
        ! THE VARIBLES WILL NOT BE MODIFIED FURTHERMORE.
        call read_arg(1)
        write(*,'(a, f16.5)')"t_epoch = ",tepoch
        write(*,'(a, f16.5)')"t_start = ",tstart
        write(*,'(a, f16.5)')"t_int   = ",tint
        flush(6)

        n_bodies = NB ! needed to be used by python wrapper ... to check if I can avoid it
        if (allocated(e_bounds)) deallocate (e_bounds)
        allocate (e_bounds(2, NB))
        e_bounds(1, :) = TOL_dp
        e_bounds(2, :) = 1.0_dp-TOL_dp

        ! IT READS THE FILES AND THE NAMES OF THE BODIES AND DETERMINES THE PARAMETERS TO BE FITTED
        call read_list(1)
        flush(6)

        ! IT READS RV DATA
        nRV = 0
        call read_RVobs(1)
        flush(6)
        nRV = obsData%obsRV%nRV

        ! IT READS T0 DATA
        if (.not. allocated(obsData%obsT0)) allocate (obsData%obsT0(NB-1))
        call read_T0obs(1)
        if (idtra .ne. 0) then
            write (*, '(a,1000(i5,1x))') " T0 DATA: nT0 = ",&
                &obsData%obsT0(:)%nT0
            if (durcheck .eq. 1) write (*, '(a,1000(i5,1x))') " DUR DATA: nT0 = ",&
                &obsData%obsT0(:)%nT0
            write (*, '(a,1000(l2,1x))') ' DO TRANSIT FLAG: ', do_transit
        end if
        flush(6)

        ! IT SETS THE LINEAR EPHEMERIS FROM T0 DATA
        if (obsData%nTTs .gt. 0) then
            call set_ephem()
            call compute_oc(obsData%obsT0)
        end if
        flush(6)

        ! IT DETERMINES THE NDATA
        nTTs = obsData%nTTs
        nDurs = obsData%nDurs

        obsData%ndata = nRV+nTTs+nDurs
        obsData%nfree = 0 ! gamma now fitted

        ndata = obsData%ndata
        nfree = obsData%nfree

        write (*, '(a,a)') " NUMBER OF ORBITAL PARAMETERS TO FIT: nfit = ", trim(adjustl(string(nfit)))
        nfit = nfit+rv_trend_order+nRVset+nRVset ! 2 x nRVset: gamma + jitter for each RV dataset
        if (rv_trend_order .gt. 0) then
            write (*, '(a,i2)') " RV trend of order ", rv_trend_order
        end if
        if (nRVset .gt. 0) then
            write (*, '(a,i2)') " number RV dataset ", nRVset
        end if
        write (*, '(a,a)') " NUMBER OF PARAMETERS TO FIT: nfit = ", trim(adjustl(string(nfit)))
        obsData%dof = (obsData%ndata-nfit-obsData%nfree)
        dof = obsData%dof
        flush(6)

        if (dof .le. 0) then
            write (*, '(a,a)') ' FOUND dof <= 0 SO IT IS FORCED TO 1 IN CASE',&
            &' THE USER WANT TO SIMULATE/INTEGRATE AND NOT CHECK THE FIT.'
            obsData%dof = 1
            dof = obsData%dof
        end if
        flush(6)

        obsData%inv_dof = one/real(obsData%dof, dp)
        inv_dof = obsData%inv_dof

        ! IT DEFINES THE ID OF THE PARAMETERS TO BE FITTED
        call idpar()
        flush(6)

        ! IT READS BOUNDARIES OF THE KEPLERIAN ELEMENTS
        call read_fullpar(1, mass, radius, period, sma, ecc, argp, meanA, inc, longN, system_parameters)
        flush(6)

        ! IT DETERMINES THE LN_ERR_CONST TO COMPUTE LOGLIKELIHOOD
        ln_err_const = get_lnec(obsData)

        ! IT SETS THE LIST OF THE PARAMETERS TO FIT
        call set_parid_list()
        ! IT SETS FITNESS PARAMETERS
        if (nRV .ne. 0 .and. nTTs .ne. 0) then
            if (durcheck .eq. 0) then
                allocate (nset(2))
                nset(1) = nRV
                nset(2) = nTTs
            else
                allocate (nset(3))
                nset(1) = nRV
                nset(2) = nTTs
                nset(3) = nDurs
            end if
        else if (nRV .ne. 0 .and. nTTs .eq. 0) then
            allocate (nset(1))
            nset(1) = nRV
        else if (nRV .eq. 0 .and. nTTs .ne. 0) then
            if (durcheck .eq. 0) then
                allocate (nset(1))
                nset(1) = nTTs
            else
                allocate (nset(2))
                nset(1) = nTTs
                nset(2) = nDurs
            end if
        else
            allocate (nset(1))
            nset(1) = 1
        end if
        deallocate (nset)
        flush(6)
        

        write (*, '(a)') ''
        write (*, '(a)') 'Initial-input Keplerian orbital elements: val'
        write (*, '(a, 1000(1x,es23.16))') "mass     [Msun] = ", mass
        write (*, '(a, 1000(1x,es23.16))') "radius   [Rsun] = ", radius
        write (*, '(a, 1000(1x,es23.16))') "period   [days] = ", period
        write (*, '(a, 1000(1x,es23.16))') "sma      [au]   = ", sma
        write (*, '(a, 1000(1x,es23.16))') "ecc             = ", ecc
        write (*, '(a, 1000(1x,es23.16))') "argp     [deg]  = ", argp
        write (*, '(a, 1000(1x,es23.16))') "meana    [deg]  = ", meana
        write (*, '(a, 1000(1x,es23.16))') "inc      [deg]  = ", inc
        write (*, '(a, 1000(1x,es23.16))') "longn    [deg]  = ", longn
        write (*, '(a)') ''
        flush(6)

        write (*, '(a23,3(1x,a23))') 'Full System Parameters', "value", '( par_min ,', 'par_max )'
        do ipar = 1, npar
            write (*, '(a23,1x,es23.16,2(a,es23.16),a)') all_names_list(ipar), system_parameters(ipar),&
                &' ( ', par_min(ipar), ' , ', par_max(ipar), ' )'
        end do
        flush(6)

        ! read priors within fortran!
        call read_priors(1)
        
        str_len = len(parid(1))

        ! check if there are derived parameters to compute and to check
        call init_derived_parameters(1, path)

        ! deallocated variables not needed anymore
        if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)

        path = trim(adjustl(path_in))//trim(adjustl(sub_folder))

        return
    end subroutine initialize_trades
    ! ============================================================================

    subroutine get_system_parameters_and_minmax(npar, sys_params, sys_par_min, sys_par_max)
        ! Input
        integer,intent(in)::npar
        ! Output
        real(dp),dimension(npar),intent(out)::sys_params, sys_par_min, sys_par_max

        sys_params = system_parameters
        sys_par_min = par_min
        sys_par_max = par_max
    
        return
    end subroutine get_system_parameters_and_minmax

    ! ============================================================================

    subroutine get_default_fitting_parameters_and_minmax(fit_params, fit_par_minmax, n_fit)
        ! Input
        integer,intent(in)::n_fit
        ! Output
        real(dp),dimension(n_fit),intent(out)::fit_params
        real(dp),dimension(n_fit, 2),intent(out)::fit_par_minmax
        ! ! Hide
        ! integer intent(hide),depend(fit_params) :: n_fit=len(fit_params)
        ! Local
        real(dp),dimension(:),allocatable::fitting_parameters

        call init_param(system_parameters, fitting_parameters)
        fit_params = fitting_parameters
        deallocate(fitting_parameters)
        fit_par_minmax(:, 1) = minpar
        fit_par_minmax(:, 2) = maxpar
    
        return
    end subroutine get_default_fitting_parameters_and_minmax

    ! ============================================================================
    subroutine get_do_transit_flag(n_bodies, transit_flag)
        ! Input
        integer,intent(in)::n_bodies
        ! Output
        logical,dimension(n_bodies),intent(out)::transit_flag

        transit_flag = do_transit

        return
    end subroutine get_do_transit_flag

    subroutine get_observed_transits_info(n_planets, n_tra_per_body, Peph, Teph)
        ! Input
        integer,intent(in)::n_planets
        ! Output
        integer,dimension(n_planets),intent(out)::n_tra_per_body
        real(dp),dimension(n_planets),intent(out):: Peph, Teph

        n_tra_per_body = obsData%obsT0(:)%nT0
        Peph = obsData%obsT0(:)%Pephem
        Teph = obsData%obsT0(:)%Tephem

        return
    end subroutine get_observed_transits_info

    ! ============================================================================

    subroutine print_obsdata()
        ! Local
        integer::n_rv
        integer::n_obj, i_obj

        write(*, '(a)') "OBSERVATIONS DATA"
        ! print RV
        n_rv = obsData%obsRV%nRV
        write(*, '(a, i5)') "Number of RV observations: ", n_rv
        write(*, '(a, i5)') "Number of RV datasets: ", nRVset
        ! print T0
        if(allocated(obsData%obsT0))then
            n_obj = size(obsData%obsT0)
            do i_obj=1,n_obj
                write(*,'(a, i2, a, i4)')"Number of T0 observations for object ", i_obj, ": ", obsData%obsT0(i_obj)%nT0
            end do
        else
            write(*,'(a)')"No T0 data"
        end if

        return
    end subroutine print_obsdata

    ! ============================================================================

    subroutine get_priors(n,  names, values)
        integer, intent(in) :: n
        character(10), dimension(n), intent(out) ::  names
        real(dp), dimension(n, 3), intent(out):: values

        names = priors_names
        values = priors_values
    
        return
    end subroutine get_priors
    ! ============================================================================

    subroutine compute_ln_priors(n_all, all_fix_par, n_fit, fit_par, n, names, values, lnpr)
        integer,intent(in)::n_all, n_fit, n
        real(dp), dimension(n_all), intent(in):: all_fix_par
        real(dp), dimension(n_fit), intent(in):: fit_par
        character(10), dimension(n), intent(in) ::  names
        real(dp), dimension(n, 3), intent(in):: values
        real(dp),intent(out)::lnpr

        call ln_priors(all_fix_par, fit_par, names, values, lnpr)

        return
    end subroutine compute_ln_priors
    

    ! ============================================================================

    ! subroutine that returns the args stored in constants

    ! subroutine that sets the args stored in constants

    ! ============================================================================

    subroutine get_parameter_names(parameter_names, n_fit, str_length)
        integer, intent(in)::n_fit, str_length
        character(1), dimension(n_fit, str_length), intent(out)::parameter_names

        integer::in, ic

        do in = 1, n_fit
            do ic = 1, str_length
                parameter_names(in, ic) = parid(in) (ic:ic)
            end do
        end do

        return
    end subroutine get_parameter_names

    subroutine get_tofit(n_par, to_fit)
        ! Input
        integer,intent(in)::n_par
        ! Output
        integer,dimension(n_par),intent(out)::to_fit
        ! Local

        to_fit = tofit

        return
    end subroutine get_tofit

    subroutine get_all_keplerian_names(all_names, n_par, str_length)
        ! Input
        integer, intent(in)::n_par, str_length
        ! Output
        character(1), dimension(n_par, str_length), intent(out)::all_names
        ! Local
        integer::in, ic

        do in = 1, n_par
            do ic = 1, str_length
                all_names(in, ic) = all_names_list(in) (ic:ic)
            end do
        end do

        return
    end subroutine get_all_keplerian_names

    ! ============================================================================

    subroutine init_fit_parameters(all_parameters, n_par, fit_parameters, n_fit)
        use parameters_conversion, only: init_param
        integer, intent(in)::n_par
        real(dp), dimension(n_par), intent(in)::all_parameters
        integer, intent(in)::n_fit
        real(dp), dimension(n_fit), intent(out)::fit_parameters

        real(dp), dimension(:), allocatable::temp_fit

        call init_param(all_parameters, temp_fit)
        fit_parameters = temp_fit
        if (allocated(temp_fit)) deallocate (temp_fit)

        return
    end subroutine init_fit_parameters

    ! ============================================================================

    subroutine set_one_fit_par_boundaries(ifit, min_val, max_val, n_fit, fit_par_minmax)
        ! Input
        integer,intent(in)::ifit, n_fit
        real(dp),intent(in)::min_val, max_val
        ! Input/Output
        real(dp),dimension(n_fit, 2),intent(inout)::fit_par_minmax

        minpar(ifit) = min_val
        fit_par_minmax(ifit, 1) = minpar(ifit)
        maxpar(ifit) = max_val
        fit_par_minmax(ifit, 2) = maxpar(ifit)

        return
    end subroutine set_one_fit_par_boundaries

    subroutine reset_all_fit_boundaries(n_fit, fit_par_minmax)
        ! Input
        integer,intent(in)::n_fit
        ! Output
        real(dp),dimension(n_fit, 2),intent(out)::fit_par_minmax

        minpar = minpar_bck
        fit_par_minmax(:, 1) = minpar
        maxpar = maxpar_bck
        fit_par_minmax(:, 2) = maxpar

        return
    end subroutine reset_all_fit_boundaries

    ! ============================================================================
    subroutine init_pso(cpuid, path_in)
        ! Input
        integer, intent(in)::cpuid
        character(*), intent(in)::path_in
        ! Local variables
        character(512)::path_temp

        path_temp = trim(adjustl(path))
        path = trim(adjustl(path_in))
        call read_pso_opt(cpuid)
        n_global = nGlobal
        path = trim(adjustl(path_temp))

        return
    end subroutine init_pso

    ! ============================================================================
    subroutine check_boundaries(n_all, all_fix_par, n_fit, fit_par, check)
        integer,intent(in)::n_all, n_fit
        real(dp),dimension(n_all),intent(in)::all_fix_par
        real(dp),dimension(n_fit),intent(in)::fit_par
        logical,intent(out)::check

        check = .true.
        check = check_only_boundaries(all_fix_par, fit_par)

        return
    end subroutine

    ! ============================================================================

    subroutine trades_rchisq(fit_parameters, rchisq, n_fit)
        ! Input
        integer, intent(in)::n_fit
        real(dp), dimension(n_fit), intent(in)::fit_parameters
        ! Output
        real(dp), intent(out)::rchisq
        ! local variables
        real(dp)::chi_square, lnL, lnc, lnprior, bic

        rchisq = zero

        ! rchisq = base_fitness_function(system_parameters, fit_parameters)
        call base_fitness_function(system_parameters, fit_parameters,&
            &chi_square, rchisq, lnL, lnprior, lnc, bic)

        return
    end subroutine trades_rchisq

    ! ============================================================================

    subroutine fortran_fitness_function(fit_parameters, chi_square,&
         &reduced_chi_square, lgllhd, lnprior, ln_const, bic, check, n_fit)
        ! Input
        integer, intent(in)::n_fit
        real(dp), dimension(n_fit), intent(in)::fit_parameters
        ! Output
        real(dp), intent(out)::chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic
        logical, intent(out)::check
        ! Local variables

        check = .true.
        chi_square = zero
        lgllhd = zero

        call bound_fitness_function(system_parameters, fit_parameters,&
            &chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic)
        if (chi_square .ge. resmax) check = .false.

        return
    end subroutine fortran_fitness_function

    ! ============================================================================

    subroutine fortran_loglikelihood(fit_parameters, lgllhd, check, n_fit)
        ! Input
        integer, intent(in)::n_fit
        real(dp), dimension(n_fit), intent(in)::fit_parameters
        ! Output
        real(dp), intent(out)::lgllhd
        logical, intent(out)::check
        ! Local variables
        ! integer::ns, ne
        ! real(dp)::fitness, lnconst
        ! real(dp), dimension(:), allocatable::jitter
        real(dp)::chi_square, reduced_chi_square, lnprior, ln_const, bic

        check = .true.
        chi_square = zero
        lgllhd = zero

        ! write(*,*)"DEBUG: bound_fitness_function"
        call bound_fitness_function(system_parameters, fit_parameters,&
            &chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic)
        if (chi_square .ge. resmax) check = .false.
        ! write(*,*)"DEBUG: chi_square, reduced_chi_square, lgllhd, ln_const, bic"
        ! write(*,*)chi_square, reduced_chi_square, lgllhd, ln_const, bic

        return
    end subroutine fortran_loglikelihood

    ! ============================================================================

    subroutine fortran_logprob(fit_parameters, lgllhd, lnp, check, n_fit)
        ! Input
        integer, intent(in)::n_fit
        real(dp), dimension(n_fit), intent(in)::fit_parameters
        ! Output
        real(dp), intent(out)::lgllhd, lnp
        logical, intent(out)::check
        ! Local variables
        ! integer::ns, ne
        ! real(dp)::fitness, lnconst
        ! real(dp), dimension(:), allocatable::jitter
        real(dp)::chi_square, reduced_chi_square, ln_const, bic

        check = .true.
        chi_square = zero
        lgllhd = zero
        lnp = zero

        ! write(*,*)"DEBUG: bound_fitness_function"
        call bound_fitness_function(system_parameters, fit_parameters,&
            &chi_square, reduced_chi_square, lgllhd, lnp, ln_const, bic)
        if (chi_square .ge. resmax) check = .false.

        return
    end subroutine fortran_logprob

    ! ============================================================================

    subroutine write_summary_files(write_number, parameters_values,&
        &chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic, check, n_fit)
        use driver, only: write_summary_nosigma
        ! Input
        integer, intent(in)::n_fit
        integer, intent(in)::write_number
        real(dp), dimension(n_fit), intent(in)::parameters_values
        ! Output
        real(dp), intent(out)::chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic
        logical, intent(out)::check
        ! local variables
        real(dp), dimension(:), allocatable::run_all_par
        logical::check_status
        ! integer::ns, ne
!     logical::wrt_info=.true.

        check = .true.
        check_status = .true.

        check = check_only_boundaries(system_parameters, parameters_values)
        ! write (*, *) "DEBUG check_only_boundaries : ", check
        lnprior = zero

        allocate (run_all_par(npar))
        run_all_par = system_parameters
        if (check_derived) check_status = check_derived_parameters(parameters_values)
        if (fix_derived) call fix_derived_parameters(parameters_values, run_all_par, check_status)
        call write_summary_nosigma(1, write_number, 0, run_all_par, parameters_values,&
            &chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic) ! add .true. at the end for full debug
        deallocate (run_all_par)
        ! write (*, *) "DEBUG check = ", check, " check_status = ", check_status
        ! write (*, *) "DEBUG after write_summary_nosigma: reduced_chi_square = ", reduced_chi_square
        ! flush (6)
        if (.not. check .or. .not. check_status) chi_square = resmax
        call set_fitness_values(parameters_values,&
            &chi_square, reduced_chi_square, lgllhd, ln_const, bic)

        if (.not. check) then
            write (*, '(a)') '*******'
            write (*, '(a)') 'WARNING'
            write (*, '(a)') 'WARNING'
            write (*, '(a)') 'FITTED PARAMETERS COULD NOT BE PHYSICAL!'
            write (*, '(a)') 'BE VERY CAREFUL WITH THIS PARAMETER SET!'
            write (*, '(a)') 'WARNING'
            write (*, '(a)') 'WARNING'
            write (*, '(a)') '*******'
        end if

        return
    end subroutine write_summary_files

    ! ============================================================================

    subroutine write_summary_files_long(write_number, all_parameters, n_par, parameters_values, n_fit,&
        &chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic, check)
        use driver, only: write_summary_nosigma
!     use ode_run,only:ode_out
!     use output_files,only:write_parameters
        integer, intent(in)::n_par
        real(dp), dimension(n_par), intent(in)::all_parameters
        integer, intent(in)::n_fit
        integer, intent(in)::write_number
        real(dp), dimension(n_fit), intent(in)::parameters_values

        real(dp), intent(out)::chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic
        logical, intent(out)::check

        real(dp), dimension(:), allocatable::run_all_par
        logical::check_status
        ! real(dp),dimension(:),allocatable::jitter
        ! integer::ns,ne
        ! real(dp)::ln_const

        check = .true.
        check_status = .true.

!     check=check_fit_boundaries(parameters_values)
        check = check_only_boundaries(system_parameters, parameters_values)
        lnprior = zero

        if (check) then
            allocate (run_all_par(n_par))
            run_all_par = all_parameters
            if (check_derived) check_status = check_derived_parameters(parameters_values)
            if (fix_derived) call fix_derived_parameters(parameters_values, run_all_par, check_status)
            if (check_status) then
                call write_summary_nosigma(1, write_number, 0, run_all_par, parameters_values,&
                &chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic)
            else
                chi_square = resmax ! set it to resmax
                call set_fitness_values(parameters_values,&
                & chi_square, reduced_chi_square, lgllhd, ln_const, bic)
                check = .false.
            end if
            deallocate (run_all_par)

        else
            chi_square = resmax ! set it to resmax
            call set_fitness_values(parameters_values,&
                & chi_square, reduced_chi_square, lgllhd, ln_const, bic)
            write (*, '(a)') '*******'
            write (*, '(a)') 'WARNING'
            write (*, '(a)') 'WARNING'
            write (*, '(a)') 'FITTED PARAMETERS COULD NOT BE PHYSICAL!'
            write (*, '(a)') 'BE VERY CAREFUL WITH THIS PARAMETER SET!'
            write (*, '(a)') 'WARNING'
            write (*, '(a)') 'WARNING'
            write (*, '(a)') '*******'
        end if

        ! call fortran_loglikelihood(parameters_values, lgllhd, check, n_fit)
        ! if (chi_square .ge. resmax) check = .false.

        return
    end subroutine write_summary_files_long

    ! ============================================================================

    subroutine set_fitness(nrv, nt0, ndur, err_rv, err_t0, err_dur, dof, chi_square_in, chi_square, &
        &reduced_chi_square, lnLikelihood, ln_const, bic)
        integer,intent(in)::nrv, nt0, ndur, dof
        real(dp),dimension(nrv),intent(in)::err_rv
        real(dp),dimension(nt0),intent(in)::err_t0
        real(dp),dimension(ndur),intent(in)::err_dur
        real(dp),intent(in)::chi_square_in
        real(dp), intent(out)::chi_square,reduced_chi_square, lnLikelihood, ln_const, bic

        ! Local variables
        real(dp)::inv_dof, nd, ln_rv=zero, ln_t0=zero, ln_dur=zero

        ! write(*,*)"in set_fitness"
        chi_square = chi_square_in
        ! write(*,*) "chi_square_in -> chi_square"

        if (chi_square .ge. resmax) then
            chi_square = resmax
        end if
        ! write(*,*)"chi_square = ",chi_square

        inv_dof =one/real(dof,dp)
        reduced_chi_square = chi_square*inv_dof
        ! write(*,*)"reduced_chi_square = ",reduced_chi_square

        ln_const = zero
        if (nrv .gt. 0)then
            ln_rv = sum(log(pack(err_rv, err_rv /= zero)*pack(err_rv, err_rv /= zero)))
        end if
        if (nt0 .gt. 0)then
            ln_t0 = sum(log(pack(err_t0, err_t0 /= zero)*pack(err_t0, err_t0 /= zero)))
        end if
        if (ndur .gt. 0)then
            ln_dur = sum(log(pack(err_dur, err_dur /= zero)*pack(err_dur, err_dur /= zero)))
        end if
        nd = real(nrv+nt0+ndur, dp)
        ln_const = ln_const -half*nd*log(dpi)
        ln_const = ln_const -half*(ln_rv + ln_t0 + ln_dur)

        lnLikelihood = -half*chi_square+ln_const
        ! write(*,*)"lnLikelihood = ",lnLikelihood

        ! bic = -two*lnLikelihood + real(nfit,dp)*log(real(ndata,dp))
        bic = -two*lnLikelihood+bic_const ! bic_const global variable
        ! write(*,*)"bic = ",bic


        return
    end subroutine set_fitness

    ! ============================================================================

    ! pso
    subroutine pyrun_pso(n_fit, i_global, best_parameters, best_fitness)
        use opti_pso, only: pso_driver, evaluate_pso
        integer, intent(in)::n_fit
        integer, intent(in)::i_global
        real(dp), dimension(n_fit), intent(out)::best_parameters
        real(dp), intent(out)::best_fitness

        path = trim(adjustl(path))
        best_parameters = zero
        best_fitness = zero
        call pso_driver(i_global, evaluate_pso, n_fit, system_parameters, minpar, maxpar,&
          &best_parameters, best_fitness) ! PSO DRIVER

        return
    end subroutine pyrun_pso

    ! ============================================================================

    ! init both cases for derived parameters
    ! 1)
    ! check if there are derived parameters to compute and to check
    subroutine init_check_parameters(cpuid, path_in)
        integer, intent(in)::cpuid ! cpu number: use 1
        character(*), intent(in)::path_in ! path of the folder with derived_boundaries.dat

        call init_check_derived_parameters(cpuid, path_in)

        return
    end subroutine init_check_parameters

    ! ============================================================================

    ! 2)
    ! check if there are derived parameters to compute and to check
    subroutine init_fix_parameters(n_derived_in, in_names, in_parameters)
        integer, intent(in)::n_derived_in
        character(15), dimension(n_derived_in), intent(in)::in_names
        real(dp), dimension(n_derived_in), intent(in)::in_parameters

        call init_fix_derived_parameters(n_derived_in, in_names, in_parameters)

        return
    end subroutine init_fix_parameters

    ! ============================================================================

    subroutine deallocate_variables()

        call deallocate_all() ! from 'parameters' module

        return
    end subroutine deallocate_variables

    ! ============================================================================
    ! +++***+++
    ! CHECKS THE NEXT!
    ! +++***+++

    subroutine set_time_epoch(t_epoch)
        ! Input
        ! t_epoch        == reference time epoch
        real(dp), intent(in)::t_epoch
        ! Output
        ! NONE
        tepoch = t_epoch

        return
    end subroutine set_time_epoch

    subroutine set_time_start(t_start)
        ! Input
        ! t_start        == start of the integration
        real(dp), intent(in)::t_start
        ! Output
        ! NONE
        tstart = t_start

        return
    end subroutine set_time_start

    subroutine set_time_int(t_int)
        ! Input
        ! t_int          == total integration time in days
        real(dp), intent(in)::t_int
        ! Output
        ! NONE
        tint = t_int

        return
    end subroutine set_time_int

    subroutine set_hill_check(hill_check)
        ! Input
        logical, intent(in)::hill_check
        ! Output
        ! NONE
        do_hill_check = hill_check

        return
    end subroutine set_hill_check

    subroutine set_amd_hill_check(ah_check)
        ! Input
        logical, intent(in)::ah_check
        ! Output
        ! NONE
        amd_hill_check = ah_check

        return
    end subroutine set_amd_hill_check

    subroutine set_rv_res_gls(res_gls)
        ! Input
        logical, intent(in)::res_gls
        ! Output
        ! NONE
        rv_res_gls = res_gls

        return
    end subroutine set_rv_res_gls

    subroutine set_durcheck(duration_check)
        ! Input
        integer,intent(in)::duration_check

        durcheck = duration_check
    
        return
    end subroutine set_durcheck

    ! SUBROUTINE TO INITIALISE TRADES WITHOUT READING FILES
    subroutine args_init(n_body, duration_check)
        ! INPUT
        ! n_body         == number of bodies (take into account the star)
        ! duration_check == check or not the duration, 0 not check, 1 to check
        ! t_int          == total integration time in days
        ! OUTPUT
        ! None ==> some variables set globally

        ! Input
        integer, intent(in)::n_body
        integer, intent(in)::duration_check

        NB = n_body
        NBDIM = n_body*6

        if (allocated(e_bounds)) deallocate (e_bounds)
        allocate (e_bounds(2, n_body))
        e_bounds(1, :) = zero
        e_bounds(2, :) = one-TOL_dp

        call deallocate_dataObs(obsData)
        allocate (obsData%obsT0(n_body-1))

        call set_durcheck(duration_check) ! updating durcheck global variable

        amin = TOL_dp
        amax = 1.0e4_dp

        return
    end subroutine args_init

    ! ============================================================================
    subroutine set_rv_dataset(time_rv, obs_rv, obs_erv, rv_setid, n_rvset, n_rv)
        ! Input
        integer, intent(in)::n_rv
        real(dp), dimension(n_rv), intent(in)::time_rv, obs_rv, obs_erv
        integer, dimension(n_rv), intent(in):: rv_setid
        integer, intent(in)::n_rvset
        ! Local
        integer::n

        call deallocate_dataRV(obsData%obsRV)

        n = size(time_rv)
        nRV = n
        if (n .gt. 0) then
            rvcheck = 1
            call init_dataRV(n, obsData%obsRV)
            obsData%obsRV%jd = time_rv
            obsData%obsRV%RV = obs_rv
            obsData%obsRV%eRV = obs_erv
            obsData%obsRV%RVsetID = rv_setid
            obsData%obsRV%nRVset = n_rvset
        else
            rvcheck = 0
        end if

        return
    end subroutine

    subroutine deallocate_rv_dataset()

        call deallocate_dataRV(obsData%obsRV)
        nRV = 0

        return
    end subroutine deallocate_rv_dataset

    subroutine set_t0_dataset(body_id, epo, obs_t0, obs_et0, sources_id, n_t0)
        ! Input
        integer, intent(in)::body_id
        integer, intent(in)::n_t0
        integer, dimension(n_t0), intent(in)::epo
        real(dp), dimension(n_t0), intent(in)::obs_t0, obs_et0
        integer, dimension(n_t0), intent(in)::sources_id
        ! Local
        integer::n, i_body

        n = size(epo)
        i_body = body_id-1
        if (n .gt. 0) then
            call init_dataT0(n, obsData%obsT0(i_body), durcheck)
            obsData%obsT0(i_body)%epo = epo
            obsData%obsT0(i_body)%T0 = obs_t0
            obsData%obsT0(i_body)%eT0 = obs_et0
            obsData%obsT0(i_body)%source_id = sources_id
            obsData%nTTs = sum(obsData%obsT0(:)%nT0)
            nTTs = obsData%nTTs
            call set_ephem()
            call compute_oc_one_planet(obsData%obsT0(i_body))
        end if

        return
    end subroutine

    subroutine set_t14_dataset(body_id, t14, et14, n_t14)
        ! Input
        integer, intent(in)::body_id
        integer, intent(in)::n_t14
        real(dp), dimension(n_t14), intent(in)::t14, et14
        ! Local
        integer::i_body

        i_body = body_id-1
        obsData%obsT0(i_body)%dur = t14
        obsData%obsT0(i_body)%edur = et14
        obsData%obsT0(i_body)%nDur = n_t14
        obsData%nDurs = sum(obsData%obsT0(:)%nDur)
        nDurs = obsData%nDurs

        return
    end subroutine

    subroutine deallocate_t0_dataset(body_id)
        !Input
        integer, intent(in)::body_id

        call deallocate_dataT0(obsData%obsT0(body_id-1))
        obsData%nTTs = sum(obsData%obsT0(:)%nT0)
        nTTs=obsData%nTTs
        if (durcheck .eq. 1) then
            obsData%nDurs = sum(obsData%obsT0(:)%nDur)
        end if
        nDurs=obsData%nDurs

        return
    end subroutine deallocate_t0_dataset

    subroutine deallocate_all_data()

        call deallocate_dataObs(obsData)
        nRV=0
        nTTS=0
        nDurs=0

        return
    end subroutine deallocate_all_data
    ! ============================================================================

    ! ============================================================================
  !!! SUBROUTINE TO RUN TRADES INTEGRATION FROM PyORBIT
    subroutine kelements_to_rv_and_t0s(t_start, t_epoch, t_int,&
        &m_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg,&
        &transit_flag,& ! input
        &rv_sim,&
        &body_t0_sim, epo_sim, t0_sim, t14_sim, lambda_rm_sim, kel_sim, stable,& ! output
        &n_body, n_rv, n_tt, n_kep) ! dimensions to be not provided

        ! INPUT
        ! t_start      == start of the integration
        ! t_epoch      == reference time epoch
        ! t_int        == total integration time in days

        ! m_msun       == masses of all the bodies in Msun m_sun(n_body)
        ! R_rsun       == radii of all the bodies in Rsun r_rsun(n_body)
        ! P_day        == periods of all the bodies in days p_day(n_body); p_day(0) = 0
        ! ecc          == eccentricities of all the bodies ecc(n_body); ecc(0) = 0
        ! argp_deg     == argument of pericentre of all the bodies argp_deg(n_body); argp_deg(0) = 0
        ! mA_deg       == mean anomaly of all the bodies mA_deg(n_body); mA_deg(0) = 0
        ! inc_deg      == inclination of all the bodies inc_deg(n_body); inc_deg(0) = 0
        ! lN_deg       == longitude of node of all the bodies lN_deg(n_body); lN_deg(0) = 0

        ! transit_flag == logical/boolean vector with which bodies should transit (.true.) or not (.false.) transit_flag(n_body); transit_flag(0) = False

        ! OUTPUT
        ! rv_sim       == rv simulated in m/s, same dimension of t_rv
        ! body_t0_sim  == id of the body (2...NB) corresponding to the transit
        ! epo_sim      == epoch of each all T0s based on the linear ephemeris of the body
        ! t0_sim       == t0 simulated in days, same dimension of body_t0_sim
        ! t14_sim      == Total Duration = T4 - T1 simulated in minutes, same dimension of body_T0_sim
        ! lambda_rm_sim ==  lambda of the Rossiter-Mclaughing effect computed at each t0_sim
        ! kel_sim      == period, sma, ecc, inc, meana, argp, truea, longn computed at each t0_sim

        ! DIMENSIONS
        ! n_body       == number of bodies (take into account the star)
        ! n_rv, n_tt, n_kep == number of RVs and all T0s

        ! Input
        integer, intent(in)::n_body, n_rv, n_tt, n_kep
        real(dp), intent(in)::t_start, t_epoch, t_int
        real(dp), dimension(n_body), intent(in)::m_msun, R_rsun, P_day
        real(dp), dimension(n_body), intent(in)::ecc, argp_deg, mA_deg, inc_deg, lN_deg
        logical, dimension(n_body), intent(in)::transit_flag
        ! Output
        real(dp), dimension(n_rv), intent(out)::rv_sim
        integer, dimension(n_tt), intent(out)::body_t0_sim, epo_sim
        real(dp), dimension(n_tt), intent(out)::t0_sim, t14_sim, lambda_rm_sim
        real(dp), dimension(n_tt, n_kep), intent(out)::kel_sim
        logical, intent(out)::stable
        ! Local
        integer::id_transit_body = 1 ! needed to be == 1
        real(dp), dimension(:), allocatable::rv_temp
        integer, dimension(:), allocatable::body_t0_temp, epo_temp
        real(dp), dimension(:), allocatable::t0_temp, t14_temp, lbd_rm_temp
        real(dp), dimension(:, :), allocatable::kel_temp
        real(dp)::step

        step = minval(P_day(2:NB))/ten

        call integration_info_to_data(t_start, t_epoch, step, t_int,&
            &m_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg,&
            &id_transit_body, transit_flag, durcheck,&
            &rv_temp,&
            &body_t0_temp, epo_temp,&
            &t0_temp, t14_temp, lbd_rm_temp, kel_temp,&
            &stable)

        if (allocated(rv_temp))then
            rv_sim = rv_temp
        end if

        if (allocated(body_t0_temp)) then
            body_t0_sim = body_t0_temp
            epo_sim = epo_temp
            t0_sim = t0_temp
            t14_sim = t14_temp
            lambda_rm_sim = lbd_rm_temp
            kel_sim = kel_temp
        end if

        if (allocated(rv_temp)) deallocate (rv_temp)
        if (allocated(body_T0_temp)) deallocate (body_T0_temp, epo_temp)
        if (allocated(t0_temp)) deallocate (t0_temp, t14_temp, lbd_rm_temp, kel_temp)

        return
    end subroutine kelements_to_rv_and_t0s

    subroutine kelements_to_rv(t_start, t_epoch, t_int,&
        &m_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg,&
        &rv_sim, stable,&
        &n_body, n_rv) ! dimensions to be not provided

        ! INPUT
        ! t_start      == start of the integration
        ! t_epoch      == reference time epoch
        ! t_int        == total integration time in days

        ! m_msun       == masses of all the bodies in Msun m_sun(n_body)
        ! R_rsun       == radii of all the bodies in Rsun r_rsun(n_body)
        ! P_day        == periods of all the bodies in days p_day(n_body); p_day(0) = 0
        ! ecc          == eccentricities of all the bodies ecc(n_body); ecc(0) = 0
        ! argp_deg     == argument of pericentre of all the bodies argp_deg(n_body); argp_deg(0) = 0
        ! mA_deg       == mean anomaly of all the bodies mA_deg(n_body); mA_deg(0) = 0
        ! inc_deg      == inclination of all the bodies inc_deg(n_body); inc_deg(0) = 0
        ! lN_deg       == longitude of node of all the bodies lN_deg(n_body); lN_deg(0) = 0

        ! OUTPUT
        ! rv_sim       == rv simulated in m/s, same dimension of t_rv

        ! DIMENSIONS
        ! n_body       == number of bodies (take into account the star)
        ! n_rv         == number of RVs

        ! Input
        integer, intent(in)::n_body, n_rv
        real(dp), intent(in)::t_start, t_epoch, t_int
        real(dp), dimension(n_body), intent(in)::m_msun, R_rsun, P_day
        real(dp), dimension(n_body), intent(in)::ecc, argp_deg, mA_deg, inc_deg, lN_deg
        ! Output
        real(dp), dimension(n_rv), intent(out)::rv_sim
        logical, intent(out)::stable
        ! Local
        integer::id_transit_body = 0 ! needed to be == 0, so no needs to check for transits
        real(dp), dimension(:), allocatable::rv_temp
        logical, dimension(n_body)::transit_flag
        integer, dimension(:), allocatable::body_t0_temp, epo_temp
        real(dp), dimension(:), allocatable::t0_temp, t14_temp, lbd_rm_temp
        real(dp), dimension(:, :), allocatable::kel_temp
        real(dp)::step

        ! transit_flag(1) = .false.
        ! transit_flag(2:NB) = .true.
        transit_flag = .false.
        step = minval(P_day(2:NB))/ten

        call integration_info_to_data(t_start, t_epoch, step, t_int,&
            &m_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg,&
            &id_transit_body, transit_flag, durcheck,&
            &rv_temp,&
            &body_t0_temp, epo_temp,&
            &t0_temp, t14_temp, lbd_rm_temp, kel_temp, stable)

        rv_sim = rv_temp

        if (allocated(rv_temp)) deallocate (rv_temp)
        if (allocated(body_T0_temp)) deallocate (body_T0_temp, epo_temp)
        if (allocated(t0_temp)) deallocate (t0_temp, t14_temp, lbd_rm_temp, kel_temp)

        return
    end subroutine kelements_to_rv

    subroutine kelements_to_t0s(t_start, t_epoch, t_int,&
        &m_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg,&
        &transit_flag,& ! input
        &body_t0_sim, epo_sim, t0_sim, t14_sim, lambda_rm_sim, kel_sim, stable,& ! output
        &n_body, n_tt, n_kep) ! dimensions to be not provided

        ! INPUT
        ! t_start      == start of the integration
        ! t_epoch      == reference time epoch
        ! t_int        == total integration time in days

        ! m_msun       == masses of all the bodies in Msun m_sun(n_body)
        ! R_rsun       == radii of all the bodies in Rsun r_rsun(n_body)
        ! P_day        == periods of all the bodies in days p_day(n_body); p_day(0) = 0
        ! ecc          == eccentricities of all the bodies ecc(n_body); ecc(0) = 0
        ! argp_deg     == argument of pericentre of all the bodies argp_deg(n_body); argp_deg(0) = 0
        ! mA_deg       == mean anomaly of all the bodies mA_deg(n_body); mA_deg(0) = 0
        ! inc_deg      == inclination of all the bodies inc_deg(n_body); inc_deg(0) = 0
        ! lN_deg       == longitude of node of all the bodies lN_deg(n_body); lN_deg(0) = 0

        ! transit_flag == logical/boolean vector with which bodies should transit (.true.) or not (.false.) transit_flag(n_body); transit_flag(0) = False

        ! OUTPUT
        ! body_t0_sim  == id of the body (2...NB) corresponding to the transit
        ! epo_sim      == epoch of each all T0s based on the linear ephemeris of the body
        ! t0_sim       == t0 simulated in days, same dimension of body_t0_sim
        ! t14_sim      == Total Duration = T4 - T1 simulated in minutes, same dimension of body_T0_sim
        ! lambda_rm_sim == lambda of the Rossiter-Mclaughin effect computed at each t0_sim
        ! kel_sim      == period, sma, ecc, inc, meana, argp, truea, longn computed at each t0_sim

        ! DIMENSIONS
        ! n_body       == number of bodies (take into account the star)
        ! n_tt, n_kep  == number of all T0s

        ! Input
        integer, intent(in)::n_body, n_tt, n_kep
        real(dp), intent(in)::t_start, t_epoch, t_int
        real(dp), dimension(n_body), intent(in)::m_msun, R_rsun, P_day
        real(dp), dimension(n_body), intent(in)::ecc, argp_deg, mA_deg, inc_deg, lN_deg
        logical, dimension(n_body), intent(in)::transit_flag
        ! Output
        integer, dimension(n_tt), intent(out)::body_t0_sim, epo_sim
        real(dp), dimension(n_tt), intent(out)::t0_sim, t14_sim, lambda_rm_sim
        real(dp), dimension(n_tt, n_kep), intent(out)::kel_sim
        logical, intent(out)::stable
        ! Local
        integer::id_transit_body = 1 ! needed to be == 1
        real(dp), dimension(:), allocatable::rv_temp
        integer, dimension(:), allocatable::body_t0_temp, epo_temp
        real(dp), dimension(:), allocatable::t0_temp, t14_temp, lbd_rm_temp
        real(dp), dimension(:, :), allocatable::kel_temp
        real(dp)::step

        step = minval(P_day(2:NB))/ten

        call integration_info_to_data(t_start, t_epoch, step, t_int,&
            &m_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg,&
            &id_transit_body, transit_flag, durcheck,&
            &rv_temp,&
            &body_t0_temp, epo_temp,&
            &t0_temp, t14_temp, lbd_rm_temp, kel_temp,stable)

        ! rv_sim      = rv_temp

        body_t0_sim = body_t0_temp
        epo_sim = epo_temp
        t0_sim = t0_temp
        t14_sim = t14_temp
        lambda_rm_sim = lbd_rm_temp
        kel_sim = kel_temp

        if (allocated(rv_temp)) deallocate (rv_temp)
        if (allocated(body_T0_temp)) deallocate (body_T0_temp, epo_temp)
        if (allocated(t0_temp)) deallocate (t0_temp, t14_temp, lbd_rm_temp, kel_temp)

        return
    end subroutine kelements_to_t0s
    ! ============================================================================

    subroutine kelements_to_orbits(n_steps, n_body, nb_dim, time_steps, mass, radius, period, ecc, argp, meanA, inc, longN, orbits, check)
        ! Input
        integer, intent(in)::n_steps, n_body, nb_dim
        real(dp), dimension(n_steps), intent(in)::time_steps
        real(dp), dimension(n_body), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN
        ! Output
        real(dp), dimension(n_steps, nb_dim), intent(out)::orbits
        logical,intent(out)::check

        call args_init(n_body, 1)
        call ode_keplerian_elements_to_orbits(time_steps, mass, radius, period, ecc, argp, meanA, inc, longN, orbits, check)

        return
    end subroutine kelements_to_orbits

    subroutine orbits_to_rvs(n_steps, n_body, nb_dim, mass, orbits, rvs)
        ! Input
        integer, intent(in)::n_body, n_steps, nb_dim
        real(dp), dimension(n_body), intent(in)::mass
        real(dp), dimension(n_steps, nb_dim), intent(in):: orbits
        ! Output
        real(dp), dimension(n_steps), intent(out):: rvs
        ! Locals
        integer::i_steps

        do i_steps = 1, n_steps
            call get_RV(mass, orbits(i_steps, :), rvs(i_steps))
        end do

        return
    end subroutine orbits_to_rvs

    subroutine orbits_to_transits(n_steps, n_body, nb_dim, n_all_transits,&
        &time_steps, mass, radius, orbits,&
        &transiting_body, transits, durations, lambda_rm, kep_elem, body_flag)
        ! Input
        integer, intent(in)::n_steps, n_body, nb_dim, n_all_transits
        real(dp), dimension(n_steps), intent(in)::time_steps
        real(dp), dimension(n_body), intent(in)::mass, radius
        real(dp), dimension(n_steps, nb_dim), intent(in)::orbits
        integer, intent(in)::transiting_body
        ! Output
        real(dp), dimension(n_all_transits), intent(out)::transits, durations, lambda_rm
        real(dp), dimension(n_all_transits, 8), intent(out)::kep_elem
        integer, dimension(n_all_transits), intent(out)::body_flag
        ! Local
        integer, dimension(:), allocatable::X, Y, Z
        integer, dimension(:), allocatable::VX, VY, VZ
        real(dp)::A, B, AB
        logical::ABflag, Zflag
        logical::do_transit_check
        integer::body_transiting_start, body_transiting_end
        integer::i_steps, i_body, i_tra
        real(dp)::trun1, trun2, integration_step
        real(dp), dimension(:), allocatable::r1, r2, rtra
        logical, dimension(:), allocatable::do_transit_flag
        logical::check_tra
        real(dp)::tra, dur, lbd_rm, dummy

        call set_checking_coordinates(n_body, X, Y, Z, VX, VY, VZ)
        call set_transiting_bodies(transiting_body, do_transit_check, body_transiting_start, body_transiting_end)
        allocate (do_transit_flag(n_body))
        do_transit_flag = .false.
        check_tra = .true.

        body_flag = 0 ! set all flags to 0

        allocate (r1(nb_dim), r2(nb_dim), rtra(nb_dim))
        trun1 = time_steps(1)
        r1 = orbits(1, :)
        i_tra = 0
        steps: do i_steps = 2, n_steps
            trun2 = time_steps(i_steps)
            integration_step = trun2-trun1
            r2 = orbits(i_steps, :)

            if (do_transit_check) then
                body: do i_body = body_transiting_start, body_transiting_end

                    call transit_conditions(r1(X(i_body):VZ(i_body)), r2(X(i_body):VZ(i_body)),&
                        &A, B, AB, ABflag, Zflag)

                    do_transit_flag(i_body) = .true.
                    if (ABflag .and. Zflag) then
                        call compute_transit_time(zero, i_body, mass, radius, r1, r2, trun1, integration_step, do_transit_flag,&
                            &tra, dur, lbd_rm, rtra, check_tra)
                        if (check_tra) then
                            i_tra = i_tra+1
                            transits(i_tra) = tra
                            durations(i_tra) = dur
                            lambda_rm(i_tra) = lbd_rm
                            body_flag(i_tra) = i_body
                            call elements_one_body(i_body, mass, rtra,&
                                &kep_elem(i_tra, 1), kep_elem(i_tra, 2), kep_elem(i_tra, 3),& ! period, sma, ecc,
                                &kep_elem(i_tra, 4), kep_elem(i_tra, 5), kep_elem(i_tra, 6),& ! inc, meanA, argp,
                                &kep_elem(i_tra, 7), kep_elem(i_tra, 8),& ! longN, trueA
                                &dummy) ! dttau not used
                            ! call spin_orbit_misalignment(i_body, mass, radius, rtra, tra, lambda_rm(i_tra))
                        end if
                    end if

                end do body
            end if
            r1 = r2
            trun1 = trun2

        end do steps

        deallocate (X, Y, Z)
        deallocate (VX, VY, VZ)
        deallocate (r1, r2, rtra)

        return
    end subroutine orbits_to_transits
    ! ============================================================================
    
    subroutine orbits_to_elements(n_steps, n_body, nb_dim, mass, orbits, period, sma, ecc, inc, meana, argp, longn, truea, dttau)
        integer, intent(in)::n_steps, n_body, nb_dim
        real(dp), dimension(n_body), intent(in)::mass
        real(dp), dimension(n_steps, nb_dim), intent(in)::orbits
        real(dp), dimension(n_steps, n_body), intent(out)::period, sma, ecc, inc, meana, argp, longn, truea, dttau

        integer::i_steps

        do i_steps=1,n_steps
            call elements(mass, orbits(i_steps, :),&
                &period(i_steps, :), sma(i_steps, :), ecc(i_steps, :),&
                &inc(i_steps, :), meana(i_steps, :), argp(i_steps, :),&
                &truea(i_steps, :), longn(i_steps, :), dttau(i_steps, :))
        end do

        return
    end subroutine


    ! ============================================================================

    ! create wrappers to init the grid of a perturber body (a)
    ! and set properly the parameters to use in TRADES (b)
    ! the output will be all the transit times within the t_int
    ! for each parameter combination.

    ! (a)
    subroutine wrapper_read_grid(cpuid, ngrid, ncol_grid)
        integer, intent(in)::cpuid
        integer, intent(out)::ngrid, ncol_grid

        ! read input file <-> idpert
        call read_parameters_grid(cpuid, perturber_parameters_grid)
        ! perturber_parameters_grid in 'parameters' module

        ! set properly the fields of the perturber_parameters_grid variable
        call set_parameters_grid(perturber_parameters_grid, ngrid)

        ncol_grid = 10 ! fixed in module grid_search -> build_grid

        return
    end subroutine wrapper_read_grid

    subroutine wrapper_grid_init(cpuid, ngrid, ncol_grid, perturber_grid)
        integer, intent(in)::cpuid
        integer, intent(in)::ngrid, ncol_grid
        real(dp), dimension(ngrid, ncol_grid), intent(out)::perturber_grid

!     !f2py    integer,intent(hide),depend(perturber_grid)::ngrid=shape(perturber_grid,0), ncol_grid=shape(perturber_grid,1)

        real(dp), dimension(:, :), allocatable::temp_grid, fitness_grid

        ! create/build the full grid, with all the combination of the parameters of perturber body
        call build_grid(MR_star(1, 1), perturber_parameters_grid, temp_grid, fitness_grid, ngrid)
        ! perturber_parameters_grid in 'parameters' module

        perturber_grid = temp_grid

        deallocate (temp_grid, fitness_grid) ! temporary -> useless here

        ! needed to update the max allowed value of the semi-major axis
        amax = 5.0_dp*maxval(perturber_grid(:, 4))

        return
    end subroutine wrapper_grid_init

    ! ============================================================================

    subroutine get_max_nt0_nrv(P, n_bds, nt0_full, nrv_nmax)
        ! Input
        integer, intent(in)::n_bds
        real(dp), dimension(n_bds), intent(in)::P
        ! Output
        integer, intent(out)::nt0_full, nrv_nmax
        ! Local
        integer, dimension(:), allocatable::nT0_perbody

        ! compute number maximum of all transit times of all the planets from the
        ! integration time (global variable tint)
        allocate (nT0_perbody(n_bds-1))
        nT0_perbody = int((1.5_dp*tint)/P(2:n_bds))+2 ! compute the number of transit for each planet ... avoid P(1)<->star
        nt0_full = sum(nT0_perbody) ! sum them all
        deallocate (nT0_perbody)

        ! compute the number of maximum rv from integration and write time
        ! (global variables tint, wrttime
        nrv_nmax = int(tint/wrttime)+2

        return
    end subroutine get_max_nt0_nrv

    ! ============================================================================

    subroutine get_ntt_nrv(par_fit, n_fit, ntt, nrv_max)
        ! Input
        integer, intent(in)::n_fit
        real(dp), dimension(n_fit), intent(in)::par_fit
        ! Output
        integer, intent(out)::ntt, nrv_max
        ! Local
        real(dp), dimension(NB)::mass, radius
        real(dp), dimension(NB)::period, sma, ecc, argp, meanA, inc, longN
        logical::tempcheck

        ! parameter to take into account that sometimes the real number of transits
        ! is greater than ntra = int(integration_time/period)
        ! integer, parameter::safe_counter=3

        call par2kel_fit(system_parameters, par_fit, mass, radius,&
          &period, sma, ecc, argp, meanA, inc, longN, tempcheck)

        ! ntt = sum(int(tint/period(2:NB)))+(NB-1)
        ! ntt = sum(int(tint/period(2:NB)))+((NB-1)*safe_counter)
        ntt = int(sum(int(tint/period(2:NB)))*1.25_dp)

        ! nrv_max = int(tint/wrttime)+2
        nrv_max = int(1.25_dp*tint/wrttime)

        return
    end subroutine get_ntt_nrv

    ! ============================================================================

    subroutine wrapper_set_grid_orbelem_fullgrid(sim_id, id_pert, perturber_grid, ngrid, ncol,&
      &m, R, P, sma, ecc, w, mA, inc, lN, n_bds, nt0_full, nrv_nmax)
        integer, intent(in)::sim_id, id_pert
        integer, intent(in)::ngrid, ncol
        real(dp), dimension(ngrid, ncol), intent(in)::perturber_grid
        integer, intent(in)::n_bds
        real(dp), dimension(n_bds), intent(out)::m, R, P, sma, ecc, w, mA, inc, lN
        integer, intent(out)::nt0_full, nrv_nmax

!f2py    integer,intent(hide),depend(perturber_grid)::ngrid=shape(perturber_grid,0), ncol=shape(perturber_grid,1)

! !f2py    integer,intent(hide),depend(m)::n_bds=len(m)

        real(dp), dimension(:), allocatable::sim_all_parameters, sim_fit_parameters
        logical::tempcheck

        integer, dimension(:), allocatable::nT0_perbody

        tempcheck = .true.
        allocate (sim_all_parameters(npar), sim_fit_parameters(nfit))
        sim_all_parameters = system_parameters
        sim_fit_parameters = zero

        ! 1. select proper set of parameters to use: perturber_grid to sim_all_parameters
        call perturber_grid2parameters(sim_id, id_pert, perturber_grid, sim_all_parameters)
        ! 2. update from sim_all_parameters to sim_fit_parameters
        call init_param(sim_all_parameters, sim_fit_parameters)

        call par2kel_fit(sim_all_parameters, sim_fit_parameters, m, R, P, sma, ecc, w, mA, inc, lN, tempcheck)
        deallocate (sim_all_parameters, sim_fit_parameters)

!     ! compute number maximum of all transit times of all the planets from the
!     ! integration time (global variable tint)
!     allocate(nT0_perbody(n_bds-1))
!     nT0_perbody = int((1.5_dp*tint)/P(2:n_bds))+2 ! compute the number of transit for each planet ... avoid P(1)<->star
!     nt0_full=sum(nT0_perbody) ! sum them all
!     deallocate(nT0_perbody)
!
!     ! comput the number of maximum rv from integration and write time
!     ! (global variables tint, wrttime
!     nrv_nmax=int(tint/wrttime)+2

        call get_max_nt0_nrv(P, n_bds, nt0_full, nrv_nmax)

        return
    end subroutine wrapper_set_grid_orbelem_fullgrid

    ! ============================================================================

    subroutine wrapper_set_grid_orbelem(sim_id, id_pert, perturber_par, ncol,&
      &m, R, P, sma, ecc, w, mA, inc, lN, n_bds, nt0_full, nrv_nmax)
        integer, intent(in)::sim_id, id_pert
        integer, intent(in)::ncol
        real(dp), dimension(ncol), intent(in)::perturber_par
        integer, intent(in)::n_bds
        real(dp), dimension(n_bds), intent(out)::m, R, P, sma, ecc, w, mA, inc, lN
        integer, intent(out)::nt0_full, nrv_nmax

!f2py    integer,intent(hide),depend(perturber_par)::ncol=len(perturber_par)
! !f2py    integer,intent(hide),depend(m)::n_bds=len(m)

        real(dp), dimension(:), allocatable::sim_all_parameters, sim_fit_parameters
        logical::tempcheck

!     integer,dimension(:),allocatable::nT0_perbody

        tempcheck = .true.
        allocate (sim_all_parameters(npar), sim_fit_parameters(nfit))
        sim_all_parameters = system_parameters
        sim_fit_parameters = zero

        ! 1. select proper set of parameters to use: perturber_grid to sim_all_parameters
!     call perturber_grid2parameters(sim_id,id_pert,perturber_grid,sim_all_parameters)
        call perturber_par2all_par(id_pert, perturber_par, sim_all_parameters)
        ! 2. update from sim_all_parameters to sim_fit_parameters
        call init_param(sim_all_parameters, sim_fit_parameters)

        call par2kel_fit(sim_all_parameters, sim_fit_parameters, m, R, P, sma, ecc, w, mA, inc, lN, tempcheck)
        deallocate (sim_all_parameters, sim_fit_parameters)

!     ! compute number maximum of all transit times of all the planets from the
!     ! integration time (global variable tint)
!     allocate(nT0_perbody(n_bds-1))
!     nT0_perbody = int((1.5_dp*tint)/P(2:n_bds))+2 ! compute the number of transit for each planet ... avoid P(1)<->star
!     nt0_full=sum(nT0_perbody) ! sum them all
!     deallocate(nT0_perbody)
!
!     ! compute the number of maximum rv from integration and write time
!     ! (global variables tint, wrttime
!     nrv_nmax=int(tint/wrttime)+2
        call get_max_nt0_nrv(P, n_bds, nt0_full, nrv_nmax)

        return
    end subroutine wrapper_set_grid_orbelem

    ! ============================================================================

    subroutine wrapper_run_grid_combination(m, R, P, sma, ecc, w, mA, inc, lN, n_bds,&
      &ttra_full, lambda_rm_full, id_ttra_full, stats_ttra, nt0_full,&
      &time_rv_nmax, rv_nmax, stats_rv, nrv_nmax)
        integer, intent(in)::n_bds
        real(dp), dimension(n_bds), intent(in)::m, R, P, sma, ecc, w, mA, inc, lN
        integer, intent(in)::nt0_full

        real(dp), dimension(nt0_full), intent(out)::ttra_full, lambda_rm_full
        integer, dimension(nt0_full), intent(out)::id_ttra_full
        logical, dimension(nt0_full), intent(out)::stats_ttra
        integer, intent(in)::nrv_nmax
        real(dp), dimension(nrv_nmax), intent(out)::time_rv_nmax, rv_nmax
        logical, dimension(nrv_nmax), intent(out)::stats_rv

        real(dp), dimension(nt0_full)::dur_full

        call ode_all_ttra_rv(wrttime, m, R, P, sma, ecc, w, mA, inc, lN,&
          &ttra_full, dur_full, lambda_rm_full, id_ttra_full, stats_ttra,&
          &time_rv_nmax, rv_nmax, stats_rv)

        return
    end subroutine wrapper_run_grid_combination

    ! ============================================================================

    subroutine fit_par_to_ttra_rv(fit_parameters, n_fit,&
      &ttra_full, dur_full, lambda_rm_full, id_ttra_full, stats_ttra, nt0_full,&
      &time_rv_nmax, rv_nmax, stats_rv, nrv_nmax)
        ! Input
        integer, intent(in)::n_fit
        real(dp), dimension(n_fit), intent(in)::fit_parameters
        integer, intent(in)::nt0_full
        integer, intent(in)::nrv_nmax
        ! Output
        real(dp), dimension(nt0_full), intent(out)::ttra_full, dur_full, lambda_rm_full
        integer, dimension(nt0_full), intent(out)::id_ttra_full
        logical, dimension(nt0_full), intent(out)::stats_ttra
        real(dp), dimension(nrv_nmax), intent(out)::time_rv_nmax, rv_nmax
        logical, dimension(nrv_nmax), intent(out)::stats_rv
        ! Local
        real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN
        logical::checkpar = .true.

        integer::ns, ne
        real(dp), dimension(:), allocatable::rv_trend

        allocate (mass(NB), radius(NB), period(NB), sma(NB),&
            &ecc(NB), argp(NB), meanA(NB), inc(NB), longN(NB))
        call only_convert_parameters(system_parameters, fit_parameters,&
            &mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)

        ! write(*,*)" === DEBUG fit_par_to_ttra_rv: max nRV = ",nrv_nmax," nT0s = ",nt0_full
        ! flush(6)
        call ode_all_ttra_rv(wrttime, mass, radius, period, sma, ecc, argp, meanA, inc, longN,&
            &ttra_full, dur_full, lambda_rm_full, id_ttra_full, stats_ttra,&
            &time_rv_nmax, rv_nmax, stats_rv)

        deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)

        ! set proper values of the RV: trend
        if (nRVset .gt. 0) then
            ! rv trend
            if (rv_trend_order .gt. 0) then
                ns = nkel+2*nrvset+1
                ne = n_fit
                allocate (rv_trend(nrv_nmax))
                call compute_rv_trend(time_rv_nmax, nrv_nmax, fit_parameters(ns:ne), rv_trend_order, rv_trend)
                rv_nmax = rv_nmax+rv_trend
                deallocate (rv_trend)
            end if
        end if

        return
    end subroutine fit_par_to_ttra_rv

    ! ============================================================================

    subroutine convert_trades_par_to_kepelem(all_par, n_par, fit_par, n_fit,&
      &n_bds, mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
        integer, intent(in)::n_par, n_fit
        real(dp), dimension(n_par), intent(in)::all_par
        real(dp), dimension(n_fit), intent(in)::fit_par
        integer, intent(in)::n_bds

        real(dp), dimension(n_bds), intent(out)::mass, radius, period, sma, ecc, argp, meanA, inc, longN
        logical, intent(out)::checkpar

        checkpar = .true.
        ! call convert_parameters(all_par, fit_par, mass, radius, period, sma,&
        !     &ecc, argp, meanA, inc, longN, checkpar)
        call only_convert_parameters(all_par, fit_par, mass, radius, period, sma,&
            &ecc, argp, meanA, inc, longN, checkpar)

        return
    end subroutine convert_trades_par_to_kepelem

    ! ============================================================================

    ! adding RV trend if needed --> re-call addRVtrend subroutine...
    subroutine compute_rv_trend(time, ntime, coeff, ncoeff, rv_trend)
        integer, intent(in)::ntime, ncoeff
        real(dp), dimension(ntime), intent(in)::time
        real(dp), dimension(ncoeff), intent(in)::coeff

        real(dp), dimension(ntime), intent(out)::rv_trend

        ! call addRVtrend(time,tepoch,coeff,rv_trend)
        call addRVtrend(time, zero, coeff, rv_trend)

        return
    end subroutine compute_rv_trend

    ! ============================================================================

    ! subroutine to update the system (all_par) and fitting (par) parameters
    ! from the Keplerian elements
    subroutine update_parameters_from_keplerian(mass, radius, period, ecc, argp, meanA, inc, longN, nbodies,&
        &allpar, nallpar, par, npar)
        integer, intent(in)::nbodies, nallpar, npar
        real(dp), dimension(nbodies), intent(in)::mass, radius
        real(dp), dimension(nbodies), intent(in)::period, ecc, argp, meanA, inc, longN
        real(dp), dimension(nallpar), intent(out)::allpar
        real(dp), dimension(npar), intent(out)::par

        real(dp), dimension(:), allocatable:: allpar_tmp, par_tmp

        call set_all_param(mass, radius, period, ecc, argp, meanA, inc, longN, allpar_tmp)
        allpar = allpar_tmp
        call init_param(allpar_tmp, par_tmp)
        par = par_tmp

        deallocate (allpar_tmp, par_tmp)

        return
    end subroutine update_parameters_from_keplerian

    ! ============================================================================
    subroutine astrocentric_to_barycentric_orbits(mass, astro_orbits, bary_orbits, barycenter_position, nbodies, n_steps)
        ! Input
        integer, intent(in)::nbodies, n_steps
        real(dp), dimension(nbodies), intent(in)::mass
        real(dp), dimension(n_steps, nbodies*6), intent(in)::astro_orbits
        ! Output
        real(dp), dimension(n_steps, nbodies*6), intent(out)::bary_orbits
        real(dp), dimension(n_steps, 6), intent(out)::barycenter_position
        ! Locals
        integer::i_step

        do i_step = 1, n_steps
            call barycenter(mass, astro_orbits(i_step, :), barycenter_position(i_step, :), bary_orbits(i_step, :))
        end do

        return
    end subroutine astrocentric_to_barycentric_orbits

    ! ============================================================================

    subroutine linear_fit_no_errors(nx, x, y, m, err_m, q, err_q)
        ! Input
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in)::x, y
        ! Output
        real(dp), intent(out):: m, err_m, q, err_q
        ! Locals

        call linfit(x, y, m, err_m, q, err_q)

        return
    end subroutine linear_fit_no_errors

    subroutine linear_fit_errors(nx, x, y, ey, m, err_m, q, err_q)
        ! Input
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in)::x, y, ey
        ! Output
        real(dp), intent(out):: m, err_m, q, err_q
        ! Locals

        call linfit(x, y, ey, m, err_m, q, err_q)

        return
    end subroutine linear_fit_errors

    ! ============================================================================

    ! add celestial mechanics subroutines to be used from python
    subroutine f90_period_to_sma(mass_star, mass_planet, period, sma)
        ! Input
        real(dp), intent(in)::mass_star, mass_planet, period
        ! Output
        real(dp), intent(out)::sma

        call period_to_sma(mass_star, mass_planet, period, sma)

        return
    end subroutine f90_period_to_sma

    subroutine f90_period_to_sma_vec(n, mass_star, mass_planet, period, sma)
        ! Input
        integer, intent(in)::n
        real(dp), dimension(n), intent(in)::mass_star, mass_planet, period
        ! Output
        real(dp), dimension(n), intent(out)::sma
        ! Local
        integer::i

        do i = 1, n
            call period_to_sma(mass_star(i), mass_planet(i), period(i), sma(i))
        end do

        return
    end subroutine f90_period_to_sma_vec

    subroutine f90_sma_to_period(mass_star, mass_planet, sma, period)
        ! Input
        real(dp), intent(in)::mass_star, mass_planet, sma
        ! Output
        real(dp), intent(out)::period

        call sma_to_period(mass_star, mass_planet, sma, period)

        return
    end subroutine f90_sma_to_period

    subroutine f90_sma_to_period_vec(n, mass_star, mass_planet, sma, period)
        ! Input
        integer, intent(in)::n
        real(dp), dimension(n), intent(in)::mass_star, mass_planet, sma
        ! Output
        real(dp), dimension(n), intent(out)::period
        ! Local
        integer::i

        do i=1,n
            call sma_to_period(mass_star(i), mass_planet(i), sma(i), period(i))
        end do

        return
    end subroutine f90_sma_to_period_vec

    ! ============================================================================

    subroutine angular_momentum_deficit(n_body, masses, eccs, incs, periods, lambda_bodies, amd_bodies, amd)
        ! Input
        integer,intent(in)::n_body
        real(dp), dimension(n_body), intent(in)::masses, eccs, incs, periods
        ! Output
        real(dp), dimension(n_body), intent(out)::lambda_bodies, amd_bodies
        real(dp), intent(out)::amd

        call amd_system(masses, eccs, incs, periods, lambda_bodies, amd_bodies, amd)

        return
    end subroutine angular_momentum_deficit

    subroutine angular_momentum_deficit_fit_parameters(n_bodies, n_fit, fit_pars, n_all, all_pars, lambda_bodies, amd_bodies, amd, n_pairs, amd_r_pairs, amd_h_pairs, stable)
        ! Input
        integer,intent(in)::n_bodies, n_fit, n_all, n_pairs
        real(dp), dimension(n_fit), intent(in)::fit_pars
        real(dp), dimension(n_all), intent(in)::all_pars
        ! Output
        real(dp), dimension(n_bodies), intent(out)::lambda_bodies, amd_bodies
        real(dp), intent(out)::amd
        real(dp), dimension(n_pairs), intent(out)::amd_r_pairs, amd_h_pairs
        logical, intent(out)::stable
        ! Local
        real(dp), dimension(n_bodies)::mass, radius, period, sma, ecc, argp, meana, inc, longn
        logical::checkpar

        call convert_trades_par_to_kepelem(all_pars, n_all, fit_pars, n_fit,&
            &n_bodies, mass, radius, period, sma, ecc, argp, meana, inc, longn,&
            &checkpar)

        call amd_system(mass, ecc, inc, period, lambda_bodies, amd_bodies, amd)
        call amd_hill_stability(mass, ecc, inc, period, amd_r_pairs, amd_h_pairs, stable)

        return
    end subroutine angular_momentum_deficit_fit_parameters

    subroutine angular_momentum_deficit_posterior(n_bodies, n_post, n_fit, post_fit, n_all, all_pars, lambda_bodies, amd_bodies, amd, n_pairs, amd_r_pairs, amd_h_pairs, stable)
        ! Input
        integer,intent(in)::n_bodies, n_post, n_fit, n_all, n_pairs
        real(dp), dimension(n_post, n_fit), intent(in)::post_fit
        real(dp), dimension(n_all), intent(in)::all_pars
        ! Output
        real(dp), dimension(n_post, n_bodies), intent(out)::lambda_bodies, amd_bodies
        real(dp), dimension(n_post), intent(out)::amd
        real(dp), dimension(n_post, n_pairs), intent(out)::amd_r_pairs, amd_h_pairs
        logical, dimension(n_post), intent(out)::stable
        ! Local
        real(dp), dimension(n_fit)::fit_pars
        real(dp), dimension(n_bodies)::lbd, amd_bd
        real(dp), dimension(n_pairs)::amd_r, amd_h
        integer::i_post

        do i_post=1,n_post
            fit_pars = post_fit(i_post, :)
            call angular_momentum_deficit_fit_parameters(n_bodies, n_fit, fit_pars, n_all, all_pars,&
                &lbd, amd_bd, amd(i_post),&
                &n_pairs, amd_r, amd_h, stable(i_post))
            
            lbd = lambda_bodies(i_post, :)
            amd_bd =  amd_bodies(i_post, :)
            amd_r_pairs(i_post, :) = amd_r
            amd_h_pairs(i_post, :) = amd_h
        end do

        return
    end subroutine angular_momentum_deficit_posterior

    ! ============================================================================
    subroutine check_rv_res_periodogram(n, jd, rv, erv, nbody, periods, gls_check)
        ! Input
        integer, intent(in)::n, nbody
        real(dp), dimension(n), intent(in)::jd, rv, erv
        real(dp), dimension(nbody), intent(in)::periods
        ! Output
        logical, intent(out)::gls_check
        ! Local

        call check_periodogram(jd, rv, erv, periods, gls_check)

        return
    end subroutine check_rv_res_periodogram

    ! ============================================================================
    ! ============================================================================

end module f90trades
