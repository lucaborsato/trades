module parameters
    use constants
    use custom_type
    implicit none

    ! character(512)::path_0, path
    character(512)::path_0=".", path="."

    ! integration argument with default values
    integer::progtype = 0, NB = 2
    integer::wrtorb = 1, wrtconst = 1, wrtel = 1
    integer::idtra = 1, durcheck = 0, idpert = 2
    integer::rvcheck = 0, rv_trend_order = 0
    integer::lmon = 0
    real(dp)::tstart = zero, tepoch = zero, tint = 365.25_dp,&
      &step_0 = 1.0e-3_dp, wrttime = 0.04167_dp, tol_int = 1.0e-13_dp
    integer::ncpu_in = 1 ! default value of the max number of cpus
    ! for Bootstap
    integer::nboot = 0
    logical::bootstrap_scaling = .false.
    ! stability-like
    logical::close_encounter_check=.true., do_hill_check = .false., amd_hill_check = .false.
    ! check if it has been added a signal in rv close to planetary periods, if so, discard
    logical::rv_res_gls = .false.
    ! FITNESS PARAMETERS
    ! integer::oc_fit = 0

    integer::NBDIM
    integer, parameter::DIMMAX = 500000

    ! names and files of the bodies
    character(128), dimension(:), allocatable::bnames, bfiles

    ! fitting variable
    integer, dimension(:), allocatable::tofit, idall, id
    character(10), dimension(:), allocatable::parid
    character(512)::paridlist, sig_paridlist
    character(10), dimension(:), allocatable::all_names_list
    character(1024)::all_names_str
    ! ---
    character(10), dimension(:), allocatable::mass_id, radius_id, period_id, sma_id
    character(10), dimension(:), allocatable::ecc_id, argp_id, meana_id, inc_id, longn_id

    ! ---
    integer::npar, nkel, nfit
    integer::nRVset = 0

    real(dp)::bic_const = zero

    ! complete list of all the parameters for the whole system: M1,R1,M2,P2,a2,e2,...and so on
    real(dp), dimension(:), allocatable::system_parameters ! needed by ode_lm -> dimension == npar

    ! single maximum values for the weighted residuals
    real(dp), parameter::resmax = 1.0e300_dp

    ! save initial parameters of the star
    real(dp), dimension(2, 2)::MR_star

    ! for LM
    integer::maxfev, nprint
    real(dp)::epsfcn, ftol, gtol, xtol
    real(dp), dimension(:, :), allocatable::lmtols

    type(dataObs)::obsData

    logical, dimension(:), allocatable::do_transit ! it needs to determine if a planet should transit or not

    ! to read excluded observations 
    integer::n_excluded = 0
    integer, dimension(:), allocatable::excluded_body
    real(dp), dimension(:,:), allocatable::excluded_time_ranges

    ! loglikelihood constant: - 1/2 * dof * ln(2 * pi) - 1/2 sum(ln(sigma**2)
    real(dp)::ln_err_const

    ! dynamical parameters
    real(dp)::amin, amax

    ! grid
    type(parameter_grid), dimension(10)::perturber_parameters_grid ! global grid parameter: M R P a e w mA tau i lN

    ! for PIKAIA
    real(dp)::pik_ctrl(12)
    integer::seed_pik, npop, ngen

    ! for PSO
    integer::seed_pso, np_pso, nit_pso, wrt_pso
    real(dp)::inertia_in = 0.9_dp, self_in = two, swarm_in = two
    real(dp)::randsearch_in = 1.0e-5_dp, vmax_in = half, vrand_in = 0.07_dp

    ! boundaries
    real(dp), dimension(:), allocatable::par_min, par_max ! dimension: system_parameters
    real(dp), dimension(:), allocatable::minpar, maxpar ! dimension: fitting parameters
    real(dp), dimension(:), allocatable::minpar_bck, maxpar_bck ! backup


    ! for PIK/PSO
    integer::wrtAll, nGlobal

    real(dp), dimension(:, :, :), allocatable::population
    real(dp), dimension(:, :), allocatable::population_fitness
    real(dp), dimension(:, :), allocatable::pso_best_evolution

    ! other boundaries
    real(dp), dimension(:, :), allocatable::e_bounds

    ! priors
    integer::n_priors = 0
    character(10), dimension(:), allocatable::priors_names
    real(dp), dimension(:, :), allocatable::priors_values

    ! derived parameters
    integer::n_derived = 0
    integer::secondary_parameters = 0
    logical::check_derived = .false.
    logical::fix_derived = .false.
    character(10), dimension(:), allocatable::derived_names
    real(dp), dimension(:, :), allocatable::derived_boundaries

contains

! ==============================================================================

    ! deallocate all variables in 'parameters' module
    subroutine deallocate_all()

        if (allocated(bnames)) deallocate (bnames)
        if (allocated(bfiles)) deallocate (bfiles)

        if (allocated(tofit)) deallocate (tofit)
        if (allocated(id)) deallocate (id)
        if (allocated(idall)) deallocate (idall)
        if (allocated(parid)) deallocate (parid, all_names_list)
        if (allocated(all_names_list)) deallocate (all_names_list)

        if (allocated(mass_id)) deallocate (mass_id, radius_id, period_id, sma_id)
        if (allocated(ecc_id)) deallocate (ecc_id, argp_id, meana_id, inc_id, longn_id)

        if (allocated(system_parameters)) deallocate (system_parameters)

        if (allocated(lmtols)) deallocate (lmtols)

        call deallocate_dataObs(obsData)
        if (allocated(do_transit)) deallocate (do_transit)

        ! to check if needed to deallocate something within type perturber_parameters_grid

        if (allocated(par_min)) deallocate (par_min, par_max)
        if (allocated(minpar)) deallocate (minpar, maxpar)
        if (allocated(minpar_bck)) deallocate (minpar_bck, maxpar_bck)

        if (allocated(population)) deallocate (population)
        if (allocated(population_fitness)) deallocate (population_fitness)
        if (allocated(pso_best_evolution)) deallocate (pso_best_evolution)

        if (allocated(e_bounds)) deallocate (e_bounds)

        if (allocated(priors_names)) deallocate (priors_names)
        if (allocated(priors_values)) deallocate (priors_values)

        if (allocated(derived_names)) deallocate (derived_names)
        if (allocated(derived_boundaries)) deallocate (derived_boundaries)

        return
    end subroutine deallocate_all

end module parameters
