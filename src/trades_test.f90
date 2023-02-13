! +++++++++++++++++++++++++++
! trades_int_lm_bootstrap.f90
! +++++++++++++++++++++++++++

! initialise TRADES
! read parameters:
! (e,w) if progtype == 0
! (ecosw,esinw) if progtype == 2
! integration of the orbits and computes loglikelihood
! repeat n_sim times and print the mean time!

program main_test
    use constants
    use parameters
    use init_trades
    ! use ode_run, only: ode_integrates
    use pytrades
    integer::cpuid
    real(dp), dimension(:), allocatable::mass, radius
    real(dp), dimension(:), allocatable::period, sma, ecc, argp, meanA, inc, longN ! Keplerian orbital elements
    real(dp), dimension(:), allocatable::fit_par
    integer::sim_id
    logical::check
    real(dp), dimension(:), allocatable::elapsed_times
    real(dp)::lnL, lnp, mean_ela, std_ela, elapsed_start, elapsed_end

    ! ------------------
    ! initialises TRADES
    ! ------------------

    ! IT INITIALIZES THE CPU AND NUMBER OF CPU ... USEFUL FOR OPENMP COMPUTATION
    cpuid = 1
    ncpu = 1
    call initu(nfiles, ncpu) ! prepares writing file units
    ! nfiles,ncpu in module 'parameters'

    ! IT READS THE COMMAND ARGUMENT THAT DEFINE THE PATH OF THE FOLDER WITH THE FILES
    call read_com_arg(path)
    ! path in module 'parameters'

    ! IT CALLS ALL THE SUBROUTINES TO READ ALL PARAMETERS AND DATA TO STORE IN COMMON MODULE PARAMETERS
    call read_first(cpuid, mass, radius, period, sma, ecc, argp, meanA, inc, longN) ! new version already set system_parameters and boundaries (physical and user-provided)
    ! reset initu with proper ncpu_in argument from arg.in file

    ! IT SETS TWO VECTOR WITH PARAMETERS (ALL THE PARAMETERS AND ONLY THE FITTED ONES)
    call init_param(system_parameters, fit_par)

    ! check if there are derived parameters to compute and to check
    call init_derived_parameters(cpuid, path)

    ! ---------------------------------------
    ! run integration and write summary files
    ! ---------------------------------------
    sim_id = 1
    check = .true.

    write (*, *) " TESTING INTEGRATES ORBITS AND COMPUTED OBSDATA ... "
    n_sim = 1
    allocate (elapsed_times(n_sim))
    elapsed_times = zero
    mean_ela = zero
    std_ela = zero

    write (*, *) " CHECKING TIMING FOR n_sim = ", n_sim

    write (*, *)
    write (*, *) " fortran_loglikelihood"
    do i_sim = 1, n_sim
        write (*, *)
        call cpu_time(elapsed_start)
        call fortran_loglikelihood(fit_par, lnL, check, nfit)
        call cpu_time(elapsed_end)
        elapsed_times(i_sim) = elapsed_end-elapsed_start
        write (*, "(a,es23.16, a, es23.16)") "lnL = ", lnL, " lnp = ", lnp
        write (*, '(a,i4,a,f6.2,a)') " Elapsed time for sim number ", i_sim, " : ", elapsed_times(i_sim), " seconds"
        write (*, *)
    end do
    write (*, '(a,i4,a,f6.2,a)') " Elapsed time for ", n_sim, " simulations: ", sum(elapsed_times), " seconds"
    mean_ela = sum(elapsed_times)/real(n_sim, dp)
    if (n_sim .gt. 1) then
        std_ela = sqrt(sum((elapsed_times-mean_ela)**2)/real(n_sim-1, dp))/sqrt(real(n_sim, dp))
    else
        std_ela = zero
    end if
    write (*, '(a,f6.2,a,f6.2,a)') " Mean Elapsed time for 1 simulation:     ", mean_ela, "+/-", std_ela, " seconds"
    write (*, *)

    write (*, *)
    write (*, *) " fortran_logprob"
    do i_sim = 1, n_sim
        write (*, *)
        call cpu_time(elapsed_start)
        call fortran_logprob(fit_par, lnL, lnp, check, nfit)
        call cpu_time(elapsed_end)
        elapsed_times(i_sim) = elapsed_end-elapsed_start
        write (*, "(a,es23.16, a, es23.16)") "lnL = ", lnL, " lnp = ", lnp
        write (*, '(a,i4,a,f6.2,a)') " Elapsed time for sim number ", i_sim, " : ", elapsed_times(i_sim), " seconds"
        write (*, *)
    end do
    write (*, '(a,i4,a,f6.2,a)') " Elapsed time for ", n_sim, " simulations: ", sum(elapsed_times), " seconds"
    mean_ela = sum(elapsed_times)/real(n_sim, dp)
    if (n_sim .gt. 1) then
        std_ela = sqrt(sum((elapsed_times-mean_ela)**2)/real(n_sim-1, dp))/sqrt(real(n_sim, dp))
    else
        std_ela = zero
    end if
    write (*, '(a,f6.2,a,f6.2,a)') " Mean Elapsed time for 1 simulation:     ", mean_ela, "+/-", std_ela, " seconds"
    write (*, *)

    write (*, *)
    write (*, *) " ln_priors"
    deallocate (elapsed_times)
    n_sim = 100000
    allocate (elapsed_times(n_sim))
    elapsed_times = zero
    mean_ela = zero
    std_ela = zero

    do i_sim = 1, n_sim
        ! write (*,*)
        call cpu_time(elapsed_start)
        call ln_priors(system_parameters, fit_par, priors_names, priors_values, lnp)
        call cpu_time(elapsed_end)
        elapsed_times(i_sim) = elapsed_end-elapsed_start
        ! write(*, "(a,es23.16, a, es23.16)")"lnL = ",lnL, " lnp = ", lnp
        ! write (*, '(a,i4,a,f6.2,a)') " Elapsed time for sim number ", i_sim, " : ", elapsed_times(i_sim), " seconds"
        ! write (*,*)
    end do
    write (*, '(a,i8,a,f6.2,a)') " Elapsed time for ", n_sim, " simulations: ", sum(elapsed_times), " seconds"
    mean_ela = sum(elapsed_times)/real(n_sim, dp)
    if (n_sim .gt. 1) then
        std_ela = sqrt(sum((elapsed_times-mean_ela)**2)/real(n_sim-1, dp))/sqrt(real(n_sim, dp))
    else
        std_ela = zero
    end if
    write (*, '(a,es23.16,a,es23.16,a)') " Mean Elapsed time for 1 simulation:     ", mean_ela, "+/-", std_ela, " seconds"
    write (*, *)
    write (*, *) " last lnp = ", lnp
    write (*, *)

    ! deallocate all the variables

    deallocate (elapsed_times)
    if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)
    if (allocated(fit_par)) deallocate (fit_par)
    call deallocate_all()

end program
