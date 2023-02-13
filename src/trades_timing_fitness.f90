! +++++++++++++++++++++++++++
! trades_timing_fitness.f90
! +++++++++++++++++++++++++++

program main
    use constants
    use parameters
    use parameters_conversion
    use init_trades
!   use bootstrap
    use linear_ephem
    use derived_parameters_mod
    use fitness_module
    use driver
    integer::cpuid
    real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN ! Keplerian orbital elements
    real(dp), dimension(:), allocatable::fitting_parameters
    real(dp)::chi_square, reduced_chi_square, lnL, lnprior, ln_const, bic
    real(dp)::elapsed_start, elapsed_end, elapsed_time
    integer::n_sims, i_sim

!   integer::j

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

!   write(*,'(a)')"SYSTEM PARAMETERS 1"
!   write(*,'(1000(1x,es23.16))')system_parameters

    ! IT SETS TWO VECTOR WITH PARAMETERS (ALL THE PARAMETERS AND ONLY THE FITTED ONES)
!   call set_par(m,R,P,a,e,w,mA,inc,lN,system_parameters,fitting_parameters)
    call init_param(system_parameters, fitting_parameters)
!   write(*,'(a)')"SYSTEM PARAMETERS 2"
!   write(*,'(1000(1x,es23.16))')system_parameters
!   write(*,'(a)')"FITTING PARAMETERS"
!   write(*,'(1000(1x,es23.16))')fitting_parameters

    ! check if there are derived parameters to compute and to check
    call init_derived_parameters(cpuid, path)

    ! write (*, '(a)') " ID Parameters to fit:"
    ! write (*, '(a)') trim(paridlist)
    ! write (*, '(a)') trim(sig_paridlist)
    ! write (*, *)
    ! flush (6)

    ! ---------------------------------------
    ! run integration and fitness
    ! ---------------------------------------
    n_sims = 1

    write (*, *)
    write (*, *) " ===================================================="
    write (*, *) " RUNNING ", n_sims, " BASE_FITNESS_FUNCTION"
    write (*, *)
    flush (6)

    call cpu_time(elapsed_start)

    do i_sim = 1, n_sims
        call base_fitness_function(system_parameters, fitting_parameters,&
                &chi_square, reduced_chi_square, lnL, lnprior, ln_const, bic)
    end do
    call cpu_time(elapsed_end)

    write (*, *)
    flush (6)
    elapsed_time = elapsed_end-elapsed_start
    write (*, '(a,i3,a,f6.2,a)') " Elapsed time for ", n_sims, " sim: ", elapsed_time, " seconds"
    write (*, '(a,i3,a,f6.2,a)') " Mean Elapsed time for ", 1, " sim: ", elapsed_time/real(n_sims, dp), " seconds"
    write (*, *) " ===================================================="
    write (*, *)
    flush (6)

    ! deallocate all the variables

    if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)
    if (allocated(fitting_parameters)) deallocate (fitting_parameters)
    call deallocate_all()

end program
