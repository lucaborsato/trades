! +++++++++++++++++++++++++++
! trades_int_lm_bootstrap.f90
! +++++++++++++++++++++++++++

! initialise TRADES
! read parameters:
! (e,w) if progtype == 0
! (ecosw,esinw) if progtype == 2
! integration of the parameters and write summary
! if lmon = 1 run the LM on the parameters
! if nboot > 0 run bootstrap

program main_integrates
    use constants
    use parameters
    use init_trades
    use ode_run, only: ode_integrates
    integer::cpuid
    real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN ! Keplerian orbital elements
    integer::sim_id

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

    ! ---------------------------------------
    ! run integration and write summary files
    ! ---------------------------------------
    sim_id = 1

    ! IT INTEGRATES THE ORBITS AND WRITES RESULTS TO SCREEN AND FILES
    call ode_integrates(cpuid, sim_id, 0, mass, radius, period, sma, ecc, argp, meanA, inc, longN)

    ! deallocate all the variables
    if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)
    call deallocate_all()

end program
