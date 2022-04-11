! +++++++++++++++++++++++++++
! trades_grid.f90
! +++++++++++++++++++++++++++

! initialise TRADES
! read parameters:
! (e,w)
! create the grid
! integration of the parameters and write summary for each combination of the grid
! if lmon = 1 run the LM on the parameters
! if nboot > 0 run bootstrap

program main_grid
    use constants
    use parameters
    use parameters_conversion
    use init_trades
!   use bootstrap
    use linear_ephem
    use derived_parameters_mod
    use driver
    integer::cpuid
    real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meana, inc, longn ! Keplerian orbital elements
    real(dp), dimension(:), allocatable::fitting_parameters
!   integer::sim_id

    ! ------------------
    ! initialises TRADES
    ! ------------------

    ! IT INITIALIZES THE CPU AND NUMBER OF CPU ... USEFUL FOR OPENMP COMPUTATION
    cpuid = 1
    ncpu = 1
    call initu(nfiles, ncpu) ! prepares writing file units
    ! nfiles,ncpu in module 'parameters' and 'constants'

    ! IT READS THE COMMAND ARGUMENT THAT DEFINE THE PATH OF THE FOLDER WITH THE FILES
    call read_com_arg(path)
    ! path in module 'parameters'

    ! IT CALLS ALL THE SUBROUTINES TO READ ALL PARAMETERS AND DATA TO STORE IN COMMON MODULE PARAMETERS
    call read_first(cpuid, mass, radius, period, sma, ecc, argp, meana, inc, longn) ! new version already set system_parameters and boundaries (physical and user-provided)
    ! reset initu with proper ncpu_in argument from arg.in file

    ! IT SETS TWO VECTOR WITH PARAMETERS (ALL THE PARAMETERS AND ONLY THE FITTED ONES)
!   call set_par(m,R,P,a,e,w,mA,i,lN,system_parameters,fitting_parameters)
    call init_param(system_parameters, fitting_parameters)
    ! system_parameters in module 'parameters'

    ! check if there are derived parameters to compute and to check
    call init_derived_parameters(cpuid, path)

    write (*, '(a)') " ID Parameters to fit:"
    write (*, '(a)') trim(paridlist)
    write (*, '(a)') trim(sig_paridlist)
    write (*, *)

    flush (6)

    ! ---------------------------------------
    ! run the grid driver
    ! ---------------------------------------
!   call run_grid(m(1),system_parameters)
    call run_grid(system_parameters)

    ! deallocate all the variables

    if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meana, inc, longn)
    if (allocated(fitting_parameters)) deallocate (fitting_parameters)
    call deallocate_all()

end program
