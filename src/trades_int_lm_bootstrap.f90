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


program main_int_lm_bootstrap
  use constants
  use parameters
  use parameters_conversion
  use init_trades
!   use bootstrap
  use transits,only:set_ephem
  use derived_parameters_mod
  use driver
  integer::cpuid
  real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN ! Keplerian orbital elements
  real(dp),dimension(:),allocatable::fitting_parameters
  real(dp)::fitness
  integer::sim_id
  
!   integer::j
  
  ! ------------------
  ! initialises TRADES
  ! ------------------
  
  ! IT INITIALIZES THE CPU AND NUMBER OF CPU ... USEFUL FOR OPENMP COMPUTATION
  cpuid=1
  ncpu=1
  call initu(nfiles,ncpu) ! prepares writing file units
  ! nfiles,ncpu in module 'parameters'
  
  ! IT READS THE COMMAND ARGUMENT THAT DEFINE THE PATH OF THE FOLDER WITH THE FILES
  call read_com_arg(path)
  ! path in module 'parameters'
  
  ! IT CALLS ALL THE SUBROUTINES TO READ ALL PARAMETERS AND DATA TO STORE IN COMMON MODULE PARAMETERS
  call read_first(cpuid,m,R,P,a,e,w,mA,inc,lN) ! new versione already set system_parameters and boundaries (physical and user-provided)
  
!   write(*,'(a)')"SYSTEM PARAMETERS 1"
!   write(*,'(1000(1x,es23.16))')system_parameters
  
  ! IT SETS TWO VECTOR WITH PARAMETERS (ALL THE PARAMETERS AND ONLY THE FITTED ONES)
!   call set_par(m,R,P,a,e,w,mA,inc,lN,system_parameters,fitting_parameters)
  call init_param(system_parameters,fitting_parameters)
!   write(*,'(a)')"SYSTEM PARAMETERS 2"
!   write(*,'(1000(1x,es23.16))')system_parameters
!   write(*,'(a)')"FITTING PARAMETERS"
!   write(*,'(1000(1x,es23.16))')fitting_parameters
  
  ! IT SETS THE LINEAR EPHEMERIS FROM T0 DATA
  call set_ephem()
  
  ! check if there are derived parameters to compute and to check
  call init_derived_parameters(cpuid,path)
  
  write(*,'(a)')" ID Parameters to fit:"
  write(*,'(a)')trim(paridlist)
  write(*,'(a)')trim(sig_paridlist)
  write(*,*)
  flush(6)
  
  
  ! ---------------------------------------
  ! run integration and write summary files
  ! ---------------------------------------
  sim_id=1
  
  ! IT INTEGRATES THE ORBITS AND WRITES RESULTS TO SCREEN AND FILES
  call write_summary_nosigma(cpuid,sim_id,0,system_parameters,fitting_parameters,fitness)
  ! if lmon equal 1, it fits with LM and then it writes summary files
  ! fitting_parameters will be updated
  if(lmon.eq.1)then
    call run_and_write_lm(cpuid,sim_id,system_parameters,fitting_parameters)
  end if
  
  
  ! --------------------------
  ! run bootstrap if nboot > 0
  ! --------------------------
  if(nboot.gt.0)then
    call run_bootstrap(sim_id,system_parameters,fitting_parameters)
  end if

! deallocate all the variables

  if(allocated(m)) deallocate(m,R,P,a,e,w,mA,inc,lN)
  if(allocated(fitting_parameters)) deallocate(fitting_parameters)
  call deallocate_all()

end program
