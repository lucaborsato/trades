! +++++++++++++++++++++++++++
! trades_polychord.f90
! +++++++++++++++++++++++++++

! initialise TRADES
! read parameters:
! (ecosw,esinw)
! PolyChord if progtype == 5
! writes posterior and other files (see PolyChord documentation)

program main_polychord
  use constants
  use parameters
  use init_trades
  use bootstrap
  use transits,only:set_ephem
  use derived_parameters_mod
  use driver
  use PolyChord_driver,only:PC_driver
  integer::cpuid
  real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,i,lN ! Keplerian orbital elements
  real(dp),dimension(:),allocatable::fitting_parameters
  real(dp),dimension(5)::polychord_info
  
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
  call read_first(cpuid,m,R,P,a,e,w,mA,i,lN)
  
  ! IT SETS TWO VECTOR WITH PARAMETERS (ALL THE PARAMETERS AND ONLY THE FITTED ONES)
  call set_par(m,R,P,a,e,w,mA,i,lN,system_parameters,fitting_parameters)
  ! system_parameters in module 'parameters'
  
  ! IT SETS THE LINEAR EPHEMERIS FROM T0 DATA
  call set_ephem()
  
  ! check if there are derived parameters to compute and to check
  call init_derived_parameters(cpuid,path)
  
  flush(6)
  
  ! --------------------------------------------
  ! progtype == 5 => PolyChord
  ! run PolyChord and write posterior into files
  ! --------------------------------------------

  call PC_driver(polychord_info)

! deallocate all the variables

  if(allocated(m)) deallocate(m,R,P,a,e,w,mA,i,lN)
  if(allocated(fitting_parameters)) deallocate(fitting_parameters)
  call deallocate_all()
!   if(allocated(bnames))   deallocate(bnames,bfiles)
!   if(allocated(tofit))    deallocate(tofit)
!   if(allocated(e_bounds)) deallocate(e_bounds)
!   if(allocated(jdRV))     deallocate(jdRV,RVobs,eRVobs)
!   if(allocated(epoT0obs)) deallocate(epoT0obs,T0obs,eT0obs)
!   if(allocated(id))       deallocate(id,idall,parid)
!   if(allocated(system_parameters)) deallocate(system_parameters)
!   if(allocated(par_min))  deallocate(par_min,par_max)
!   if(allocated(minpar))   deallocate(minpar,maxpar)
!   if(allocated(k_b))      deallocate(k_b)
!   if(allocated(population)) deallocate(population,population_fitness)
!   if(allocated(derived_names)) deallocate(derived_names,derived_boundaries)
!   if(allocated(pso_best_evolution)) deallocate(pso_best_evolution)

end program

