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
  use init_trades
!   use bootstrap
  use transits,only:set_ephem
  use derived_parameters_mod
  use driver
  integer::cpuid
  real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,i,lN ! Keplerian orbital elements
  real(dp),dimension(:),allocatable::fitting_parameters
  integer::sim_id
  
  
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
  
  
  ! ---------------------------------------
  ! run integration and write summary files
  ! ---------------------------------------
  sim_id=1
  
  ! IT INTEGRATES THE ORBITS AND WRITES RESULTS TO SCREEN AND FILES
  call write_summary_nosigma(cpuid,sim_id,0,system_parameters,fitting_parameters)
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
