module parameters
!   use constants,only:dp,TOLERANCE,zero,one,parameter_grid
  use constants
!   use settings_module,only:program_settings ! from POLYCHORD OLD VERSION
!   use priors_module,only:prior
  implicit none

  character(512)::path_0,path
  
  ! integration argument with default values
  integer::progtype=0,NB=2
  integer::wrtorb=1,wrtconst=1,wrtel=1
  integer::idtra=1,durcheck=0,rvcheck=0,idpert=2,lmon=0
  real(dp)::tstart=zero,tepoch=zero,tint=365.25_dp,&
    &step_0=1.e-3_dp,wrttime=0.04167_dp,tol_int=1.e-13_dp
  integer::ncpu_in=1 ! default value of the max number of cpus
  ! for Bootstap
  integer::nboot=0
  logical::bootstrap_scaling=.false.,do_hill_check=.false.
  ! FITNESS PARAMETERS
!   logical::oc_fit=.false.
  integer::oc_fit=0
  real(dp)::k_chi2r=one, k_chi2wr=zero
  real(dp)::k_a
  real(dp),dimension(:),allocatable::k_b
  
  integer::NBDIM
  integer,parameter::DIMMAX=500000

  ! names and files of the bodies
  character(128),dimension(:),allocatable::bnames,bfiles

  ! fitting variable
  integer,dimension(:),allocatable::tofit,idall,id
  character(10),dimension(:),allocatable::parid
  character(128)::paridlist,sig_paridlist
  character(10),dimension(:),allocatable::all_names_list
  character(1024)::all_names_str
  integer::ndata,npar,nfit,nfree,dof
  real(dp)::inv_dof

  ! complete list of all the parameters for the whole system: M1,R1,M2,P2,a2,e2,...and so on
  real(dp),dimension(:),allocatable::system_parameters ! needed by ode_lm -> dimension == npar
  
  ! single maximum values for the weighted residuals
!   real(dp),parameter::resmax=1.e20_dp
  real(dp),parameter::resmax=1.e10_dp

  ! save initial parameters of the star
  real(dp),dimension(2,2)::MR_star
  
  ! for LM
  integer::maxfev,nprint
  real(dp)::epsfcn,ftol,gtol,xtol
  real(dp),dimension(:,:),allocatable::lmtols

  ! RV data
  integer::nRV,nRVset
  real(dp),dimension(:),allocatable::jdRV,RVobs,eRVobs
  integer,dimension(:),allocatable::RVsetID,nRVsingle

  ! T0 data
  integer,dimension(:),allocatable::nT0
  real(dp),dimension(:,:),allocatable::T0obs,eT0obs
  integer,dimension(:,:),allocatable::epoT0obs
  logical,dimension(:),allocatable::do_transit ! it needs to determine if a planet should transit or not

  ! loglikelihood constant: - 1/2 * dof * ln(2 * pi) - 1/2 sum(ln(sigma**2)
  real(dp)::ln_err_const
  
  ! ephemeris: Tephem, Pephem
  real(dp),dimension(:),allocatable::Tephem,eTephem,Pephem,ePephem

  ! dynamical parameters
  real(dp)::amin,amax

  ! grid
  type (parameter_grid),dimension(10)::perturber_parameters_grid ! global grid parameter: M R P a e w mA tau i lN
  
  ! for PIKAIA
  real(dp)::ctrl(12)
  integer::seed_pik,npop,ngen

  ! for PSO
  integer::seed_pso,np_pso,nit_pso,wrt_pso
  real(dp)::inertia_in=0.9_dp,self_in=two,swarm_in=two
  real(dp)::randsearch_in=1.e-5_dp,vmax_in=half,vrand_in=0.07_dp
  
  ! boundaries
  real(dp),dimension(:),allocatable::par_min,par_max ! dimension: system_parameters
  real(dp),dimension(:),allocatable::minpar,maxpar ! dimension: fitting parameters

  ! for PIK/PSO
  integer::wrtAll,nGlobal
  
  real(dp),dimension(:,:,:),allocatable::population
  real(dp),dimension(:,:),allocatable::population_fitness
  real(dp),dimension(:,:),allocatable::pso_best_evolution

  
!   ! for PolyChord
!   type(program_settings)::settings
  
  ! other boundaries
  real(dp),dimension(:,:),allocatable::e_bounds
  
  ! derived parameters
  integer::n_derived=0
  integer::secondary_parameters=0
  logical::check_derived=.false.
  logical::fix_derived=.false.
  character(10),dimension(:),allocatable::derived_names
  real(dp),dimension(:,:),allocatable::derived_boundaries

  contains

  ! deallocate all variables in 'parameters' module
  subroutine deallocate_all()
  
    if(allocated(bnames))   deallocate(bnames,bfiles)
    if(allocated(tofit))    deallocate(tofit)
    if(allocated(e_bounds)) deallocate(e_bounds)
    if(allocated(jdRV))     deallocate(jdRV,RVobs,eRVobs)
    if(allocated(epoT0obs)) deallocate(epoT0obs,T0obs,eT0obs,do_transit)
    if(allocated(nT0))      deallocate(nT0)
    if(allocated(id))       deallocate(id,idall,parid,all_names_list)
    if(allocated(system_parameters)) deallocate(system_parameters)
    if(allocated(par_min))  deallocate(par_min,par_max)
    if(allocated(minpar))   deallocate(minpar,maxpar)
    if(allocated(k_b))      deallocate(k_b)
    if(allocated(population)) deallocate(population,population_fitness)
    if(allocated(derived_names)) deallocate(derived_names,derived_boundaries)
    if(allocated(pso_best_evolution)) deallocate(pso_best_evolution)
  
    return
  end subroutine deallocate_all
  
end module parameters
