module parameters
  use constants,only:dp,TOLERANCE,zero,one,parameter_grid
  use settings_module,only:program_settings
  use priors_module,only:prior
  implicit none

  character(512)::path
  
  ! integration argument with default values
  integer::progtype=0,NB=2
  integer::wrtorb=1,wrtconst=1,wrtel=1
  integer::idtra=1,durcheck=0,rvcheck=0,idpert=2,lmon=0
  real(dp)::tstart=zero,tepoch=zero,tint=365.25_dp,&
    &step_0=1.e-3_dp,wrttime=0.04167_dp,tol_int=1.e-13_dp
  ! for Bootstap
  integer::nboot=0
  logical::bootstrap_scaling=.false.
  ! FITNESS PARAMETERS
  real(dp)::k_chi2r = one, k_chi2wr = zero
  real(dp)::k_a
  real(dp),dimension(:),allocatable::k_b
  
  integer::NBDIM
  integer,parameter::DIMMAX=100000

  ! names and files of the bodies
  character(128),dimension(:),allocatable::bnames,bfiles

  ! fitting variable
  integer,dimension(:),allocatable::tofit,idall,id
  character(10),dimension(:),allocatable::parid
  character(128)::paridlist,sig_paridlist
  character(10),dimension(:),allocatable::all_names_list
  character(1024)::all_names_str
  integer::ndata,npar,nfit,dof
  real(dp)::inv_dof

  ! complete list of all the parameters for the whole system: M1,R1,M2,P2,a2,e2,...and so on
  real(dp),dimension(:),allocatable::system_parameters ! needed by ode_lm -> dimension == npar
  
  ! single maximum values for the weighted residuals
!   real(dp),parameter::resmax=1.e20_dp
  real(dp),parameter::resmax=1.e10_dp

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

  ! ephemeris: Tephem, Pephem
  real(dp),dimension(:),allocatable::Tephem,eTephem,Pephem,ePephem

  ! dynamical parameters
  real(dp)::amin,amax

  ! grid
  type (parameter_grid),dimension(10)::perturber_parameters_grid ! global grid parameter: M R P a e w mA tau i lN
  
  ! for PIKAIA
  real(dp)::ctrl(12)
  integer::seed_pik,npop,ngen
  real(dp),dimension(:),allocatable::par_min,par_max

  ! for PSO
  integer::seed_pso,np_pso,nit_pso,wrt_pso
  real(dp),dimension(:),allocatable::minpar,maxpar

  ! for PIK/PSO
  integer::wrtAll,nGlobal
  
  real(dp),dimension(:,:,:),allocatable::population
  real(dp),dimension(:,:),allocatable::population_fitness
  real(dp),dimension(:,:),allocatable::pso_best_evolution

  
  ! for PolyChord
  type(program_settings)::settings
  
  ! other boundaries
  real(dp),dimension(:,:),allocatable::e_bounds
  
  ! derived parameters
  integer::n_derived=0
  integer::secondary_parameters=0
  logical::check_derived=.false.
  logical::fix_derived=.false.
  character(10),dimension(:),allocatable::derived_names
  real(dp),dimension(:,:),allocatable::derived_boundaries

  
  
  interface norm2par
    module procedure norm2par_1,norm2par_2
  end interface norm2par
  
  contains

  ! ------------------------------------------------------------------ !
  ! given the computed parameters with the L-M it adjusts some parameters
  ! i.e. setting all the angles between 0 and 360 (no < 0 angles) etc.
  subroutine param_adj(par,sigpar)
    real(dp),dimension(:),intent(inout)::par,sigpar
    integer::j1,j2
    real(dp),parameter::circ=360._dp,hcirc=180._dp

    do j1=1,nfit
      j2=id(j1)
!       if( (j2.eq.7) .or. (j2.eq.8) .or. (j2.eq.10) )then
      if(j2.eq.10)then
        par(j1)=mod(par(j1),circ)
        if(par(j1).lt.0._dp) par(j1)=par(j1)+circ
        sigpar(j1)=mod(sigpar(j1),circ)
      else if(j2.eq.9)then
        par(j1)=mod(par(j1),hcirc)
        if(par(j1).lt.0._dp) par(j1)=par(j1)+hcirc
        sigpar(j1)=mod(sigpar(j1),hcirc)
      end if
    end do

    return
  end subroutine param_adj
  ! ------------------------------------------------------------------ !
  
  ! ------------------------------------------------------------------ !
  ! conversion from parameter boundaries [0, 1] --> [ParMin, ParMax]
  subroutine norm2par_1(norm,par)
    real(dp),dimension(:),intent(in)::norm
    real(dp),dimension(:),intent(out)::par
    real(dp)::dpar
    integer::j,j1

    j1=0
    do j=1,npar
      if(tofit(j).eq.1)then
        j1=j1+1
        dpar=abs(par_max(j)-par_min(j))
        par(j1)=par_min(j)+dpar*norm(j1)
      end if
    end do

    return
  end subroutine norm2par_1

  ! conversion from parameter boundaries [0, 1] --> [ParMin, ParMax]
  subroutine norm2par_2(norm,par,allpar)
    real(dp),dimension(:),intent(in)::norm
    real(dp),dimension(:),intent(out)::par
!     real(dp),dimension(:),intent(out)::allpar
    real(dp),dimension(:)::allpar
    real(dp)::dpar
    integer::j,j1

    j1=0
    do j=1,npar
      if(tofit(j).eq.1)then
        j1=j1+1
        dpar=abs(par_max(j)-par_min(j))
        par(j1)=par_min(j)+dpar*norm(j1)
        allpar(j)=par(j1)
      end if
    end do

    return
  end subroutine norm2par_2

  ! deallocate all variables in 'parameters' module
  subroutine deallocate_all()
  
    if(allocated(bnames))   deallocate(bnames,bfiles)
    if(allocated(tofit))    deallocate(tofit)
    if(allocated(e_bounds)) deallocate(e_bounds)
    if(allocated(jdRV))     deallocate(jdRV,RVobs,eRVobs)
    if(allocated(epoT0obs)) deallocate(epoT0obs,T0obs,eT0obs)
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
