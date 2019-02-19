! TRADES fortran module to be called fron python

module pytrades
  !$ use omp_lib
  use constants
  use parameters
  use parameters_conversion
  use convert_type,only:string
  use init_trades
  use derived_parameters_mod
  use linear_ephem
  use fitness_module
  use ode_run,only:orbits_to_data,ode_all_ttra_rv
  use grid_search
  
  implicit none
  ! exposing variables in parameters to trades_lib
  !f2py integer,parameter::dp=selected_real_kind(8)
  !f2py character(512)::path
  
  !f2py integer::npar,nfit
  
  !f2py real(dp),parameter::resmax=1.e10_dp
  !f2py real(dp)::tepoch,tint
  
  !f2py real(dp),dimension(2,2)::MR_star
  !f2py real(dp),dimension(:),allocatable::system_parameters
  !f2py real(dp),dimension(:),allocatable::par_min,par_max ! dimension: system_parameters
  
  
  !f2py real(dp),dimension(:,:,:),allocatable::population
  !f2py real(dp),dimension(:,:),allocatable::population_fitness
  !f2py real(dp),dimension(:,:),allocatable::pso_best_evolution
  !f2py integer::seed_pso,np_pso,nit_pso,wrt_pso

  !f2py integer::ncpu_in
  
  ! needed because TRADES now uses these variables in data type
  integer::ndata,nfree,dof
  integer::nRV,nRVset
  integer,dimension(:),allocatable::nT0
  real(dp),dimension(:),allocatable::Pephem
  integer::nTTs,nDurs
  
!   real(dp)::ln_err_const,inv_dof
  real(dp)::inv_dof
  !f2py real(dp)::ln_err_const
  
  ! variables:  parameters to fit
  real(dp),dimension(:),allocatable::fitting_parameters
  real(dp),dimension(:,:),allocatable::parameters_minmax
!   character(10),dimension(:),allocatable::parameter_names
  integer::str_len
  integer::n_global,n_bodies
  
  contains
  
  ! ============================================================================
  
  subroutine initialize_trades(path_in, sub_folder, n_threads_in)
    !f2py real(dp),dimension(:),allocatable::eRVobs
    !f2py real(dp),dimension(:,:),allocatable::eT0obs
    character*(*),intent(in)::path_in, sub_folder
    integer,intent(in)::n_threads_in
    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
    integer,dimension(:),allocatable::nset
!     character(80)::fmt

    integer::inb
    
    ! subroutine: initu -> init_trades
    ! variables:  nfiles -> constants (default = 90)
    !             ncpu = 1
    call initu(nfiles, 1)
    
    ! IT READS THE COMMAND ARGUMENT THAT DEFINE THE PATH OF THE FOLDER WITH THE FILES
    ! subroutine: read_com_arg -> init_trades
    ! variables:  path -> parameters (default = empty)
    !call read_com_arg(path)
    path_0=trim(adjustl(path_in))
    path=trim(adjustl(path_in))
!     write(*,'(a,a)')" READING IN PATH = ",trim(path)
    
    ! IT DEFINES THE STRING TO WRITE THE REAL WITH RIGHT DECIMAL: PRECISION
    ! variables:  sprec -> constants
!     sprec="g27.15"
    sprec='es23.16'
    
    ! ---
    ! IT CALLS ALL THE SUBROUTINES TO READ ALL PARAMETERS AND DATA TO STORE IN COMMON MODULE PARAMETERS
    ! subroutine: read_first -> init_trades
    ! by scratch below:
    ! ---
    ! IT READS THE ARGUMENTS OF INTEGRATION AND STORE IN COMMON MODULE PARAMETERS.
    ! THE VARIBLES WILL NOT BE MODIFIED FURTHERMORE.
    ! subroutine: read_arg -> init_trades
    ! variables:  cpuid = 1
    call read_arg(1)
    
!     ! if executed in parallel it uses the same cpus defined in the python script
!     !$ ncpu_in=n_threads_in

    !$ call omp_set_num_threads(ncpu_in)
    !$ call initu(nfiles, ncpu_in)
    
!     write(*,*)' in pytrades_lib ncpu_in = ',ncpu_in
    
!     write(*,'(a,a,a)')" READ ",trim(path)//"arg.in"
!     progtype=6 ! needed for other subroutines
    
    ! variables:  e_bounds -> parameters
    !             NB -> parameters (default =2, updated with read_arg)
    n_bodies=NB ! needed to be used by python wrapper ... to check if I can avoid it
    allocate(e_bounds(2,NB))
    e_bounds(1,:)=TOL_dp
    e_bounds(2,:)=1._dp-TOL_dp
    ! IT READS THE FILES AND THE NAMES OF THE BODIES AND DETERMINES THE PARAMETERS TO BE FITTED
    ! subroutine: read_list -> init_trades
    ! variables:  cpuid = 1
    call read_list(1)
!     write(*,'(a,a,a)')" READ ",trim(path)//"bodies.lst"
    ! IT DEFINES THE ID OF THE PARAMETERS TO BE FITTED
    ! subroutine: idpar/idpar_fit -> init_trades
!     if(progtype.le.1)then
!       call idpar() ! IT DEFINES THE ID OF THE PARAMETERS TO BE FITTED
!     else
!       call idpar_fit()
!     end if
    call idpar()
    
    ! IT READS THE PARAMETERS FROM THE FILES
    ! subroutine: read_par -> init_trades
    ! variables:  cpuid=1
    !             m,R,P,a,e,w,mA,i,lN
!     call read_par(1,m,R,P,a,e,w,mA,i,lN)
!     ! IT READS BOUNDARIES OF THE KEPLERIAN ELEMENTS
!     ! subroutine: read_par_boundaries -> init_trades
!     ! variables:  cpuid=1
!     !             m,R
!     call read_par_boundaries(1,m,R) ! it sets minpar(nfit) and maxpar(nfit)
    call read_fullpar(1,m,R,P,a,e,w,mA,inc,lN,system_parameters)
    
    ! IT READS RV DATA
    ! variables:  nRV -> parameters
    nRV=0
    ! subroutine: read_RVobs -> init_trades
    ! variables:  cpuid=1
    call read_RVobs(1)
    ! function: string -> 
    ! variables:  rvcheck -> parameters
    !             nRVset -> parameters
    !             nRVsingle -> parameters
!     fmt=adjustl("(a,i4,a,i4,a,"//trim(string(nRVset))//"i4))")
!     if(rvcheck.eq.1) write(*,trim(fmt))" RV DATA: nRV = ",nRV,&
!       &" in ",nRVset," set of RV: ",nRVsingle
    nRV=obsData%obsRV%nRV
    nRVset=obsData%obsRV%nRVset

    ! IT READS T0 DATA
    ! allocate number of derived dataT0 for each body in the dataObs type
    allocate(obsData%obsT0(NB-1),nT0(NB-1))
    nT0=0
    ! subroutine: cpuid=1
    call read_T0obs(1)
    ! variables:  idtra -> parameters
    !             nT0 -> parameters
!     if(idtra.ne.0) write(*,'(a,1000(i5,1x))') " T0 DATA: nT0 = ",nT0(2:)

    nT0=obsData%obsT0(:)%nT0
    
    ! IT SETS THE LINEAR EPHEMERIS FROM T0 DATA
    if(obsData%nTTs.gt.0)then
      call set_ephem()
      call compute_oc(obsData%obsT0)
      allocate(Pephem(NB))
      Pephem=zero
      do inb=2,NB
        if(obsData%obsT0(inb-1)%nT0.gt.0) Pephem(inb)=obsData%obsT0(inb-1)%Pephem
      end do
    end if
      
    ! IT DETERMINES THE NDATA
    ! variables:  ndata -> parameters
    !             dof -> parametershttp://www.r-bloggers.com/wilcoxon-signed-rank-test/
    !             nfit -> parameters
    !             inv_dof -> parameters
    !             one -> constants
    nTTs=obsData%nTTs
    nDurs = obsData%nDurs
    
    obsData%ndata=nRV+nTTs+nDurs
    obsData%nfree=obsData%obsRV%nRVset
    obsData%dof=(obsData%ndata-nfit-obsData%nfree)
    
    ndata=obsData%ndata
    nfree=obsData%nfree
    dof=obsData%dof

    if(dof.le.0)then
      write(*,'(a,a)')' FOUND dof <= 0 SO IT IS FORCED TO 1 IN CASE',&
         &' THE USER WANT TO SIMULATE/INTEGRATE AND NOT CHECK THE FIT.'
      obsData%dof=1
      dof=obsData%dof
    end if
    
    obsData%inv_dof = one / real(obsData%dof,dp)
    inv_dof=obsData%inv_dof

    ! write(*,*)' ln_err_const'
    ! IT DETERMINES THE LN_ERR_CONST TO COMPUTE LOGLIKELIHOOD
!     ln_err_const = get_ln_err_const(eRVobs,eT0obs)
    ln_err_const = get_lnec(obsData)
    
    ! write(*,*)' set id of parameters'
    ! IT SETS THE LIST OF THE PARAMETERS TO FIT
    ! subroutine: set_parid_list -> init_trades
    call set_parid_list()
    ! IT SETS FITNESS PARAMETERS
    ! variables:  nset -> parameters
    !             k_a -> parameters
    !             k_b -> parameters
    if(nRV.ne.0.and.nTTs.ne.0)then
      if(durcheck.eq.0) then
        allocate(nset(2),k_b(2))
        nset(1)=nRV
        nset(2)=nTTs
      else
        allocate(nset(3),k_b(3))
        nset(1)=nRV
        nset(2)=nTTs
        nset(3)=nDurs
      end if
    else if(nRV.ne.0.and.nTTs.eq.0)then
      allocate(nset(1),k_b(1))
      nset(1)=nRV
    else if(nRV.eq.0.and.nTTs.ne.0)then
      if(durcheck.eq.0) then
        allocate(nset(1),k_b(1))
        nset(1)=nTTs
      else
        allocate(nset(2),k_b(2))
        nset(1)=nTTs
        nset(2)=nDurs
      end if
    else
      allocate(nset(1),k_b(1))
      nset(1)=1
    end if
    ! PARAMETER TO PROPERLY SCALE THE RESIDUALS FOR THE CHI2
    k_a = sqrt(k_chi2r*inv_dof)
    ! PARAMETER TO PROPERLY SCALE THE RESIDUALS FOR THE CHI2_WEIGHTED
!     k_b = sqrt(k_chi2wr/real(dof,dp))*(real(ndata,dp)/real(nset,dp))
    k_b = sqrt((k_chi2wr*inv_dof)*(real(ndata,dp)/real(nset,dp)))
!     write(*,'(2(a,f7.4))')" k_chi2r = ",k_chi2r," k_chi2wr = ",k_chi2wr
!     if(size(nset).eq.2)then
!       write(*,'(a,i5,i5,a,f16.12,a,f16.12,f16.12)')" nset = ",nset," k_a = ",k_a," k_b = ",k_b
!     else if(size(nset).eq.1)then
!       write(*,'(a,i5,a,f16.12,a,f16.12)')" nset = ",nset," k_a = ",k_a," k_b = ",k_b
!     end if
    deallocate(nset)
    
    ! ---
    ! IT SETS THE VARIABLES system_parameters and par with fitting parameters
!     allocate(system_parameters(npar), fitting_parameters(nfit), parameters_minmax(nfit,2), parameter_names(nfit))
!     allocate(parameters_minmax(nfit,2), parameter_names(nfit))
    allocate(parameters_minmax(nfit,2))
    
    ! subroutine: set_par -> init_trades
!     call set_par(m,R,P,a,e,w,mA,i,lN,system_parameters,fitting_parameters)
! !     call set_par(m,R,P,a,e,w,mA,i,lN,mc_allpar,par)
!     ! subroutine: fix_allpar -> init_trades
! !     call fix_system_parameters()
!     call fix_all_parameters(system_parameters)
!   
    ! write(*,*)' init_param'
    call init_param(system_parameters,fitting_parameters)

    parameters_minmax(:,1)=minpar
    parameters_minmax(:,2)=maxpar
    str_len=len(parid(1))
    
    ! write(*,*)' init_derived_parameters'
    ! check if there are derived parameters to compute and to check
    call init_derived_parameters(1,path)
    
    ! deallocated variables not needed anymore
    if(allocated(m)) deallocate(m,R,P,a,e,w,mA,inc,lN)
  
    ! write(*,*)' set new path'
    path=trim(adjustl(path_in))//trim(adjustl(sub_folder))
!     write(*,'(a,a)')" RUNNING IN PATH = ",trim(path)
  
    ! stop 'STOP end of initialize_trades'

    return
  end subroutine initialize_trades
  
  ! ============================================================================
  
  subroutine get_parameter_names(parameter_names, nfit, str_len)
    integer,intent(in)::nfit, str_len
    character(1),dimension(nfit,str_len),intent(out)::parameter_names
    
    integer::in,ic
    
    do in=1,nfit
      do ic=1,str_len
        parameter_names(in,ic) = parid(in)(ic:ic)
      end do
    end do
  
    return
  end subroutine get_parameter_names
  
  ! ============================================================================
  
  subroutine init_fit_parameters(all_parameters,n_par,fit_parameters,n_fit)
    use parameters_conversion,only:init_param
    integer,intent(in)::n_par
    real(dp),dimension(n_par),intent(in)::all_parameters
    integer,intent(in)::n_fit
    real(dp),dimension(n_fit),intent(out)::fit_parameters
    
    real(dp),dimension(:),allocatable::temp_fit
  
    call init_param(all_parameters,temp_fit)
    fit_parameters = temp_fit
    if(allocated(temp_fit)) deallocate(temp_fit)
  
    return
  end subroutine init_fit_parameters
  
  ! ============================================================================
  
  subroutine init_pso(cpuid, path_in)
    integer,intent(in)::cpuid
    character(512),intent(in)::path_in
    character(512)::path_temp
    
    path_temp=trim(adjustl(path))
    path=trim(adjustl(path_in))
    call read_pso_opt(cpuid)
    n_global = nGlobal
    path=trim(adjustl(path_temp))

    return
  end subroutine init_pso
  
  ! ============================================================================
  
  subroutine fortran_loglikelihood(fit_parameters, lgllhd, check, nfit)
    integer,intent(in)::nfit
    real(dp),dimension(nfit),intent(in)::fit_parameters
    real(dp),intent(out)::lgllhd
    logical,intent(out)::check
    real(dp)::fitness
    
    check=.true.
    fitness=zero
    
    fitness=bound_fitness_function(system_parameters,fit_parameters)
    lgllhd=-half*fitness*real(obsData%dof,dp) ! lgllh = - chi2 / 2 || fitness =~ chi2 / dof
!     lgllhd=-half*fitness
    
    if(fitness.ge.resmax)check=.false.

    return
  end subroutine fortran_loglikelihood

  ! ============================================================================
  
  ! subroutine that output the fitness, check for given fit_parameters and global system_parameters
  subroutine fortran_fitness_short(fit_parameters, fitness, check, n_fit)
    integer,intent(in)::n_fit
    real(dp),dimension(n_fit),intent(in)::fit_parameters
    real(dp),intent(out)::fitness
    logical,intent(out)::check
    
    check=.true.
    fitness=zero
    fitness=bound_fitness_function(system_parameters,fit_parameters)
    if(fitness.ge.resmax)check=.false.
    
    return
  end subroutine fortran_fitness_short
  
  ! ============================================================================
  
  ! subroutine that output the fitness, check for given fit_parameters and updated all_parameters (instead of system_parameters)
!   subroutine fortran_fitness_long(fit_parameters, all_parameters, fitness, check, n_fit, n_par)
  subroutine fortran_fitness_long(all_parameters, n_par, fit_parameters, n_fit, fitness, check)
    integer,intent(in)::n_par
    real(dp),dimension(n_par),intent(in)::all_parameters
    integer,intent(in)::n_fit
    real(dp),dimension(n_fit),intent(in)::fit_parameters
    real(dp),intent(out)::fitness
    logical,intent(out)::check
    
    check=.true.
    fitness=zero
    fitness=bound_fitness_function(all_parameters,fit_parameters)
!     write(*,'(a,ES23.16)')' fitness = ',fitness
    if(fitness.ge.resmax)check=.false.
    
    return
  end subroutine fortran_fitness_long
  
  ! ============================================================================
  
!   subroutine write_summary_files(write_number,parameters_values,fitness,wrt_info,lgllhd,check,nfit)
!     use driver,only:write_summary_nosigma
! !     use ode_run,only:ode_out
! !     use output_files,only:write_parameters
!     integer,intent(in)::nfit
!     integer,intent(in)::write_number
!     real(dp),dimension(nfit),intent(in)::parameters_values
!     logical,intent(in),optional::wrt_info
!     real(dp),intent(out)::fitness,lgllhd
!     logical,intent(out)::check
!     real(dp),dimension(:),allocatable::run_all_par
!     logical::check_status
!     integer::i
!     
!     check=.true.
!     check_status=.true.
!     
!     write(*,'(a,l2)') 'check begin = ', check
! !     check=check_fit_boundaries(parameters_values)
!     if(present(wrt_info))then
!       check=check_only_boundaries(system_parameters,parameters_values,wrt_info)
!     else
!       check=check_only_boundaries(system_parameters,parameters_values)
!     end if
!   
!     if(check)then
!       write(*,'(a,l2)') 'check boundaries = ', check
!       allocate(run_all_par(npar))
!       run_all_par=system_parameters
!       if(check_derived) check_status=check_derived_parameters(parameters_values)
!       if(fix_derived) call fix_derived_parameters(parameters_values,run_all_par,check_status)
!       if(check_status)then
!         call write_summary_nosigma(1,write_number,0,run_all_par,parameters_values,fitness)
!       else
!         fitness=resmax ! set it to resmax
!         check=.false.
!       end if
!       deallocate(run_all_par)
!       
!     else
!       write(*,'(a,l2)') 'check boundaries = ', check
!       fitness=resmax
!     end if
!     lgllhd=-0.5_dp*fitness*real(dof,dp)
!     if(fitness.ge.resmax)check=.false.
!     
!     return
!   end subroutine write_summary_files

  ! ============================================================================
  
  subroutine write_summary_files(write_number,parameters_values,fitness,lgllhd,check,nfit)
    use driver,only:write_summary_nosigma
    integer,intent(in)::nfit
    integer,intent(in)::write_number
    real(dp),dimension(nfit),intent(in)::parameters_values
    real(dp),intent(out)::fitness,lgllhd
    logical,intent(out)::check
    real(dp),dimension(:),allocatable::run_all_par
    logical::check_status
!     logical::wrt_info=.true.
        
    check=.true.
    check_status=.true.
    
!     write(*,'(a,l2)') 'check begin = ', check
!     if(present(wrt_info).and.wrt_info)then
!     check=check_only_boundaries(system_parameters,parameters_values,wrt_info)
!     else
!       check=check_only_boundaries(system_parameters,parameters_values)
!     end if
    check=check_only_boundaries(system_parameters,parameters_values)
  
!     write(*,'(a,l2)') 'check boundaries = ', check
    allocate(run_all_par(npar))
    run_all_par=system_parameters
    if(check_derived) check_status=check_derived_parameters(parameters_values)
    if(fix_derived) call fix_derived_parameters(parameters_values,run_all_par,check_status)
    call write_summary_nosigma(1,write_number,0,run_all_par,parameters_values,fitness)
    deallocate(run_all_par)
    if(.not.check.or..not.check_status)then
      fitness=resmax
      write(*,'(a)')'*******'
      write(*,'(a)')'WARNING'
      write(*,'(a)')'WARNING'
      write(*,'(a)')'FITTED PARAMETERS COULD NOT BE PHYSICAL!'
      write(*,'(a)')'BE VERY CAREFUL WITH THIS PARAMETER SET!'
      write(*,'(a)')'WARNING'
      write(*,'(a)')'WARNING'
      write(*,'(a)')'*******'
    end if
    lgllhd=-half*fitness*real(obsData%dof,dp)+ln_err_const
!     lgllhd=-half*fitness+ln_err_const
!     if(fitness.ge.resmax)check=.false.

        
    return
  end subroutine write_summary_files
  
  ! ============================================================================
  
  subroutine write_summary_files_long(write_number,all_parameters,npar,parameters_values,nfit,fitness,lgllhd,check)
    use driver,only:write_summary_nosigma
!     use ode_run,only:ode_out
!     use output_files,only:write_parameters
    integer,intent(in)::npar
    real(dp),dimension(npar),intent(in)::all_parameters
    integer,intent(in)::nfit
    integer,intent(in)::write_number
    real(dp),dimension(nfit),intent(in)::parameters_values
    
    real(dp),intent(out)::fitness,lgllhd
    logical,intent(out)::check
    
    real(dp),dimension(:),allocatable::run_all_par
    logical::check_status
    
    check=.true.
    check_status=.true.
    
!     check=check_fit_boundaries(parameters_values)
    check=check_only_boundaries(system_parameters,parameters_values)
    
    if(check)then
      allocate(run_all_par(npar))
      run_all_par=all_parameters
      if(check_derived) check_status=check_derived_parameters(parameters_values)
      if(fix_derived) call fix_derived_parameters(parameters_values,run_all_par,check_status)
      if(check_status)then
        call write_summary_nosigma(1,write_number,0,run_all_par,parameters_values,fitness)
      else
        fitness=resmax ! set it to resmax
        check=.false.
      end if
      deallocate(run_all_par)
      
    else
      fitness=resmax
    end if
    lgllhd=-half*fitness*real(obsData%dof,dp)+ln_err_const
!     lgllhd=-half*fitness+ln_err_const
    if(fitness.ge.resmax)check=.false.
    
    return
  end subroutine write_summary_files_long
  
  ! ============================================================================
  
  ! pso
  subroutine pyrun_pso(nfit,i_global,best_parameters,best_fitness)
    use opti_pso,only:pso_driver,evaluate_pso
    integer,intent(in)::nfit
    integer,intent(in)::i_global
    real(dp),dimension(nfit),intent(out)::best_parameters
    real(dp),intent(out)::best_fitness
    real(dp)::best_inv_fitness

    path=trim(adjustl(path))
    best_parameters=zero
    best_inv_fitness=one
    call pso_driver(i_global,evaluate_pso,nfit,system_parameters,minpar,maxpar,&
      &best_parameters,best_inv_fitness) ! PSO DRIVER
    best_fitness=one/best_inv_fitness

    return
  end subroutine pyrun_pso
  
  ! ============================================================================
  
  ! subroutine useful to modify the working path fo TRADES from python
  subroutine path_change(new_path)
    character(512),intent(in)::new_path
    
    path=trim(adjustl(new_path))
!     write(*,*)trim(adjustl(path))
    
    return
  end subroutine path_change

  ! ============================================================================
  
  ! init both cases for derived parameters
  ! 1)
  ! check if there are derived parameters to compute and to check
  subroutine init_check_parameters(cpuid,path_in)
    integer,intent(in)::cpuid ! cpu number: use 1
    character(512),intent(in)::path_in ! path of the folder with derived_boundaries.dat
    
    call init_check_derived_parameters(cpuid,path_in)
    
    return
  end subroutine init_check_parameters
  
  ! ============================================================================
  
  ! 2)
  ! check if there are derived parameters to compute and to check
  subroutine init_fix_parameters(n_derived_in,in_names,in_parameters)
    integer,intent(in)::n_derived_in
    character(15),dimension(n_derived_in),intent(in)::in_names
    real(dp),dimension(n_derived_in),intent(in)::in_parameters
    
    call init_fix_derived_parameters(n_derived_in,in_names,in_parameters)
    
    return
  end subroutine init_fix_parameters

  ! ============================================================================
  
  subroutine deallocate_variables()
  
    call deallocate_all() ! from 'parameters' module
    
    return
  end subroutine deallocate_variables
  
  ! ============================================================================
  ! +++***+++
  ! CHECKS THE NEXT!
  ! +++***+++
  
  ! SUBROUTINE TO INITIALISE TRADES WITHOUT READING FILES
  subroutine args_init(t_start,t_epoch,t_int,n_body,&
    &n_t0,t0_num,t0_obs,et0_obs,&
    &n_max_t0,n_col)
    
    ! INPUT
    ! t_start      == start of the integration
    ! t_epoch      == reference time epoch
    ! t_int        == total integration time in days
    ! n_body       == number of bodies (take into account the star)
    ! n_t0         == number of transits per each body n_t0(n_body); n_t0(0) = 0
    ! t0_num       == epochs/transit number for each body t0_num(n_max_t0,n_body); t0_num(:,0) = 0
    ! t0_obs       == transit times for each body t0_obs(n_max_t0,n_body); t0_obs(:,0) = 0
    ! et0_obs      == errors on the transit times in days for each body et0_obs(n_max_t0,n_body); et0_obs(:,0) = 0
    
    ! OUTPUT
    ! None ==> some variables set globally
    
    ! DIMENSIONS: do not provide it
    ! n_max_t0     == maxval(n_t0) == maxval of transits available
    
    
    integer::n_max_t0,n_col
    real(dp),intent(in)::t_start,t_epoch,t_int
    integer,intent(in)::n_body
    integer,dimension(n_col),intent(in)::n_t0
    integer,dimension(n_max_t0,n_col),intent(in)::t0_num
    real(dp),dimension(n_max_t0,n_col),intent(in)::t0_obs
    real(dp),dimension(n_max_t0,n_col),intent(in)::et0_obs
    
!f2py integer intent(hide),depend(t0_num,t0_obs,et0_obs)::n_max_t0=shape(t0_num,0),n_col=shape(t0_num,1)
    
    integer::ibd,nx
    
    tstart=t_start
    tepoch=t_epoch
    tint=t_int
    NB=n_body
    NBDIM=n_body*6
    
    allocate(e_bounds(2,n_body))
    e_bounds(1,:)=zero
    e_bounds(2,:)=one-TOL_dp
    
!     call set_ephem(n_body,n_t0,t0_num,t0_obs,et0_obs)
    obsData%nTTs=sum(n_t0)
    allocate(obsData%obsT0(n_body-1))
    do ibd=1,n_body-1
      nx=n_t0(ibd+1)
      obsData%obsT0(ibd)%nT0=nx
      allocate(obsData%obsT0(ibd)%epo(nx),obsData%obsT0(ibd)%T0(nx),&
        &obsData%obsT0(ibd)%eT0(nx))
      obsData%obsT0(ibd)%epo=t0_num(1:nx,ibd+1)
      obsData%obsT0(ibd)%T0=t0_obs(1:nx,ibd+1)
      obsData%obsT0(ibd)%eT0=et0_obs(1:nx,ibd+1)
    end do
    call set_ephem()
    call compute_oc(obsData%obsT0)
    
    rvcheck=1
    durcheck=0
    amin=TOL_dp
    amax=1.e4_dp
    
    return
  end subroutine args_init
  
  ! ============================================================================
  
  !!! SUBROUTINE TO RUN TRADES INTEGRATION AND RETURN RV_SIM AND T0_SIM
!   subroutine kelements_to_data(t_start,t_epoch,step_in,t_int,&
!     &m_msun,R_rsun,P_day,ecc,argp_deg,mA_deg,inc_deg,lN_deg,&
!     &t_rv,transit_flag,n_t0,t0_num,& ! input
!     &rv_sim,t0_sim,& ! output
!     &n_body,n_rv,n_max_t0) ! dimensions
  subroutine kelements_to_data(t_start,t_epoch,step_in,t_int,&
    &m_msun,R_rsun,P_day,ecc,argp_deg,mA_deg,inc_deg,lN_deg,&
    &t_rv,transit_flag,& ! input
    &rv_sim,t0_sim,& ! output
    &n_max_t0,& ! input dimension to be provided
    &n_body,n_rv) ! dimensions to be not provided

    ! INPUT
    ! t_start      == start of the integration
    ! t_epoch      == reference time epoch
    ! step_in      == initial step size of the integration
    ! t_int        == total integration time in days
    
    ! m_msun       == masses of all the bodies in Msun m_sun(n_body)
    ! R_rsun       == radii of all the bodies in Rsun r_rsun(n_body)
    ! P_day        == periods of all the bodies in days p_day(n_body); p_day(0) = 0
    ! ecc          == eccentricities of all the bodies ecc(n_body); ecc(0) = 0
    ! argp_deg     == argument of pericentre of all the bodies argp_deg(n_body); argp_deg(0) = 0
    ! mA_deg       == mean anomaly of all the bodies mA_deg(n_body); mA_deg(0) = 0
    ! inc_deg      == inclination of all the bodies inc_deg(n_body); inc_deg(0) = 0
    ! lN_deg       == longitude of node of all the bodies lN_deg(n_body); lN_deg(0) = 0

    ! t_rv         == time of the RV datapoints t_rv(n_rv)
    ! transit_flag == logical/boolean vector with which bodies should transit (.true.) or not (.false.) transit_flag(n_body); transit_flag(0) = False
    ! n_max_t0     == maxval(n_t0) == maxval of transits available
    
    ! OUTPUT
    ! rv_sim       == rv simulated in m/s, same dimension of t_rv
    ! t0_sim       == t0 simulated in days, same dimension of t0_num
    
    ! DIMENSIONS
    ! n_body       == number of bodies (take into account the star)
    ! n_rv         == number of radial velocities datapoints
    integer::n_body,n_rv
    
    real(dp),intent(in)::t_start,t_epoch,step_in,t_int
    real(dp),dimension(n_body),intent(in)::m_msun,R_rsun,P_day
    real(dp),dimension(n_body),intent(in)::ecc,argp_deg,mA_deg,inc_deg,lN_deg
    real(dp),dimension(n_RV),intent(in)::t_rv
    logical,dimension(n_body),intent(in)::transit_flag
    integer,intent(in)::n_max_t0

!     integer,dimension(n_body),intent(in)::n_t0
!     integer,dimension(n_max_t0,n_body),intent(in)::t0_num
    
    real(dp),dimension(n_RV),intent(out)::rv_sim
    real(dp),dimension(n_max_t0,n_body),intent(out)::t0_sim
    

!f2py    integer,intent(hide),depend(t_rv)::n_rv=len(t_rv)


! !f2py    integer,intent(hide),depend(t0_num)::n_max_t0=shape(t0_num,0), n_body=shape(t0_num,1)

! !f2py    integer,intent(hide),depend(n_t0)::n_body=len(m_msun)
    
    integer::id_transit_body=1 ! needed to be == 1
    real(dp),dimension(:),allocatable::rv_temp
    real(dp),dimension(:,:),allocatable::t0_temp
    
    ! DATA
    ! set obsData from here!!
    ! RV
    ! method 1 ...
!     obsData%obsRV%nRV=n_RV
!     if(.not.allocated(obsData%obsRV%jd)) allocate(obsData%obsRV%jd(n_RV))
!     obsData%obsRV%jd = t_rv
    ! method 2 ...
    if(allocated(obsData%obsRV%jd)) deallocate(obsData%obsRV%jd)
    obsData%obsRV%nRV=n_RV
    allocate(obsData%obsRV%jd(n_RV))
    obsData%obsRV%jd = t_rv
    ! TT in args_init
    
    call orbits_to_data(t_start,t_epoch,step_in,t_int,&
      &m_msun,R_rsun,P_day,ecc,argp_deg,mA_deg,inc_deg,lN_deg,&
      &t_rv,rv_temp,&
      &id_transit_body,transit_flag,durcheck,t0_temp)
    
    rv_sim=rv_temp
    t0_sim=t0_temp
    if(allocated(rv_temp)) deallocate(rv_temp)
    if(allocated(t0_temp)) deallocate(t0_temp)
  
    return
  end subroutine kelements_to_data

  ! ============================================================================
  
  ! create wrappers to init the grid of a perturber body (a)
  ! and set properly the parameters to use in TRADES (b)
  ! the output will be all the transit times within the t_int
  ! for each parameter combination.
  
  ! (a)
  subroutine wrapper_read_grid(cpuid,ngrid,ncol_grid)
    integer,intent(in)::cpuid
    integer,intent(out)::ngrid,ncol_grid
    
    ! read input file <-> idpert
    call read_parameters_grid(cpuid,perturber_parameters_grid)
    ! perturber_parameters_grid in 'parameters' module
    
    ! set properly the fields of the perturber_parameters_grid variable
    call set_parameters_grid(perturber_parameters_grid,ngrid)
    
    ncol_grid=10 ! fixed in module grid_search -> build_grid
  
    return
  end subroutine wrapper_read_grid
    
  subroutine wrapper_grid_init(cpuid,ngrid,ncol_grid,perturber_grid)
    integer,intent(in)::cpuid
    integer,intent(in)::ngrid,ncol_grid
    real(dp),dimension(ngrid,ncol_grid),intent(out)::perturber_grid

!     !f2py    integer,intent(hide),depend(perturber_grid)::ngrid=shape(perturber_grid,0), ncol_grid=shape(perturber_grid,1)
    
    real(dp),dimension(:,:),allocatable::temp_grid,fitness_grid
    
    ! create/build the full grid, with all the combination of the parameters of perturber body
    call build_grid(MR_star(1,1),perturber_parameters_grid,temp_grid,fitness_grid,ngrid)
    ! perturber_parameters_grid in 'parameters' module
    
    perturber_grid=temp_grid
    
    deallocate(temp_grid,fitness_grid) ! temporary -> useless here
    
    ! needed to update the max allowed value of the semi-major axis
    amax=5._dp*maxval(perturber_grid(:,4))
    
    return
  end subroutine wrapper_grid_init
  
  ! ============================================================================
  
  ! (b) - TOBESPLITTED
!   subroutine wrapper_set_grid_orbelem(sim_id,idpert,perturber_grid,ngrid,ncol,&
!     &sim_all_parameters,npar,&
!     &ttra_all,id_ttra_all,n_ttra_all,time_rv_all,rv_all,n_rv_all)
!     integer,intent(in)::sim_id,idpert
!     integer::ngrid,ncol
!     real(dp),dimension(ngrid,ncol),intent(in)::perturber_grid
!     integer::npar
!     real(dp),dimension(npar),intent(out)::sim_all_parameters
!     integer::n_ttra_all
!     real(dp),dimension(n_ttra_all),intent(out)::ttra_all
!     integer,dimension(n_ttra_all),intent(out)::id_ttra_all
!     integer::n_rv_all
!     real(dp),dimension(n_rv_all),intent(out)::time_rv_all
!     real(dp),dimension(n_rv_all),intent(out)::rv_all
!     
! !f2py    integer,intent(hide),depend(perturber_grid)::ngrid=shape(perturber_grid,0), ncol=shape(perturber_grid,1)
! !f2py    integer,intent(hide),depend(sim_all_parameters)::npar=len(sim_all_parameters)
! !f2py    integer,intent(hide),depend(ttra_all)::n_ttra_all=len(ttra_all)
! !f2py    integer,intent(hide),depend(time_rv_all)::n_rv_all=len(time_rv_all)
! 
!     real(dp),dimension(:),allocatable::sim_fit_parameters
!     real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
!     real(dp),dimension(:),allocatable::temp_ttra,temp_time_rv,temp_rv
!     integer,dimension(:),allocatable::temp_id_ttra
!     logical::tempcheck
!     
!     allocate(sim_fit_parameters(nfit))
!     sim_all_parameters=system_parameters
!     sim_fit_parameters=zero
!     
!     ! 1. select proper set of parameters to use: perturber_grid to sim_all_parameters
!     call perturber_grid2parameters(sim_id,idpert,perturber_grid,sim_all_parameters)
!     ! 2. update from sim_all_parameters to sim_fit_parameters
!     call init_param(sim_all_parameters,sim_fit_parameters)
!     
!     allocate(m(NB),R(NB),P(NB),a(NB),e(NB),w(NB),mA(NB),inc(NB),lN(NB))
!     call par2kel_fit(sim_all_parameters,sim_fit_parameters,m,R,P,a,e,w,mA,inc,lN,tempcheck)
!     deallocate(sim_fit_parameters)
! 
!     call ode_all_ttra_rv(wrttime,m,R,P,a,e,w,mA,inc,lN,&
!       &temp_ttra,temp_id_ttra,temp_time_rv,temp_rv)
!     ttra_all=temp_ttra
!     id_ttra_all=temp_id_ttra
!     time_rv_all=temp_time_rv
!     rv_all=temp_rv
!     deallocate(temp_time_rv,temp_rv,temp_ttra,temp_id_ttra)
!     
!     return
!   end subroutine wrapper_run_grid_parameters
  
  ! ============================================================================
  
    subroutine get_max_nt0_nrv(P,n_bodies,nt0_full,nrv_nmax)
    integer,intent(in)::n_bodies
    real(dp),dimension(n_bodies),intent(in)::P
    integer,intent(out)::nt0_full,nrv_nmax
    
    integer,dimension(:),allocatable::nT0_perbody
  
!     write(*,*)'f90: n_bodies = ',n_bodies
!     write(*,*)'f90: size(P) = ',size(P),' => P = ',P
  
  ! compute number maximum of all transit times of all the planets from the
    ! integration time (global variable tint)
    allocate(nT0_perbody(n_bodies-1))
    nT0_perbody = int((1.5_dp*tint)/P(2:n_bodies))+2 ! compute the number of transit for each planet ... avoid P(1)<->star
    nt0_full=sum(nT0_perbody) ! sum them all
    deallocate(nT0_perbody)

    ! compute the number of maximum rv from integration and write time
    ! (global variables tint, wrttime
    nrv_nmax=int(tint/wrttime)+2
  
    return
  end subroutine get_max_nt0_nrv
  
  ! ============================================================================

  
  subroutine wrapper_set_grid_orbelem_fullgrid(sim_id,idpert,perturber_grid,ngrid,ncol,&
    &m,R,P,sma,ecc,w,mA,inc,lN,n_bodies,nt0_full,nrv_nmax)
    integer,intent(in)::sim_id,idpert
    integer,intent(in)::ngrid,ncol
    real(dp),dimension(ngrid,ncol),intent(in)::perturber_grid
    integer,intent(in)::n_bodies
    real(dp),dimension(n_bodies),intent(out)::m,R,P,sma,ecc,w,mA,inc,lN
    integer,intent(out)::nt0_full,nrv_nmax
    
!f2py    integer,intent(hide),depend(perturber_grid)::ngrid=shape(perturber_grid,0), ncol=shape(perturber_grid,1)

! !f2py    integer,intent(hide),depend(m)::n_bodies=len(m)

    real(dp),dimension(:),allocatable::sim_all_parameters,sim_fit_parameters
    logical::tempcheck
    
    integer,dimension(:),allocatable::nT0_perbody
    
    tempcheck=.true.
    allocate(sim_all_parameters(npar),sim_fit_parameters(nfit))
    sim_all_parameters=system_parameters
    sim_fit_parameters=zero
    
    ! 1. select proper set of parameters to use: perturber_grid to sim_all_parameters
    call perturber_grid2parameters(sim_id,idpert,perturber_grid,sim_all_parameters)
    ! 2. update from sim_all_parameters to sim_fit_parameters
    call init_param(sim_all_parameters,sim_fit_parameters)
    
    call par2kel_fit(sim_all_parameters,sim_fit_parameters,m,R,P,sma,ecc,w,mA,inc,lN,tempcheck)
    deallocate(sim_all_parameters,sim_fit_parameters)
    
!     ! compute number maximum of all transit times of all the planets from the
!     ! integration time (global variable tint)
!     allocate(nT0_perbody(n_bodies-1))
!     nT0_perbody = int((1.5_dp*tint)/P(2:n_bodies))+2 ! compute the number of transit for each planet ... avoid P(1)<->star
!     nt0_full=sum(nT0_perbody) ! sum them all
!     deallocate(nT0_perbody)
! 
!     ! comput the number of maximum rv from integration and write time
!     ! (global variables tint, wrttime
!     nrv_nmax=int(tint/wrttime)+2
    
    call get_max_nt0_nrv(P,n_bodies,nt0_full,nrv_nmax)
    
    return
  end subroutine wrapper_set_grid_orbelem_fullgrid
  
  ! ============================================================================

  subroutine wrapper_set_grid_orbelem(sim_id,idpert,perturber_par,ncol,&
    &m,R,P,sma,ecc,w,mA,inc,lN,n_bodies,nt0_full,nrv_nmax)
    integer,intent(in)::sim_id,idpert
    integer,intent(in)::ncol
    real(dp),dimension(ncol),intent(in)::perturber_par
    integer,intent(in)::n_bodies
    real(dp),dimension(n_bodies),intent(out)::m,R,P,sma,ecc,w,mA,inc,lN
    integer,intent(out)::nt0_full,nrv_nmax
    
!f2py    integer,intent(hide),depend(perturber_par)::ncol=len(perturber_par)
! !f2py    integer,intent(hide),depend(m)::n_bodies=len(m)

    real(dp),dimension(:),allocatable::sim_all_parameters,sim_fit_parameters
    logical::tempcheck
    
!     integer,dimension(:),allocatable::nT0_perbody
    
    tempcheck=.true.
    allocate(sim_all_parameters(npar),sim_fit_parameters(nfit))
    sim_all_parameters=system_parameters
    sim_fit_parameters=zero
    
    ! 1. select proper set of parameters to use: perturber_grid to sim_all_parameters
!     call perturber_grid2parameters(sim_id,idpert,perturber_grid,sim_all_parameters)
    call perturber_par2all_par(idpert,perturber_par,sim_all_parameters)
    ! 2. update from sim_all_parameters to sim_fit_parameters
    call init_param(sim_all_parameters,sim_fit_parameters)
    
    call par2kel_fit(sim_all_parameters,sim_fit_parameters,m,R,P,sma,ecc,w,mA,inc,lN,tempcheck)
    deallocate(sim_all_parameters,sim_fit_parameters)
    
!     ! compute number maximum of all transit times of all the planets from the
!     ! integration time (global variable tint)
!     allocate(nT0_perbody(n_bodies-1))
!     nT0_perbody = int((1.5_dp*tint)/P(2:n_bodies))+2 ! compute the number of transit for each planet ... avoid P(1)<->star
!     nt0_full=sum(nT0_perbody) ! sum them all
!     deallocate(nT0_perbody)
! 
!     ! compute the number of maximum rv from integration and write time
!     ! (global variables tint, wrttime
!     nrv_nmax=int(tint/wrttime)+2
    call get_max_nt0_nrv(P,n_bodies,nt0_full,nrv_nmax)
    
    return
  end subroutine wrapper_set_grid_orbelem
  
  ! ============================================================================
  
  subroutine wrapper_run_grid_combination(m,R,P,sma,ecc,w,mA,inc,lN,n_bodies,&
    &ttra_full,id_ttra_full,stats_ttra,nt0_full,&
    &time_rv_nmax,rv_nmax,stats_rv,nrv_nmax)
    integer,intent(in)::n_bodies
    real(dp),dimension(n_bodies),intent(in)::m,R,P,sma,ecc,w,mA,inc,lN
    integer,intent(in)::nt0_full
    
    real(dp),dimension(nt0_full),intent(out)::ttra_full
    integer,dimension(nt0_full),intent(out)::id_ttra_full
    logical,dimension(nt0_full),intent(out)::stats_ttra
    integer,intent(in)::nrv_nmax
    real(dp),dimension(nrv_nmax),intent(out)::time_rv_nmax,rv_nmax
    logical,dimension(nrv_nmax),intent(out)::stats_rv
    
    real(dp),dimension(nt0_full)::dur_full
    
    call ode_all_ttra_rv(wrttime,m,R,P,sma,ecc,w,mA,inc,lN,&
      &ttra_full,dur_full,id_ttra_full,stats_ttra,&
      &time_rv_nmax,rv_nmax,stats_rv)
  
    return
  end subroutine wrapper_run_grid_combination
  
  ! ============================================================================
  
  subroutine fit_par_to_ttra_rv(fit_parameters,nfit,&
    &ttra_full,dur_full,id_ttra_full,stats_ttra,nt0_full,&
    &time_rv_nmax,rv_nmax,stats_rv,nrv_nmax)
    integer,intent(in)::nfit
    real(dp),dimension(nfit),intent(in)::fit_parameters
    integer,intent(in)::nt0_full
    
    real(dp),dimension(nt0_full),intent(out)::ttra_full,dur_full
    integer,dimension(nt0_full),intent(out)::id_ttra_full
    logical,dimension(nt0_full),intent(out)::stats_ttra
    
    integer,intent(in)::nrv_nmax
    real(dp),dimension(nrv_nmax),intent(out)::time_rv_nmax,rv_nmax
    logical,dimension(nrv_nmax),intent(out)::stats_rv
    
    real(dp),dimension(:),allocatable::m,R,P,sma,ecc,w,mA,inc,lN
    logical::checkpar=.true.
    
    allocate(m(NB),R(NB),P(NB),sma(NB),ecc(NB),w(NB),mA(NB),inc(NB),lN(NB))
    call convert_parameters(system_parameters,fit_parameters,&
      &m,R,P,sma,ecc,w,mA,inc,lN,checkpar)
      
    call ode_all_ttra_rv(wrttime,m,R,P,sma,ecc,w,mA,inc,lN,&
      &ttra_full,dur_full,id_ttra_full,stats_ttra,&
      &time_rv_nmax,rv_nmax,stats_rv)
    
    deallocate(m,R,P,sma,ecc,w,mA,inc,lN)
    
    return
  end subroutine fit_par_to_ttra_rv

  ! ============================================================================
  
  subroutine convert_trades_par_to_kepelem(all_par,npar,fit_par,nfit,&
    &m,R,P,sma,ecc,w,mA,inc,lN,n_bodies)
    integer,intent(in)::npar,nfit
    real(dp),dimension(npar),intent(in)::all_par
    real(dp),dimension(nfit),intent(in)::fit_par
    
    integer,intent(in)::n_bodies
    real(dp),dimension(n_bodies),intent(out)::m,R,P,sma,ecc,w,mA,inc,lN
  
    logical::checkpar
  
    call convert_parameters(all_par,fit_par,m,R,P,sma,ecc,w,mA,inc,lN,checkpar)
      
   return
  end subroutine convert_trades_par_to_kepelem
  
  ! ============================================================================
  
end module pytrades
