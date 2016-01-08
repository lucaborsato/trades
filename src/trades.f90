! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
! MAIN PROGRAM -- TRADES -- Luca Borsato 2014
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
program main_trades
  use constants
  use parameters
  use init_trades
  use timing,only:timer
  use transits,only:set_ephem
  use output_files,only:write_lm_inputpar,write_par,write_simlst
  use ode_run,only:ode_out,ode_lm
  use grid_search,only:grid_build_2,set_grid_par_2
  use Levenberg_Marquardt,only:lm_driver
  use Genetic_Algorithm,only:ga_driver,fpik
  use opti_pso,only:pso_driver,evfpso
  use bootstrap
  use mpi_module
  use PolyChord_driver,only:PC_driver
  !$ use omp_lib
  implicit none
  real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,i,lN ! Keplerian orbital elements
  integer::cpuid ! CPU number
  real(dp),dimension(:),allocatable::allpar,par ! vectors with all the parameters and only the parameters to be fitted, needed to LM, GA, and PSO

  real(dp),dimension(:),allocatable::resw,copar,sigpar ! vectors for weighte residuals, covariance matrix values of fitted parameters, sigma of fitted parameters
  integer,dimension(:),allocatable::iwa ! working vector for LM
  integer::info ! returned by LM. It depends on the LM convergence

  integer::Ngrid,jgrid ! Number of grid combinations, iterator in grid search
  !real(dp),dimension(:),allocatable::Pgrid,agrid,egrid,wgrid ! parameter in grid search
  real(dp),dimension(:),allocatable::mw,Pw,aw,ew,ww,xpar ! working Keplerian elements in the grid search
  real(dp),dimension(:,:),allocatable::grid ! array with the orbital combinations for the gread search, plus Chi2 and Chi2r
  real(dp)::inv_fitness,fitness ! inverse Fitness and Fitness: Fitness = Chi2r*k_chi2r + Chi2wr*k_chi2wr
  integer,dimension(:),allocatable::infos ! LM convergence type for each grid search combination

  character(80)::fmt ! write file format string

  real(dp)::start1,end1,sec1 ! timing serial
  integer::hour1,minute1
  real(dp)::starta,enda,seca
  integer::houra,minutea
  !$ real(dp)::ostart,oend ! timing openMP
  integer,dimension(8)::date_values

  !   integer::nseed ! number of seeds to generate for random numbers

  ! PolyChord
  real(dp),dimension(5)::PolyChord_info

  !   integer::my_rank,mpi_nthreads

  call date_and_time(VALUES=date_values)
  starta=zero
  enda=zero
  call cpu_time(starta)
  !$ ostart=zero
  !$ oend=zero
  !$ ostart=omp_get_wtime()

  write(*,'(a)')""
  write(*,'(a)')" ======================================================"
  write(*,'(a)')" TRADES v. 2.5.0 -- Luca Borsato 2016"
  write(*,'(a)')" ======================================================"
  write(*,'(a)')""
  write(*,'(a,i4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2)')" START TIME: ",&
  &date_values(1),"-",date_values(2),"-",date_values(3)," -- ",&
  &date_values(5),":",date_values(6),":",date_values(7)
  write(*,'(a)')""

  ! IT INITIALIZES THE CPU AND NUMBER OF CPU ... USEFUL FOR OPENMP COMPUTATION
  cpuid=1
  ncpu=1
  call initu(nfiles,ncpu) ! prepares writing file units

  ! IT READS THE COMMAND ARGUMENT THAT DEFINE THE PATH OF THE FOLDER WITH THE FILES
  call read_com_arg(path)

  ! IT DEFINES THE STRING TO WRITE THE REAL WITH RIGHT DECIMAL: dpECISION
  sprec=trim(adjustl("g"//&
  &trim(adjustl(string(2*prec)))//&
  &"."//&
  &trim(adjustl(string(prec)))))

  write(*,'(a)') " -----------------------------------------------------"
  write(*,'(a)') " SUMMARY INFO "
  write(*,'(a,a,a,a)')" SELECTED PATH: ",trim(path),&
  &" FORMAT EDIT DESCRIPTOR: ",trim(sprec)

  nRV=0 ! set number of RV points to zero (warning -Wunitialized ...)
  ! IT CALLS ALL THE SUBROUTINES TO READ ALL PARAMETERS AND DATA TO STORE IN COMMON MODULE PARAMETERS
  call read_first(cpuid,m,R,P,a,e,w,mA,i,lN)
  ! IT SETS THE LINEAR EPHEMERIS FROM T0 DATA
  call set_ephem()

  write(*,'(a)') " -----------------------------------------------------"
  write(*,'(a)') " ====================================================="
  write(*,'(a)')""
  write(*,'(4(a))')" PROGTYPE ",trim(adjustl(string(progtype))),&
  &" IDPERT ",trim(adjustl(string(idpert)))
  write(*,'(a)')""

  flush(6)

  ! --- GRID SEARCH ---
  if(progtype.eq.1)then
    !$omp parallel default(shared) &
    !$omp private(cpuid,ncpu)
    !$ cpuid=omp_get_thread_num()+1
    !$ ncpu=omp_get_num_threads()
    !$omp master
    !$ call initu(nfiles,ncpu)
    ! IT CALLS THE SUBROUTINE TO BUILD THE GRID PARAMETERS
    call grid_build_2(cpuid,m,P,a,e,w,Ngrid,grid)
    !write(*,'(a,2i5)')" shape(grid) = ",shape(grid)
    write(*,'(a)') " -----------------------------------------------------"
    allocate(infos(Ngrid))
    infos=0
    if(.not.allocated(mw)) allocate(mw(NB),Pw(NB),aw(NB),ew(NB),ww(NB))
    !$omp end master
    !$ WRITE(*,'(a,i4,a)')" CPU ",cpuid," REACHED BARRIER"
    !$omp barrier
    !$omp do private(jgrid,mw,Pw,aw,ew,ww,allpar,par,resw,iwa,copar,&
    !$omp& sigpar,info,start1,end1,hour1,minute1,sec1) &
    !$omp& schedule(dynamic)

    looppar: do jgrid=1,Ngrid
      write(*,'(a,i4,a,i8,a,i8)')" CPU ",cpuid,&
        &" START GRID SIMULATION NUMBER ",jgrid," / ",Ngrid
      write(*,'(a)')""
      call set_grid_par_2(jgrid,m,P,a,e,w,grid,mw,Pw,aw,ew,ww)
      ! IT SETS THE PARAMETERS LONG/SHORT VECTORs TO BE USED IN LM-dif
      call set_par(mw,R,Pw,aw,ew,ww,mA,i,lN,allpar,par)
      write(*,'(a)') " -----------------------------------------------------"
      if(.not.allocated(resw)) allocate(resw(ndata),iwa(nfit),copar(nfit),sigpar(nfit))
      resw=zero
      copar=zero
      sigpar=zero
      info=0
      write(*,*)""

      ! LM fit if idpert is greater than 0
      if(lmon.gt.0)then
        if(idpert.gt.0)then
          fmt=adjustl("(a,i6,a,1000"//trim(sprec)//")")
          write(*,*)""
          write(*,'(a,i4,a)')" CPU ",cpuid," CALLING LM with INPUT PARAMETERS:"
          write(*,*)""
          call write_lm_inputpar(cpuid,par) ! write to screen parameters to fit
          write(*,*)""
          write(*,'(a)') " -----------------------------------------------------"
          call cpu_time(start1)
          call lm_driver(cpuid,jgrid,ode_lm,allpar,ndata,nfit,par,&
            &resw,copar,sigpar,info,iwa) ! IT CALLS L-M
          call cpu_time(end1)
          call timer(start1,end1,hour1,minute1,sec1)
          write(*,'(a,i2)')" info = ",info
          infos(jgrid)=info
          write(*,'(a)') " -----------------------------------------------------"
          write(*,*)""
          write(*,'(a,i4,a,i4,a,i3,a,f5.2,a)')&
            &" CPU ",cpuid," DONE L-M in ",&
            &hour1," h : ",minute1," m : ",sec1," s"
          write(*,*)""
          write(*,'(a,i4,a)')" CPU ",cpuid," OUTPUT PARAMETERS:"
          call write_lm_inputpar(cpuid,par)
          write(*,*)""
          call param_adj(par,sigpar) ! adjust parameters and sigma (i.e., angle [0-360]deg
          call write_par(cpuid,par,sigpar,resw) ! write to screen parameters
          call write_par(cpuid,jgrid,lmon,par,sigpar,resw) ! write to files parameters
          grid(6,jgrid)=sum(resw*resw)
        end if
      end if

      ! IT WRITES FIT RESULTS TO SCREEN AND FILES
      if(.not.allocated(resw)) allocate(resw(ndata))
      resw=zero
      call cpu_time(start1)
      call ode_out(cpuid,jgrid,lmon,allpar,par,resw) ! IT CALLS SUBROUTINE TO CALCULATE ORBIT, CONST. of MOTION, KEP. ELEMENTS, ALL THE TRANSITS, RV. IT WRITES THEM INTO FILES

      if(idpert.eq.0 .or. lmon.eq.0) grid(6,jgrid)=sum(resw*resw)
      call cpu_time(end1)
      call timer(start1,end1,hour1,minute1,sec1)
      write(*,*)""
      write(*,'(a,i4,a,a,a,i4,a,i3,a,f5.2,a)')&
        &" CPU ",cpuid," DONE INTEGRATION ",&
        &trim(adjustl(string(jgrid)))," AND WRITE in ",&
        &hour1," h : ",minute1," m : ",sec1," s"
      write(*,*)""
      !        fmt=adjustl("(2(a,"//trim(sprec)//"))")
      !        if(ndata.gt.0)then
      !           write(*,trim(fmt))" Fitness(Chi2r*k_chi2r + Chi2wr*k_chi2wr) = ",&
      !             &sum(resw*resw),&
      !             &" => Fitness*dof = ",(sum(resw*resw)*real(dof,dp))
      !           write(*,*)""
      !        end if
      deallocate(resw)
      if(allocated(iwa)) deallocate(iwa,copar,sigpar)

      ! BOOTSTRAP ANALYSIS for each grid combination
      if(nboot.gt.0)then
        write(*,'(a)') " -----------------------------------------------------"
        call cpu_time(start1)
        call strap_driver(jgrid,allpar,par) ! IT CALLS THE BOOTSTRAP-DRIVER
        call cpu_time(end1)
        call timer(start1,end1,hour1,minute1,sec1)
        write(*,'(a)') " -----------------------------------------------------"
        write(*,*)""
        write(*,'(a,i4,a,i4,a,i3,a,f5.2,a)')&
          &" CPU ",cpuid," DONE BOOSTRAP in ",&
          &hour1," h : ",minute1," m : ",sec1," s"
        write(*,*)""
      end if

    end do looppar

    !$omp end do
    !$ write(*,'(a,i4,a)')" CPU ",cpuid," ENDED ALL OWN SIMULATIONS"
    !$omp barrier
    !$omp master
    if(allocated(Pw)) deallocate(mw,Pw,aw,ew,ww)
    call write_simlst(cpuid,Ngrid,grid,infos) ! IT WRITES THE ORIGINAL PARAMETERS OF THE GRID WITH THE FINAL FITNESS
    write(*,'(a)')""
    !$omp end master
    !$omp barrier
    !$omp end parallel

  ! === OTHER METHODS (PIK-PSO-LM-BOOTSTRAP-ONLY INTEGRATION) === 
  else if((progtype.eq.0).or.(progtype.gt.1))then
    if(progtype.lt.5)then
      write(*,'(a,i3)')" PROGTYPE: ",progtype
      if(progtype.gt.1)then
        write(*,*)
        write(*,'(a,i3)')" GLOBAL SEARCH -- NUMBER OF GLOBAL SIMULATIONS = ",nGlobal
        write(*,*)
      end if

      doGlobal: do jgrid=1,nGlobal
        call set_par(m,R,P,a,e,w,mA,i,lN,allpar,par)
        cpuid = 1
        ! --- PIKAIA/GENETIC ---
        if(progtype.eq.3)then
          write(*,'(a)')" +++++++++++++++++++++++++++++++++++++++++++++++++++ "
          write(*,'(a)')" +++++++++++++++++++++++++++++++++++++++++++++++++++ "
          write(*,'(a,i4,a,i3,a,i3)')" CPU ",cpuid,&
            &" START PIKAIA SIMULATION number ",jgrid," / ",nGlobal
          allocate(xpar(nfit))
          xpar=zero
          !call random_number(xpar)
          inv_fitness=zero
          fitness=zero
          call cpu_time(start1)
          call ga_driver(jgrid,fpik,nfit,allpar,xpar,inv_fitness) ! GA DRIVER
          fitness=one/inv_fitness
          call cpu_time(end1)
          call timer(start1,end1,hour1,minute1,sec1)
          write(*,*)""
          write(*,'(a,i4,a,i4,a,i3,a,f5.2,a)')&
            &" CPU ",cpuid," DONE PIKAIA in ",&
            &hour1," h : ",minute1," m : ",sec1," s"
          write(*,*)""
          write(*,'(a,i4,a,a,a)')" CPU ",cpuid," PIKAIA ",&
            &trim(adjustl(string(jgrid)))," OUTPUT PARAMETERS"
          call write_lm_inputpar(cpuid,xpar)
          write(*,*)""
          write(*,'(2(a,g25.14),a,a,a,g25.14)')&
            &" inv_fitness = ",inv_fitness," => fitness = ",fitness," => fitness*",&
            &trim(adjustl(string(dof)))," = ",fitness*real(dof,dp)
          write(*,*)""
          write(*,'(a)') " -----------------------------------------------------"
          call norm2par(xpar,par,allpar) ! from [0-1] parameter values to physical values
          deallocate(xpar)

        ! --- PARTICLE SWARM OPTIMIZATION ---
        else if(progtype.eq.4)then
          write(*,'(a)')" +++++++++++++++++++++++++++++++++++++++++++++++++++ "
          write(*,'(a)')" +++++++++++++++++++++++++++++++++++++++++++++++++++ "
          write(*,'(a,i4,a,i3,a,i3)')" CPU ",cpuid,&
            &" START PSO SIMULATION number ",jgrid," / ",nGlobal
          allocate(xpar(nfit))
          xpar=par
          fitness=zero
          call cpu_time(start1)
          call pso_driver(jgrid,evfpso,nfit,allpar,minpar,maxpar,xpar,inv_fitness) ! PSO DRIVER
          fitness=one/inv_fitness
          par=xpar
          call cpu_time(end1)
          call timer(start1,end1,hour1,minute1,sec1)
          write(*,*)""
          write(*,'(a,i4,a,i4,a,i3,a,f5.2,a)')&
            &" CPU ",cpuid," DONE PSO in ",&
            &hour1," h : ",minute1," m : ",sec1," s"
          write(*,*)""
          write(*,'(a,i4,a,a,a)')" CPU ",cpuid," PSO ",&
            &trim(adjustl(string(jgrid)))," OUTPUT PARAMETERS"
          call write_lm_inputpar(cpuid,xpar)
          write(*,*)""
          write(*,'(2(a,g25.14),a,a,a,g25.14)')&
            &" inv_fitness = ",inv_fitness," => fitness = ",fitness," => fitness*",&
            &trim(adjustl(string(dof)))," = ",fitness*real(dof,dp)
          write(*,*)""
          write(*,'(a)') " -----------------------------------------------------"
          deallocate(xpar)

        end if

        ! write files with parameters from global search, without LM fit
        if(.not.allocated(resw)) allocate(resw(ndata),iwa(nfit),copar(nfit),sigpar(nfit))
        resw=zero
        copar=zero
        sigpar=zero
        info=0

        ! IT WRITES RESULTS TO SCREEN AND FILES
        write(*,*)
        write(*,'(a)') " -----------------------------------------------------"
        call cpu_time(start1)
        call ode_out(cpuid,jgrid,0,allpar,par,resw)
        call cpu_time(end1)
        call timer(start1,end1,hour1,minute1,sec1)
        write(*,'(a)') " -----------------------------------------------------"
        write(*,*)""
        write(*,'(a,i4,a,a,a,i4,a,i3,a,f5.2,a)')&
          &" CPU ",cpuid," DONE INTEGRATION ",&
          &trim(adjustl(string(jgrid)))," AND WRITE in ",&
          &hour1," h : ",minute1," m : ",sec1," s"
        write(*,*)""
        write(*,'(a)') " -----------------------------------------------------"
        if(ndata.gt.0)then
          call write_par(cpuid,par,sigpar,resw)
          call write_par(cpuid,jgrid,0,par,sigpar,resw)
          !             fmt=adjustl("(2(a,"//trim(sprec)//"))")
          !             write(*,trim(fmt))" Chi Square = ",sum(resw*resw),&
          !                  &" => Reduced Chi Square = ",(sum(resw*resw)/real(dof,dp))
          write(*,'(a)') " -----------------------------------------------------"
          write(*,*)""
        end if

        ! -- LM --- alone or after PIK/PSO
        if((progtype.eq.2).or.(lmon.eq.1))then
          write(*,'(a,i4,a)')" CPU ",cpuid," CALLING LM with INPUT PARAMETERS:"
          write(*,*)""
          call write_lm_inputpar(cpuid,par)
          write(*,*)""
          write(*,'(a)') " -----------------------------------------------------"
          call cpu_time(start1)
          call lm_driver(jgrid,ode_lm,allpar,ndata,nfit,par,&
            &resw,copar,sigpar,info,iwa) ! IT CALLS L-M
          call cpu_time(end1)
          call timer(start1,end1,hour1,minute1,sec1)
          write(*,'(a)') " -----------------------------------------------------"
          write(*,*)""
          write(*,'(a,i4,a,i4,a,i3,a,f5.2,a)')&
            &" CPU ",cpuid," DONE L-M in ",&
            &hour1," h : ",minute1," m : ",sec1," s"
          write(*,*)""
          write(*,'(a)') " -----------------------------------------------------"
          write(*,*)""
          call write_par(cpuid,par,sigpar,resw)
          call write_par(cpuid,jgrid,1,par,sigpar,resw)
          write(*,*)""
          write(*,'(a)') " -----------------------------------------------------"
          write(*,*)""
          flush(6)

          ! IT WRITES FIT RESULTS TO SCREEN AND FILES
          if(.not.allocated(resw)) allocate(resw(ndata))
          resw=zero
          write(*,*)
          write(*,'(a)') " -----------------------------------------------------"
          call cpu_time(start1)
          call ode_out(cpuid,jgrid,1,allpar,par,resw)
          call cpu_time(end1)
          call timer(start1,end1,hour1,minute1,sec1)
          write(*,'(a)') " -----------------------------------------------------"
          write(*,*)""
          write(*,'(a,i4,a,a,a,i4,a,i3,a,f5.2,a)')&
            &" CPU ",cpuid," DONE INTEGRATION ",&
            &trim(adjustl(string(jgrid)))," AND WRITE in ",&
            &hour1," h : ",minute1," m : ",sec1," s"
          write(*,*)""
          write(*,'(a)') " -----------------------------------------------------"
          !             if(ndata.gt.0)then
          !                fmt=adjustl("(2(a,"//trim(sprec)//"))")
          !                write(*,trim(fmt))" Chi Square = ",sum(resw*resw),&
          !                     &" => Reduced Chi Square = ",(sum(resw*resw)/real(dof,dp))
          !                write(*,'(a)') " -----------------------------------------------------"
          !               write(*,*)""
          !             end if

        end if

        deallocate(resw)
        if(allocated(iwa)) deallocate(iwa,copar,sigpar)

        ! BOOTSTRAP ANALYSIS
        if((progtype.lt.5).and.(nboot.gt.0))then
          write(*,'(a)') " -----------------------------------------------------"
          call cpu_time(start1)
          call strap_driver(jgrid,allpar,par)
          call cpu_time(end1)
          call timer(start1,end1,hour1,minute1,sec1)
          write(*,'(a)') " -----------------------------------------------------"
          write(*,*)""
          write(*,'(a,i4,a,i4,a,i3,a,f5.2,a)')&
            &" CPU ",cpuid," DONE BOOSTRAP in ",&
            &hour1," h : ",minute1," m : ",sec1," s"
          write(*,*)""
        end if

      end do doGlobal
    
    end if

    ! PolyChord - testing
    if(progtype.eq.5)then
      call set_par(m,R,P,a,e,w,mA,i,lN,PC_allpar,par) ! set common variable PC_allpar needed by loglikelihood
      call PC_driver(PolyChord_info)

    end if

  end if

  flush(6)

  !$ oend=omp_get_wtime()
  call cpu_time(enda)
  call timer(starta,enda,houra,minutea,seca)
  write(*,*)""
  write(*,'(a)') " -----------------------------------------------------"
  write(*,'(a)') " -----------------------------------------------------"
  write(*,'(a,i4,a,i3,a,f5.2,a)')" DONE TRADES IN ",&
    &houra," h : ",minutea," m : ",seca," s"
  write(*,'(a)') " -----------------------------------------------------"
  write(*,'(a)') " -----------------------------------------------------"
  !$ write(*,'(a)') " Execution time by omp_get_wtime() for openMP: "
  !$ call timer(ostart,oend,houra,minutea,seca)
  !$ write(*,'(a,i4,a,i3,a,f5.2,a)')" DONE PARALLEL TRADES IN ",&
  !$     &houra," h : ",minutea," m : ",seca," s"
  !$ write(*,'(a)') " -----------------------------------------------------"
  !$ write(*,'(a)') " -----------------------------------------------------"
  !$ write(*,*)""
  call date_and_time(VALUES=date_values)
  write(*,'(a,i4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2)')" END TIME: ",&
    &date_values(1),"-",date_values(2),"-",date_values(3)," -- ",&
    &date_values(5),":",date_values(6),":",date_values(7)
  write(*,'(a)')""

  if(allocated(bnames))   deallocate(bnames,bfiles)
  if(allocated(tofit))    deallocate(tofit)
  if(allocated(m))        deallocate(m,R,P,a,e,w,mA,i,lN,e_bounds)
  if(allocated(jdRV))     deallocate(jdRV,RVobs,eRVobs)
  if(allocated(epoT0obs)) deallocate(epoT0obs,T0obs,eT0obs)
  if(allocated(allpar))   deallocate(allpar,par)
  if(allocated(id))       deallocate(id,idall,parid)
  if(allocated(PC_allpar))deallocate(PC_allpar)
  if(allocated(par_min))  deallocate(par_min,par_max)
  if(allocated(minpar))   deallocate(minpar,maxpar)
  if(allocated(k_b))      deallocate(k_b)
end program




