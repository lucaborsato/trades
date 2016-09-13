module grid_search
  use constants,only:dp,zero,TOLERANCE,sprec,Mjups,Rjups
  use parameters
  use parameters_conversion
  use convert_type,only:string
  use init_trades,only:get_unit
  use random_trades
  use sorting,only:sort
  use celestial_mechanics,only:semax,period,semax_vec,period_vec,tau2mA,mA2tau,tau2mA_vec,mA2tau_vec
  implicit none

  contains

  ! ------------------------------------------------------------------ !
  ! set properly a parameter in the range
  subroutine grid_1par(Xmin,Xmax,Xdummy,step_type,NX,Xgrid)
    real(dp),intent(in)::Xmin,Xmax,Xdummy
    character(2),intent(in)::step_type
    integer,intent(out)::NX
    real(dp),dimension(:),allocatable,intent(out)::Xgrid
    real(dp)::step

    !if((Xdummy.eq.zero).or.(Xmax.le.Xmin))then
    if((abs(Xdummy).le.TOLERANCE)&
        &.or.(Xmax.le.Xmin))then
      ! when Xdummy = 0 set NX = 1, Xgrid(1) = Xmin
      NX=1
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      Xgrid=Xmin
      return
    end if

    ! Xdummy not 0
    if(step_type.eq.'rn')then
      ! Xdummy = random number of steps
      ! Xmin + Xmax + NX random steps
      NX = int(abs(Xdummy))+2
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      if(Xdummy.gt.zero)then
        ! linear scale
        call grid_rand(Xmin,Xmax,NX,Xgrid)
      else if(Xdummy.lt.zero)then
        ! logarithmic scale
        call grid_rand(log10(Xmin),log10(Xmax),NX,Xgrid)
        Xgrid=10._dp**Xgrid
      end if
      
    else if(step_type.eq.'ss')then
      ! Xdummy = step size
      ! NX ... + 1 so it uses at least the min parameter
!       NX=int(((Xmax-Xmin)/abs(Xdummy))+0.5_dp)+1
      NX=nint((Xmax-Xmin)/abs(Xdummy))+1
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      if(Xdummy.gt.zero)then
        ! linear scale
        call grid_set(Xmin,Xdummy,NX,Xgrid)
      else if(Xdummy.lt.zero)then
        ! logarithmic scale
        call grid_set(log10(Xmin),log10(abs(Xdummy)),NX,Xgrid)
        Xgrid=10._dp**Xgrid
      end if
      
    else if(step_type.eq.'sn')then
      ! Xdummy = number of steps
      NX=int(abs(Xdummy))
      step=(Xmax-Xmin)/NX
      ! NX ... + 1 so it uses at least the min and max parameters
      NX=NX+1
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      if(Xdummy.gt.zero)then
        ! linear scale
        call grid_set(Xmin,step,NX,Xgrid)
      else if(Xdummy.lt.zero)then
        ! logarithmic scale
        step=(log10(xmax)-log10(xmin))/real((NX-1),dp)
        call grid_set(log10(Xmin),step,NX,Xgrid)
        Xgrid=10._dp**Xgrid
      end if
    else
      ! if nothing in step_type set NX = 1, Xgrid(1) = Xmin
      NX=1
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      Xgrid=Xmin
    end if

    return
  end subroutine grid_1par

  ! set the single parameter with NX random steps
  ! between Xmin and Xman
  subroutine grid_rand(Xmin,Xmax,NX,Xgrid)
    real(dp),intent(in)::Xmin,Xmax
    integer,intent(in)::NX
    real(dp),dimension(:),intent(out)::Xgrid
    real(dp)::dX
    real(dp),dimension(NX)::Xrand

!     write(*,'(a)')" Creating a random grid with seeds:"
    call init_random_seed_clock(NX)
    dX=Xmax-Xmin
    call random_number(Xrand)
    Xgrid=Xmin+dX*Xrand
    Xgrid(1)=Xmin
    Xgrid(NX)=Xmax
    call sort(Xgrid)

    return
  end subroutine grid_rand

  ! set single parameter with NX steps between Xmin e Xman
  subroutine grid_set(Xmin,Xstep,NX,Xgrid)
    real(dp),intent(in)::Xmin,Xstep
    integer,intent(in)::NX
    real(dp),dimension(:),intent(out)::Xgrid
    integer::j

    do j=1,NX
      Xgrid(j)=Xmin+(j-1)*Xstep
    end do

    return
  end subroutine grid_set

  ! reads the parameters from the perturber body file: bfiles(idpert)
  subroutine grid_read_1(cpuid,&
      &Pmin,Pmax,deltaP,Ptype,aXmin,aXmax,da,atype,&
      &emin,emax,de,etype,wmin,wmax,dw,wtype)
    integer,intent(in)::cpuid
    real(dp),intent(out)::Pmin,Pmax,deltaP,aXmin,aXmax,da,emin,emax,de,wmin,wmax,dw
    character(2),intent(out)::Ptype,atype,etype,wtype
    integer::uread
    character(512)::fpert
    logical::fstat

    fpert = trim(adjustl(trim(path)//trim(bfiles(idpert))))
    inquire(file=trim(fpert),exist=fstat)
    if(fstat)then
      ! read parameters
      uread=get_unit(cpuid)
      open(uread,file=trim(fpert),status='OLD')
      read(uread,*)
      read(uread,*)
      read(uread,*) Pmin,Pmax,deltaP,Ptype
      read(uread,*) aXmin,aXmax,da,atype
      read(uread,*) emin,emax,de,etype
      read(uread,*) wmin,wmax,dw,wtype
      close(uread)
    else
      write(*,'(a,a,a)')" ERROR in grid_read: file ",&
          &trim(adjustl(fpert))," does NOT exist"
      stop
    end if

    return
  end subroutine grid_read_1

  ! reads the parameters from the perturber body file: bfiles(idpert)
  ! with this it will read the mass!!
  ! I will change it as soon as possible in order to provide
  ! the possibility to read all the parameters in the grid way
  subroutine grid_read_2(cpuid,&
      &Mmin,Mmax,dM,Mtype,&
      &Pmin,Pmax,deltaP,Ptype,aXmin,aXmax,da,atype,&
      &emin,emax,de,etype,wmin,wmax,dw,wtype)
    integer,intent(in)::cpuid
    real(dp),intent(out)::Mmin,Mmax,dM,Pmin,Pmax,deltaP,aXmin,aXmax,da,&
        &emin,emax,de,wmin,wmax,dw
    character(2),intent(out)::Mtype,Ptype,atype,etype,wtype
    integer::uread
    character(512)::fpert
    logical::fstat

    fpert = trim(adjustl(trim(path)//trim(bfiles(idpert))))
    inquire(file=trim(fpert),exist=fstat)
    if(fstat)then
      ! read parameters
      uread=get_unit(cpuid)
      open(uread,file=trim(fpert),status='OLD')
      read(uread,*) Mmin,Mmax,dM,Mtype
      read(uread,*)
      read(uread,*) Pmin,Pmax,deltaP,Ptype
      read(uread,*) aXmin,aXmax,da,atype
      read(uread,*) emin,emax,de,etype
      read(uread,*) wmin,wmax,dw,wtype
      close(uread)
      Mmin=Mmin*Mjups
      Mmax=Mmax*Mjups
      if(Mtype.eq.'ss') dM=dM*Mjups
    else
      write(*,'(a,a,a)')" ERROR in grid_read: file ",&
          &trim(adjustl(fpert))," does NOT exist"
      stop
    end if

    return
  end subroutine grid_read_2
  

  ! reads the parameter files and defines al the parameter grids
  subroutine grid_def_1(cpuid,m,NPgrid,Nagrid,Negrid,Nwgrid,Pgrid,agrid,egrid,wgrid)
    integer,intent(in)::cpuid
    real(dp),dimension(:),intent(in)::m
    integer,intent(out)::NPgrid,Nagrid,Negrid,Nwgrid
    real(dp),dimension(:),allocatable,intent(out)::Pgrid,egrid,wgrid,agrid
    real(dp)::Pmin,Pmax,deltaP,aXmin,aXmax,da,emin,emax,de,wmin,wmax,dw
    integer::j
    character(2)::Ptype,atype,etype,wtype

    !call grid_read(cpuid,Pmin,Pmax,deltaP,Ptype,aXmin,aXmax,da,atype,&
    !     &emin,emax,de,etype,wmin,wmax,dw,wtype)
    call grid_read_1(cpuid,Pmin,Pmax,deltaP,Ptype,&
        aXmin,aXmax,da,atype,emin,emax,de,etype,wmin,wmax,dw,wtype)
    if(Pmin.ge.9.e6_dp)then
      call grid_1par(aXmin,aXmax,da,atype,Nagrid,agrid)
      NPgrid=Nagrid
      if(.not.allocated(Pgrid)) allocate(Pgrid(NPgrid))
      a2P: do j=1,Nagrid
        Pgrid(j)=period(m(1),m(idpert),agrid(j))
      end do a2P
    else if(aXmin.ge.999.0_dp)then
      call grid_1par(Pmin,Pmax,deltaP,Ptype,NPgrid,Pgrid)
      Nagrid=NPgrid
      if(.not.allocated(agrid)) allocate(agrid(Nagrid))
      P2a: do j=1,NPgrid
        agrid(j)=semax(m(1),m(idpert),Pgrid(j))
      end do P2a
    end if
    call grid_1par(emin,emax,de,etype,Negrid,egrid)
    call grid_1par(wmin,wmax,dw,wtype,Nwgrid,wgrid)

    return
  end subroutine grid_def_1

  ! reads the parameter files and defines al the parameter grids
  ! version with grid Mass
  subroutine grid_def_2(cpuid,Nmgrid,Negrid,Nwgrid,Mgrid,egrid,wgrid,&
      &Pmin,Pmax,deltaP,Ptype,aXmin,aXmax,da,atype)
    integer,intent(in)::cpuid
    integer,intent(out)::Nmgrid,Negrid,Nwgrid
    real(dp),intent(out)::Pmin,Pmax,deltaP,aXmin,aXmax,da
    real(dp),dimension(:),allocatable,intent(out)::Mgrid,egrid,&
        &wgrid
    real(dp)::Mmin,Mmax,dM,emin,emax,de,wmin,wmax,dw
    character(2)::Mtype,Ptype,atype,etype,wtype

    !call grid_read(cpuid,Pmin,Pmax,deltaP,Ptype,aXmin,aXmax,da,atype,&
    !     &emin,emax,de,etype,wmin,wmax,dw,wtype)
    call grid_read_2(cpuid,Mmin,Mmax,dM,Mtype,Pmin,Pmax,deltaP,Ptype,&
        &aXmin,aXmax,da,atype,emin,emax,de,etype,wmin,wmax,dw,wtype)

    ! create vector for M, e, and w 
    call grid_1par(Mmin,Mmax,dM,Mtype,Nmgrid,Mgrid)
    call grid_1par(emin,emax,de,etype,Negrid,egrid)
    call grid_1par(wmin,wmax,dw,wtype,Nwgrid,wgrid)

    return
  end subroutine grid_def_2

  ! it creates the final grid array from the 'small' grid vectors
  ! - mass not included here -
  subroutine grid_build_1(cpuid,m,P,e,w,NPgrid,Nagrid,Negrid,Nwgrid,Ngrid,grid)
    implicit none
    integer,intent(in)::cpuid
    real(dp),dimension(:),intent(in)::m,P,e,w
    integer,intent(out)::NPgrid,Nagrid,Negrid,Nwgrid,Ngrid
    real(dp),dimension(:,:),intent(out),allocatable::grid
    real(dp),dimension(:),allocatable::Pgrid,agrid,egrid,wgrid
    integer::ijkgrid,ii,jj,kk
    character(80)::fmt

    ! IN CASE OF A PERTURBER IT DEFINES THE VECTORS OF PARAMETERS
    if((idpert.ge.2).and.(idpert.le.NB))then
      write(*,'(" CALLING GRID_INIT with: ",a)') bfiles(idpert)
      call grid_def_1(cpuid,m,NPgrid,Nagrid,Negrid,Nwgrid,&
          &Pgrid,agrid,egrid,wgrid)
    else
      write(*,'(" ID perturber not in range 2-nbody")')
      write(*,'(" Set perturber to nbody = ",I3)') NB
      NPgrid=1
      Nagrid=1
      Negrid=1
      Nwgrid=1
      if(.not.allocated(Pgrid)) &
          &allocate(Pgrid(NPgrid),agrid(Nagrid),&
          &egrid(Negrid),wgrid(Nwgrid))
      Pgrid=P(NB)
      agrid(1)=semax(m(1),m(NB),Pgrid(1))
      egrid=e(NB)
      wgrid=w(NB)
    end if
    ! ALL THE POSSIBLE COMBINATION OF THE PARAMETERS GRID ARRAY
    Ngrid=NPgrid*Negrid*Nwgrid
    ! IT ALLOCATES THE SPACE FOR GRID PARAMETERS (P,a,e,w)
    if(.not.allocated(grid)) allocate(grid(Ngrid,5))
    write(*,'(4(a,I5))')" NPgrid = ",NPgrid," (= Nagrid) Negrid = ",&
        &Negrid," Nwgrid = ",Nwgrid," => Ngrid = ",Ngrid
    fmt=adjustl("(a,1000(1x,"//trim(sprec)//"))")
    write(*,trim(fmt)) " Pgrid = ",Pgrid
    write(*,trim(fmt)) " agrid = ",agrid
    write(*,trim(fmt)) " egrid = ",egrid
    write(*,trim(fmt)) " wgrid = ",wgrid
    ijkgrid=0
    do ii=1,NPgrid,1
      do jj=1,Negrid,1
        do kk=1,Nwgrid,1
          ijkgrid=ijkgrid+1
          grid(ijkgrid,1)=Pgrid(ii)
          grid(ijkgrid,2)=agrid(ii)
          grid(ijkgrid,3)=egrid(jj)
          grid(ijkgrid,4)=wgrid(kk)
        end do
      end do
    end do
    !grid(:,5)=huge(0._dp)
    grid(:,5)=1.e6_dp
    deallocate(Pgrid,agrid,egrid,wgrid)

    return
  end subroutine grid_build_1

  ! it creates the final grid array from the 'small' grid vectors
  ! - mass included -
  subroutine grid_build_2(cpuid,m,P,a,e,w,Ngrid,grid)
    implicit none
    integer,intent(in)::cpuid
    real(dp),dimension(:),intent(in)::m,P,a,e,w
    integer,intent(out)::Ngrid
    real(dp),dimension(:,:),intent(out),allocatable::grid
    real(dp),dimension(:),allocatable::xgrid,Mgrid,egrid,wgrid
    integer::NMgrid,Nxgrid,Negrid,Nwgrid
    real(dp)::Pmin,Pmax,deltaP,aXmin,aXmax,da
    character(2)::Ptype,atype
    integer::ijkgrid,ii,jj,kk,ll

    ! IN CASE OF A PERTURBER IT DEFINES THE VECTORS OF PARAMETERS
    if((idpert.ge.2).and.(idpert.le.NB))then
      write(*,'(" CALLING GRID_INIT with: ",a)') bfiles(idpert)
      call grid_def_2(cpuid,NMgrid,Negrid,Nwgrid,Mgrid,egrid,wgrid,&
          &Pmin,Pmax,deltaP,Ptype,aXmin,aXmax,da,atype)
      !write(*,'(2(3g25.15,a))')Pmin,Pmax,deltaP,Ptype,aXmin,aXmax,da,atype
      if(aXmin.ge.999._dp)then
        call grid_1par(Pmin,Pmax,deltaP,Ptype,Nxgrid,xgrid)
      else
        call grid_1par(aXmin,aXmax,da,atype,Nxgrid,xgrid)
      end if

    else
      ! no grid, everything set to 1 dimension
      write(*,'(" ID perturber not in range 2-nbody")')
      write(*,'(" Set perturber to nbody = ",I3)') NB
      ! only one combination
      NMgrid=1
      Nxgrid=1
      Negrid=1
      Nwgrid=1
      ! allocate vectors needed for grid
      if(.not.allocated(xgrid)) &
          &allocate(Mgrid(NMgrid),xgrid(Nxgrid),&
          &egrid(Negrid),wgrid(Nwgrid))
      ! assign default value
      Mgrid=m(NB) 
      if(a(NB).ge.999._dp)then
        xgrid=P(NB)
      else
        xgrid=a(NB)
      end if
      egrid=e(NB)
      wgrid=w(NB)
    end if

    ! ALL THE POSSIBLE COMBINATION OF THE PARAMETERS GRID ARRAY
    write(*,'(a)')" Number of single grid parameters"
    write(*,'(10(a,i5))')" NMgrid = ",NMgrid," Nxgrid = ",Nxgrid,&
        &" Negrid = ",Negrid," Nwgrid = ",Nwgrid
    Ngrid=NMgrid*Nxgrid*Negrid*Nwgrid
    write(*,'(a,i5)')" Total combination of parameters: Ngrid = ",Ngrid
    ! IT ALLOCATES THE SPACE FOR GRID PARAMETERS (M,P,a,e,w,Chi^2)
    if(.not.allocated(grid)) allocate(grid(6,Ngrid))
    grid=zero

    ijkgrid=0
    do ll=1,Nmgrid,1
      do ii=1,Nxgrid,1
        do jj=1,Negrid,1
          do kk=1,Nwgrid,1
            ijkgrid=ijkgrid+1
            grid(1,ijkgrid)=Mgrid(ll)
            if(aXmin.ge.999._dp)then
              ! computes the semi-major axis from the Period
              grid(2,ijkgrid)=xgrid(ii)
              grid(3,ijkgrid)=semax(m(1),Mgrid(ll),xgrid(ii))
            else
              ! else computes the Period from the semi-major axis
              grid(3,ijkgrid)=xgrid(ii)
              grid(2,ijkgrid)=period(m(1),Mgrid(ll),xgrid(ii))
            end if
            grid(4,ijkgrid)=egrid(jj)
            grid(5,ijkgrid)=wgrid(kk)
          end do
        end do
      end do
    end do
    !grid(6,:)=huge(zero)
    grid(6,:)=1.e6_dp

    if(allocated(Mgrid)) deallocate(xgrid,Mgrid,egrid,wgrid)

    return
  end subroutine grid_build_2

  ! from the grid array to working parameter vectors
  ! - no mass -
  subroutine set_grid_par_1(jgrid,P,a,e,w,grid,Pw,aw,ew,ww)
    integer,intent(in)::jgrid
    real(dp),dimension(:),intent(in)::P,a,e,w
    real(dp),dimension(:,:),intent(in)::grid
    real(dp),dimension(:),intent(out)::Pw,aw,ew,ww
    integer::sel

    sel=idpert
    if( (idpert.lt.1).or.(idpert.gt.NB) ) sel=NB

    ! IT ASSIGNES THE RIGHT VALUE TO THE PARAMETERS
    Pw=P
    aw=a
    ew=e
    ww=w
    Pw(sel)=grid(jgrid,1)
    aw(sel)=grid(jgrid,2)
    ew(sel)=grid(jgrid,3)
    ww(sel)=grid(jgrid,4)

    return
  end subroutine set_grid_par_1

  ! from the grid array to working parameter vectors
  ! - mass -
  subroutine set_grid_par_2(jgrid,m,P,a,e,w,grid,mw,Pw,aw,ew,ww)
    integer,intent(in)::jgrid
    real(dp),dimension(:),intent(in)::m,P,a,e,w
    real(dp),dimension(:,:),intent(in)::grid
    real(dp),dimension(:),intent(out)::mw,Pw,aw,ew,ww
    integer::sel

    sel=idpert
    if( (idpert.lt.1).or.(idpert.gt.NB) ) sel=NB

    ! IT ASSIGNES THE RIGHT VALUE TO THE PARAMETERS
    mw=m
    Pw=P
    aw=a
    ew=e
    ww=w
    mw(sel)=grid(1,jgrid)
    Pw(sel)=grid(2,jgrid)
    aw(sel)=grid(3,jgrid)
    ew(sel)=grid(4,jgrid)
    ww(sel)=grid(5,jgrid)

    write(*,'(a,i5,a)')" PARAMETERS FOR THE GRID SEARCH ",jgrid," ITERATION "
    write(*,'(a,i5,a,g18.8)')" ITERATION ",jgrid," mw = ",mw(sel)
    write(*,'(a,i5,a,g18.8)')" ITERATION ",jgrid," Pw = ",Pw(sel)
    write(*,'(a,i5,a,g18.8)')" ITERATION ",jgrid," aw = ",aw(sel)
    write(*,'(a,i5,a,g18.8)')" ITERATION ",jgrid," ew = ",ew(sel)
    write(*,'(a,i5,a,g18.8)')" ITERATION ",jgrid," ww = ",ww(sel)

    return
  end subroutine set_grid_par_2

  
  ! -------------------------------------
  ! -------------------------------------
  ! NEW GRID SUBROUTINES - 2016-03-22
  ! -------------------------------------
  ! -------------------------------------
  
  
  ! -------------------------------------
  ! new subroutine to read input file and
  ! create grid for each parameter
  ! -------------------------------------
  subroutine read_parameters_grid(cpuid,parameters_grid)
    integer,intent(in)::cpuid
    type (parameter_grid),dimension(:),intent(out)::parameters_grid
    
    character(512)::perturber_file
    logical::file_stat
    integer::uread,i_read
    
    perturber_file=trim(adjustl(trim(path)//trim(bfiles(idpert))))
    ! path, bfiles, idpert in 'parameters' module
    inquire(file=trim(perturber_file),exist=file_stat)
    if(file_stat)then
    
      uread=get_unit(cpuid)
      open(uread,file=trim(perturber_file),status='OLD')
      do i_read=1,10
        read(uread,*) parameters_grid(i_read)%input_values,parameters_grid(i_read)%step_type
      end do
      close(uread)

      ! fix mass unit from M_Jup to M_sun
      parameters_grid(1)%input_values(1:2)=parameters_grid(1)%input_values(1:2)*Mjups
      if(parameters_grid(1)%step_type.eq.'ss') parameters_grid(1)%input_values(3)=parameters_grid(1)%input_values(3)*Mjups
      
      ! fix radius unit from R_Jup to R_sun
      parameters_grid(2)%input_values(1:2)=parameters_grid(2)%input_values(1:2)*Rjups
      if(parameters_grid(2)%step_type.eq.'ss') parameters_grid(2)%input_values(3)=parameters_grid(2)%input_values(3)*Rjups
      
      
    else

      write(*,'(a,a,a)')" ERROR in grid_read: file ",&
        &trim(adjustl(perturber_file))," does NOT exist"
      stop
      
    end if
    
    return
  end subroutine read_parameters_grid
  

  ! driver for grid_1par to be used with type parameter_grid
  subroutine values_to_type_grid(single_grid)
    type (parameter_grid)::single_grid
  
  ! call grid_1par(Xmin,Xmax,Xdummy,step_type,NX,Xgrid)
    call grid_1par(single_grid%input_values(1),& ! min
                  &single_grid%input_values(2),& ! max
                  &single_grid%input_values(3),& ! input step
                  &single_grid%step_type,&       ! step size, if zero no grid or random number 'rn'
                  &single_grid%n_steps,&         ! number of steps
                  &single_grid%grid_values)      ! grid for the parameter

    ! n_steps > 1
    if(single_grid%n_steps.gt.1)then
      
      ! random steps: single_grid%step_type = 'rn'
      if(single_grid%step_type.eq.'rn')then
        single_grid%step_grid=zero
      
      ! single_grid%step_type = 'ss/sn'
      else
        
        ! linear scale: single_grid%input_values(3) > 0
        if(single_grid%input_values(3).gt.zero)then
          single_grid%step_grid=single_grid%grid_values(2)-single_grid%grid_values(1)
        
        ! log scale: single_grid%input_values(3) < 0
        else
        
          single_grid%step_grid=log10(single_grid%grid_values(2))-log10(single_grid%grid_values(1))
        
        end if
      end if
      
    ! n_steps <= 1
    else
      single_grid%step_grid=zero
    end if
    
                  
    return
  end subroutine values_to_type_grid
  
  ! set the parameter_grid to zero/default: i.e., sma not used in the grid
  subroutine zero_default_grid(single_grid)
    type (parameter_grid)::single_grid
    integer::n_one
  
    single_grid%input_values=9.e8_dp
    single_grid%step_type='xx'
    n_one=1
    single_grid%n_steps=n_one
    allocate(single_grid%grid_values(n_one))
    single_grid%grid_values=zero
    single_grid%step_grid=zero
  
    return
  end subroutine zero_default_grid
  
  ! -----------------------------------------------------
  ! new subrotine to set properly the values for each
  ! parameter, i.e. name, n_steps, step_grid, grid_values
  ! -----------------------------------------------------
  subroutine set_parameters_grid(parameters_grid)
    type (parameter_grid),dimension(:)::parameters_grid
    
    integer::jpar,ipar
    
    ipar=(idpert-2)*8
    
    ! 1 == Mass
    parameters_grid(1)%name='m'//trim(string(idpert))
    call values_to_type_grid(parameters_grid(1))
    if(parameters_grid(1)%n_steps.gt.1)then
      par_min(ipar+3)=parameters_grid(1)%input_values(1)
      par_max(ipar+3)=parameters_grid(1)%input_values(2)
    end if
    
    ! 2 == Radius
    parameters_grid(2)%name='R'//trim(string(idpert))
    call values_to_type_grid(parameters_grid(2))
    if(parameters_grid(2)%n_steps.gt.1)then
      par_min(ipar+4)=parameters_grid(2)%input_values(1)
      par_max(ipar+4)=parameters_grid(2)%input_values(2)
    end if
    
    ! 3 == Period / 4 == semi-major axis [a]
    parameters_grid(3)%name='P'//trim(string(idpert))
    parameters_grid(4)%name='a'//trim(string(idpert))
    ! Pmin < 9.e6_dp --> create Period grid and then keep unset semi-major axis: n_steps = 1, step_type = 'xx'
    if(parameters_grid(3)%input_values(1).lt.9.e6_dp)then
      call values_to_type_grid(parameters_grid(3))
      call zero_default_grid(parameters_grid(4))
      if(parameters_grid(3)%n_steps.gt.1)then
        par_min(ipar+5)=parameters_grid(3)%input_values(1)
        par_max(ipar+5)=parameters_grid(3)%input_values(2)
      end if
      
    ! Pmin >= 9.e6_dp --> no Period ... check semi-major axis
    else
    
      ! sma-min < 999. --> create semi-major axis grid and then keep unset period: n_steps = 1, step_type = 'xx'
      if(parameters_grid(4)%input_values(1).lt.999._dp)then
        call values_to_type_grid(parameters_grid(4))
        call zero_default_grid(parameters_grid(3))
        if(parameters_grid(4)%n_steps.gt.1)then
          par_min(ipar+5)=period(MR_star(1,1),parameters_grid(1)%input_values(2),parameters_grid(4)%input_values(1))
          par_max(ipar+5)=period(MR_star(1,1),parameters_grid(1)%input_values(1),parameters_grid(4)%input_values(2))
        end if
      
      ! Pmin >= 9.e6_dp & sma-min >= 999.: STOP
      else
        write(*,'(a,i3)')' WARNING: SET PERIOD >= 9.e6 AND SEMI-MAJOR AXIS >= 999 FOR PLANET WITH ID = ',idpert
        stop
      end if
    end if
    
    ! 5 = Eccentricity
    parameters_grid(5)%name='e'//trim(string(idpert))
    call values_to_type_grid(parameters_grid(5))
    if(parameters_grid(5)%n_steps.gt.1)then
      par_min(ipar+6)=parameters_grid(5)%input_values(1)
      par_max(ipar+6)=parameters_grid(5)%input_values(2)
    end if
    
    ! 6 = omega / Argument of Pericenter / w
    parameters_grid(6)%name='w'//trim(string(idpert))
    call values_to_type_grid(parameters_grid(6))
    if(parameters_grid(6)%n_steps.gt.1)then
      par_min(ipar+7)=parameters_grid(6)%input_values(1)
      par_max(ipar+7)=parameters_grid(6)%input_values(2)
    end if
    
    ! 7 = Mean Anomaly / nu / 8 = Time pericenter / tau
    parameters_grid(7)%name='mA'//trim(string(idpert))
    parameters_grid(8)%name='tau'//trim(string(idpert))
    ! mAmin < 999. --> create Mean Anomaly grid, unset tau: n_steps = 1
    if(parameters_grid(7)%input_values(1).lt.999._dp)then
      call values_to_type_grid(parameters_grid(7))
      call zero_default_grid(parameters_grid(8))
      if(parameters_grid(6)%n_steps.gt.1)then
        par_min(ipar+8)=parameters_grid(7)%input_values(1)
        par_max(ipar+8)=parameters_grid(7)%input_values(2)
      end if
    
    ! mAmin >= 999. --> create Time Pericenter grid, unset mA: n_steps = 1
    else
    
      if(parameters_grid(8)%input_values(1).lt.9.e8_dp)then
        call values_to_type_grid(parameters_grid(8))
        call zero_default_grid(parameters_grid(7))
        if(parameters_grid(8)%n_steps.gt.1)then
          par_min(ipar+8)=tau2mA(parameters_grid(8)%input_values(1),tepoch,par_min(ipar+5))
          par_max(ipar+8)=tau2mA(parameters_grid(8)%input_values(2),tepoch,par_max(ipar+5))
        end if
      ! mAmin >= 999. & taumin >= 9.e8 : STOP
      else
        write(*,'(a,i3)')' WARNING: SET MEAN ANOMALY >= 999 AND TIME PERICENTER >= 9e8 FOR PLANET WITH ID = ',idpert
        stop
      end if
    
    end if
    
    ! 9 = Inclination
    parameters_grid(9)%name='i'//trim(string(idpert))
    call values_to_type_grid(parameters_grid(9))
    if(parameters_grid(9)%n_steps.gt.1)then
      par_min(ipar+9)=parameters_grid(9)%input_values(1)
      par_max(ipar+9)=parameters_grid(9)%input_values(2)
    end if
      
    ! 10 = Omega / Longitude of the ascending Node
    parameters_grid(10)%name='lN'//trim(string(idpert))
    call values_to_type_grid(parameters_grid(10))
    if(parameters_grid(10)%n_steps.gt.1)then
      par_min(ipar+10)=parameters_grid(10)%input_values(1)
      par_max(ipar+10)=parameters_grid(10)%input_values(2)
    end if
      
    call set_minmax()
    
    write(*,'(a)')' SET PARAMETERS GRID'
    write(*,'(a15,2(1x,a23))')'name','min','max'
    do jpar=3,10
      write(*,'(a15,2(1x,es23.16))')all_names_list(ipar+jpar),par_min(ipar+jpar),par_max(ipar+jpar)
    end do
    
    write(*,'(a15,2(1x,a23))')'name','fitmin','fitmax'
    do jpar=1,nfit
      write(*,'(a15,2(1x,es23.16))')parid(jpar),minpar(jpar),maxpar(jpar)
    end do
    write(*,*)
    flush(6)
    
    return
  end subroutine set_parameters_grid
  
  
  ! create the whole grid for the perturber
  subroutine build_grid(Mstar,parameters_grid,perturber_grid,fitness_grid,n_grid)
    real(dp),intent(in)::Mstar
    type (parameter_grid),dimension(:),intent(in)::parameters_grid
    integer,intent(out)::n_grid
    real(dp),dimension(:,:),allocatable,intent(out)::perturber_grid
    real(dp),dimension(:,:),allocatable,intent(out)::fitness_grid
  
    integer::i_col
    integer::n_spread,n_in,n_out
    
    n_grid=product(parameters_grid(:)%n_steps)
    allocate(perturber_grid(n_grid,10),fitness_grid(n_grid,2))
    fitness_grid=resmax ! resmax in 'parameters' module
    
    ! old bad way
!     do i_col=1,10
!       n_spread=n_grid/parameters_grid(i_col)%n_steps
!       perturber_grid(:,i_col)=reshape(spread(parameters_grid(i_col)%grid_values, 1, n_spread),(/n_grid/))
!     end do
!   
    ! first parameter M
    n_spread=n_grid/parameters_grid(1)%n_steps
    perturber_grid(:,1) = reshape(spread(parameters_grid(1)%grid_values,1,n_spread), (/n_grid/))
    
    ! other parameters: it already takes into account non-grid parameters such as P/a or mA/tau
    do i_col=2,9
      n_in=product(parameters_grid(i_col+1:)%n_steps)
      n_out=product(parameters_grid(1:i_col-1)%n_steps)
      perturber_grid(:,i_col)=reshape(&
        &spread(&
        &reshape(&
        &spread(parameters_grid(i_col)%grid_values,1,n_in),&
        &(/n_in*parameters_grid(i_col)%n_steps/)), 2, n_out),&
        &(/n_grid/))
    end do
    
    ! last parameter lN
    n_spread=n_grid/parameters_grid(10)%n_steps
    perturber_grid(:,10) = reshape(spread(parameters_grid(10)%grid_values,2,n_spread), (/n_grid/))

    ! set properly P or a, and mA or tau
    if(parameters_grid(4)%step_type.eq.'xx')then
      ! calculate semi-major axis from periods
      call semax_vec(Mstar,perturber_grid(:,1),perturber_grid(:,3),perturber_grid(:,4))
    else
      ! calculates periods from semi-major axis
      call period_vec(Mstar,perturber_grid(:,1),perturber_grid(:,4),perturber_grid(:,3))
    end if
    
    if(parameters_grid(8)%step_type.eq.'xx')then
      ! calculates tau from mean anomaly
      call mA2tau_vec(perturber_grid(:,7),tepoch,perturber_grid(:,3),perturber_grid(:,8))
    else
      ! calculates mean anomaly from tau
      call tau2mA_vec(perturber_grid(:,8),tepoch,perturber_grid(:,3),perturber_grid(:,7))
    end if

    return
  end subroutine build_grid
  
  subroutine perturber_grid2parameters(sim_id,id_perturber,perturber_grid,all_parameters)
    integer,intent(in)::sim_id,id_perturber
    real(dp),dimension(:,:),intent(in)::perturber_grid
    real(dp),dimension(:)::all_parameters ! no intent, because it will update only a part of it
    
    real(dp),dimension(8)::perturber_parameters
    integer,dimension(8)::list_id=(/1,2,3,5,6,7,9,10/)
    
    integer::id_start,id_end
    
    perturber_parameters=perturber_grid(sim_id,list_id)
    id_start=3+(id_perturber-2)*8
    id_end=id_start+7
!     write(*,*)' npar = ',npar
!     write(*,*)' id_perturber = ',id_perturber
!     write(*,*)' id_start = ',id_start
!     write(*,*)' id_end = ',id_end
    all_parameters(id_start:id_end)=perturber_parameters
  
    return
  end subroutine perturber_grid2parameters
  ! -------------------------------------
  ! -------------------------------------
  
end module grid_search
