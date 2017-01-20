module grid
  use constants,only:PR,zero,TOLERANCE,sprec,Mjups
  use parameters
  use init_trades,only:get_unit,init_random_seed
  use sorting,only:sort
  use celestial_mechanics,only:semax,period
  implicit none

contains

  ! ------------------------------------------------------------------ !
  ! set properly a parameter in the range
  subroutine grid_1par(Xmin,Xmax,Xdummy,type,NX,Xgrid)
    real(PR),intent(in)::Xmin,Xmax,Xdummy
    character(2),intent(in)::type
    integer,intent(out)::NX
    real(PR),dimension(:),allocatable,intent(out)::Xgrid
    real(PR)::step

    !if((Xdummy.eq.0._PR).or.(Xmax.le.Xmin))then
    if((abs(Xdummy).le.TOLERANCE)&
        &.or.(Xmax.le.Xmin))then
      ! when Xdummy = 0 set NX = 1, Xgrid(1) = Xmin
      NX=1
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      Xgrid=Xmin
      return
    end if

    ! Xdummy not 0
    if(type.eq.'rn')then
      ! Xdummy = random number of steps
      ! Xmin + Xmax + NX random steps
      NX = int(abs(Xdummy))+2
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      if(Xdummy.gt.0._PR)then
        ! linear scale
        call grid_rand(Xmin,Xmax,NX,Xgrid)
      else if(Xdummy.lt.0._PR)then
        ! logarithmic scale
        call grid_rand(log10(Xmin),log10(Xmax),NX,Xgrid)
        Xgrid=10._PR**Xgrid
      end if
    else if(type.eq.'ss')then
      ! Xdummy = step size
      ! NX ... + 1 so it uses at least the min parameter
      NX=int(((Xmax-Xmin)/abs(Xdummy))+0.5_PR)+1
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      if(Xdummy.gt.0._PR)then
        ! linear scale
        call grid_set(Xmin,Xdummy,NX,Xgrid)
      else if(Xdummy.lt.0._PR)then
        ! logarithmic scale
        call grid_set(log10(Xmin),log10(abs(Xdummy)),NX,Xgrid)
        Xgrid=10._PR**Xgrid
      end if
    else if(type.eq.'sn')then
      ! Xdummy = number of steps
      NX=int(abs(Xdummy))
      step=(Xmax-Xmin)/NX
      ! NX ... + 1 so it uses at least the min and max parameters
      NX=NX+1
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      if(Xdummy.gt.0._PR)then
        ! linear scale
        call grid_set(Xmin,step,NX,Xgrid)
      else if(Xdummy.lt.0._PR)then
        ! logarithmic scale
        step=(log10(xmax)-log10(xmin))/real((NX-1),PR)
        call grid_set(log10(Xmin),step,NX,Xgrid)
        Xgrid=10._PR**Xgrid
      end if
    else
      ! if nothing in type set NX = 1, Xgrid(1) = Xmin
      NX=1
      if(.not.allocated(Xgrid)) allocate(Xgrid(NX))
      Xgrid=Xmin
    end if

    return
  end subroutine grid_1par

  ! set the single parameter with NX random steps
  ! between Xmin and Xman
  subroutine grid_rand(Xmin,Xmax,NX,Xgrid)
    real(PR),intent(in)::Xmin,Xmax
    integer,intent(in)::NX
    real(PR),dimension(:),intent(out)::Xgrid
    real(PR)::dX
    real(PR),dimension(NX)::Xrand

    write(*,'(a)')" Creating a random grid with seeds:"
    call init_random_seed(NX)
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
    real(PR),intent(in)::Xmin,Xstep
    integer,intent(in)::NX
    real(PR),dimension(:),intent(out)::Xgrid
    integer::j

    do j=1,NX
      Xgrid(j)=Xmin+(j-1)*Xstep
    end do

    return
  end subroutine grid_set

  ! reads the parameters from the perturber body file: bfiles(idpert)
  subroutine grid_read_1(cpuid,&
      &Pmin,Pmax,dP,Ptype,aXmin,aXmax,da,atype,&
      &emin,emax,de,etype,wmin,wmax,dw,wtype)
    integer,intent(in)::cpuid
    real(PR),intent(out)::Pmin,Pmax,dP,aXmin,aXmax,da,emin,emax,de,wmin,wmax,dw
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
      read(uread,*) Pmin,Pmax,dP,Ptype
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
      &Pmin,Pmax,dP,Ptype,aXmin,aXmax,da,atype,&
      &emin,emax,de,etype,wmin,wmax,dw,wtype)
    integer,intent(in)::cpuid
    real(PR),intent(out)::Mmin,Mmax,dM,Pmin,Pmax,dP,aXmin,aXmax,da,&
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
      read(uread,*) Pmin,Pmax,dP,Ptype
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
    real(PR),dimension(:),intent(in)::m
    integer,intent(out)::NPgrid,Nagrid,Negrid,Nwgrid
    real(PR),dimension(:),allocatable,intent(out)::Pgrid,egrid,wgrid,agrid
    real(PR)::Pmin,Pmax,dP,aXmin,aXmax,da,emin,emax,de,wmin,wmax,dw
    integer::j
    character(2)::Ptype,atype,etype,wtype

    !call grid_read(cpuid,Pmin,Pmax,dP,Ptype,aXmin,aXmax,da,atype,&
    !     &emin,emax,de,etype,wmin,wmax,dw,wtype)
    call grid_read_1(cpuid,Pmin,Pmax,dP,Ptype,&
        aXmin,aXmax,da,atype,emin,emax,de,etype,wmin,wmax,dw,wtype)
    if(Pmin.ge.9.e6_PR)then
      call grid_1par(aXmin,aXmax,da,atype,Nagrid,agrid)
      NPgrid=Nagrid
      if(.not.allocated(Pgrid)) allocate(Pgrid(NPgrid))
      a2P: do j=1,Nagrid
        Pgrid(j)=period(m(1),m(idpert),agrid(j))
      end do a2P
    else if(aXmin.ge.999.0_PR)then
      call grid_1par(Pmin,Pmax,dP,Ptype,NPgrid,Pgrid)
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
      &Pmin,Pmax,dP,Ptype,aXmin,aXmax,da,atype)
    integer,intent(in)::cpuid
    integer,intent(out)::Nmgrid,Negrid,Nwgrid
    real(PR),intent(out)::Pmin,Pmax,dP,aXmin,aXmax,da
    real(PR),dimension(:),allocatable,intent(out)::Mgrid,egrid,&
        &wgrid
    real(PR)::Mmin,Mmax,dM,emin,emax,de,wmin,wmax,dw
    character(2)::Mtype,Ptype,atype,etype,wtype

    !call grid_read(cpuid,Pmin,Pmax,dP,Ptype,aXmin,aXmax,da,atype,&
    !     &emin,emax,de,etype,wmin,wmax,dw,wtype)
    call grid_read_2(cpuid,Mmin,Mmax,dM,Mtype,Pmin,Pmax,dP,Ptype,&
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
    real(PR),dimension(:),intent(in)::m,P,e,w
    integer,intent(out)::NPgrid,Nagrid,Negrid,Nwgrid,Ngrid
    real(PR),dimension(:,:),intent(out),allocatable::grid
    real(PR),dimension(:),allocatable::Pgrid,agrid,egrid,wgrid
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
    !grid(:,5)=huge(0._PR)
    grid(:,5)=1.e6_PR
    deallocate(Pgrid,agrid,egrid,wgrid)

    return
  end subroutine grid_build_1

  ! it creates the final grid array from the 'small' grid vectors
  ! - mass included -
  subroutine grid_build_2(cpuid,m,P,a,e,w,Ngrid,grid)
    implicit none
    integer,intent(in)::cpuid
    real(PR),dimension(:),intent(in)::m,P,a,e,w
    integer,intent(out)::Ngrid
    real(PR),dimension(:,:),intent(out),allocatable::grid
    real(PR),dimension(:),allocatable::xgrid,Mgrid,egrid,wgrid
    integer::NMgrid,Nxgrid,Negrid,Nwgrid
    real(PR)::Pmin,Pmax,dP,aXmin,aXmax,da
    character(2)::Ptype,atype
    integer::ijkgrid,ii,jj,kk,ll

    ! IN CASE OF A PERTURBER IT DEFINES THE VECTORS OF PARAMETERS
    if((idpert.ge.2).and.(idpert.le.NB))then
      write(*,'(" CALLING GRID_INIT with: ",a)') bfiles(idpert)
      call grid_def_2(cpuid,NMgrid,Negrid,Nwgrid,Mgrid,egrid,wgrid,&
          &Pmin,Pmax,dP,Ptype,aXmin,aXmax,da,atype)
      !write(*,'(2(3g25.15,a))')Pmin,Pmax,dP,Ptype,aXmin,aXmax,da,atype
      if(aXmin.ge.999._PR)then
        call grid_1par(Pmin,Pmax,dP,Ptype,Nxgrid,xgrid)
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
      if(a(NB).ge.999._PR)then
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
            if(aXmin.ge.999._PR)then
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
    grid(6,:)=1.e6_PR

    if(allocated(Mgrid)) deallocate(xgrid,Mgrid,egrid,wgrid)

    return
  end subroutine grid_build_2

  ! from the grid array to working parameter vectors
  ! - no mass -
  subroutine set_grid_par_1(jgrid,P,a,e,w,grid,Pw,aw,ew,ww)
    integer,intent(in)::jgrid
    real(PR),dimension(:),intent(in)::P,a,e,w
    real(PR),dimension(:,:),intent(in)::grid
    real(PR),dimension(:),intent(out)::Pw,aw,ew,ww
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
    real(PR),dimension(:),intent(in)::m,P,a,e,w
    real(PR),dimension(:,:),intent(in)::grid
    real(PR),dimension(:),intent(out)::mw,Pw,aw,ew,ww
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

end module grid
