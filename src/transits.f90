module transits
  use constants
  use parameters
  use celestial_mechanics,only:rsky,barycenter,fgfunctions,eleMD
  use lin_fit,only:linfit
  use convert_type,only:string
  implicit none

  interface assign_T0
    module procedure assign_T0_1,assign_T0_2
  end interface assign_T0

  contains

  ! ------------------------------------------------------------------ !
  ! given the T0 data it does a linear fit to the data and it finds the
  ! ephemeris T and P: tn = Tref + Pref*n
  subroutine set_ephem()
    integer,dimension(:),allocatable::x
    real(dp),dimension(:),allocatable::y,ey
    integer::j
    character(80)::fmt

    write(*,'(a)')" COMPUTING LINEAR EPHEMERIS OF: "
    if(.not.allocated(Tephem))&
        &allocate(Tephem(NB),Pephem(NB),eTephem(NB),ePephem(NB))
    Tephem=zero
    Pephem=zero
    eTephem=zero
    ePephem=zero
    do j=2,NB
      if(nT0(j).gt.0)then
        allocate(x(nT0(j)),y(nT0(j)),ey(nT0(j)))
        x=epoT0obs(1:nT0(j),j)
        y=T0obs(1:nT0(j),j)
        ey=eT0obs(1:nT0(j),j)
        call linfit(x,y,ey,Pephem(j),ePephem(j),Tephem(j),eTephem(j))
        fmt=adjustl("(a,i3,a,4("//trim(sprec)//",a))")
        write(*,trim(fmt))" body ",j,": t_N = (",&
            &Tephem(j),"+/-",eTephem(j),") + (",&
            &Pephem(j),"+/-",ePephem(j),") x N"
        deallocate(x,y,ey)
      end if
    end do

    return
  end subroutine set_ephem
  
  
  ! ------------------------------------------------------------------ !

  !
  ! Impact parameter of the transit
  !
  function impact_parameter(Rs,sma_p,inc_p,ecc_p,arg_p,R_p) result(b)
    real(dp)::b
    real(dp),intent(in)::Rs,sma_p,inc_p
    real(dp),optional,intent(in)::ecc_p,arg_p,R_p
    
    real(dp)::Rsum,aRs,rhoc
    
    Rsum=Rs*RsunAU
    if(present(R_p))Rsum=(Rs+R_p)*RsunAU
    aRs=sma_p/Rsum
    
    if(present(ecc_p).and.present(arg_p))then
    
      rhoc=(one-(ecc_p*ecc_p)) / (one+ecc_p*sin(arg_p*deg2rad))
    
    else
    
      rhoc=one
      
    end if

    b = aRs*rhoc*cos(inc_p*deg2rad)
    
    return
  end function impact_parameter
  
  
  ! move the complete state vectors (for all the bodies) of time = dt
  ! using fgfunctions subroutine
  subroutine advancefg(m,rw,dt,Hc)
    real(dp),dimension(:),intent(in)::m
    real(dp),dimension(:),intent(inout)::rw
    real(dp),intent(in)::dt
    logical,intent(inout)::Hc
!     logical::Hc
    
    integer::j1,j2,i1,i6
    real(dp)::mu

    do j1=2,NB
      j2=(j1-1)*6
      i1=1+j2
      i6=6+j2
      mu=Giau*(m(1)+m(j1))
      call fgfunctions(mu,rw(i1:i6),dt,Hc)
    end do

    return
  end subroutine advancefg
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! a bisection step, it updates a boundary and the step for the next iteration
  ! - transit -
  subroutine onetra_bis(itra,m,A,B,rw,dt,Hc)
    integer::itra
    real(dp),dimension(:),intent(in)::m
    real(dp),intent(inout)::A,B,dt
    real(dp),dimension(:),intent(inout)::rw
    logical,intent(inout)::Hc
!     logical::Hc
    
    !real(dp),dimension(:),allocatable::rin,drdt,err
    !real(dp),parameter::half=0.5_dp
    integer::i1,i2,i4,i5
    real(dp)::C

    call advancefg(m,rw,dt,Hc) ! let's use the integrator to get closer to the transit
    if(.not.Hc) return
    !allocate(rin(size(rw)),drdt(size(rw)),err(size(rw)))
    !rin=rw
    !call rkck_a(m,rin,drdt,dt,rw,err)
    !deallocate(rin,drdt,err)
    i1=1+(itra-1)*6
    i2=2+(itra-1)*6
    i4=4+(itra-1)*6
    i5=5+(itra-1)*6
    C=rw(i1)*rw(i4)+rw(i2)*rw(i5)
    if((A*C).lt.0._dp)then
      B=C
      !dt=-half*dt
      dt=-dt
    else
      A=C
      !dt=half*dt
    end if

    return
  end subroutine onetra_bis

  ! given a state vector and the itra of the body to check it computes the next step
  ! for the Newton-Raphson method
  ! - transit -
  subroutine onetra_nr(itra,m,rw,dt)
    use eq_motion,only:eqmastro
    integer,intent(in)::itra
    real(dp),dimension(:),intent(in)::m,rw
    real(dp),intent(out)::dt
    integer::ix,iy,ivx,ivy
    real(dp)::x,y,vx,vy,ax,ay,vx2,xax,vy2,yay,fk,dfk
    real(dp),dimension(:),allocatable::drw

    dt=zero
    ix=1+(itra-1)*6
    iy=2+(itra-1)*6
    ivx=4+(itra-1)*6
    ivy=5+(itra-1)*6
    x=rw(ix)
    y=rw(iy)
    vx=rw(ivx)
    vy=rw(ivy)
    ! Newton-Raphson
    ! f(k) = x*vx + y*vy
    fk=x*vx + y*vy
    allocate(drw(NBDIM))
    call eqmastro(m,rw,drw)
    ! df(k)/dt = vx^2 + x*ax + vy^2 + y*ay
    ax=drw(ivx)
    ay=drw(ivy)
    deallocate(drw)
    vx2=vx*vx
    xax=x*ax
    vy2=vy*vy
    yay=y*ay
    dfk=vx2+xax+vy2+yay
    dt=-(fk/dfk)

    return
  end subroutine onetra_nr

  ! given a state vector and the itra of the body to check,
  ! it computes the next step
  ! for the Newton-Raphson method
  ! - contact -
  subroutine onecont_nr(itra,Rcheck,rw,dt)
    integer,intent(in)::itra
    real(dp),intent(in)::Rcheck
    real(dp),dimension(:),intent(in)::rw
    real(dp),intent(out)::dt
    integer::ix,iy,ivx,ivy,jtra
    real(dp)::r_sky,rs2,hk,dhk,xvx,yvy

    dt=zero

    jtra=(itra-1)*6
    ix=1+jtra
    iy=2+jtra
    ivx=4+jtra
    ivy=5+jtra

    ! Newton-Raphson
    r_sky=rsky(rw(ix:iy))
    rs2=r_sky*r_sky
    hk=rs2-Rcheck

    xvx=rw(ix)*rw(ivx)
    yvy=rw(iy)*rw(ivy)
    dhk=2._dp*(xvx+yvy)
!     if (abs(dhk).le.TOLERANCE) write(*,*)" dhk <= TOLERANCE: ",dhk 
    if (abs(dhk).le.TOLERANCE)dhk=TOLERANCE
    dt=-(hk/dhk)

    return
  end subroutine onecont_nr

  ! a bisection step, it updates a boundary and the step for the next iteration
  ! - contact -
  subroutine onecont_bis(itra,Rcheck,m,A,B,rw,dt,Hc)
    integer::itra
    real(dp),intent(in)::Rcheck
    real(dp),dimension(:),intent(in)::m
    real(dp),intent(inout)::A,B,dt
    real(dp),dimension(:),intent(inout)::rw
    logical,intent(inout)::Hc
!     logical::Hc
    
    !real(dp),dimension(:),allocatable::rin,drdt,err
    real(dp),parameter::half=0.5_dp
    integer::i1,i2,i4,i5
    real(dp)::C

    call advancefg(m,rw,dt,Hc)
    if(.not.Hc) return
    !allocate(rin(size(rw)),drdt(size(rw)),err(size(rw)))
    !rin=rw
    !call rkck_a(m,rin,drdt,dt,rw,err)
    !deallocate(rin,drdt,err)
    i1=1+(itra-1)*6
    i2=2+(itra-1)*6
    i4=4+(itra-1)*6
    i5=5+(itra-1)*6
    C=rsky(rw(i1:i2))**2-Rcheck
    if((A*C).lt.0._dp)then
      B=C
      dt=-half*dt
    else
      A=C
      dt=half*dt
    end if

    return
  end subroutine onecont_bis

  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! seeks the root of rsky = 0 using an hybrid seeker: Bisection + Newton-Raphson
  subroutine find_transit(itra,m,r1,r2,itime,hok,tmidtra,lte,ro,Hc)
    integer,intent(in)::itra
    real(dp),dimension(:),intent(in)::m,r1,r2
    real(dp),intent(in)::itime,hok
    real(dp),intent(out)::tmidtra,lte
    real(dp),dimension(:),intent(inout)::ro
    logical,intent(inout)::Hc
!     logical::Hc
    
    real(dp),dimension(:),allocatable::rw,rwbar
    real(dp),dimension(6)::bar
    real(dp)::A,B,dt1,dt2
    integer::ix,iy,ivx,ivy
    integer::loop,many_iter

    if(.not.Hc)return
    
    ix=1+(itra-1)*6
    iy=2+(itra-1)*6
    ivx=4+(itra-1)*6
    ivy=5+(itra-1)*6
    ! A
    A=r1(ix)*r1(ivx) + r1(iy)*r1(ivy)
    ! B
    B=r2(ix)*r2(ivx) + r2(iy)*r2(ivy)
    dt1=hok
    dt2=zero
    tmidtra=zero
    allocate(rw(NBDIM),rwbar(NBDIM))
    ro=zero
    rw=r1
    loop=0
    
    many_iter=1000

    traloop: do
      loop=loop+1
      if(abs(dt1).le.TOLERANCE) exit traloop
      if(loop.ge.many_iter)then
        write(*,'(a,a,a)')" Reached ",trim(adjustl(string(many_iter))),&
        &"-th iteration in one_contact"
        exit traloop
      end if
      call onetra_nr(itra,m,rw,dt2)
      if(abs(dt2).le.abs(dt1))then
        ! continue with N-R
        dt1=dt2
        call advancefg(m,rw,dt1,Hc)
        if(.not.Hc) return
        tmidtra=tmidtra+dt1
      else
        ! dt2 >= dt1 so use Bisection
        dt1=dt1*0.5_dp
        tmidtra=tmidtra+dt1
        call onetra_bis(itra,m,A,B,rw,dt1,Hc)
        if(.not.Hc) return
      end if
      !tmidtra=tmidtra+dt1
    end do traloop

    call barycenter(m,rw,bar,rwbar)
    lte=-bar(3)/speedaud
    ro=rw
    tmidtra=tmidtra+tepoch+itime+lte
    deallocate(rw,rwbar)

    return
  end subroutine find_transit

  ! IT DEFINES THE RIGHT DIMENSION FOR THE RADIUS CHECK IN TRANSIT AND CONTACT TIMES
  subroutine Rbounds(itra,R,Rs,Rp,Rmin,Rmax)
    integer,intent(in)::itra
    real(dp),dimension(:),intent(in)::R
    real(dp),intent(out)::Rs,Rp,Rmin,Rmax

    Rs=R(1)*RsunAU
    Rp=R(itra)*RsunAU
    Rmax=Rs+Rp
    Rmin=Rs-Rp

    return
  end subroutine Rbounds

  ! IT CALCULATES A CONTACT TIME
  subroutine one_contact(icon,itra,m,R,rtra,ttra,tcont)
    integer,intent(in)::icon,itra
    real(dp),dimension(:),intent(in)::m,R,rtra
    real(dp),intent(in)::ttra
    real(dp),intent(out)::tcont
    real(dp),dimension(:),allocatable::rw,rwbar
    real(dp),dimension(6)::bar
    real(dp)::vmid,dt1,dt2,A,B
    real(dp)::Rp,Rs,Rmin,Rmax,Rcheck,tmid,tt,lte
    integer::jtra,ix,iy,ivx,ivy
    logical::Hc
    integer::loop

    jtra=(itra-1)*6
    ix=1+jtra
    iy=2+jtra
    ivx=4+jtra
    ivy=5+jtra

    call Rbounds(itra,R,Rs,Rp,Rmin,Rmax)
    vmid=rsky(rtra(ivx:ivy))
    dt1=-Rs/vmid
    dt2=zero
    if(icon.ge.3) dt1=-dt1
    Rcheck=Rmax*Rmax
    if((icon.eq.2).or.(icon.eq.3)) Rcheck=Rmin*Rmin

    allocate(rw(NBDIM),rwbar(NBDIM))
    call barycenter(m,rtra,bar,rwbar)
    lte=-bar(3)/speedaud
    tmid=ttra-lte
    bar=zero
    rwbar=zero
    tcont=zero
    tt=dt1
    rw=rtra
    
    Hc=.true.

    A=rsky(rtra(ix:iy))**2 - Rcheck
    B=A
    loop=0
    contloop: do
      loop=loop+1
      if(abs(dt1).le.TOLERANCE) exit contloop
      if(loop.ge.500)then
        write(*,'(a)')" Reached 500-th iteration in one_contact"
        exit contloop
      end if
      call onecont_nr(itra,Rcheck,rw,dt2)
      if(abs(dt2).le.abs(dt1))then
        dt1=dt2
        call advancefg(m,rw,dt1,Hc)
        if(.not.Hc) return          
      else
        call onecont_bis(itra,Rcheck,m,A,B,rw,dt1,Hc)
        if(.not.Hc) return
      end if
      tt=tt+dt1
    end do contloop
    call barycenter(m,rw,bar,rwbar)
    lte=-bar(3)/speedaud
    deallocate(rw,rwbar)
    tcont=tmid+tt+lte

    return
  end subroutine one_contact

  ! IT DETERMINES ALL CONTACT TIMES (IF THEY EXIST) OF TRANSIT
  subroutine find_contacts(itra,m,R,rtra,ttra,tcont)
    integer,intent(in)::itra
    real(dp),dimension(:),intent(in)::m,R,rtra
    real(dp),intent(in)::ttra
    real(dp),dimension(4),intent(out)::tcont
    real(dp)::r_sky,Rs,Rp,Rmin,Rmax
    integer::jcont,jtra,ix,iy,step

    jtra=(itra-1)*6
    ix=1+jtra
    iy=2+jtra
    r_sky=rsky(rtra(ix:iy))
    call Rbounds(itra,R,Rs,Rp,Rmin,Rmax)
    tcont=zero
    if(r_sky.le.Rmax)then
      step=3
      if(r_sky.lt.Rmin) step=1
      !write(*,'(a,i2)')" contact step ",step
      do jcont=1,4,step
        !write(*,'(a,i2)',advance='no')" contact id ",jcont
        call one_contact(jcont,itra,m,R,rtra,ttra,tcont(jcont))
        !write(*,'(a,g25.15)')" tcont = ",tcont(jcont)
      end do
      if(r_sky.ge.Rmin)then
        !write(*,'(2(a,g18.15))',advance='no')" r_sky ",r_sky," >= Rmin ",Rmin
        tcont(2)=tcont(1)
        tcont(3)=tcont(4)
        !write(*,'(a)',advance='no')" swapping tcont: 2-->1, 3-->4 "
      end if
    end if

    return
  end subroutine find_contacts

  ! call find_transit to compute the transit time (TT) and it assigns the right place
  ! of the TT comparing with the observations
  !! TO DO THE FIT OF THE DURATION OF TRANSIT, SO WE NEED TO FIND CONTACTS AND THEN
  !! DECIDE IF DURATION IS D_TOT = t4-t1, D_FULL = t3-t1 or D_HALF = t3.5-t1.5
  subroutine check_T0(itra,m,R,r1,r2,itime,hok,T0_stat,T0_sim,Hc)
    integer,intent(in)::itra
    real(dp),dimension(:),intent(in)::m,R,r1,r2
    real(dp),intent(in)::itime,hok
    integer,dimension(:,:),intent(inout)::T0_stat
    real(dp),dimension(:,:),intent(inout)::T0_sim
    logical,intent(inout)::Hc
!     logical::Hc ! removed intent(inout) to check
    
    real(dp)::tmidtra,lte,r_sky,Rs,Rp,Rmin,Rmax
    real(dp),dimension(:),allocatable::rtra
    real(dp),dimension(4)::tcont
    integer::nTs,jtra,ix,iy

    real(dp)::mu,P_p,sma_p,ecc_p,inc_p,mA_p,w_p,lN_p,f_p,dtau_p,b_itra
    
    allocate(rtra(NBDIM))
    rtra=one
    call find_transit(itra,m,r1,r2,itime,hok,tmidtra,lte,rtra,Hc)
    
    if(.not.Hc)then ! get out, if Hc == .false. is bad, so stop running
      deallocate(rtra)
      return
    
    else
    
      ! compute impact parameter from orbital elements at the transit:
      ! rtra -> kep. elem.
      jtra=(itra-1)*6
      mu=Giau*(m(1)+m(itra))
      call eleMD(mu,rtra(1+jtra:6+jtra),P_p,sma_p,ecc_p,inc_p,mA_p,w_p,lN_p,f_p,dtau_p)
      b_itra = abs(impact_parameter(R(1),sma_p,inc_p*rad2deg,ecc_p=ecc_p,arg_p=w_p*rad2deg))
      
      if((b_itra.lt.one).and.(.not.do_transit(itra)))then ! planet should not transit
!         write(*,'(a,i2,5(a,f23.10),a,l2)')' FOUND for planet ',itra,&
!           &' : i = ',inc_p*rad2deg,' a = ',sma_p,&
!           &' e = ',ecc_p,' w = ',w_p*rad2deg,&
!           &' ==> b = ',b_itra,&
!           &' < 1 && do_transit is ',do_transit(itra)
!         flush(6)
        Hc=.false.
        deallocate(rtra)
        return
      
      else
      
        call Rbounds(itra,R,Rs,Rp,Rmin,Rmax)
        
        ix=1+jtra
        iy=2+jtra
        r_sky=rsky(rtra(ix:iy))
        if(r_sky.gt.Rmax)then
          deallocate(rtra)
          return ! it is not transiting
        end if
        nTs=nT0(itra)

        if(nTs.gt.0) call assign_T0(itra,nTs,epoT0obs(1:nTs,itra),tmidtra,T0_stat,T0_sim)

    !     if(icont.eq.1) call find_contacts(itra,m,R,rtra,tmidtra,tcont)
        if(durcheck.eq.1) call find_contacts(itra,m,R,rtra,tmidtra,tcont)
        
        
        deallocate(rtra)
    
      end if
      
    end if

    return
  end subroutine check_T0

  ! IT DETERMINES WHICH IS THE RIGHT T_0,obs TO BE ASSOCIATED WITH
  ! THE SIMULATED T_0,sim = tmidtra
  ! v1
  subroutine assign_T0_1(itra,nTs,Tobs,tmidtra,T0_stat,T0_sim)
    integer,intent(in)::itra,nTs
    real(dp),dimension(:),intent(in)::Tobs
    real(dp),intent(in)::tmidtra
    integer,dimension(:,:)::T0_stat
    real(dp),dimension(:,:)::T0_sim
    integer,dimension(:),allocatable::intdt
    real(dp),dimension(:),allocatable::dttemp
    integer,dimension(1)::pos
    real(dp)::dtsel,Pchk
    real(dp),parameter::limP=1._dp/3._dp

    allocate(intdt(nTs),dttemp(nTs))
    dttemp=Tobs-tmidtra
    intdt=abs(int(dttemp))
    pos=minloc(intdt)
    dtsel=abs(dttemp(pos(1)))
    Pchk=Pephem(itra)*limP
    if(dtsel.le.Pchk)then
      if(T0_stat(pos(1),itra).eq.0)then
        T0_sim(pos(1),itra)=tmidtra
        T0_stat(pos(1),itra)=1
      end if
    end if
    deallocate(intdt,dttemp)

    return
  end subroutine assign_T0_1

  ! IT DETERMINES WHICH IS THE RIGHT T_0,obs TO BE ASSOCIATED WITH
  ! THE SIMULATED T_0,sim = tmidtra
  ! v2
  subroutine assign_T0_2(itra,nTs,epoTobs,tmidtra,T0_stat,T0_sim)
    integer,intent(in)::itra,nTs
    integer,dimension(:),intent(in)::epoTobs
    real(dp),intent(in)::tmidtra
    integer,dimension(:,:),intent(inout)::T0_stat
    real(dp),dimension(:,:),intent(inout)::T0_sim
    real(dp)::dT,dTP
    integer::ntmid,in

    dT=tmidtra-Tephem(itra)
    dTP=dT/Pephem(itra)
    if(dT.ge.0._dp)then
      dTP=dTP+0.5_dp
    else
      dTP=dTP-0.5_dp
    end if
    ntmid=int(dTP)
    do in=1,nTs
      if(ntmid.eq.epoTobs(in))then
        T0_sim(in,itra)=tmidtra
        T0_stat(in,itra)=1
      end if
    end do

    return
  end subroutine assign_T0_2

  ! it finds all transits of the selected planet (body id is itra) and store them
  ! in storetra variable, ready to be write into file
  subroutine all_transits(pos,itra,m,R,r1,r2,itime,hok,stat_tra,storetra)
    integer,intent(in)::pos,itra
    real(dp),dimension(:),intent(in)::m,R,r1,r2
    real(dp),intent(in)::itime,hok
    integer,dimension(:,:),intent(out)::stat_tra
    real(dp),dimension(:,:),intent(out)::storetra
    real(dp)::ttra,lte,Rs,Rp,Rmin,Rmax,r_sky
    real(dp),dimension(:),allocatable::rtra
    real(dp),dimension(4)::tcont
    logical::Hc
    integer::jtra
    
    real(dp)::mu,P_p,sma_p,ecc_p,inc_p,mA_p,w_p,lN_p,f_p,dtau_p,b_itra
    
    Hc=.true.
    allocate(rtra(NBDIM))
    rtra=one
    call find_transit(itra,m,r1,r2,itime,hok,ttra,lte,rtra,Hc)
    
    if(.not.Hc)then
      ttra=zero
      lte=zero
      tcont=zero
      rtra=zero
      
    else
    
      jtra=(itra-1)*6
      mu=Giau*(m(1)+m(itra))
      call eleMD(mu,rtra(1+jtra:6+jtra),P_p,sma_p,ecc_p,inc_p,mA_p,w_p,lN_p,f_p,dtau_p)
      b_itra = abs(impact_parameter(R(1),sma_p,inc_p*rad2deg,ecc_p=ecc_p,arg_p=w_p*rad2deg))
      
      if((b_itra.lt.one).and.(.not.do_transit(itra)))then ! planet should not transit
!         write(*,'(a,i2,5(a,f23.10),a,l2)')' FOUND for planet ',itra,&
!           &' : i = ',inc_p*rad2deg,' a = ',sma_p,&
!           &' e = ',ecc_p,' w = ',w_p*rad2deg,&
!           &' ==> b = ',b_itra,&
!           &' < 1 && do_transit is ',do_transit(itra)
!         flush(6)
        Hc=.false.
        deallocate(rtra)
        return
      
      else
      
        call Rbounds(itra,R,Rs,Rp,Rmin,Rmax)
        r_sky=rsky(rtra(1+jtra:2+jtra))
        if(r_sky.gt.Rmax) then ! it is not transiting
          ttra=zero
          lte=zero
          tcont=zero
          rtra=zero
        
        else
        
          call find_contacts(itra,m,R,rtra,ttra,tcont)

        end if

      end if
            
    end if
    
    stat_tra(itra,pos)=1
    storetra(1,pos)=ttra
    storetra(2,pos)=lte
    storetra(3:6,pos)=tcont
    storetra(7:,pos)=rtra
    deallocate(rtra)

    return
  end subroutine all_transits
  ! ------------------------------------------------------------------ !

end module transits
