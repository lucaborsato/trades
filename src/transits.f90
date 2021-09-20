module transits
  use constants
  use custom_type
  use parameters
  use celestial_mechanics,only:rsky,barycenter,fgfunctions,eleMD=>elem_mer
  use linear_ephem
  use convert_type,only:string
  use numerical_integrator,only:integrates_rk
  implicit none


!   interface assign_T0
!     module procedure assign_T0_byTime,assign_T0_byNumber
!   end interface assign_T0

  interface check_T0
    module procedure check_T0_1,check_T0_2
  end interface check_T0

  contains

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

  ! ------------------------------------------------------------------ !
  ! move the whole state vectors (for all the bodies) of time = dt
  ! using fgfunctions subroutine
  subroutine advancefg(m,rw,dt,Hc)
    real(dp),dimension(:),intent(in)::m
    real(dp),dimension(:),intent(inout)::rw
    real(dp),intent(in)::dt
    logical,intent(inout)::Hc

    integer::j1,j2,i1,i6
    real(dp)::mu

!     ! TEST: SPLIT IN 3 SUB-STEPS
!     real(dp)::dtx
!     real(dp),dimension(6)::rtemp

    integer::n_body
    n_body=size(m)

!     ! TEST
!     dtx=dt*onethird

!     do j1=2,NB
    do j1=2,n_body
      j2=(j1-1)*6
      i1=1+j2
      i6=6+j2
      mu=Giau*(m(1)+m(j1))
      call fgfunctions(mu,rw(i1:i6),dt,Hc)
!       ! TEST
!       rtemp=rw(i1:i6)
!       call fgfunctions(mu,rtemp,dtx,Hc)
!       call fgfunctions(mu,rtemp,dtx,Hc)
!       call fgfunctions(mu,rtemp,(dt-2*dtx),Hc)
!       rw(i1:i6)=rtemp

      ! 2018-05-16 try to not exit the loop, but move all the bodies
      if(.not.Hc) exit
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

!     real(dp),dimension(:),allocatable::rin,drdt,err
    integer::i1,i2,i4,i5
!     integer::dim_rin
    real(dp)::C

    ! let's use the F&G functions to get closer to the transit
    call advancefg(m,rw,dt,Hc)
!     if(.not.Hc) return

    ! let's use the integrator to get closer to the transit
!     dim_rin=size(rw)
!     allocate(rin(dim_rin),drdt(dim_rin),err(dim_rin))
!     rin=rw
!     call rkck_a(m,rin,drdt,dt,rw,err)
!     deallocate(rin,drdt,err)

    i1=1+(itra-1)*6
    i2=2+(itra-1)*6
    i4=4+(itra-1)*6
    i5=5+(itra-1)*6
    C=rw(i1)*rw(i4)+rw(i2)*rw(i5)
    if((A*C).lt.zero)then
      B=C
!       dt=-half*dt
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
    vx2=vx*vx
    xax=x*ax
    vy2=vy*vy
    yay=y*ay
    dfk=vx2+xax+vy2+yay
    dt=-(fk/dfk)
    deallocate(drw)

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
!     real(dp)::r_sky,rs2,hk,dhk,xvx,yvy
    real(dp)::rs2,hk,dhk,xvx,yvy

    dt=zero

    jtra=(itra-1)*6
    ix=1+jtra
    iy=2+jtra
    ivx=4+jtra
    ivy=5+jtra

    ! Newton-Raphson
!     r_sky=rsky(rw(ix:iy))
!     rs2=r_sky*r_sky
    rs2 = rw(ix)*rw(ix) + rw(iy)*rw(iy)
    hk=rs2-Rcheck

    xvx=rw(ix)*rw(ivx)
    yvy=rw(iy)*rw(ivy)
    dhk=two*(xvx+yvy)
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

!     real(dp),dimension(:),allocatable::rin
!     real(dp),dimension(:),allocatable::drdt,err
    integer::i1,i2,i4,i5
!     integer::dim_rin
    real(dp)::C

    ! let's use the F&G functions to get closer to the transit
    call advancefg(m,rw,dt,Hc)
!     if(.not.Hc) return

    ! let's use the integrator to get closer to the transit
!     dim_rin=size(rw)
!     allocate(rin(dim_rin),drdt(dim_rin),err(dim_rin))
!     rin=rw
!     call rkck_a(m,rin,drdt,dt,rw,err)
!     deallocate(rin,drdt,err)

    i1=1+(itra-1)*6
    i2=2+(itra-1)*6
    i4=4+(itra-1)*6
    i5=5+(itra-1)*6
    C=(rsky(rw(i1:i2))**2)-Rcheck
    if((A*C).lt.zero)then
      B=C
      dt=-half*dt
    else
      A=C
      dt=half*dt
    end if

    return
  end subroutine onecont_bis

  ! ------------------------------------------------------------------ !

  ! IT DEFINES THE RIGHT DIMENSION FOR THE RADIUS CHECK IN TRANSIT AND CONTACT TIMES
  subroutine Rbounds(itra,radii,Rs,Rp,Rmin,Rmax)
    integer,intent(in)::itra
    real(dp),dimension(:),intent(in)::radii
    real(dp),intent(out)::Rs,Rp,Rmin,Rmax

    Rs=radii(1)*RsunAU
    Rp=radii(itra)*RsunAU
    Rmax=Rs+Rp
    Rmin=Rs-Rp

    return
  end subroutine Rbounds

  ! ------------------------------------------------------------------ !
  ! seeks the root of rsky = 0 using an hybrid seeker: Bisection + Newton-Raphson
  subroutine find_transit(itra,m,r1,r2,itime,hok,tmidtra,lte,ro,Hc)
    integer,intent(in)::itra
    real(dp),dimension(:),intent(in)::m,r1,r2
    real(dp),intent(in)::itime,hok
    real(dp),intent(out)::tmidtra,lte
    real(dp),dimension(:),intent(inout)::ro
    logical,intent(inout)::Hc

    real(dp),dimension(:),allocatable::rw,rwbar
!     real(dp),dimension(:),allocatable::rwx,drdt,err
    real(dp),dimension(6)::bar
    real(dp)::A,B,dt1,dt2
    integer::ix,iy,ivx,ivy
    integer::loop,many_iter
    integer::n_body,nb_dim

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
!     dt1=hok*half
!     dt2=dt1
!     tmidtra=zero
    tmidtra=dt1
    
    n_body=size(m)
    nb_dim=n_body*6
    allocate(rw(nb_dim),rwbar(nb_dim))
!     allocate(rwx(nb_dim),drdt(nb_dim),err(nb_dim))
    ro=r1
    rw=r2

    ! 2018-05-17 add a first bisection iteration ...
!     tmidtra=tmidtra+dt1
!     call onetra_bis(itra,m,A,B,rw,dt1,Hc)

    loop=0

    many_iter=1000

    traloop: do
      loop=loop+1
      if(abs(dt1).le.TOLERANCE) exit traloop
      if(loop.ge.many_iter)then
        ! write(*,'(a,a,a)')" Reached ",trim(adjustl(string(many_iter))),&
        ! &"-th iteration in find_transit"
        exit traloop
      end if
      call onetra_nr(itra,m,rw,dt2)
      if(abs(dt2).le.abs(dt1))then
        ! continue with N-R
        dt1=dt2
!         call advancefg(m,rw,dt1,Hc)

        ! -- Luca 2018-06-18
        call integrates_rk(m,rw,dt1,ro)
        rw=ro
        ! --
        
! !         if(.not.Hc) write(*,'(a)')' find_transit: advancefg ==> Hc == False'
        ! 2017-11-22 use the integrator...
!         rwx=zero
!         drdt=zero
!         err=zero
!         call rkck_a(m,rw,drdt,dt1,rwx,err)
!         rw=rwx

        tmidtra=tmidtra+dt1
        if(.not.Hc) exit traloop
      else
        ! dt2 >= dt1 so use Bisection
!         dt1=dt1*half
!         tmidtra=tmidtra+dt1
!         call onetra_bis(itra,m,A,B,rw,dt1,Hc)

        ! -- LUCA 2018-06-18
        if(A*B.lt.zero) dt1 = -dt1*half
!         call advancefg(m,rw,dt1,Hc)

        ! -- Luca 2018-06-18
        call integrates_rk(m,rw,dt1,ro)
        rw=ro
        ! --
        
        tmidtra=tmidtra+dt1
        ! -- 
        
        if(.not.Hc) exit traloop
      end if
      ! -- LUCA 2018-06-18
      ! update B and A
      A = B
      B = rw(ix)*rw(ivx) + rw(iy)*rw(ivy)
      ! -- 
      !tmidtra=tmidtra+dt1
    end do traloop

    call barycenter(m,rw,bar,rwbar)
!     lte=-bar(3)/speedaud ! probrably wrong...
    lte = -rwbar(3)/speedaud
    ro=rw
    tmidtra=tmidtra+tepoch+itime+lte
    deallocate(rw,rwbar)
!     deallocate(rwx,drdt,err)

    return
  end subroutine find_transit

  ! IT CALCULATES A CONTACT TIME
  subroutine one_contact(icon,itra,mass,radii,rtra,ttra,tcont)
    integer,intent(in)::icon,itra
    real(dp),dimension(:),intent(in)::mass,radii,rtra
    real(dp),intent(in)::ttra
    real(dp),intent(out)::tcont
    real(dp),dimension(:),allocatable::ro,rw,rwbar
    real(dp),dimension(6)::bar
    real(dp)::vmid,dt1,dt2,A,B
    real(dp)::Rp,Rs,Rmin,Rmax,Rcheck,tmid,tt,lte
    integer::jtra,ix,iy,ivx,ivy
    logical::Hc
    integer::loop,many_iter

    jtra=(itra-1)*6
    ix=1+jtra
    iy=2+jtra
    ivx=4+jtra
    ivy=5+jtra

    ! icon = 0 ==> t_1.5
    ! icon = 1 ==> t_1
    ! icon = 2 ==> t_2
    ! icon = 3 ==> t_3
    ! icon = 4 ==> t_4
    ! icon = 5 ==> t_3.5

    call Rbounds(itra,radii,Rs,Rp,Rmin,Rmax)
    vmid=rsky(rtra(ivx:ivy))
    dt1=-Rs/vmid ! t_1.5, t_1, t_2: - sign
    dt2=zero
    if(icon.ge.3) dt1=-dt1 ! t_3, t_4, t_3.5: + sign

    if((icon.eq.2).or.(icon.eq.3))then
      Rcheck=Rmin*Rmin
    else if((icon.eq.0).or.(icon.eq.5))then
      Rcheck=Rs*Rs
    else
      Rcheck=Rmax*Rmax
    end if

    allocate(ro(NBDIM),rw(NBDIM),rwbar(NBDIM))
    call barycenter(mass,rtra,bar,rwbar)
!     lte=-bar(3)/speedaud
    lte=-rwbar(3)/speedaud
    tmid=ttra-lte
    bar=zero
    rwbar=zero
    tcont=zero
    tt=dt1
    ro=rtra
    rw=rtra

    Hc=.true.

    ! move the state vector of dt1 as suggested by Fabricky 2010
!     call advancefg(mass,rw,dt1,Hc)
    ! -- Luca 2018-06-18
    call integrates_rk(mass,rtra,dt1,rw)
    ! --

!     A=(rsky(rtra(ix:iy))**2)-Rcheck
    A=(rtra(ix)*rtra(ix)+rtra(iy)*rtra(iy))-Rcheck
    B=(rw(ix)*rw(ix)+rw(iy)*rw(iy))-Rcheck
!     A=(rw(ix)*rw(ix)+rw(iy)*rw(iy))-Rcheck
!     B=A

    many_iter = 1000
    loop=0
    contloop: do
      loop=loop+1
      if(abs(dt1).le.TOLERANCE) exit contloop
      if(loop.ge.many_iter)then
        ! write(*,'(a)')" Reached max iterations in one_contact"
        exit contloop
      end if
      call onecont_nr(itra,Rcheck,rw,dt2)
      if(abs(dt2).le.abs(dt1))then
        dt1=dt2
!         call advancefg(mass,rw,dt1,Hc)
        ! -- Luca 2018-06-18
        call integrates_rk(mass,rw,dt1,ro)
        rw=ro
        ! --
        tt=tt+dt1
        ! if(.not.Hc) return
        if(.not.Hc) exit contloop
      else ! bisection
!         tt=tt+dt1
!         call onecont_bis(itra,Rcheck,mass,A,B,rw,dt1,Hc)
        ! if(.not.Hc) return
        ! -- LUCA 2018-06-18
        if(A*B.lt.zero) dt1 = -dt1*half
!         call advancefg(mass,rw,dt1,Hc)
        ! -- Luca 2018-06-18
        call integrates_rk(mass,rw,dt1,ro)
        rw=ro
        ! --
        tt=tt+dt1
        ! -- 
        if(.not.Hc) exit contloop
      end if
!       tt=tt+dt1
!       if(.not.Hc) return
        ! -- LUCA 2018-06-18
        A = B
        B = (rw(ix)*rw(ix)+rw(iy)*rw(iy))-Rcheck
        ! -- 
    end do contloop
    
    call barycenter(mass,rw,bar,rwbar)
!     lte=-bar(3)/speedaud
    lte=-rwbar(3)/speedaud
    deallocate(rw,ro,rwbar)
    tcont=tmid+tt+lte

    return
  end subroutine one_contact

  ! IT DETERMINES ALL CONTACT TIMES (IF THEY EXIST) OF TRANSIT
  subroutine find_contacts(itra,mass,radii,rtra,ttra,tcont)
    integer,intent(in)::itra
    real(dp),dimension(:),intent(in)::mass,radii,rtra
    real(dp),intent(in)::ttra
    real(dp),dimension(4),intent(out)::tcont
    real(dp)::r_sky,Rs,Rp,Rmin,Rmax
    integer::jcont,jtra,ix,iy,step

    jtra=(itra-1)*6
    ix=1+jtra
    iy=2+jtra
    r_sky=rsky(rtra(ix:iy))
    call Rbounds(itra,radii,Rs,Rp,Rmin,Rmax)
    tcont=zero
    if(r_sky.le.Rmax)then
      step=3
      if(r_sky.lt.Rmin) step=1
      ! write(*,'(a,i2)')" contact step ",step
      do jcont=1,4,step
        ! write(*,'(a,i2)',advance='no')" contact id ",jcont
        call one_contact(jcont,itra,mass,radii,rtra,ttra,tcont(jcont))
        ! write(*,'(a,g25.15)')" tcont = ",tcont(jcont)
      end do
      if(r_sky.ge.Rmin)then
        ! write(*,'(2(a,g18.15))',advance='no')" r_sky ",r_sky," >= Rmin ",Rmin
        tcont(2)=tcont(1)
        tcont(3)=tcont(4)
        ! write(*,'(a)',advance='no')" swapping tcont: 2-->1, 3-->4 "
      end if
    end if

    return
  end subroutine find_contacts

  ! computes transit duration as:
  ! duration = t_3.5 - t_1.5
  ! where
  ! t_1.5 = time between contact time 1 and 2,
  ! t_3.5 = time between contact time 3 and 4,
  ! when project distance of the centre of the planet is on the edge of the star:
  ! rsky == Rstar (ingress <-> t_1.5, egress <-> t_3.5)
  subroutine compute_transit_duration_c2c(id_body,mass,radii,rtra,ttra,duration)
    integer,intent(in)::id_body
    real(dp),dimension(:),intent(in)::mass,radii,rtra
    real(dp),intent(in)::ttra
    real(dp),intent(out)::duration

    integer::sel_r
    real(dp)::r_sky,Rs,Rp,Rmin,Rmax
    real(dp)::t_hing,t_hegr

    t_hing=zero
    t_hegr=zero

    call Rbounds(id_body,radii,Rs,Rp,Rmin,Rmax)
    sel_r=(id_body-1)*6
    r_sky=rsky(rtra(1+sel_r:2+sel_r))

    if(r_sky.le.Rmax)then
      ! computes the t_1.5 == t_hing = planet on the edge of the star
      call one_contact(0,id_body,mass,radii,rtra,ttra,t_hing)
      call one_contact(5,id_body,mass,radii,rtra,ttra,t_hegr)

      ! in case of grazing? rsky > Rstar!
      ! computing duration = t_4 - t_1
!       call one_contact(1,id_body,mass,radii,rtra,ttra,t_hing)
!       call one_contact(4,id_body,mass,radii,rtra,ttra,t_hegr)

      duration = t_hegr - t_hing

    end if

   return
  end subroutine compute_transit_duration_c2c

  ! computes transit duration as in Kipping 2010, eq. (15):
  ! duration = T1 = (P/pi) * (rhoc^2 / sqrt(1-ecc^2)) * arcsin( sqrt(1 - (a/Rs)^2 * rhoc^2 * cos(inc)^2) / (a/Rs) * rhoc * sin(inc) )
  ! where rhoc = (1-ecc^2) / (1 +/- ecc*sin(argp)) with + <-> transit, - <-> occultation (W11_eq7-8, K10 rhoc)
  subroutine compute_transit_duration_K10_15(id_body,mass,radii,rtra,duration)
    integer,intent(in)::id_body
    real(dp),dimension(:),intent(in)::mass,radii,rtra
    real(dp),intent(out)::duration

    integer::sel_r
    real(dp)::mu,P_p,sma_p,ecc_p,inc_p,mA_p,w_p,lN_p,f_p,dtau_p
    real(dp)::ome2,rhoc,aRs,num1,den1,asin_num,asin_den


    mu=Giau*(mass(1)+mass(id_body))
    sel_r=(id_body-1)*6
    call eleMD(mu,rtra(1+sel_r:6+sel_r),P_p,sma_p,ecc_p,inc_p,mA_p,w_p,lN_p,f_p,dtau_p)

    ! 1 - ecc^2
    ome2=one-(ecc_p*ecc_p)
    rhoc=ome2/(one+ecc_p*sin(w_p*deg2rad)) ! in some cases the + would be change to - for occultation
    ! sma / Rstar
    aRs=sma_p/(radii(1)*RsunAU)
    ! T1 Kipping 2010 eq. 15: One-term expression
    num1=P_p*rhoc*rhoc
    den1=dpi*sqrt(ome2)
    asin_num=sqrt(one-aRs*aRs*rhoc*rhoc*cos(inc_p*deg2rad))
    asin_den=aRs*rhoc*sin(inc_p*deg2rad)
    duration = asin(asin_num/asin_den)*num1/den1

   return
  end subroutine compute_transit_duration_K10_15


  ! call find_transit to compute the transit time (TT) and it assigns the right place
  ! of the TT comparing with the observations
  ! Computes duration of the transit (does not assign it now) as in Kipping 2010 eq.15
!   subroutine check_T0_1(id_body,mass,radii,r1,r2,itime,hok,T0_stat,T0_sim,Hc)
subroutine check_T0_1(id_body,mass,radii,r1,r2,itime,hok,simT0,Hc)
    integer,intent(in)::id_body
    real(dp),dimension(:),intent(in)::mass,radii,r1,r2
    real(dp),intent(in)::itime,hok
!     integer,dimension(:,:),intent(inout)::T0_stat
!     real(dp),dimension(:,:),intent(inout)::T0_sim
    type(dataT0),dimension(:),intent(inout)::simT0
    logical,intent(inout)::Hc

    real(dp)::tmidtra,ttra_temp,lte,duration,r_sky,Rs,Rp,Rmin,Rmax
    real(dp),dimension(:),allocatable::rtra
    real(dp),dimension(4)::tcont
!     integer::nTs,sel_r
    integer::sel_r

!     write(*,'(a)')'check_T0_1'
!     flush(6)

    tmidtra=zero
    duration=zero
    tcont=zero
    allocate(rtra(NBDIM))
!     rtra=one
    rtra=r1

    call find_transit(id_body,mass,r1,r2,itime,hok,ttra_temp,lte,rtra,Hc)

    if(Hc)then

      call Rbounds(id_body,radii,Rs,Rp,Rmin,Rmax)
      sel_r=(id_body-1)*6
      r_sky=rsky(rtra(1+sel_r:2+sel_r))

      if(r_sky.le.Rmax)then ! planet transits the star (b <= 1)

        if(do_transit(id_body))then ! planet has to transit!

          tmidtra=ttra_temp+lte
          ! old
!           nTs=nT0(id_body)
!           if(durcheck.eq.1) call compute_transit_duration_K10_15(id_body,mass,radii,rtra,duration)
!           if(nTs.gt.0) call assign_T0(id_body,nTs,epoT0obs(1:nTs,id_body),tmidtra,T0_stat,T0_sim)
          if(durcheck.eq.1)then
            call find_contacts(id_body,mass,radii,rtra,ttra_temp,tcont) !tcont take into account LTE
            duration=(tcont(4)-tcont(1))*1440.0_dp
            ! write(*,*)"tcont: ",tcont
            ! write(*,*)"tmidtra: ",tmidtra," duration: ",duration
          end if
          if(obsData%obsT0(id_body-1)%nT0.gt.0) &
            &call assign_T0_byNumber(id_body,obsData,tmidtra,duration,simT0)

        else ! ... but it should not!

          Hc=.false.

        end if ! do_transit

      end if ! r_sky

    end if ! Hc

    deallocate(rtra)

    return
  end subroutine check_T0_1

    ! same as check_T0_1, but the epochs of the transit have to be provided
  ! as the array T0_num
!   subroutine check_T0_2(id_body,mass,radii,r1,r2,itime,hok,transit_flag,dur_check,n_T0,T0_num,T0_stat,T0_sim,Hc)
  subroutine check_T0_2(id_body,mass,radii,r1,r2,itime,hok,transit_flag,dur_check,oDataIn,simT0,Hc)
    integer,intent(in)::id_body
    real(dp),dimension(:),intent(in)::mass,radii,r1,r2
    real(dp),intent(in)::itime,hok
    logical,dimension(:),intent(in)::transit_flag
    integer,intent(in)::dur_check
!     integer,dimension(:),intent(in)::n_T0
!     integer,dimension(:,:),intent(in)::T0_num
!     integer,dimension(:,:),intent(inout)::T0_stat
!     real(dp),dimension(:,:),intent(inout)::T0_sim
    type(dataObs),intent(in)::oDataIn
    type(dataT0),dimension(:),intent(inout)::simT0
    logical,intent(inout)::Hc

    real(dp)::tmidtra,ttra_temp,lte,duration,r_sky,Rs,Rp,Rmin,Rmax
    real(dp),dimension(:),allocatable::rtra
    real(dp),dimension(4)::tcont
!     integer::nTs,sel_r
    integer::sel_r

    tmidtra=zero
    duration=zero
    tcont=zero
    allocate(rtra(NBDIM))
!     rtra=one
    rtra=r1

    call find_transit(id_body,mass,r1,r2,itime,hok,ttra_temp,lte,rtra,Hc)

    if(Hc)then

      call Rbounds(id_body,radii,Rs,Rp,Rmin,Rmax)
      sel_r=(id_body-1)*6
      r_sky=rsky(rtra(1+sel_r:2+sel_r))

      if(r_sky.le.Rmax)then ! planet transits the star (b <= 1)

        if(transit_flag(id_body))then ! planet has to transit!

          tmidtra=ttra_temp+lte
          ! old
!           nTs=n_T0(id_body)
!           if(dur_check.eq.1) call compute_transit_duration_K10_15(id_body,mass,radii,rtra,duration)
!           if(nTs.gt.0) call assign_T0(id_body,nTs,T0_num(1:nTs,id_body),tmidtra,T0_stat,T0_sim)
          if(dur_check.eq.1)then
            call find_contacts(id_body,mass,radii,rtra,ttra_temp,tcont) !tcont take into account LTE
            duration=(tcont(4)-tcont(1))*1440.0_dp
          end if
          if(obsData%obsT0(id_body-1)%nT0.gt.0) &
            &call assign_T0_byNumber(id_body,oDataIn,tmidtra,duration,simT0)

        else ! ... but it should not!

          Hc=.false.

        end if ! transit_flag

      end if ! r_sky

    end if ! Hc

    deallocate(rtra)

    return
  end subroutine check_T0_2


  ! IT DETERMINES WHICH IS THE RIGHT T_0,obs TO BE ASSOCIATED WITH
  ! THE SIMULATED T_0,sim = tmidtra
  ! v1
!   subroutine assign_T0_byTime(id_body,nTs,Tobs,tmidtra,T0_stat,T0_sim)
  subroutine assign_T0_byTime(id_body,oDataIn,tmidtra,duration,simT0)
    integer,intent(in)::id_body
    type(dataObs),intent(in)::oDataIn
    real(dp),intent(in)::tmidtra,duration
    type(dataT0),dimension(:),intent(inout)::simT0

    integer,dimension(:),allocatable::intdt
    real(dp),dimension(:),allocatable::dttemp
!     integer,dimension(1)::pos
    integer::pos
    real(dp)::dtsel,Pchk
    integer::nTs,ibd
!     real(dp),parameter::limP=1._dp/3._dp

    ibd=id_body-1
    nTs=oDataIn%obsT0(ibd)%nT0
    allocate(intdt(nTs),dttemp(nTs))
    dttemp=oDataIn%obsT0(ibd)%T0-tmidtra
    intdt=abs(int(dttemp))
!     pos=minloc(intdt)
    pos=minloc(intdt,dim=1)
    dtsel=abs(dttemp(pos))
!     Pchk=Pephem(id_body)*onethird
    Pchk=oDataIn%obsT0(ibd)%Pephem*onethird
    if(dtsel.le.Pchk)then
      if(simT0(ibd)%T0_stat(pos).eq.0)then
        simT0(ibd)%epo=oDataIn%obsT0(ibd)%epo(pos)
        simT0(ibd)%T0(pos)=tmidtra
        simT0(ibd)%T0_stat(pos)=1
        simT0(ibd)%nT0=simT0(ibd)%nT0+1
        if(durcheck.eq.1)then
          simT0(ibd)%dur(pos)=duration
          simT0(ibd)%dur_stat(pos)=1
          simT0(ibd)%nDur=simT0(ibd)%nDur+1
        end if
      end if
    end if
    deallocate(intdt,dttemp)

    return
  end subroutine assign_T0_byTime

  ! IT DETERMINES WHICH IS THE RIGHT T_0,obs TO BE ASSOCIATED WITH
  ! THE SIMULATED T_0,sim = tmidtra
  ! v2
  subroutine assign_T0_byNumber(id_body,oDataIn,tmidtra,duration,simT0)
    integer,intent(in)::id_body
    type(dataObs),intent(in)::oDataIn
    real(dp),intent(in)::tmidtra,duration
    type(dataT0),dimension(:),intent(inout)::simT0

    real(dp)::dT,dTP,Tref,Pref
    integer::ntmid,in,nTs,ibd,epo

    ibd=id_body-1
    nTs=oDataIn%obsT0(ibd)%nT0
    Tref=oDataIn%obsT0(ibd)%Tephem
    Pref=oDataIn%obsT0(ibd)%Pephem

    dT=tmidtra-Tref
    dTP=dT/Pref
    ntmid = nint(dTP)

    do in=1,nTs
      epo=oDataIn%obsT0(ibd)%epo(in)
      if(ntmid.eq.epo)then
        simT0(ibd)%epo(in)=ntmid
        simT0(ibd)%T0(in)=tmidtra
        simT0(ibd)%T0_stat(in)=1
        simT0(ibd)%nT0=simT0(ibd)%nT0+1
        if(durcheck.eq.1)then
          simT0(ibd)%dur(in)=duration
          simT0(ibd)%dur_stat(in)=1
          simT0(ibd)%nDur=simT0(ibd)%nDur+1
        end if
      end if
    end do

    return
  end subroutine assign_T0_byNumber

  ! it finds all transits of the selected planet (id_body) and store them
  ! in storetra variable, ready to be written into file
  subroutine all_transits(pos,id_body,mass,radii,r1,r2,itime,hok,stat_tra,storetra)
    integer,intent(in)::pos,id_body
    real(dp),dimension(:),intent(in)::mass,radii,r1,r2
    real(dp),intent(in)::itime,hok
    integer,dimension(:,:),intent(out)::stat_tra
    real(dp),dimension(:,:),intent(out)::storetra

    real(dp)::ttra,lte,Rs,Rp,Rmin,Rmax,r_sky
    real(dp),dimension(:),allocatable::rtra
    real(dp),dimension(4)::tcont
    logical::Hc
    integer::sel_r

!     real(dp)::mu,P_p,sma_p,ecc_p,inc_p,mA_p,w_p,lN_p,f_p,dtau_p,b_id_body

    ttra=zero
    lte=zero
    tcont=zero
    Hc=.true.
    allocate(rtra(NBDIM))
!     rtra=one
    rtra=r1
    call find_transit(id_body,mass,r1,r2,itime,hok,ttra,lte,rtra,Hc)

    if(Hc)then

      call Rbounds(id_body,radii,Rs,Rp,Rmin,Rmax)
      sel_r=(id_body-1)*6
      r_sky=rsky(rtra(1+sel_r:2+sel_r))

      if(r_sky.le.Rmax)then ! planet transits the star (b <= 1)

        if(do_transit(id_body))then ! planet has to transit!

          call find_contacts(id_body,mass,radii,rtra,ttra,tcont)
          stat_tra(id_body,pos)=1
          storetra(1,pos)=ttra
          storetra(2,pos)=lte
          storetra(3:6,pos)=tcont
          storetra(7:,pos)=rtra

        else

          Hc=.false.

        end if ! do_transit/transit_flag

      end if ! r_sky

    end if ! Hc

    deallocate(rtra)

    return
  end subroutine all_transits
  ! ------------------------------------------------------------------ !

! ==============================================================================
! compute the transit time and proper duration from K10 eq. 15
  subroutine transit_time(id_body,mass,radii,r1,r2,iter_time,step_ok,ttra,dur_tra,check_ttra)
    integer,intent(in)::id_body ! value: 2 to NB
    real(dp),dimension(:),intent(in)::mass,radii,r1,r2
    real(dp),intent(in)::iter_time,step_ok
    real(dp),intent(out)::ttra,dur_tra
    logical,intent(out)::check_ttra

    real(dp)::ttra_temp,lte
    real(dp),dimension(4)::tcont
    real(dp),dimension(:),allocatable::rtra

    real(dp)::Rs,Rp,Rmin,Rmax,r_sky
!     real(dp)::mu,P_p,sma_p,ecc_p,inc_p,mA_p,w_p,lN_p,f_p,dtau_p,b_p
    integer::sel_r

    check_ttra=.true.
    ttra=-9.e10_dp
    dur_tra=-9.e10_dp

    allocate(rtra(NBDIM))
    rtra=one
    call find_transit(id_body,mass,r1,r2,iter_time,step_ok,ttra_temp,lte,rtra,check_ttra)

    if(check_ttra)then

      call Rbounds(id_body,radii,Rs,Rp,Rmin,Rmax)
      sel_r=(id_body-1)*6
      r_sky=rsky(rtra(1+sel_r:2+sel_r))

      if(r_sky.le.Rmax)then ! planet transits the star (b <= 1)

        if(do_transit(id_body))then ! planet has to transit!

          ttra=ttra_temp+lte
!           call compute_transit_duration_K10_15(id_body,mass,radii,rtra,dur_tra)
          call find_contacts(id_body,mass,radii,rtra,ttra_temp,tcont) !tcont take into account LTE
          dur_tra=(tcont(4)-tcont(1))*1440.0_dp

        else ! ... but it should not!

          check_ttra=.false.

        end if ! do_transit

      end if ! r_sky

    end if ! check_ttra

    deallocate(rtra)

    return
  end subroutine transit_time


! ==============================================================================

end module transits
