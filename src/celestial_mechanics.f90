module celestial_mechanics
  use constants
  use parameters
  implicit none

  interface dist
    module procedure dist_1,dist_2
  end interface dist

  interface rHill
    module procedure rHill_1,rHill_2
  end interface rHill
  
  contains

!   function to compute the semi-major axis of an orbit from Period
  function semax(ms,mp,P) result(sma)
    real(dp)::sma
    real(dp),intent(in)::ms,mp,P
    real(dp)::mu,P2
    real(dp)::dpi2=dpi*dpi

    mu=Giau*(ms+mp)
    P2=P*P
    sma=((mu*P2)/dpi2)**(onethird)

    return
  end function semax
  
!   function to compute the period given the semi-major axis
  function period(ms,mp,a) result(a2P)
    real(dp)::a2P
    real(dp),intent(IN)::ms,mp,a
    real(dp)::mu,a3

    mu=Giau*(ms+mp)
    a3=a**3
    a2P=dpi*sqrt(a3/mu)

    return
  end function period

  ! vector version
  !   function to compute the semi-major axis of an orbit from Period
  subroutine semax_vec(ms,mp,P,sma)
    real(dp),intent(IN)::ms
    real(dp),dimension(:),intent(in)::mp,P
    real(dp),dimension(:),intent(out)::sma
    real(dp),dimension(:),allocatable::mu,P2
    real(dp)::dpi2=dpi*dpi

    allocate(mu(size(mp)),P2(size(P)))
    mu=Giau*(ms+mp)
    P2=P*P
    sma=((mu*P2)/dpi2)**(onethird)
    deallocate(mu,P2)
    
    return
  end subroutine semax_vec
  
!   function to compute the period given the semi-major axis
  subroutine period_vec(ms,mp,a,a2P)
    real(dp),intent(IN)::ms
    real(dp),dimension(:),intent(IN)::mp,a
    real(dp),dimension(:),intent(out)::a2P
    real(dp),dimension(:),allocatable::mu,a3

    allocate(mu(size(mp)),a3(size(a)))
    mu=Giau*(ms+mp)
    a3=a**3
    a2P=dpi*sqrt(a3/mu)
    deallocate(mu,a3)

    return
  end subroutine period_vec

  
  
  ! time pericentre to mean anomaly
  function tau2mA(tau, t_ref, per) result(mA)
    real(dp)::mA
    real(dp),intent(in)::tau,t_ref,per
    
    mA = mod( ((360._dp/per)*(t_ref-tau)) , 360._dp)
    
    return
  end function tau2mA
  
  ! mean anomaly to time pericentre
  function mA2tau(mA, t_ref, per) result(tau)
    real(dp)::tau
    real(dp),intent(in)::mA,t_ref,per
    
    tau = t_ref - (mA*per/360._dp)
    
    return
  end function mA2tau
  
  ! vector version
  ! time pericentre to mean anomaly
  subroutine tau2mA_vec(tau, t_ref, per, mA)
    real(dp),dimension(:),intent(in)::tau,per
    real(dp),intent(in)::t_ref
    real(dp),dimension(:),intent(out)::mA
    
    mA = mod( ((360._dp/per)*(t_ref-tau)) , 360._dp)
    
    return
  end subroutine tau2mA_vec
  
  ! mean anomaly to time pericentre
  subroutine mA2tau_vec(mA, t_ref, per, tau)
    real(dp),dimension(:),intent(in)::mA,per
    real(dp),intent(in)::t_ref
    real(dp),dimension(:),intent(out)::tau
    
    tau = t_ref - (mA*per/360._dp)
    
    return
  end subroutine mA2tau_vec

  
  
!   calculates Eccentric anomaly from meanAnomaly [deg] and eccentricity
  function EAnom(mA,ecc) result(EA)
    real(dp)::EA
    real(dp),intent(IN)::mA,ecc
    real(dp)::mArad,E,fE,dfE
    real(dp),parameter::tol=1.e-9_dp
    integer::i

    EA=zero
    mArad=mod(mA*deg2rad,dpi)
    if(mArad.lt.zero) mArad=mArad+dpi
    i=0
    E=mArad
    if(ecc.gt.TOLERANCE)then
      if(ecc.gt.0.6_dp)E=pi
      loopo: do
        i=i+1
        fE=E-ecc*sin(E)-mArad
        dfE=one-ecc*cos(E)
        EA=E-(fE/dfE)
        if(abs(E-EA).le.tol) exit loopo
        E=EA
      end do loopo
      EA=EA*rad2deg
    else
      EA=mA
    end if

    return
  end function EAnom

  function calculate_true_anomaly(mean_anomaly, eccentricity) result(true_anomaly)
    real(dp)::true_anomaly
    
    real(dp),intent(in)::mean_anomaly,eccentricity
    
    real(dp)::EA,tan_EA,ecc_coeff
    
    EA=EAnom(mean_anomaly,eccentricity)
    tan_EA=tan(EA*half)
    ecc_coeff=sqrt((one+eccentricity)/(one-eccentricity))
    true_anomaly = two * atan(ecc_coeff*tan_EA)
  
    return
  end function calculate_true_anomaly

  
!   calculate the module of a 3-D vector
  function dist_1(r) result(out)
    real(dp)::out
    real(dp),intent(in)::r(3)
    real(dp)::x,y,z

    x=r(1)
    y=r(2)
    z=r(3)

    out=sqrt(x*x + y*y + z*z)

    return
  end function dist_1

!   function to compute the distance between to vectors
  function dist_2(r1,r2) result(out)
    real(dp)::out
    real(dp),intent(in)::r1(3),r2(3)
    real(dp)::dx,dy,dz

    dx=r1(1)-r2(1)
    dy=r1(2)-r2(2)
    dz=r1(3)-r2(3)
    out=sqrt(dx*dx + dy*dy + dz*dz)

    return
  end function dist_2

!   module of a 2-D vector
  function rsky(r) result(out)
    real(dp)::out
    real(dp),intent(in)::r(2)
    real(dp)::x,y

    x=r(1)
    y=r(2)
    out=sqrt(x*x + y*y)

    return
  end function rsky

!   calculate module of Rvec x Vvec 
  function rvprod(r,v) result(out)
    real(dp)::out
    real(dp),dimension(3),intent(in)::r,v
    real(dp)::x,y,z,vx,vy,vz

    x=r(1)
    y=r(2)
    z=r(3)
    vx=v(1)
    vy=v(2)
    vz=v(3)

    out=x*vx + y*vy + z*vz

    return
  end function rvprod

! -----------------------
! bad Hill radius check
  
!   function to compute the Hill radius w/o eccentricity
  function rHill_1(ms,mp,sma) result(rH)
    real(dp)::rH
    real(dp),intent(in)::ms,mp,sma

    rH=sma*((onethird*(mp/ms))**onethird)

    return
  end function rHill_1

!   function to compute the Hill radius w/ eccentricity
  function rHill_2(ms,mp,sma,ecc) result(rH)
    real(dp)::rH
    real(dp),intent(in)::ms,mp,sma,ecc
    real(dp)::e1

    e1=one-ecc
    rH=sma*e1*((onethird*(mp/ms))**onethird)

    return
  end function rHill_2

!   given the state vector it check the Hill condition
  function Hillcheck(m,rin) result(rH)
    logical::rH
    real(dp),dimension(:),intent(in)::m,rin
    real(dp),dimension(3)::rvec1,vvec1,rvec2,vvec2
    real(dp)::r1,v1,r2,v2,rv1,rv2,a1,a2,dr,mu1,mu2,n1,es,ec,e1,rH1,rH2
    integer::j,j1,j2
    character(72)::fmt

    rH=.true.
    ! rin -> a -> RHill -> 1=overlap,bad, 0=ok
    do j=2,(NB-1)

      j1=(j-1)*6
      j2=j*6

      rvec1=rin(1+j1:3+j1)
      vvec1=rin(4+j1:6+j1)
      !rv1=rvec1(1)*vvec1(1)+rvec1(2)*vvec1(2)+rvec1(3)*vvec1(3)
      rv1=rvprod(rvec1,vvec1)
      r1=dist(rvec1)
      v1=dist(vvec1)
      mu1=Giau*(m(1)+m(j))
      a1=(two/r1)-((v1**2)/mu1) ! 1/a
      if(a1.le.zero)then
        fmt=adjustl('(a,1x,'//trim(sprec)//')')
        write(*,trim(fmt),advance='no')" 1/a1 < 0 :",a1
        rH=.false.
        return
      end if
      a1=one/a1
      n1=sqrt(mu1/a1**3)
      es=rv1/(n1*a1*a1)
      ec=((r1*v1*v1)/mu1)-one
      e1=sqrt(ec*ec+es*es)

      rvec2=rin(1+j2:3+j2)
      vvec2=rin(4+j2:6+j2)
      !rv2=rvec2(1)*vvec2(1)+rvec2(2)*vvec2(2)+rvec2(3)*vvec2(3)
      rv2=rvprod(rvec2,vvec2)
      r2=dist(rvec2)
      v2=dist(vvec2)
      mu2=Giau*(m(1)+m(j+1))
      a2=(two/r2)-((v2**2)/mu2) ! 1/a
      if(a2.le.zero)then
        fmt=adjustl('(a,1x,'//trim(sprec)//')')
        write(*,trim(fmt),advance='no')" 1/a2 < 0 :",a2
        rH=.false.
        return
      end if
      a2=one/a2
      n1=sqrt(mu1/a2**3)
      es=rv2/(n1*a2*a2)
      ec=((r2*v2*v2)/mu2)-one
      e1=sqrt(ec*ec+es*es)
      if((a1.le.amin).or.(a1.ge.amax).or.(a2.le.amin).or.(a2.ge.amax))then
        !write(*,'(a)')" Bound semi-major axis limits overcome "
        rH=.false.
        return
      end if
      rH1=rHill(m(1),m(j),a1,e1)
      rH2=rHill(m(1),m(j+1),a2,e1)
      dr=dist(rvec1,rvec2)
      if(dr.le.(rH1+rH2))then
        rH=.false.
        return
      end if
    end do

    return
  end function Hillcheck
  
! -----------------------

! -----------------------
! mutual Hill radius

  function mutual_Hill_radius(ms,mp1,sma1,mp2,sma2) result(rH)
    real(dp)::rH
    real(dp),intent(in)::ms,mp1,sma1,mp2,sma2
    real(dp)::sma_mean,mass_ratio
    
    real(dp)::min_ratio
    
    rH=zero
    sma_mean=half*(sma1+sma2)
    mass_ratio=(mp1+mp2)/ms
    min_ratio = TOLERANCE**onethird
    
    
    if(mass_ratio.le.TOLERANCE)then
    
      rH=sma_mean*min_ratio
      write(*,'(2(a,es23.16),a)')'mass_ratio = ',mass_ratio,' <= ',TOLERANCE,' = TOLERANCE ==> USING TOLERANCE'
      write(*,'(3(a,es23.16))')'Mstar [Msun] = ',ms,' Mi [Msun]= ',mp1,' Mj [Msun]= ',mp2
      flush(6)
    
    else if(ms.le.TOLERANCE)then
    
      write(*,'(2(a,es23.16))')' Mstar = ',ms,' <= ',TOLERANCE,' = TOLERANCE'
      write(*,'(3(a,es23.16))')'Mstar [Msun] = ',ms,' Mi [Msun]= ',mp1,' Mj [Msun]= ',mp2
      flush(6)
      stop
      
    else if(sma_mean.le.TOLERANCE)then
    
      write(*,'(2(a,es23.16))')' sma_mean = ',sma_mean,' <= ',TOLERANCE,' = TOLERANCE'
      write(*,'(3(a,es23.16))')'Mstar [Msun] = ',ms,' Mi [Msun]= ',mp1,' Mj [Msun]= ',mp2
      flush(6)
      stop
      
    else if(ms.le.(mp1+mp2))then
    
      write(*,'(2(a,es23.16))')' Mstar [Msun] = ',ms,' <= ',mp1+mp2,' = Mi + Mj  [Msun]'
      flush(6)
      stop
      
    else
    
      rH=sma_mean*((onethird*mass_ratio)**onethird)
    
    end if
      
    return
  end function mutual_Hill_radius

  subroutine rrdot_to_invsma(body_id,m,rin,inv_sma,check)
    integer,intent(in)::body_id
    real(dp),dimension(:),intent(in)::m,rin
    real(dp),intent(out)::inv_sma
    logical,intent(out)::check
    integer::ii
    real(dp),dimension(3)::rvec,vvec
    real(dp)::rv,rr,vv,mu
    
    ii=(body_id-1)*6
    rvec=rin(1+ii:3+ii)
    vvec=rin(4+ii:6+ii)
    rv=rvprod(rvec,vvec)
    rr=dist(rvec)
    vv=dist(vvec)
    mu=Giau*(m(1)+m(body_id))
    inv_sma=(two/rr)-((vv**2)/mu) ! 1/a
    check=.true.
    if(inv_sma.le.zero) check=.false. ! if false is bad

    return
  end subroutine rrdot_to_invsma
  
    
  function mutual_Hill_check_0(m,rin) result(hill_check)
    logical::hill_check
    real(dp),dimension(:),intent(in)::m,rin
    real(dp)::sma_i, sma_j
    real(dp)::Hill_radius_ij,delta_ij,stability_criterion
    integer::i,j
    
    hill_check=.true.
    
    do i=2,(NB-1)
      call rrdot_to_invsma(i,m,rin,sma_i,hill_check)
      if(.not.hill_check) return
      sma_i=one/sma_i
      
      do j=i+1,NB
        call rrdot_to_invsma(j,m,rin,sma_j,hill_check)
        if(.not.hill_check) return
        sma_j=one/sma_j
        
        Hill_radius_ij = mutual_Hill_radius(m(1),m(i),sma_i,m(j),sma_j)
        delta_ij = abs(sma_j - sma_i)
        stability_criterion = sqrt_12 * Hill_radius_ij
        if(delta_ij.lt.stability_criterion)then
          hill_check = .false. ! if false is bad/unstable
          return
        end if
      end do
    
    end do
  
    return
  end function mutual_Hill_check_0
  
  function mutual_Hill_check(m,rin) result(hill_check)
    logical::hill_check
    real(dp),dimension(:),intent(in)::m,rin
    real(dp)::sma_i, sma_j
    real(dp)::Hill_radius_ij,delta_ij,stability_criterion
    integer::i,j
    
    hill_check=.true.
    
    if(NB.gt.2)then
      do i=2,(NB-1)
        call rrdot_to_invsma(i,m,rin,sma_i,hill_check)
        if(.not.hill_check) return
        sma_i=one/sma_i

        j=i+1
        call rrdot_to_invsma(j,m,rin,sma_j,hill_check)
        if(.not.hill_check) return
        sma_j=one/sma_j

        Hill_radius_ij = mutual_Hill_radius(m(1),m(i),sma_i,m(j),sma_j)
        delta_ij = abs(sma_j - sma_i)
        stability_criterion = sqrt_12 * Hill_radius_ij
        if(delta_ij.lt.stability_criterion)then
          hill_check = .false. ! if false is bad/unstable
          return
        end if
      end do
    end if
        
    return
  end function mutual_Hill_check
  
  ! -------------------------------
  
  ! ------------------------------------------------------------------ !
  ! computes the initial state vector in the orbital reference frame
  subroutine initial_state(P,a,e,mA,output)
    real(dp),dimension(:),intent(in)::P,a,e,mA
    real(dp),dimension(:),intent(out)::output
    real(dp),parameter::circ=360._dp
    real(dp)::EA,dEA,n,cosE,sinE,e2,rad1e2,mA2
    integer::i,nci

    output=zero

    do i=2,NB
      !computation of the anomalies of the bodies
      mA2=mod(mA(i),circ)+circ
      n=dpi/P(i)
      EA=EAnom(mA2,e(i))*deg2rad
      e2=e(i)*e(i)
      rad1e2=sqrt(one-e2)
      cosE=cos(EA)
      sinE=sin(EA)
      dEA=n/(one-e(i)*cosE)
      !calculation of the radius vector and velocity of the bodies
      nci=(i-1)*6
      output(1+nci)=a(i)*(cosE-e(i))  !! x
      output(2+nci)=a(i)*rad1e2*sinE  !! y
      output(4+nci)=-a(i)*sinE*dEA  !! vx
      output(5+nci)=a(i)*rad1e2*cosE*dEA  !! vy
    end do

    return
  end subroutine initial_state

  ! ------------------------------------------------------------------ !
  ! function to determine if angle is in some predefined range
  function checklN(lN) result(clN)
    integer::clN
    real(dp),intent(in)::lN
    real(dp)::lN1

    clN=0
    lN1=mod(lN,360._dp)
    if(lN1.lt.zero) lN1=lN1+360._dp
    ! IT DEFINES THE COORDINATES FOR ECLIPSE CONDITION DUE TO THE LONGITUDE OF NODES
    if((lN1.ge.zero.and.lN1.le.45._dp)&
        &.or.&
        &(lN1.gt.135._dp.and.lN1.le.225._dp)&
        &.or.&
        &(lN1.ge.315._dp.and.lN1.le.360._dp))then
      clN=0
    else
      clN=1
    end if

    return
  end function checklN
  
  ! sets the angle check to define the coordinate to check for the transits
  subroutine lNset(lN,clN)
    real(dp),dimension(:),intent(in)::lN
    integer,dimension(:),intent(out),allocatable::clN
    integer::j

    if(.not.allocated(clN)) allocate(clN(NB))
    clN=0
    do j=2,NB
      clN(j)=checklN(lN(j))
    end do

    return
  end subroutine lNset
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! computes the barycenter of the system
  subroutine barycenter(m,ri,rbar,ro)
    real(dp),dimension(:),intent(in)::m,ri
    real(dp),dimension(:),intent(out)::rbar,ro
    real(dp)::mtot
    integer::j,nj

    mtot=sum(m)
    rbar=zero
    ro=zero
    do j=1,NB
      nj=(j-1)*6
      rbar=rbar+ri(1+nj:6+nj)*m(j)
    end do
    rbar=rbar/mtot
    !compute the position of the star and other bodies respect to the barycenter
    do j=1,NB
      nj=(j-1)*6
      ro(1+nj:6+nj)=ri(1+nj:6+nj)-rbar
      !end if
    end do

    return
  end subroutine barycenter
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! move a body coordinates (x,y,z,vx,vy,vz) of time = dt using the
  ! f and g functions (Murray and Dermont)
  subroutine fgfunctions(mu,rout,dt,Hc)
    real(dp),intent(in)::mu,dt
    real(dp),dimension(:),intent(inout)::rout
    logical,intent(inout)::Hc
    
    real(dp),dimension(3)::rtemp,vtemp
    real(dp)::r0,v0,reca,a,n,sinE,cosE,ecc
    real(dp)::Ek1,MAk1,MAk2,Ek2,dE
    real(dp)::ar0,art
    real(dp)::Ffunc,dFfunc,Gfunc,dGfunc
    real(dp)::h2,x,y,z,vx,vy,vz,hx,hy,hz
    real(dp)::e2

    r0=dist(rout(1:3))
    v0=dist(rout(4:6))
    reca=(two/r0)-((v0**2)/mu)
    if(reca.lt.zero)then
!       write(*,*)" In fgfunctions reca < 0"
      Hc=.false.
      return
    end if
    a=one/reca
    if((a.le.amin).or.(a.ge.amax))then
!       write(*,*)" In fgfunctions a < amin or a > amax"
      Hc=.false.
      return
    end if
    n=sqrt(mu/(a**3))
    sinE=(sum(rout(1:3)*rout(4:6)))/(n*a**2)
    cosE=((r0*v0**2)/mu)-one

    !hx = y*vz-z*vy
    !hy = z*vx-x*vz
    !hz = x*vy-y*vx
    !h2=hx*hx+hy*hy+hz*hz
    !ecc=sqrt( 1 - h2/(a*mu))
    x=rout(1)
    y=rout(2)
    z=rout(3)
    vx=rout(4)
    vy=rout(5)
    vz=rout(6)
    hx = y*vz-z*vy
    hy = z*vx-x*vz
    hz = x*vy-y*vx
    h2=hx*hx+hy*hy+hz*hz
    e2=one-h2/(a*mu)
    if((abs(e2).le.TOLERANCE).or.(e2.lt.zero)) e2=zero
!     write(*,*)"e2 = ",e2
    ecc=sqrt(e2)

    !ecc=sqrt(sinE**2+cosE**2)
    !write(*,'(10000(a,g25.15))')" ecc = ",ecc," a=",a," P = ",dpi/n
    !write(*,'(10000(a,g25.15))')" h2 = ",h2," a*mu = ",a*mu,&
    !     &" (one - h2/(a*mu)) = ",e2
    !read(*,*)
    Ek1=mod(atan2(sinE,cosE),dpi)
    !write(*,'(a,g25.14,a,g25.14)',advance='no')" Ek1 1 r= ",Ek1," d= ",Ek1*rad2deg
    if(Ek1.lt.zero) Ek1=Ek1+dpi

    !write(*,'(a,g25.14,a,g25.14)')" Ek1 2 r= ",Ek1," d= ",Ek1*rad2deg

    MAk1=Ek1-sinE
    !write(*,'(a,g25.14,a,g25.14)')" MAk1 1 r= ",MAk1," d= ",MAk1*rad2deg
    MAk1=mod(MAk1,dpi)
    !write(*,'(a,g25.14,a,g25.14)')" MAk1 2 r= ",MAk1," d= ",MAk1*rad2deg
    if(MAk1.lt.zero) MAk1=MAk1+dpi

    !write(*,'(a,g25.14,a,g25.14)',advance='no')" MAk1 3 r= ",MAk1," d= ",MAk1*rad2deg
    !write(*,'(a,g25.14,a,g25.14,a,g25.14)',advance='no')" n= ",n," dt= ",dt," n*dt= ",n*dt

    MAk2=(MAk1+n*dt)*rad2deg
    !write(*,'(a,g25.14)')" MAk2 1 d= ",MAk2
    MAk2=mod(MAk2,360._dp)
    !write(*,'(a,g25.14)',advance='no')" MAk2 2 d= ",MAk2

    if(MAk2.lt.zero) MAk2=MAk2+360._dp
    !write(*,'(a,g25.14)')" MAk2 3 d= ",MAk2

    Ek2=EAnom(MAk2,ecc)*deg2rad
    !write(*,'(a,g25.14,a,g25.14)',advance='no')" Ek2 r= ",Ek2," d= ",Ek2*rad2deg
    dE=Ek2-Ek1
    !write(*,'(a,g25.14,a,g25.14)')" dE 1 r= ",dE," d= ",dE*rad2deg
    if(dE.ge.pi)then
      dE=dE-dpi
    else if(dE.le.-pi)then
      dE=dE+dpi
    end if
    !write(*,'(a,g25.14,a,g25.14)')" dE 2 r= ",dE," d= ",dE*rad2deg
    ar0=a/r0
    Ffunc=ar0*(cos(dE)-one)+one
    Gfunc=dt+(sin(dE)-dE)/n
    rtemp=Ffunc*rout(1:3)+Gfunc*rout(4:6)
    art=a/dist(rtemp)
    dFfunc=-ar0*art*n*sin(dE)
    dGfunc=art*(cos(dE)-one)+one
    vtemp=dFfunc*rout(1:3)+dGfunc*rout(4:6)

    rout(1:3)=rtemp
    rout(4:6)=vtemp

    return
  end subroutine fgfunctions
  
    ! computing the singular specific angular momentum
  subroutine angmom(rvin,Lvec)
    real(dp),intent(in)::rvin(6)
    real(dp),intent(out)::Lvec(3)
    real(dp)::x,y,z,vx,vy,vz

    x=rvin(1)
    y=rvin(2)
    z=rvin(3)
    vx=rvin(4)
    vy=rvin(5)
    vz=rvin(6)

    Lvec(1)=y*vz-z*vy
    Lvec(2)=z*vx-x*vz
    Lvec(3)=x*vy-y*vx

    return
  end subroutine angmom

  ! computes the constant of the motion: Energy and Angular Momentum
  subroutine const_motion(m,rin,Etot,htot)
    real(dp),dimension(:),intent(in)::m
    real(dp),dimension(:),intent(in)::rin
    real(dp),intent(out)::Etot,htot
    real(dp),dimension(3)::hvec,htvec,rj1,rj2,r0
    real(dp),dimension(6)::rvj1
    real(dp)::Ekin,Epot,temp,rj10,vj1,rj21
    integer j1,j2,nj1,nj2

    htot=zero
    Etot=zero
    Ekin=zero
    Epot=zero
    hvec=zero
    htvec=zero
    r0=rin(1:3)

    ! Total Angular Momentum & kinetic Energy partial
    j1do: do j1=1,NB
      nj1=(j1-1)*6
      temp=zero
      rvj1=rin(1+nj1:6+nj1)
      rj1=rvj1(1:3)
      vj1=dist(rvj1(4:6))
      call angmom(rvj1,hvec)
      htvec=htvec+m(j1)*hvec
      Ekin=Ekin+m(j1)*(vj1*vj1)
      if(j1.gt.1)then
        rj10=dist(rj1,r0)
        j2do: do j2=j1+1,NB
          nj2=(j2-1)*6
          rj2=rin(1+nj2:3+nj2)
          rj21=dist(rj2,rj1)
          temp=temp+m(j2)/rj21
        end do j2do
        Epot=Epot-Giau*m(j1)*(temp+m(1)/rj10)
      end if
    end do j1do
    htot=dist(htvec)
    Etot=half*Ekin+Epot

    return
  end subroutine const_motion

  ! small subroutine to call sequentially the subroutines
  ! needed to update the E and h constant of motion
  subroutine compute_con(m,rin,Etot,htot)
    real(dp),dimension(:),intent(in)::m,rin
    real(dp),intent(inout)::Etot,htot
    real(dp),dimension(6)::bar
    real(dp),dimension(:),allocatable::rbar

    allocate(rbar(NBDIM))
    call barycenter(m,rin,bar,rbar)
    call const_motion(m,rbar,Etot,htot)
    deallocate(rbar)

    return
  end subroutine compute_con
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! from state vector (cartesian coordinates) to keplerian orbital elements
  ! for one body
  subroutine eleMD(mu,svec,P,a,e,inc,mA,w,lN,f,dt)
    real(dp),intent(in)::mu
    real(dp),intent(in),dimension(:)::svec
    real(dp),intent(out)::P,a,e,inc,mA,w,lN,f,dt
    real(dp)::x,y,z,vx,vy,vz,R,R2,V,V2,rv,rd
    real(dp)::hx,hy,hz,h2,h
    real(dp)::Energy
    real(dp)::cosi,sini
    real(dp)::coslN,sinlN,cpslN,cmslN
    real(dp)::wf,sinwf,coswf
    real(dp)::Ea,mmotion
    real(dp)::sinf,cosf
    real(dp)::arg1,arg2

    x=svec(1)
    y=svec(2)
    z=svec(3)
    R2=x*x+y*y+z*z
    R=sqrt(R2)

    vx=svec(4)
    vy=svec(5)
    vz=svec(6)
    V2=vx*vx+vy*vy+vz*vz
    V=sqrt(V2)

    rv=x*vx+y*vy+z*vz

    hx = y*vz-z*vy
    hy = z*vx-x*vz
    hz = x*vy-y*vx
    h2=hx*hx+hy*hy+hz*hz
    h=sqrt(h2)

    arg1 = h2/R2
    if(V2.lt.arg1)then
      rd = zero
    else
      rd=sqrt(V2-(h2/R2))
    end if
    rd=sign(rd,rv)

    Energy=(half*V2)-(mu/R)
    a=-half*mu/Energy

    P=dpi*sqrt((a**3)/mu)

    mmotion=sqrt(mu/(a**3))

    lN=mod((atan2(hx,-hy)+dpi),dpi)
!     if(hx.lt.zero) lN=lN+dpi
    coslN=cos(lN)
    sinlN=sin(lN)
    
    cosi=hz/h
!     inc=acos(cosi)
!     sini=sin(inc)
        
    sini=sqrt( ((hx*hx)+(hy*hy)) / h2 )
!     sini=hx/(h*sinlN)
    
    inc=atan2(sini,cosi)
!     write(*,*)'cosi = ',cosi,' sini = ',sini,' inc = ',inc
!     flush(6)

!     sinwf=z/(R*sini)
    !if(sini.eq.zero)then
    if(abs(sini-zero).le.TOLERANCE)then
      e=zero
      w=zero
      f=zero
      lN=zero
      mA=zero
      dt=zero
      return
    else
      sinwf=z/(R*sini)
    end if

    arg1=h2/(a*mu)
    arg2=one-arg1
    if(arg2.le.TOLERANCE)then
      ! circular orbit
      e=zero
      w=pi/two
      cpslN=coslN+sinlN
      cmslN=coslN-sinlN
      cosf=sinwf
      sinf=(one/cmslN)*( ((y-x)/R)-(cosi*cosf*cpslN))
      f=atan2(sinf,cosf)
      if(sinf.lt.zero) f=f+dpi
      mA=f
    else
      ! elliptical orbit
      e=sqrt(arg2)
      !if((abs(coslN).le.TOLERANCE).or.(coslN.eq.zero))then
      if((abs(coslN).le.TOLERANCE)&
          &.or.(abs(coslN-zero).le.TOLERANCE))then
        coswf=(one/sinlN)*((y/R)-coslN*cosi*sinwf)
      else
        coswf=(one/coslN)*((x/R)+sinlN*cosi*sinwf)
      end if
      wf=atan2(sinwf,coswf)
      if(isnan(coswf))then
        write(*,'(a,es23.16)')" coswf ",coswf
        write(*,'(a,es23.16)')" sinwf = ",sinwf
        write(*,'(a,es23.16)')" sinlN = ",sinlN
        write(*,'(a,es23.16)')" coslN = ",coslN
        write(*,'(a,es23.16)')" y = ",y
        write(*,'(a,es23.16)')" x = ",x
        write(*,'(a,es23.16)')" R = ",R
        write(*,'(a,es23.16)')" cosi = ",cosi
        write(*,'(a,es23.16)')" sini = ",sini
!         stop
        return
      end if
      if(sinwf.lt.zero) wf=wf+dpi

      arg1=a*(one-e*e)
      !sinf=rd*arg1/(h*e)
      sinf=rd*arg1/h
      if(isnan(sinf))then
        write(*,*)" WARNING: sinf is NaN"
        write(*,'( "mu = ",es23.16)')mu
        write(*,*)" a = ",a," e = ",e," h = ",h," h2 = ",h2
        write(*,*)" rd = ",rd," rv = ",rv," V2 = ",V2," R2 = ",R2
        write(*,*)" h2/R2 = ",(h2/R2)," V2 -(h2/R2) = ",V2-(h2/R2)
        write(*,*)" rd[before sign(rd,rv)] = ",sqrt(V2-(h2/R2))
!         stop
        return
      end if
      !cosf=(one/e)*((arg1/R)-one)
      cosf=(arg1/R)-one
      if(isnan(cosf))then
        write(*,*)" WARNING: cosf is NaN"
        write(*,*)" a = ",a," e = ",e
!         stop
        return
      end if

      f=atan2(sinf,cosf)
      if(sinf.lt.zero) f=f+dpi

      w=wf-f
      if(w.lt.zero) w=w+dpi

      arg1=sqrt((one-e)/(one+e))
      arg2=tan(half*f)
      Ea=two*atan(arg1*arg2)
      if(Ea.lt.zero) Ea=Ea+dpi

      mA=Ea-e*sin(Ea)
    end if

    dt=mA/mmotion
    if(mA.lt.zero) mA=mA+dpi

    return
  end subroutine eleMD

  ! from state vector (cartesian coordinates) to keplerian orbital elements
  ! for all the bodies
  subroutine elements(m,rin,P,a,e,inc,mA,w,f,lN,dttau)
    real(dp),dimension(:),intent(in)::m,rin
    real(dp),dimension(:),intent(out)::P,a,e,inc,mA,w,f,lN,dttau
    real(dp),dimension(6)::svec
    real(dp)::mu
    integer::j,ncj

    P=zero
    a=zero
    e=zero
    inc=zero
    mA=zero
    w=zero
    f=zero
    lN=zero
    dttau=zero

    cicle: do j=2,NB
      ncj=(j-1)*6
      mu=Giau*(m(1)+m(j))
      ! state vector j-th body
      svec=rin(1+ncj:6+ncj)
      call eleMD(mu,svec,P(j),a(j),e(j),inc(j),mA(j),w(j),lN(j),f(j),dttau(j))
      inc(j)=inc(j)*rad2deg
      mA(j)=mA(j)*rad2deg
      lN(j)=lN(j)*rad2deg
      w(j)=w(j)*rad2deg
      f(j)=f(j)*rad2deg
    end do cicle

    return
  end subroutine elements
  ! ------------------------------------------------------------------ !
  
end module celestial_mechanics
