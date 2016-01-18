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
    real(dp)::dpi2=dpi*dpi,athird=1._dp/3._dp

    mu=Giau*(ms+mp)
    P2=P*P
    sma=((mu*P2)/dpi2)**(athird)

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
!     real(dp),parameter::csqrt=1._dp/3._dp

    rH=sma*(mp/(3._dp*ms))**sqrt_1_3

    return
  end function rHill_1

!   function to compute the Hill radius w/ eccentricity
  function rHill_2(ms,mp,sma,ecc) result(rH)
    real(dp)::rH
    real(dp),intent(in)::ms,mp,sma,ecc
!     real(dp),parameter::csqrt=1._dp/3._dp
    real(dp)::e1

    e1=1._dp-ecc
    rH=sma*e1*(mp/(3._dp*ms))**sqrt_1_3

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

    rH=.false.
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
      a1=(2._dp/r1)-((v1**2)/mu1) ! 1/a
      if(a1.le.zero)then
        fmt=adjustl('(a,1x,'//trim(sprec)//')')
        write(*,trim(fmt),advance='no')" 1/a1 < 0 :",a1
        rH=.true.
        return
      end if
      a1=1._dp/a1
      n1=sqrt(mu1/a1**3)
      es=rv1/(n1*a1*a1)
      ec=((r1*v1*v1)/mu1)-1._dp
      e1=sqrt(ec*ec+es*es)

      rvec2=rin(1+j2:3+j2)
      vvec2=rin(4+j2:6+j2)
      !rv2=rvec2(1)*vvec2(1)+rvec2(2)*vvec2(2)+rvec2(3)*vvec2(3)
      rv2=rvprod(rvec2,vvec2)
      r2=dist(rvec2)
      v2=dist(vvec2)
      mu2=Giau*(m(1)+m(j+1))
      a2=(2._dp/r2)-((v2**2)/mu2) ! 1/a
      if(a2.le.zero)then
        fmt=adjustl('(a,1x,'//trim(sprec)//')')
        write(*,trim(fmt),advance='no')" 1/a2 < 0 :",a2
        rH=.true.
        return
      end if
      a2=1._dp/a2
      n1=sqrt(mu1/a2**3)
      es=rv2/(n1*a2*a2)
      ec=((r2*v2*v2)/mu2)-1._dp
      e1=sqrt(ec*ec+es*es)
      if((a1.le.amin).or.(a1.ge.amax).or.(a2.le.amin).or.(a2.ge.amax))then
        !write(*,'(a)')" Bound semi-major axis limits overcome "
        rH=.true.
        return
      end if
      rH1=rHill(m(1),m(j),a1,e1)
      rH2=rHill(m(1),m(j+1),a2,e1)
      dr=dist(rvec1,rvec2)
      if(dr.le.(rH1+rH2))then
        rH=.true.
        return
      end if
    end do

    return
  end function Hillcheck
  
! -----------------------

! -----------------------
! mutual Hill radius

  function mutual_Hill_radius(ms,mp1,sma1,mp2, sma2) result(rH)
    real(dp)::rH
    real(dp),intent(in)::ms,mp1,sma1,mp2,sma2
    real(dp)::sma_mean,mass_ratio
    
    sma_mean = 0.5_dp*(sma1+sma2)
    mass_ratio = (mp1+mp2)/(3._dp*ms)
    rH=sma_mean*mass_ratio**sqrt_1_3

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
    inv_sma=(2._dp/rr)-((vv**2)/mu) ! 1/a
    check=.false.
    if(inv_sma.le.zero) check=.true. ! if true is bad

    return
  end subroutine rrdot_to_invsma
  
    
  function mutual_Hill_check_0(m,rin) result(hill_check)
    logical::hill_check
    real(dp),dimension(:),intent(in)::m,rin
    real(dp)::sma_i, sma_j
    real(dp)::Hill_radius_ij,delta_ij,stability_criterion
    integer::i,ii,j
    
    do i=2,(NB-1)
      call rrdot_to_invsma(i,m,rin,sma_i,hill_check)
      if(hill_check) return
      sma_i=one/sma_i
      
      do j=i+1,NB
        call rrdot_to_invsma(j,m,rin,sma_j,hill_check)
        if(hill_check) return
        sma_j=one/sma_j
        
        Hill_radius_ij = mutual_Hill_radius(m(1),m(i),sma_i,m(j),sma_j)
        delta_ij = abs(sma_j - sma_i)
        stability_criterion = sqrt_12 * Hill_radius_ij
        if(delta_ij.lt.stability_criterion)then
          hill_check = .true. ! if true is bad/unstable
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
    integer::i,ii,j
    
    hill_check=.false.
    if(NB.gt.2)then
      do i=2,(NB-1)
        call rrdot_to_invsma(i,m,rin,sma_i,hill_check)
        if(hill_check) return
        sma_i=one/sma_i

        j=i+1
        call rrdot_to_invsma(j,m,rin,sma_j,hill_check)
        if(hill_check) return
        sma_j=one/sma_j

        Hill_radius_ij = mutual_Hill_radius(m(1),m(i),sma_i,m(j),sma_j)
        delta_ij = abs(sma_j - sma_i)
        stability_criterion = sqrt_12 * Hill_radius_ij
        if(delta_ij.lt.stability_criterion)then
          hill_check = .true. ! if true is bad/unstable
          return
        end if
      end do
    end if
        
    return
  end function mutual_Hill_check
  
! -------------------------------
  
!   check the physical bounds of the Keplerian elements, e.g., Mass cannot be lower than 0 ...etc
  function checkbounds(theta) result(check)
    real(dp),dimension(:),intent(in)::theta
    logical::check
    integer::j,body

    !all thetas must be positive
    check = .true.

    ! Mass (id=3), Radius (id=4), Period (id=5), eccentricity (id=6), arg.peric. (id=7), 
    ! meanAnom. (id=8), inclination (id=9), long.Node (id=10)
    do j=1,nfit
      if(theta(j).lt.zero)then ! check if parameter is negative
        if(id(j).eq.3)then ! mass
          check=.false.
          return
        else if(id(j).eq.4)then ! radius
          check=.false.
          return
        else if(id(j).eq.5)then ! period
          check=.false.
          return
        else if(id(j).eq.6)then ! eccentricity
          check=.false.
          return
        else if(id(j).eq.9)then ! inclination
          check=.false.
          return
        end if
      else ! parameter is zero or positive
        ! eccentricity cannot be equal or greater than one
        if(id(j).eq.6)then
          body=int((idall(j)-3)/8)+2
          if((theta(j).lt.e_bounds(1,body)).or.(theta(j).gt.e_bounds(2,body)))then
            check=.false.
            return
          end if
        end if
        ! inclination cannot be greater than 180degrees
        if((id(j).eq.9).and.(theta(j).ge.180._dp))then
          check=.false.
          return
        end if
      end if
    end do

    return
  end function checkbounds

!   check the physical bounds of the Keplerian elements, e.g., Mass cannot be lower than 0 ...etc
  function checkbounds_fit(theta) result(check)
    real(dp),dimension(:),intent(in)::theta
    logical::check
    integer::j,body
    real(dp)::temp_ecc

    !all thetas must be positive
    check = .true.

    ! Mass (id=3), Radius (id=4), Period (id=5), ecosw (id=6), esinw (id=7), 
    ! meanAnom. (id=8), inclination (id=9), long.Node (id=10)
    do j=1,nfit
        if(id(j).eq.3)then ! mass
          if(theta(j).lt.zero)then ! check if parameter is negative
            check=.false.
            return
          end if
        else if(id(j).eq.4)then ! radius
          if(theta(j).lt.zero)then ! check if parameter is negative
            check=.false.
            return
          end if
        else if(id(j).eq.5)then ! period
          if(theta(j).lt.zero)then ! check if parameter is negative
            check=.false.
            return
          end if
        else if(id(j).eq.6)then ! ecosw, esinw
          body=int((idall(j)-3)/8)+2
          temp_ecc = sqrt(theta(j)*theta(j) + theta(j+1)*theta(j+1)) ! ecc
          if((temp_ecc.lt.e_bounds(1,body)).or.(temp_ecc.gt.e_bounds(2,body)))then
            check=.false.
            return
          end if
        else if(id(j).eq.9)then ! inclination
          if((theta(j).le.TOLERANCE).or.(theta(j).ge.180._dp))then
            check=.false.
            return
          end if
        end if
      
      end do

    return
  end function checkbounds_fit

  
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
      rad1e2=sqrt(1._dp-e2)
      cosE=cos(EA)
      sinE=sin(EA)
      dEA=n/(1._dp-e(i)*cosE)
      !calculation of the radius vector and velocity of the bodies
      nci=(i-1)*6
      output(1+nci)=a(i)*(cosE-e(i))  !! x
      output(2+nci)=a(i)*rad1e2*sinE  !! y
      output(4+nci)=-a(i)*sinE*dEA  !! vx
      output(5+nci)=a(i)*rad1e2*cosE*dEA  !! vy
    end do

    return
  end subroutine initial_state

  ! from allpar and par to keplerian elements
  subroutine par2kel(allpar,par,m,R,P,a,e,w,mA,i,lN)
    real(dp),dimension(:),intent(in)::allpar,par
    real(dp),dimension(:),intent(out)::m,R,P,a,e,w,mA,i,lN
    real(dp),dimension(:),allocatable::atemp
    integer::j1,cnt

    m=zero
    R=zero
    P=zero
    a=zero
    e=zero
    w=zero
    mA=zero
    i=zero
    lN=zero
    allocate(atemp(npar))
    atemp=allpar
    cnt=0
    do j1=1,npar
      if(tofit(j1).eq.1)then
        cnt=cnt+1
        atemp(j1)=par(cnt)
      end if
    end do
    cnt=0
    m(1)=atemp(1)
    R(1)=atemp(2)
    do j1=2,NB
      cnt=(j1-2)*8
      m(j1)=atemp(3+cnt)
      R(j1)=atemp(4+cnt)
      P(j1)=atemp(5+cnt)
      a(j1)=semax(m(1),m(j1),P(j1))
      e(j1)=atemp(6+cnt)
      w(j1)=atemp(7+cnt)
      mA(j1)=atemp(8+cnt)
      i(j1)=atemp(9+cnt)
      lN(j1)=atemp(10+cnt)
      !w(j1)=mod(atemp(7+cnt),360._dp)
      !mA(j1)=mod(atemp(8+cnt),360._dp)
      !i(j1)=mod(atemp(9+cnt),180._dp)
      !lN(j1)=mod(atemp(10+cnt),360._dp)
    end do
    deallocate(atemp)

    return
  end subroutine par2kel

  ! from allpar and par to keplerian elements
  subroutine par2kel_fit(allpar,par,m,R,P,a,e,w,mA,i,lN,checkpar)
    real(dp),dimension(:),intent(in)::allpar,par
    real(dp),dimension(:),intent(out)::m,R,P,a,e,w,mA,i,lN
    real(dp),dimension(:),allocatable::atemp
    integer::j1,cnt
    logical,intent(out)::checkpar
    real(dp)::temp_ecc2

    checkpar=.true.
    m=zero
    R=zero
    P=zero
    a=zero
    e=zero
    w=zero
    mA=zero
    i=zero
    lN=zero
    allocate(atemp(npar))
    atemp=allpar
    cnt=0
    do j1=1,npar
      if(tofit(j1).eq.1)then
        cnt=cnt+1
        atemp(j1)=par(cnt)
      end if
    end do
    cnt=0
    m(1)=atemp(1)
    R(1)=atemp(2)
    do j1=2,NB
      cnt=(j1-2)*8
      m(j1)=atemp(3+cnt)
      R(j1)=atemp(4+cnt)
      P(j1)=atemp(5+cnt)
      a(j1)=semax(m(1),m(j1),P(j1))
      if(tofit(6+cnt).eq.1)then
        temp_ecc2 = atemp(6+cnt)*atemp(6+cnt) + atemp(7+cnt)*atemp(7+cnt)
        if(temp_ecc2.lt.TOLERANCE)then
          checkpar=.false.
          return
        else
          e(j1)=sqrt(temp_ecc2)
          if((e(j1).lt.e_bounds(1,j1)).or.(e(j1).gt.e_bounds(2,j1)))then
            checkpar=.false.
            return
          end if
          w(j1)=atan2(atemp(7+cnt), atemp(6+cnt))*rad2deg
          if(w(j1).lt.zero) w(j1)=w(j1)+360._dp
        end if
      else
        e(j1)=atemp(6+cnt)
        w(j1)=atemp(7+cnt)
      end if
      mA(j1)=atemp(8+cnt)
      i(j1)=atemp(9+cnt)
      lN(j1)=atemp(10+cnt)
      !w(j1)=mod(atemp(7+cnt),360._dp)
      !mA(j1)=mod(atemp(8+cnt),360._dp)
      !i(j1)=mod(atemp(9+cnt),180._dp)
      !lN(j1)=mod(atemp(10+cnt),360._dp)
    end do
    deallocate(atemp)

    return
  end subroutine par2kel_fit
  ! ------------------------------------------------------------------ !
  
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
    rbar=0._dp
    ro=0._dp
    do j=1,NB
      nj=(j-1)*6
      rbar=rbar+ri(1+nj:6+nj)*m(j)
    end do
    rbar=rbar/mtot
    !compute the position of the star and other bodies respect to the barycenter
    do j=1,NB
      nj=(j-1)*6
      !if(m(j).ne.0._dp)then
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
    reca=(2._dp/r0)-((v0**2)/mu)
    if(reca.lt.0._dp)then
      write(*,*)" In fgfunctions reca < 0"
      Hc=.true.
      return
    end if
    a=1._dp/reca
    if((a.le.amin).or.(a.ge.amax))then
      write(*,*)" In fgfunctions a < amin or a > amax"
      Hc=.true.
      return
    end if
    n=sqrt(mu/(a**3))
    sinE=(sum(rout(1:3)*rout(4:6)))/(n*a**2)
    cosE=((r0*v0**2)/mu)-1._dp

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
    e2=1._dp-h2/(a*mu)
    if((abs(e2).le.TOLERANCE).or.(e2.lt.zero)) e2=zero
!     write(*,*)"e2 = ",e2
    ecc=sqrt(e2)

    !ecc=sqrt(sinE**2+cosE**2)
    !write(*,'(10000(a,g25.15))')" ecc = ",ecc," a=",a," P = ",dpi/n
    !write(*,'(10000(a,g25.15))')" h2 = ",h2," a*mu = ",a*mu,&
    !     &" (1._dp - h2/(a*mu)) = ",e2
    !read(*,*)
    Ek1=mod(atan2(sinE,cosE),dpi)
    !write(*,'(a,g25.14,a,g25.14)',advance='no')" Ek1 1 r= ",Ek1," d= ",Ek1*rad2deg
    if(Ek1.lt.0._dp) Ek1=Ek1+dpi

    !write(*,'(a,g25.14,a,g25.14)')" Ek1 2 r= ",Ek1," d= ",Ek1*rad2deg

    MAk1=Ek1-sinE
    !write(*,'(a,g25.14,a,g25.14)')" MAk1 1 r= ",MAk1," d= ",MAk1*rad2deg
    MAk1=mod(MAk1,dpi)
    !write(*,'(a,g25.14,a,g25.14)')" MAk1 2 r= ",MAk1," d= ",MAk1*rad2deg
    if(MAk1.lt.0._dp) MAk1=MAk1+dpi

    !write(*,'(a,g25.14,a,g25.14)',advance='no')" MAk1 3 r= ",MAk1," d= ",MAk1*rad2deg
    !write(*,'(a,g25.14,a,g25.14,a,g25.14)',advance='no')" n= ",n," dt= ",dt," n*dt= ",n*dt

    MAk2=(MAk1+n*dt)*rad2deg
    !write(*,'(a,g25.14)')" MAk2 1 d= ",MAk2
    MAk2=mod(MAk2,360._dp)
    !write(*,'(a,g25.14)',advance='no')" MAk2 2 d= ",MAk2

    if(MAk2.lt.0._dp) MAk2=MAk2+360._dp
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
    Ffunc=ar0*(cos(dE)-1._dp)+1._dp
    Gfunc=dt+(sin(dE)-dE)/n
    rtemp=Ffunc*rout(1:3)+Gfunc*rout(4:6)
    art=a/dist(rtemp)
    dFfunc=-ar0*art*n*sin(dE)
    dGfunc=art*(cos(dE)-1._dp)+1._dp
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

    x   =rvin(1)
    y   =rvin(2)
    z   =rvin(3)
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
    Etot=0.5_dp*Ekin+Epot

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
  subroutine eleMD(mu,svec,P,a,e,i,mA,w,lN,f,dt)
    real(dp),intent(in)::mu
    real(dp),intent(in),dimension(:)::svec
    real(dp),intent(out)::P,a,e,i,mA,w,lN,f,dt
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

    Energy=(0.5_dp*V2)-(mu/R)
    a=-0.5_dp*mu/Energy

    P=dpi*sqrt((a**3)/mu)

    mmotion=sqrt(mu/(a**3))

    cosi=hz/h
    i=acos(cosi)
    sini=sin(i)

    lN=atan2(hx,-hy)
    if(hx.lt.zero) lN=lN+dpi
    coslN=cos(lN)
    sinlN=sin(lN)

    sinwf=z/(R*sini)
    !if(sini.eq.zero)then
    if(abs(sini-zero).le.TOLERANCE)then
      e=zero
      w=zero
      f=zero
      lN=zero
      mA=zero
      dt=zero
      return
    end if

    arg1=h2/(a*mu)
    arg2=1._dp-arg1
    if(arg2.le.TOLERANCE)then
      ! circular orbit
      e=zero
      w=pi/2._dp
      cpslN=coslN+sinlN
      cmslN=coslN-sinlN
      cosf=sinwf
      sinf=(1._dp/cmslN)*( ((y-x)/R)-(cosi*cosf*cpslN))
      f=atan2(sinf,cosf)
      if(sinf.lt.zero) f=f+dpi
      mA=f
    else
      ! elliptical orbit
      e=sqrt(arg2)
      !if((abs(coslN).le.TOLERANCE).or.(coslN.eq.zero))then
      if((abs(coslN).le.TOLERANCE)&
          &.or.(abs(coslN-zero).le.TOLERANCE))then
        coswf=(1._dp/sinlN)*((y/R)-coslN*cosi*sinwf)
      else
        coswf=(1._dp/coslN)*((x/R)+sinlN*cosi*sinwf)
      end if
      wf=atan2(sinwf,coswf)
      if(isnan(coswf))then
        write(*,'(a,g25.15)')" coswf ",coswf
        write(*,'(a,g25.15)')" sinwf = ",sinwf
        write(*,'(a,g25.15)')" sinlN = ",sinlN
        write(*,'(a,g25.15)')" coslN = ",coslN
        write(*,'(a,g25.15)')" y = ",y
        write(*,'(a,g25.15)')" x = ",x
        write(*,'(a,g25.15)')" R = ",R
        write(*,'(a,g25.15)')" cosi = ",cosi
        write(*,'(a,g25.15)')" sini = ",sini
!         stop
        return
      end if
      if(sinwf.lt.zero) wf=wf+dpi

      arg1=a*(1._dp-e*e)
      !sinf=rd*arg1/(h*e)
      sinf=rd*arg1/h
      if(isnan(sinf))then
        write(*,*)" WARNING: sinf is NaN"
        write(*,'( "mu = ",G25.14)')mu
        write(*,*)" a = ",a," e = ",e," h = ",h," h2 = ",h2
        write(*,*)" rd = ",rd," rv = ",rv," V2 = ",V2," R2 = ",R2
        write(*,*)" h2/R2 = ",(h2/R2)," V2 -(h2/R2) = ",V2-(h2/R2)
        write(*,*)" rd[before sign(rd,rv)] = ",sqrt(V2-(h2/R2))
!         stop
        return
      end if
      !cosf=(1._dp/e)*((arg1/R)-1._dp)
      cosf=(arg1/R)-1._dp
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

      arg1=sqrt((1._dp-e)/(1._dp+e))
      arg2=tan(0.5_dp*f)
      Ea=2._dp*atan(arg1*arg2)
      if(Ea.lt.zero) Ea=Ea+dpi

      mA=Ea-e*sin(Ea)
    end if

    dt=mA/mmotion
    if(mA.lt.zero) mA=mA+dpi

    return
  end subroutine eleMD

  ! from state vector (cartesian coordinates) to keplerian orbital elements
  ! for all the bodies
  subroutine elements(m,rin,P,a,e,i,mA,w,f,lN,dttau)
    real(dp),dimension(:),intent(in)::m,rin
    real(dp),dimension(:),intent(out)::P,a,e,i,mA,w,f,lN,dttau
    real(dp),dimension(6)::svec
    real(dp)::mu
    integer::j,ncj

    P=zero
    a=zero
    e=zero
    i=zero
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
      call eleMD(mu,svec,P(j),a(j),e(j),i(j),mA(j),w(j),lN(j),f(j),dttau(j))
      i(j)=i(j)*rad2deg
      mA(j)=mA(j)*rad2deg
      lN(j)=lN(j)*rad2deg
      w(j)=w(j)*rad2deg
      f(j)=f(j)*rad2deg
    end do cicle

    return
  end subroutine elements
  ! ------------------------------------------------------------------ !
  
end module celestial_mechanics
