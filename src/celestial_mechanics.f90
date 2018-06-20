module celestial_mechanics
  use constants
  use parameters
  use rotations
  implicit none

  interface dist
    module procedure dist_1,dist_2
  end interface dist

  contains

! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   function to compute the period given the semi-major axis
  subroutine period_vec(ms,mp,sma,a2P)
    real(dp),intent(IN)::ms
    real(dp),dimension(:),intent(in)::mp,sma
    real(dp),dimension(:),intent(out)::a2P
    real(dp),dimension(:),allocatable::mu

    allocate(mu(size(mp)))
    mu=Giau*(ms+mp)
    a2P=dpi*sqrt(sma**3/mu)
    deallocate(mu)

    return
  end subroutine period_vec
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  ! time pericentre to mean anomaly
  function tau2mA(tau, t_ref, per) result(mA)
    real(dp)::mA
    real(dp),intent(in)::tau,t_ref,per

    mA = mod( ((circ/per)*(t_ref-tau)) , circ)

    return
  end function tau2mA
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  ! mean anomaly to time pericentre
  function mA2tau(mA, t_ref, per) result(tau)
    real(dp)::tau
    real(dp),intent(in)::mA,t_ref,per

    tau = t_ref - (mA*per/circ)

    return
  end function mA2tau
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  ! vector version
  ! time pericentre to mean anomaly
  subroutine tau2mA_vec(tau, t_ref, per, mA)
    real(dp),dimension(:),intent(in)::tau,per
    real(dp),intent(in)::t_ref
    real(dp),dimension(:),intent(out)::mA

    mA = mod( ((circ/per)*(t_ref-tau)) , circ)

    return
  end subroutine tau2mA_vec
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  ! mean anomaly to time pericentre
  subroutine mA2tau_vec(mA, t_ref, per, tau)
    real(dp),dimension(:),intent(in)::mA,per
    real(dp),intent(in)::t_ref
    real(dp),dimension(:),intent(out)::tau

    tau = t_ref - (mA*per/circ)

    return
  end subroutine mA2tau_vec
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   calculates Eccentric anomaly from meanAnomaly [deg] and eccentricity
  function EAnom(mA,ecc) result(EA)
    real(dp)::EA
    real(dp),intent(IN)::mA,ecc
    real(dp)::mArad,E,fE,dfE
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
        if(abs(E-EA).le.TOLERANCE) exit loopo
        E=EA
      end do loopo
      EA=EA*rad2deg
    else
      EA=mA
    end if

    return
  end function EAnom
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  function calculate_true_anomaly(mean_anomaly, eccentricity) result(true_anomaly)
    real(dp)::true_anomaly

    real(dp),intent(in)::mean_anomaly,eccentricity

    real(dp)::EA,tan_EA,ecc_coeff

    if(eccentricity.le.TOLERANCE)then
      true_anomaly=mean_anomaly*deg2rad
    else
      EA=EAnom(mod(mean_anomaly+circ,circ),eccentricity)*deg2rad
      tan_EA=tan(EA*half)
      ecc_coeff=sqrt((one+eccentricity)/(one-eccentricity))
      true_anomaly = two * atan(ecc_coeff*tan_EA) ! output in rad!!!
    end if

    return
  end function calculate_true_anomaly
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  function trueAnom_ecc_to_eccAnom(trueAnom, ecc) result(eccAnom)
    real(dp)::eccAnom

    real(dp),intent(in)::trueAnom,ecc

    real(dp)::tan_htA,ecoeff

    if(ecc.le.TOLERANCE)then
      eccAnom=trueAnom
    else
      tan_htA=tan(half*trueAnom*deg2rad)

      ecoeff=sqrt((one-ecc)/(one+ecc))
      eccAnom = mod((two * atan(ecoeff*tan_htA))+dpi,dpi) ! output in rad!!!
    end if

    return
  end function trueAnom_ecc_to_eccAnom
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   calculate the module of a 3-D vector
  function dist_1(r) result(out)
    real(dp)::out
    real(dp),dimension(3),intent(in)::r
    real(dp)::x,y,z

    x=r(1)
    y=r(2)
    z=r(3)

    out=sqrt(x*x + y*y + z*z)

    return
  end function dist_1
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   function to compute the distance between to vectors
  function dist_2(r1,r2) result(out)
    real(dp)::out
    real(dp),dimension(3),intent(in)::r1,r2
    real(dp)::dx,dy,dz

    dx=r1(1)-r2(1)
    dy=r1(2)-r2(2)
    dz=r1(3)-r2(3)
    out=sqrt(dx*dx + dy*dy + dz*dz)

    return
  end function dist_2
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   module of a 2-D vector
  function rsky(r) result(out)
    real(dp)::out
    real(dp),dimension(2),intent(in)::r
    real(dp)::x,y

    x=r(1)
    y=r(2)
    out=sqrt(x*x + y*y)

    return
  end function rsky
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   function to compute the Hill radius w/o eccentricity
  function rHill_circ(ms,mp,sma) result(rH)
    real(dp)::rH
    real(dp),intent(in)::ms,mp,sma

    rH=sma*((onethird*(mp/ms))**onethird)

    return
  end function rHill_circ
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   function to compute the Hill radius w/ eccentricity
  function rHill_ecc(ms,mp,sma,ecc) result(rH)
    real(dp)::rH
    real(dp),intent(in)::ms,mp,sma,ecc
    real(dp)::e1

    e1=one-ecc
    rH=sma*e1*((onethird*(mp/ms))**onethird)

    return
  end function rHill_ecc
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! mutual Hill radius
  function mutual_Hill_radius(ms,mp1,sma1,mp2,sma2) result(rH)
    real(dp)::rH
    real(dp),intent(in)::ms,mp1,sma1,mp2,sma2
    real(dp)::sma_mean,mass_ratio

    real(dp)::min_ratio

    rH=zero
    sma_mean=half*(sma1+sma2)
    mass_ratio=(mp1+mp2)/ms
    min_ratio = TOL_dp**onethird


    if(mass_ratio.le.TOL_dp)then

      rH=sma_mean*min_ratio
      write(*,'(2(a,es23.16),a)')' mass_ratio = ',mass_ratio,' <= ',TOL_dp,' = TOL_dp ==> USING TOL_dp'
      write(*,'(3(a,es23.16))')' Mstar [Msun] = ',ms,' Mi [Msun]= ',mp1,' Mj [Msun]= ',mp2
      flush(6)
!       rH=1.e10_dp

    else if(ms.le.TOL_dp)then

      write(*,'(2(a,es23.16),a)')' Mstar = ',ms,' <= ',TOL_dp,' = TOL_dp'
      write(*,'(3(a,es23.16))')' Mstar [Msun] = ',ms,' Mi [Msun]= ',mp1,' Mj [Msun]= ',mp2
      flush(6)
      rH=1.e10_dp
!       stop

    else if(sma_mean.le.TOL_dp)then

      write(*,'(2(a,es23.16),a)')' sma_mean = ',sma_mean,' <= ',TOL_dp,' = TOL_dp'
      write(*,'(3(a,es23.16))')' Mstar [Msun] = ',ms,' Mi [Msun]= ',mp1,' Mj [Msun]= ',mp2
      flush(6)
      rH=1.e10_dp
!       stop

!    ! else if(ms.le.(mp1+mp2))then
!     else if(abs(ms-(mp1+mp2)).le.TOL_dp)then
!
!       write(*,'(2(a,es23.16),a)')' Mstar [Msun] = ',ms,' <= ',mp1+mp2,' = Mi + Mj  [Msun]'
!       flush(6)
!       rH=1.e10_dp
! !       stop

    else

      rH=sma_mean*((onethird*mass_ratio)**onethird)

    end if

    return
  end function mutual_Hill_radius
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  function separation_mutual_Hill_check(m,R,rin,do_Hr_check) result(hill_check)
    logical::hill_check
    real(dp),dimension(:),intent(in)::m,R,rin
    logical,intent(in)::do_Hr_check
    real(dp)::rij,sma_i,sma_j
    real(dp)::Hill_radius_ij,delta_ij,stability_criterion
    integer::i,nci,j,ncj

    ! temp variables
    real(dp)::mui,muj,pxx,exx,ixx,maxx,wxx,lnxx,taxx,dtxx,radiusij

    hill_check=.true.

    if(NB.gt.2)then
      do i=2,(NB-1)

        nci=(i-1)*6
        mui=Giau*(m(1)+m(i))
        call elem_mer(mui,rin(1+nci:6+nci),pxx,sma_i,exx,ixx,maxx,wxx,lnxx,taxx,dtxx)
        if(sma_i.le.TOL_dp)then
!           write(*,*)'sma_i (<=0) = ',sma_i,' p_i = ',pxx
          hill_check=.false.
          return
        else if(exx.ge.one)then
          hill_check=.false.
          return
        end if

        j=i+1
        ncj=(j-1)*6
        muj=Giau*(m(1)+m(j))
        call elem_mer(muj,rin(1+ncj:6+ncj),pxx,sma_j,exx,ixx,maxx,wxx,lnxx,taxx,dtxx)
        if(sma_j.le.TOL_dp)then
!           write(*,*)'sma_j (<=0) = ',sma_j,' p_j = ',pxx
          hill_check=.false.
          return
        else if(exx.ge.one)then
          hill_check=.false.
          return
        end if

        rij=dist(rin(1+nci:3+nci), rin(1+ncj:3+ncj))
        radiusij = (R(i)+R(j))*RsunAU
        if(abs(rij-radiusij).le.TOL_dp)then
          hill_check=.false. ! if false is bad/unstable
          write(*,*)'PLANET ',i,' AND ',j,' TOUCH EACH OTHER'
          write(*,*)'rij (au) = ',rij,' <= Ri+Rj (au) = ',radiusij
          flush(6)
          return
        end if

        if(do_Hr_check)then
          Hill_radius_ij = mutual_Hill_radius(m(1),m(i),sma_i,m(j),sma_j)
          delta_ij = abs(sma_j - sma_i)
          stability_criterion = (delta_ij/Hill_radius_ij) - sqrt_12
          if(stability_criterion.le.TOL_dp)then
            hill_check = .false. ! if false is bad/unstable
            return
          end if
        end if
      end do
    end if

    return
  end function separation_mutual_Hill_check
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! OLD
  ! computes the initial state vector in the orbital reference frame
  subroutine initial_state(P,sma,ecc,mA,output)
    real(dp),dimension(:),intent(in)::P,sma,ecc,mA
    real(dp),dimension(:),intent(out)::output
    real(dp)::EA,dEA,n,cosE,sinE,ecc2,rad1e2,mA_temp
    integer::i,nci

    output=zero

    do i=2,NB
      !computation of the anomalies of the bodies
!       mA_temp=mod(mA(i),circ)+circ
      mA_temp=mod(mA(i)+circ,circ)
!       if(mA_temp.lt.zero) mA_temp=mA_temp+circ
      n=dpi/P(i)
      EA=EAnom(mA_temp,ecc(i))*deg2rad
      ecc2=ecc(i)*ecc(i)
      rad1e2=sqrt(one-ecc2)
      cosE=cos(EA)
      sinE=sin(EA)
      dEA=n/(one-ecc(i)*cosE)
      !calculation of the radius vector and velocity of the bodies
      nci=(i-1)*6
      output(1+nci)=sma(i)*(cosE-ecc(i))  !! x
      output(2+nci)=sma(i)*rad1e2*sinE  !! y
      output(4+nci)=-sma(i)*sinE*dEA  !! vx
      output(5+nci)=sma(i)*rad1e2*cosE*dEA  !! vy
    end do

    return
  end subroutine initial_state
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  ! other subroutine to compute the initial state vector in the sky frame
  ! from the keplerian elements
  subroutine kepelements2statevector(mass,sma,ecc,meanA,argp,inc,longN,rout)
    real(dp),dimension(:),intent(in)::mass,sma,ecc,meanA,argp,inc,longN
    real(dp),dimension(:),intent(out)::rout

    integer::ibd,nci
    real(dp)::mu,trueA,cosf,sinf,oneme2,rval,vval
    real(dp),dimension(:),allocatable::rtemp
    real(dp),dimension(3,3)::Rargp,Rinc,RlN

    rout=zero
    allocate(rtemp(size(rout)))
    rtemp=zero

    do ibd=2,NB
      ! mu = G(Mstar+Mplane)
      mu=Giau*(mass(1)+mass(ibd))
      trueA = calculate_true_anomaly(meanA(ibd),ecc(ibd)) ! true anomaly -> f
      cosf = cos(trueA)
      sinf = sin(trueA)
      oneme2 = one-(ecc(ibd)*ecc(ibd))
      rval = sma(ibd)*oneme2/(one+ecc(ibd)*cosf) ! r = a (1-e^2)/ (1+e*cosf)
!       vval = sqrt((mu/sma(ibd))/oneme2) ! v = sqrt(G(ms+mp)/a/(1-e^2))
      vval = sqrt(mu/(sma(ibd)*oneme2)) ! v = sqrt(G(ms+mp)/a/(1-e^2))
!       write(*,*)"body: ",ibd
!       write(*,*)" mA     = ",meanA(ibd),"e = ",ecc(ibd)," ==> f = ",trueA*rad2deg
!       write(*,*)" cosf   = ",cosf," sinf = ",sinf
!       write(*,*)" oneme2 = ",oneme2
!       write(*,*)" sma    = ",sma(ibd)
!       write(*,*)" rval   = ",rval
!       write(*,*)" vval   = ",vval

      ! TEST ONE: x,y,z (vx,vy,vz) --> rotations --> X,Y,Z (VX,VY,VZ)
      nci=(ibd-1)*6
      rtemp(1+nci)=rval*cosf ! x = rcosf
      rtemp(2+nci)=rval*sinf ! y = rsinf
      rtemp(4+nci)=-vval*sinf ! vx = v*(-sinf)
      rtemp(5+nci)=vval*(ecc(ibd)+cosf) ! vy = v*(e+cosf)
!       write(*,*)"rtemp   = ",rtemp(1+nci:6+nci)
      Rargp=zero
      call rotmat3(argp(ibd),Rargp)
!       write(*,*)"Rargp"
!       write(*,*)Rargp(1,:)
!       write(*,*)Rargp(2,:)
!       write(*,*)Rargp(3,:)
      Rinc=zero
      call rotmat1(inc(ibd),Rinc)
!       write(*,*)"Rinc"
!       write(*,*)Rinc(1,:)
!       write(*,*)Rinc(2,:)
!       write(*,*)Rinc(3,:)
      RlN=zero
      call rotmat3(longN(ibd),RlN)
!       write(*,*)"RlN"
!       write(*,*)RlN(1,:)
!       write(*,*)RlN(2,:)
!       write(*,*)RlN(3,:)
      rout(1+nci:3+nci)=matmul(RlN,matmul(Rinc,matmul(Rargp,rtemp(1+nci:3+nci))))
      rout(4+nci:6+nci)=matmul(RlN,matmul(Rinc,matmul(Rargp,rtemp(4+nci:6+nci))))
!       write(*,*)"rot3w     = ",matmul(Rargp,rtemp(1+nci:3+nci))
!       write(*,*)"rot1i3w   = ",matmul(Rinc,matmul(Rargp,rtemp(1+nci:3+nci)))
!       write(*,*)"rot3O1i3w = ",matmul(RlN,matmul(Rinc,matmul(Rargp,rtemp(1+nci:3+nci))))
!       write(*,*)
    end do
!     call orb2obs(rtemp,longN,inc,argp,rout)
    deallocate(rtemp)

    return
  end subroutine kepelements2statevector
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  ! function to determine if angle is in some predefined range
  function checklN(lN) result(clN)
    integer::clN
    real(dp),intent(in)::lN
    real(dp)::lN1


    clN=0
    lN1=mod(lN+circ,circ)
!     if(lN1.lt.zero) lN1=lN1+circ
    ! IT DEFINES THE COORDINATES FOR ECLIPSE CONDITION DUE TO THE LONGITUDE OF NODES
    if((lN1.ge.zero.and.lN1.le.45._dp)&
        &.or.&
        &(lN1.gt.135._dp.and.lN1.le.225._dp)&
        &.or.&
        &(lN1.ge.315._dp.and.lN1.le.circ))then
      clN=0
    else
      clN=1
    end if

    return
  end function checklN
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
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
    end do

    return
  end subroutine barycenter
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
  ! move a body coordinates (x,y,z,vx,vy,vz) of time = dt using the
  ! f and g functions (Murray and Dermott 1999)
  subroutine fgfunctions(mu,rout,dt,Hc)
    real(dp),intent(in)::mu,dt
    real(dp),dimension(:),intent(inout)::rout
    logical,intent(inout)::Hc

    real(dp),dimension(3)::rtemp,vtemp
!     real(dp)::r0,v0,reca,sma,n,sinE,cosE,ecc
!     real(dp)::Ek1,MAk1,MAk2,Ek2,dE
!     real(dp)::n,sinE,cosE
!     real(dp)::h2,x,y,z,vx,vy,vz,hx,hy,hz
!     real(dp)::e2
    real(dp)::r0,v0
    real(dp)::pxx,sma,ecc,ixx,mA0,wxx,tA0,lnxx,dtxx
    real(dp)::n,EA0,mA1,EA1,dEA
    real(dp)::ar0,art
    real(dp)::Ffunc,dFfunc,Gfunc,dGfunc

    !==============
    ! OLD
!     r0=dist(rout(1:3))
!     v0=dist(rout(4:6))
!     reca=(two/r0)-((v0**2)/mu)
!     if(reca.lt.zero)then
! !       write(*,*)" In fgfunctions reca < 0"
!       Hc=.false.
!       return
!     end if
!     sma=one/reca
!     if((sma.le.amin).or.(sma.ge.amax))then
! !       write(*,*)" In fgfunctions sma < amin or sma > amax"
!       Hc=.false.
!       return
!     end if
!     n=sqrt(mu/(sma**3))
!     sinE=(sum(rout(1:3)*rout(4:6)))/(n*sma**2)
!     cosE=((r0*v0**2)/mu)-one
!
!     x=rout(1)
!     y=rout(2)
!     z=rout(3)
!     vx=rout(4)
!     vy=rout(5)
!     vz=rout(6)
!     hx = y*vz-z*vy
!     hy = z*vx-x*vz
!     hz = x*vy-y*vx
!     h2=hx*hx+hy*hy+hz*hz
!
!     e2=one-h2/(sma*mu)
!     if((abs(e2).le.TOLERANCE).or.(e2.lt.zero)) e2=zero
!     ecc=sqrt(e2)
!
!     Ek1=mod(atan2(sinE,cosE)+dpi,dpi)
!
!     MAk1=Ek1-sinE
!     MAk1=mod(MAk1+dpi,dpi)
!     MAk2=(MAk1+n*dt)*rad2deg
!     MAk2=mod(MAk2+circ,circ)
! !     if(MAk2.lt.zero) MAk2=MAk2+circ
!
!     Ek2=EAnom(MAk2,ecc)*deg2rad
!     dE=Ek2-Ek1
!     if(dE.ge.pi)then
!       dE=dE-dpi
!     else if(dE.le.-pi)then
!       dE=dE+dpi
!     end if

!     ar0=sma/r0
!     Ffunc=ar0*(cos(dE)-one)+one
!     Gfunc=dt+(sin(dE)-dE)/n
!     rtemp=Ffunc*rout(1:3)+Gfunc*rout(4:6)
!     art=sma/dist(rtemp)
!     dFfunc=-ar0*art*n*sin(dE)
!     dGfunc=art*(cos(dE)-one)+one
!     vtemp=dFfunc*rout(1:3)+dGfunc*rout(4:6)
!
!     rout(1:3)=rtemp
!     rout(4:6)=vtemp
    !==============

    call elem_mer(mu,rout,pxx,sma,ecc,ixx,mA0,wxx,lnxx,tA0,dtxx) ! in rad!
    if(sma.le.TOL_dp.or.ecc.gt.one)then
!       write(*,*)" mu = ",mu
!       write(*,*)" rout = ",rout
!       write(*,*)" sma = ",sma
!       write(*,*)" pxx = ",pxx
!       write(*,*)" ecc = ",ecc
      Hc=.false.
!       stop('ERRORRRRR')
      return
    end if
    if(ecc.le.TOLERANCE)then
      EA0=tA0
    else
!       EA0=mod(EAnom(mA0*rad2deg,ecc)+circ,circ)*deg2rad
      EA0=trueAnom_ecc_to_eccAnom(tA0*rad2deg, ecc)
    end if
    n=dpi/pxx
    mA1=mA0+dt*n ! not checking if <0 or greater than 2pi, I want to know the 'direction'
    EA1=EAnom(mA1*rad2deg,ecc)*deg2rad
    dEA=EA1-EA0

    r0=dist(rout(1:3))
    v0=dist(rout(4:6))

    ar0=sma/r0
    Ffunc=ar0*(cos(dEA)-one)+one
    Gfunc=dt+(sin(dEA)-dEA)/n
    rtemp=Ffunc*rout(1:3)+Gfunc*rout(4:6)
    art=sma/dist(rtemp)
    dFfunc=-ar0*art*n*sin(dEA)
    dGfunc=art*(cos(dEA)-one)+one
    vtemp=dFfunc*rout(1:3)+dGfunc*rout(4:6)

    rout(1:3)=rtemp
    rout(4:6)=vtemp

    return
  end subroutine fgfunctions
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! OLD, NOT USED
  ! from state vector (astrocentri cartesian coordinates)
  ! to keplerian orbital elements
  ! for one body
  subroutine eleMD(mu,svec,P,sma,ecc,inc,mA,w,lN,f,dt)
    real(dp),intent(in)::mu
    real(dp),intent(in),dimension(:)::svec
    real(dp),intent(out)::P,sma,ecc,inc,mA,w,lN,f,dt

    real(dp)::x,y,z,vx,vy,vz,R,R2,V,V2,rrd,rd2,rd
    real(dp)::hx,hy,hz,h2,h,signx,signy

    real(dp)::inv_sma,musma,ecc2
    real(dp)::cosi,sini

    real(dp)::coslN,sinlN
    real(dp)::wf,sinwf,coswf
    real(dp)::p_slr,ecosf,esinf

    real(dp)::Ea !,mmotion


    ! init kep elem: P,sma,ecc,inc,mA,w,lN,f,dt
    P=zero
    sma=zero
    ecc=one
    inc=zero
    mA=zero
    w=half*pi
    lN=pi
    f=zero
    dt=zero

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

    hx = y*vz-z*vy
    hy = z*vx-x*vz
    hz = x*vy-y*vx
    h2=hx*hx+hy*hy+hz*hz
    h=sqrt(h2)
    ! Murray & Dermott 1999
    ! if hz > 0 ==> +hx, -hy
    if(hz.gt.zero)then
      signx=one
      signy=-one
    ! if hz < 0 ==> -hx, +hy
    else
      signx=-one
      signy=one
    end if

    ! semi-major axis = sma
    ! 1/sma
!     inv_sma=(two*mu-R*V2)/(R*mu)
    inv_sma=(two/R)-(V2/mu)
    if(inv_sma.le.zero)then ! 1/sma <= 0
      write(*,*)'1/a < 0: ',inv_sma
      return
    else
     sma=one/inv_sma
    end if
    ! period = P
    P=dpi*sqrt((sma**3)/mu)

    ! eccentricity = ecc
    ! ecc2 = 1 - h^2/(mu*sma)
    musma=mu*sma
    ecc2 = one - (h2/musma)
    if(ecc2.gt.TOLERANCE)then
      ecc=sqrt(ecc2)
    else
      ecc=zero
    end if

    ! inclination = inc
    cosi=hz/h
    inc=acos(cosi) ! rad
    sini=sin(inc)

    ! longitude of node = lN
    if(abs(sini).le.TOLERANCE)then
      coslN=cos(pi)
      sinlN=sin(pi)
    else
      coslN=signy*hy/(h*sini)
      sinlN=signx*hx/(h*sini)
    end if
    lN=mod(atan2(sinlN,coslN)+dpi,dpi)

    ! argument of pericentre = w
    ! and
    ! true anomaly = f
    ! w+f
    if(abs(sini).gt.TOLERANCE)then ! sini != 0 ==> cosi cold be 0
      sinwf=z/(R*sini)
    else ! sin == 0 ==> cosi != 0
      sinwf=(y*coslN-x*sinlN)/(R*cosi)
    end if
    coswf=((x/R)+(sinlN*sinwf*cosi))/coslN
    wf=mod(atan2(sinwf,coswf)+dpi,dpi)
    ! semi-latus rectum = p_slr
    p_slr=sma*(one-ecc2)
    ! RRdot
    rrd=x*vx+y*vy+z*vz ! R Rdot
    rd2 = V2-(h2/R2)
    if(rd2.lt.zero)then
      rd = zero
    else
      rd=sqrt(rd2) ! |Rdot|
    end if
    rd=sign(rd,rrd)
    ! f
    esinf=p_slr*rd/h
    ecosf=(p_slr-R)/R ! or p_slr/R - 1 ? what is the best way to code it?
    f=mod(atan2(esinf,ecosf)+dpi,dpi)
    ! w and mean anomaly = mA
    if(ecc.le.TOLERANCE)then
      mA=f
      w=half*pi
    else
      if(ecc.eq.one)then ! e == 1
        mA=tan(half*f)+(tan(half*f)**3)/three
      else if(ecc.gt.one)then ! e > 1
        Ea=two*atanh(sqrt((ecc+one)/(ecc-one))*tan(half*f))
        mA=ecc*sinh(Ea)-Ea
      else ! e > 0 & e < 1
        Ea=two*atan(sqrt((one-ecc)/(one+ecc))*tan(half*f))
        mA=Ea-ecc*sin(Ea)
      end if
      w=mod(mod(wf-f+dpi,dpi)+dpi,dpi)
    end if
    ! t - tau
    dt=mA*P/dpi

    return
  end subroutine eleMD
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! ADAPTED FROM MERCURY
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_X2EL.FOR    (ErikSoft  20 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates Keplerian orbital elements given relative coordinates and
! velocities, and GM = G times the sum of the masses.
!
! The elements are: q = perihelion distance
!                   ecc = eccentricity
!                   inc = inclination
!                   p = longitude of perihelion (NOT argument of perihelion!!)
!                   ln = longitude of ascending node
!                   mA = mean anomaly (or mean longitude if e < 1.e-8)
!
!-------------------------------------------------------------------------------
!
  subroutine mco_x2el(mu,svec,q,ecc,inc,p,ln,mA,f)
    implicit none
    ! include 'mercury.inc'

  ! Input/Output
    real(dp),intent(in)::mu
    real(dp),intent(in),dimension(:)::svec
!     real(dp),intent(inout)::q,ecc,inc,p,ln,mA,f
    real(dp),intent(out)::q,ecc,inc,p,ln,mA,f

    ! Local
    real(dp)::x,y,z,u,v,w
    real(dp)::hx,hy,hz,h2,h,v2,r,rv,s,true
    real(dp)::ci,t_o,temp,tmp2,bige,cf,ce

    ! init all output variables
    q=zero
    ecc=zero
    inc=zero
    p=zero
    ln=zero
    mA=zero
    f=zero

    x=svec(1)
    y=svec(2)
    z=svec(3)
    u=svec(4)
    v=svec(5)
    w=svec(6)

    hx = y * w  -  z * v
    hy = z * u  -  x * w
    hz = x * v  -  y * u
    h2 = hx*hx + hy*hy + hz*hz
    v2 = u * u  +  v * v  +  w * w
    rv = x * u  +  y * v  +  z * w
    r = sqrt(x*x + y*y + z*z)
    h = sqrt(h2)
    s = h2 / mu

    ! Inclination and node
    ci = hz / h
    if (abs(ci).lt.one) then
      inc = acos(ci)
      ln = atan2(hx,-hy)
      if (ln.lt.zero) ln = ln + dpi
    else
!       if (ci.gt.zero) inc = zero
!       if (ci.lt.zero) inc = pi
      if (ci.gt.zero)then
        inc = zero
      else
        inc = pi
      end if
      ln = pi
    end if

    ! Eccentricity and perihelion distance
    ecc=zero
    temp = one  +  s * (v2 / mu  -  two / r)

!     write(*,*)
!     write(*,*)" x         = ",x
!     write(*,*)" y         = ",y
!     write(*,*)" z         = ",z
!     write(*,*)" vx        = ",u
!     write(*,*)" vy        = ",v
!     write(*,*)" vz        = ",w
!     write(*,*)" s = h2/mu = ",s
!     write(*,*)" h2        = ",h2
!     write(*,*)" mu        = ",mu
!     write(*,*)" v2        = ",v2
!     write(*,*)" v2/mu     = ",v2/mu
!     write(*,*)" r         = ",r
!     write(*,*)" 2/r       = ",two/r
!     write(*,*)" temp      = ",temp

    ! if (temp.le.zero) then
    if (temp.le.TOL_dp) then
      ecc = zero
      q = s
    else
      ecc = sqrt(temp)
      q = s / (one + ecc)
      if(ecc.gt.one) q = -q ! hyperbola added by Luca
    end if
!     write(*,*)" ecc      = ",ecc
!     write(*,*)
!     q = s / (one + ecc)

    ! True longitude
    if (hy.ne.zero) then
      t_o = -hx/hy
      temp = (one - ci) * t_o
      tmp2 = t_o * t_o
      true = mod(atan2((y*(one+tmp2*ci)-x*temp),(x*(tmp2+ci)-y*temp)) + dpi, dpi)
    else
      true = mod(atan2(y * ci, x) + dpi, dpi)
    end if
    if (ci.lt.zero) true = mod(true + pi,dpi)

    if (ecc.le.TOLERANCE) then
      ! p = zero
      p = -half*pi
      mA = true
      f = mA
    else
      ce = (v2*r - mu) / (ecc*mu)

      ! Mean anomaly for ellipse
      if (ecc.lt.one) then
        if (abs(ce).gt.one) ce = sign(one,ce)
        bige = acos(ce)
        if (rv.lt.zero) bige = dpi - bige
        mA = bige - ecc*sin(bige)
      else

        ! Mean anomaly for hyperbola
        if (ce.lt.one) ce = one
        bige = log( ce + sqrt(ce*ce-one) )
        if (rv.lt.zero) bige = - bige
        mA = ecc*sinh(bige) - bige
      end if

      ! Longitude of perihelion
      cf = (s - r) / (ecc*r)
      if (abs(cf).gt.one) cf = sign(one,cf)
      f = acos(cf)
      if (rv.lt.zero) f = dpi - f
      p = true - f
      p = mod (p + dpi + dpi, dpi)
    end if

    if (mA.lt.zero.and.ecc.lt.one) mA = mA + dpi
    if (mA.gt.dpi.and.ecc.lt.one) mA = mod (mA, dpi)

    return
  end subroutine mco_x2el
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! call the mco_x2el and 'completes' the kepler elements
  subroutine elem_mer(mu,svec,P,sma,ecc,inc,mA,w,lN,f,dt)
    real(dp),intent(in)::mu
    real(dp),intent(in),dimension(:)::svec
    real(dp),intent(out)::P,sma,ecc,inc,mA,w,lN,f,dt

    real(dp)::q,lperi
!     real(dp)::cf

    q=zero
    lperi=-half*pi
!     cf=zero

    call mco_x2el(mu,svec,q,ecc,inc,lperi,lN,mA,f)
    ! semi-major axis sma
    sma=q/(one-ecc)
    ! arg. pericentre w
    w=mod(lperi-lN+dpi,dpi)
    if(ecc.le.TOLERANCE)w=half*pi

    P=dpi*sqrt((sma**3)/mu)

    dt=mA*P/dpi

!     if(ecc.ge.one)then
!       write(*,*)'mu = ',mu
!       write(*,*)'svec = ',svec
!       write(*,*)'P,sma,ecc,inc,mA,w,lN,f,dt'
!       write(*,*)P,sma,ecc,inc,mA,w,lN,f,dt
!     end if

  end subroutine elem_mer
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  ! from state vector (astrocentric cartesian coordinates)
  ! to keplerian orbital elements
  ! for all the bodies
  subroutine elements(m,rin,P,sma,ecc,inc,mA,argp,tA,lN,dttau)
    real(dp),dimension(:),intent(in)::m,rin
    real(dp),dimension(:),intent(out)::P,sma,ecc,inc,mA,argp,tA,lN,dttau
    real(dp),dimension(6)::svec
    real(dp)::mu
    integer::j,ncj

    P=zero
    sma=zero
    ecc=zero
    inc=90._dp
    mA=zero ! mean anomaly
    argp=90._dp ! argument of pericentre
    tA=zero ! true anomaly
    lN=180._dp ! longitude of node
    dttau=zero

    cicle: do j=2,NB
      ncj=(j-1)*6
      mu=Giau*(m(1)+m(j))
      ! state vector j-th body
      svec=rin(1+ncj:6+ncj)
!       call eleMD(mu,svec,P(j),sma(j),ecc(j),inc(j),mA(j),argp(j),lN(j),tA(j),dttau(j))
      call elem_mer(mu,svec,P(j),sma(j),ecc(j),inc(j),mA(j),argp(j),lN(j),tA(j),dttau(j))
      inc(j)=inc(j)*rad2deg
      mA(j)=mA(j)*rad2deg
      lN(j)=lN(j)*rad2deg
      argp(j)=argp(j)*rad2deg
      tA(j)=tA(j)*rad2deg
    end do cicle

    return
  end subroutine elements
! ------------------------------------------------------------------------------

end module celestial_mechanics
