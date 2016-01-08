module ode_run
  use constants
  use parameters
  use celestial_mechanics
  use rotations,only:orb2obs
  use eq_motion,only:eqmastro
  use numerical_integrator,only:int_rk_a
  use transits
  use radial_velocities
  use output_files
  implicit none

  contains

  ! ------------------------------------------------------------------ !
  ! given the RV and T0 simulations it assigns the residuals value as (obs-sim)/eobs
  subroutine setval(RV_obs,RV_sim,T0_obs,T0_sim,val)
    real(dp),dimension(:),intent(in)::RV_obs,RV_sim
    real(dp),dimension(:,:),intent(in)::T0_obs,T0_sim
    real(dp),dimension(:),intent(out)::val
    real(dp)::obs,sim,eobs
    integer::j,j1,a

    val=zero
    if(nRV.gt.0) val(1:nRV)=(RV_obs-RV_sim)/eRVobs
    a=nRV
    do j=2,NB
      do j1=1,nT0(j)
        a=a+1
        obs=T0_obs(j1,j)
        sim=T0_sim(j1,j)
        eobs=eT0obs(j1,j)
        val(a)=(obs-sim)/eobs
      end do
    end do

    return
  end subroutine setval

  ! as setval but it accepts different arguments
  subroutine setval_2(RV_obs,RV_sim,e_RVobs,T0_obs,T0_sim,e_T0obs,val)
    real(dp),dimension(:),intent(in)::RV_obs,RV_sim,e_RVobs
    real(dp),dimension(:,:),intent(in)::T0_obs,T0_sim,e_T0obs
    real(dp),dimension(:),intent(out)::val
    real(dp)::obs,sim,eobs
    integer::j,j1,a

    val=zero
    if(nRV.gt.0) val(1:nRV)=(RV_obs-RV_sim)/e_RVobs
    a=nRV
    do j=2,NB
      do j1=1,nT0(j)
        a=a+1
        obs=T0_obs(j1,j)
        sim=T0_sim(j1,j)
        eobs=e_T0obs(j1,j)
        val(a)=(obs-sim)/eobs
      end do
    end do

    return
  end subroutine setval_2

  ! as setval but it accepts different arguments, i.e., gamma of RV
  subroutine setval_3(RV_obs,RV_sim,e_RVobs,gamma,T0_obs,T0_sim,e_T0obs,val)
    use statistics,only:wmean
    real(dp),dimension(:),intent(in)::RV_obs,RV_sim,e_RVobs
    real(dp),dimension(:),intent(out)::gamma
    real(dp),dimension(:,:),intent(in)::T0_obs,T0_sim,e_T0obs
    real(dp),dimension(:),intent(out)::val
    real(dp),dimension(:),allocatable::dRV,wi
    real(dp)::obs,sim,eobs
    integer::j,j1,a

    gamma=zero
    val=zero
    if(nRV.gt.0)then
      allocate(dRV(nRV),wi(nRV))
      dRV=zero
      dRV=RV_obs-RV_sim
      wi=one/(e_RVobs*e_RVobs)
      gamma(1)=wmean(dRV,wi)
      gamma(2)=sqrt(one/sum(wi))
      deallocate(dRV,wi)
      val(1:nRV)=(RV_obs-(gamma(1)+RV_sim))/e_RVobs
    end if
    a=nRV
    do j=2,NB
      do j1=1,nT0(j)
        a=a+1
        obs=T0_obs(j1,j)
        sim=T0_sim(j1,j)
        eobs=e_T0obs(j1,j)
        val(a)=(obs-sim)/eobs
      end do
    end do

    return
  end subroutine setval_3

!   as setval but it accepts different arguments, i.e., gamma of RV, but for multiple set of RV
  subroutine setval_4(RV_obs,RV_sim,e_RVobs,gamma,T0_obs,T0_sim,e_T0obs,resw)
!     use statistics,only:wmean
    real(dp),dimension(:),intent(in)::RV_obs,RV_sim,e_RVobs
    real(dp),dimension(:,:),allocatable,intent(out)::gamma
    real(dp),dimension(:,:),intent(in)::T0_obs,T0_sim,e_T0obs
    real(dp),dimension(:),intent(out)::resw
    real(dp),dimension(:),allocatable::val

    allocate(val(ndata))
    call set_RV_resw(RV_obs,RV_sim,e_RVobs,gamma,val(1:nRV))
    call set_T0_resw(T0_obs,T0_sim,e_T0obs,val(nRV+1:ndata))
    call set_fitness(val,resw)
    deallocate(val)
    
    return
  end subroutine setval_4

  ! set only RV to the weighted residuals
  subroutine set_RV_resw(RV_obs,RV_sim,e_RVobs,gamma,resw)
    use statistics,only:wmean
    real(dp),dimension(:),intent(in)::RV_obs,RV_sim,e_RVobs
    real(dp),dimension(:,:),allocatable,intent(out)::gamma
    real(dp),dimension(:)::resw
    real(dp),dimension(:),allocatable::dRV,wi
    integer::j,j1,a,aRV,bRV

    if(.not.allocated(gamma)) allocate(gamma(nRVset,2))
    gamma=zero
    resw=zero
    if(nRV.gt.0)then
      aRV=0
      bRV=0
      do j=1,nRVset
        aRV=bRV+1
        bRV=bRV+nRVsingle(j)
!         write(*,*)" aRV = ",aRV," bRV = ",bRV
        allocate(dRV(nRVsingle(j)),wi(nRVsingle(j)))
        dRV=zero
        dRV=RV_obs(aRV:bRV)-RV_sim(aRV:bRV)
        wi=one/(e_RVobs(aRV:bRV)*e_RVobs(aRV:bRV))
        gamma(j,1)=wmean(dRV,wi)
        gamma(j,2)=sqrt(one/sum(wi))
        deallocate(dRV,wi)
        resw(aRV:bRV)=(RV_obs(aRV:bRV)-(gamma(j,1)+RV_sim(aRV:bRV)))/e_RVobs(aRV:bRV)
      end do
    end if
    
    return
  end subroutine set_RV_resw
  
  ! set only T0 to the weighted residuals
  subroutine set_T0_resw(T0_obs,T0_sim,eT0_obs,resw)
    real(dp),dimension(:,:),intent(in)::T0_obs,T0_sim,eT0_obs
    real(dp),dimension(:),intent(out)::resw
    real(dp)::obs,sim,eobs
    integer::j,j1,a
    ! NB,nT0 are common variables
    resw=zero
    a=0
    do j=2,NB
      do j1=1,nT0(j)
        a=a+1
        obs=T0_obs(j1,j)
        sim=T0_sim(j1,j)
        eobs=eT0_obs(j1,j)
        resw(a)=(obs-sim)/eobs
      end do
    end do
    
    return
  end subroutine set_T0_resw
  
  subroutine set_fitness(val,resw)
    real(dp),dimension(:),intent(in)::val
    real(dp),dimension(:),intent(out)::resw
    real(dp),dimension(:),allocatable::resa1,resa2,resb1,resb2

    ! set the values scaled by the k_chi2r and k_chi2wr
    ! They are set in such a way the sum of resw^2 is:
    resw=zero
    if(k_chi2wr.eq.zero)then
      ! Chi2r
      resw=val*sqrt(inv_dof)

    else if(k_chi2r.eq.zero)then
      ! Chi2wr
      if(nRV.eq.0.or.sum(nT0).eq.0)then
        resw=val*k_b
        
      else
        resw(1:nRV)=val(1:nRV)*k_b(1)
        resw(nRV+1:ndata)=val(nRV+1:ndata)*k_b(2)
        
      end if

    else
      ! fitness = Chi2r*k_chi2r + Chi2wr*k_chi2wr
      ! NO RV: nRV=0, nT0!=0
      if(nRV.eq.0.and.sum(nT0).ne.0)then
        allocate(resa2(sum(nT0)),resb2(sum(nT0)))
        resa2=val*k_a
        resb2=val*k_b
        resw=sqrt(resa2*resa2 + resb2*resb2)
        deallocate(resa2,resb2)
      
      ! NO T0: nRV!=0, nT0=0
      else if(nRV.ne.0.and.sum(nT0).eq.0)then
        allocate(resa1(nRV),resb1(nRV))
        resa1=val*k_a
        resb1=val*k_b
        resw=sqrt(resa1*resa1 + resb1*resb1)
        deallocate(resa1,resb1)
        
      else
        allocate(resa1(nRV),resb1(nRV),resa2(sum(nT0)),resb2(sum(nT0)))
        resa1=val(1:nRV)*k_a
        resb1=val(1:nRV)*k_b(1)
        resw(1:nRV)=sqrt(resa1*resa1 + resb1*resb1)
        resa2=val(nRV+1:ndata)*k_a
        resb2=val(nRV+1:ndata)*k_b(2)
        resw(nRV+1:ndata)=sqrt(resa2*resa2 + resb2*resb2)
        deallocate(resa1,resa2,resb1,resb2)
      end if
      
    end if
  
    return
  end subroutine set_fitness
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  function set_max_residuals(ndata) result(resw_max)
    integer,intent(in)::ndata
    real(dp)::resw_max
    
    resw_max = sqrt(resmax/real(ndata,dp))
    
    return
  end function set_max_residuals
  
  subroutine check_max_residuals(resw,ndata)
    real(dp),dimension(:),intent(inout)::resw
    integer,intent(in)::ndata
    real(dp)::resw_max,resw_max_possible
    
    resw_max=maxval(resw)
    resw_max_possible=set_max_residuals(ndata)
    if(resw_max.gt.resw_max_possible) resw=resw_max_possible
  
    return
  end subroutine check_max_residuals
  
  
  ! ------------------------------------------------------------------ !
  
  
  ! ------------------------------------------------------------------ !
  ! performs the integration, integration time to be provided as argument
  subroutine ode_a(m,R,rin,time,clN,cntRV,RV_stat,RV_sim,cntT0,T0_stat,T0_sim,Hc)
    real(dp),dimension(:),intent(in)::m,R,rin
    real(dp),intent(in)::time
    integer,dimension(:),intent(in)::clN
    integer,intent(inout)::cntRV
    integer,dimension(:),intent(inout)::RV_stat
    real(dp),dimension(:),intent(inout)::RV_sim
    integer,intent(inout)::cntT0
    real(dp),dimension(:,:),intent(inout)::T0_sim
    integer,dimension(:,:),intent(inout)::T0_stat
    logical,intent(out)::Hc

    real(dp),dimension(:),allocatable::dr,r1,r2,err
    integer,dimension(:),allocatable::X,Y,Z
    real(dp),dimension(:),allocatable::cX,cY,cR,rmean
    real(dp)::hw,hok,hnext,itime
    integer::j,j1 !,nj

    !real(dp)::tmid

!     Hc=Hillcheck(m,rin)
    Hc=mutual_Hill_check(m,rin)
    if(Hc) return

    allocate(X(NB),Y(NB),Z(NB),cX(NB),cY(NB),cR(NB),rmean(NB))
    X=0
    Y=0
    Z=0
    do j=2,NB
      X(j)=1+(j-1)*6
      Y(j)=2+(j-1)*6
      Z(j)=3+(j-1)*6
    end do
    cX=1._dp
    cY=1._dp
    cR=zero
    cR(2:NB)=1.5_dp*(R(1)+R(2:NB))*RsunAU
    rmean=9.e3_dp

    allocate(dr(NBDIM),r1(NBDIM),r2(NBDIM),err(NBDIM))
    hw=step_0
    if(time.lt.zero) hw=-hw
    itime=zero
    r1=rin
    r2=zero
    err=zero

    j1=0
    integration: do
      j1=j1+1
      if(abs(itime+hw).gt.abs(time)) hw=time-itime
      call eqmastro(m,r1,dr)
      call int_rk_a(m,r1,dr,hw,hok,hnext,r2,err)

!       Hc=Hillcheck(m,r2)
      Hc=mutual_Hill_check(m,r2)
      if(Hc) return

      ! RV check
      if(cntRV.lt.nRV) call check_RV(m,r1,dr,itime,hok,cntRV,RV_stat,RV_sim)

      ! T0 check
      do j=2,NB
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))
        rmean(j)=0.5*( rsky(r1(X(j):Y(j))) + rsky(r2(X(j):Y(j))) )
      end do

      if((idtra.gt.0).and.(idtra.le.NB))then
        if(cntT0.lt.sum(nT0))then
          if(idtra.eq.1)then
            do j=2,NB
              if(rmean(j).le.cR(j))then
                if(clN(j).eq.0)then
                  if((cX(j).le.0._dp).and.(r1(Z(j)).gt.0._dp))&
                      &call check_T0(j,m,R,r1,r2,&
                      &itime,hok,T0_stat,T0_sim,Hc)
                else
                  if((cY(j).le.0._dp).and.(r1(Z(j)).gt.0._dp))&
                      &call check_T0(j,m,R,r1,r2,&
                      &itime,hok,T0_stat,T0_sim,Hc)
                end if
              end if
            end do
          else
            if(rmean(idtra).le.cR(idtra))then
              if(clN(idtra).eq.0)then
                if((cX(idtra).le.0._dp).and.(r1(Z(idtra)).gt.0._dp))&
                    &call check_T0(idtra,m,R,r1,r2,&
                    &itime,hok,T0_stat,T0_sim,Hc)
              else
                if((cY(idtra).le.0._dp).and.(r1(Z(idtra)).gt.0._dp))&
                    &call check_T0(idtra,m,R,r1,r2,&
                    &itime,hok,T0_stat,T0_sim,Hc)
              end if
            end if
          end if
          cntT0=sum(T0_stat)
        end if
      end if

      itime=itime+hok
      if(abs(itime).ge.abs(time)) exit integration
      hw=hnext
      r1=r2
    end do integration
    deallocate(X,Y,Z,cX,cY,cR,rmean)
    deallocate(dr,r1,r2,err)

    return
  end subroutine ode_a
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! subroutine that integrates orbit, computes, stores and
  ! writes orbits, orbital elements, constants of motion, transit time and RV into files
  ! it is called by ode_out, that check in which direction (in time) integrates
  subroutine ode_a_o(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
      &m,R,rin,time,clN,cntRV,RV_stat,RV_sim,cntT0,T0_stat,T0_sim,Hc)
    integer,intent(in)::uorb,ucon
    integer,dimension(:),intent(in)::uele,utra
    character(*),intent(in)::fmorb,fmcon,fmele
    real(dp),dimension(:),intent(in)::m,R,rin
    real(dp),intent(in)::time
    integer,dimension(:),intent(in)::clN
    integer,intent(inout)::cntRV
    integer,dimension(:),intent(inout)::RV_stat
    real(dp),dimension(:),intent(inout)::RV_sim
    integer,intent(inout)::cntT0
    real(dp),dimension(:,:),intent(inout)::T0_sim
    integer,dimension(:,:),intent(inout)::T0_stat
    logical,intent(out)::Hc

    real(dp),dimension(:),allocatable::dr,r1,r2,err
    integer,dimension(:),allocatable::X,Y,Z
    real(dp),dimension(:),allocatable::cX,cY,cR,rmean
    real(dp)::hw,hok,hnext,itime
    integer::j,j1,j2,jtra
    real(dp),dimension(:,:),allocatable::storeorb,storecon,storetra
    integer,dimension(:,:),allocatable::stat_tra
    integer::Norb
    real(dp)::itwrt,Etot,Eold,htot,hold

    real(dp)::ftime
    integer::iftime,compare

!     Hc=Hillcheck(m,rin)
    Hc=mutual_Hill_check(m,rin)
    if(Hc) return

    allocate(X(NB),Y(NB),Z(NB),cX(NB),cY(NB),cR(NB),rmean(NB))
    X=0
    Y=0
    Z=0
    do j=2,NB
      X(j)=1+(j-1)*6
      Y(j)=2+(j-1)*6
      Z(j)=3+(j-1)*6
    end do
    cX=1._dp
    cY=1._dp
    cR=zero
    cR(2:NB)=1.5_dp*(R(1)+R(2:NB))*RsunAU
    rmean=9.e3_dp

    Norb=NBDIM+3

    allocate(dr(NBDIM),r1(NBDIM),r2(NBDIM),err(NBDIM))
    hw=step_0
    if(time.lt.zero) hw=-hw
    itime=zero
    r1=rin
    r2=zero
    err=zero

    j1=0
    j2=1

    itwrt=wrttime
    if((wrtorb.eq.1).or.(wrtel.eq.1))then
      allocate(storeorb(Norb,DIMMAX))
      storeorb=zero
      call store_orb(j2,itime,m,r1,storeorb)
    end if
    if(wrtconst.eq.1)then
      Etot=zero
      Eold=zero
      htot=zero
      hold=zero
      call compute_con(m,r1,Eold,hold)
      allocate(storecon(5,DIMMAX))
      storecon=zero
      call store_con(j2,itime,Eold,Eold,hold,hold,storecon)
    end if
    if((idtra.ge.1).and.(idtra.le.NB))then
      allocate(storetra(NBDIM+6,DIMMAX))
      storetra=zero
      allocate(stat_tra(NB,DIMMAX))
      stat_tra=0
      jtra=0
    end if

    ftime=zero
    iftime=0
    compare=10

    integration: do
      j1=j1+1

      if(abs(itime+hw).gt.abs(time)) hw=time-itime
      call eqmastro(m,r1,dr) ! computes the eq. of motion
      call int_rk_a(m,r1,dr,hw,hok,hnext,r2,err) ! computes the next orbit step

      if(wrtconst.eq.1) call compute_con(m,r2,Etot,htot) ! if selected computed the Energy and Angular momentum

      if(abs(itime+hw).ge.(abs(itwrt)))then
        j2=j2+1
        if((wrtorb.eq.1).or.(wrtel.eq.1))then
          call store_orb(j2,itime+hok,m,r2,storeorb)
          if(j2.eq.DIMMAX)then
            if(wrtorb.eq.1) call write_file(DIMMAX,uorb,fmorb,storeorb)
            if(wrtel.eq.1) call write_elem(DIMMAX,uele,fmele,m,storeorb)
          end if
        end if
        if(wrtconst.eq.1)then
          call store_con(j2,itime+hok,Etot,Eold,htot,hold,storecon)
          if(j2.eq.DIMMAX) call write_file(DIMMAX,ucon,fmcon,storecon)
        end if
        if(j2.eq.DIMMAX) j2=0
        itwrt=itwrt+wrttime
      end if

!       Hc=Hillcheck(m,r2)
      Hc=mutual_Hill_check(m,r2)
      if(Hc) return

      ! RV check
      if(cntRV.lt.nRV) call check_RV(m,r1,dr,itime,hok,cntRV,RV_stat,RV_sim)

      ! T0 check (to compare and all)
      do j=2,NB
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))
        rmean(j)=0.5*( rsky(r1(X(j):Y(j))) + rsky(r2(X(j):Y(j))) )
      end do
      if((idtra.gt.0).and.(idtra.le.NB))then
        !if(cntT0.lt.sum(nT0))then
        if(idtra.eq.1)then
          do j=2,NB
            if(rmean(j).le.cR(j))then
              if(clN(j).eq.0)then
                if((cX(j).le.0._dp).and.(r1(Z(j)).gt.0._dp))then
                  if(nT0(j).gt.0) call check_T0(j,m,R,r1,r2,&
                      &itime,hok,T0_stat,T0_sim,Hc)
                  jtra=jtra+1
                  call all_transits(jtra,j,m,R,r1,r2,&
                      &itime,hok,stat_tra,storetra)
                end if
              else
                if((cY(j).le.0._dp).and.(r1(Z(j)).gt.0._dp))then
                  if(nT0(j).gt.0)call check_T0(j,m,R,r1,r2,&
                      &itime,hok,T0_stat,T0_sim,Hc)
                  jtra=jtra+1
                  call all_transits(jtra,j,m,R,r1,r2,&
                      &itime,hok,stat_tra,storetra)
                end if
              end if
            end if
          end do
        else
          if(rmean(idtra).le.cR(idtra))then
            if(clN(idtra).eq.0)then
              if((cX(idtra).le.0._dp).and.&
                  &(r1(Z(idtra)).gt.0._dp))then
                if(nT0(idtra).gt.0)call check_T0(idtra,m,R,r1,r2,&
                    &itime,hok,T0_stat,T0_sim,Hc)
                jtra=jtra+1
                call all_transits(jtra,idtra,m,R,r1,r2,&
                    &itime,hok,stat_tra,storetra)
              end if
            else
              if((cY(idtra).le.0._dp).and.&
                  &(r1(Z(idtra)).gt.0._dp))then
                if(nT0(idtra).gt.0)call check_T0(idtra,m,R,r1,r2,&
                    &itime,hok,T0_stat,T0_sim,Hc)
                jtra=jtra+1
                call all_transits(jtra,idtra,m,R,r1,r2,&
                    &itime,hok,stat_tra,storetra)
              end if
            end if
          end if
        end if
        if(sum(nT0).gt.0) cntT0=sum(T0_stat)

        if(jtra.eq.DIMMAX)then
          call write_tra(jtra,utra,stat_tra,storetra)
          jtra=0
          stat_tra=0
          storetra=zero
        end if

      end if

      itime=itime+hok

      ftime=(itime/time)*100._dp
      iftime=int(ftime+0.5_dp)
      if(iftime.eq.compare)then
        write(*,'(a,i3,a)')" done the ",compare," % "
        compare=compare+10
      end if

      if(abs(itime).ge.abs(time)) exit integration
      hw=hnext
      r1=r2

    end do integration

    if((idtra.ge.1).and.(idtra.le.NB)) call write_tra(jtra,utra,stat_tra,storetra)

    if((wrtorb.eq.1).or.(wrtel.eq.1))then
      if(wrtorb.eq.1) call write_file(j2,uorb,fmorb,storeorb)
      if(wrtel.eq.1) call write_elem(j2,uele,fmele,m,storeorb)
      deallocate(storeorb)
    end if
    if(wrtconst.eq.1)then
      call write_file(j2,ucon,fmcon,storecon)
      deallocate(storecon)
    end if

    deallocate(X,Y,Z,cX,cY,cR,rmean)
    deallocate(dr,r1,r2,err)

    return
  end subroutine ode_a_o
  ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
  ! subroutine called by the L-M to fit the parameters
  subroutine ode_lm(allpar,nm,nn,par,resw,iflag)
    integer,intent(in)::nm,nn
    real(dp),dimension(:),intent(in)::allpar,par
    real(dp),dimension(:),intent(out)::resw
    integer,intent(inout)::iflag

    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,i,lN
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp)::dt1,dt2
    logical::Hc
    !real(dp),parameter::resmax=huge(0._dp)
    !real(dp),parameter::resmax=1.e6_dp
    integer::cntRV
    real(dp),dimension(:),allocatable::RV_sim
    integer,dimension(:),allocatable::RV_stat
    real(dp),dimension(:,:),allocatable::gamma
    integer::cntT0,nTs
    real(dp),dimension(:,:),allocatable::T0_sim
    integer,dimension(:,:),allocatable::T0_stat
    logical::checkpar

    resw=zero
    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),w(NB),mA(NB),i(NB),lN(NB),clN(NB))

    ! from the allpar and par to keplerian elements: due to the LMdif approach/code
!     call par2kel(allpar,par,m,R,P,a,e,w,mA,i,lN)
    if(progtype.le.1)then
      call par2kel(allpar,par,m,R,P,a,e,w,mA,i,lN)
      checkpar = checkbounds(par)
    else
      call par2kel_fit(allpar,par,m,R,P,a,e,w,mA,i,lN,checkpar)
      if(.not.checkpar)then
        resw=set_max_residuals(ndata)
        return
      end if
      checkpar = checkbounds_fit(par)
    end if
    if(.not.checkpar)then
!       resw=sqrt(resmax)/real(ndata,dp)
        resw=set_max_residuals(ndata)
      return
    end if

    
    ! it is needed to define which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))
    ! IT CREATES THE INITIAL STATE VECTOR FROM KEPLERIAN ORBITAL ELEMENS IN THE ORBITAL PLANE
    call initial_state(P,a,e,mA,ra0)
    ! IT ROTATES THE STATE VECTORS FROM ORBITAL REF. SYSTEM TO THE OBSERVER REF. SYSTEM
    call orb2obs(ra0,-lN,-i,-w,ra1)

    cntRV=0
    if(nRV.gt.0)then
      allocate(RV_stat(nRV),RV_sim(nRV))
      RV_stat=0
      RV_sim=zero
    end if

    cntT0=0
    nTs=maxval(nT0)
    if(nTs.gt.0)then
      allocate(T0_stat(nTs,NB),T0_sim(nTs,NB))
      T0_stat=0
      T0_sim=zero
    end if

    dt1=tstart-tepoch
    dt2=dt1+tint
    if(dt1.lt.zero)then
      call ode_a(m,R,ra1,dt1,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
      if(.not.Hc)then
        if(abs(dt1).le.tint)then
          call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
              &cntT0,T0_stat,T0_sim,Hc)
        end if
      end if
    else
      call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
    end if

    if(Hc)then
!       resw=sqrt(resmax)/real(ndata,dp)
        resw=set_max_residuals(ndata)
    else
      !call setval(RVobs,RV_sim,T0obs,T0_sim,resw)
      !call setval_2(RVobs,RV_sim,eRVobs,T0obs,T0_sim,eT0obs,resw)
      !call setval_3(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
      call setval_4(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
      call check_max_residuals(resw,ndata)
    end if

    if(allocated(RV_sim)) deallocate(RV_stat,RV_sim)
    if(allocated(gamma)) deallocate(gamma)
    if(allocated(T0_sim)) deallocate(T0_stat,T0_sim)
    if(allocated(m))   deallocate(m,R,P,a,e,w,mA,i,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_lm
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! subroutines called by the bootstrap to calculate the new set of T_0 and RV from 
  ! the fitted parameters
  subroutine ode_parok(allpar,par,RVok,T0ok,resw,fitness)
    real(dp),dimension(:),intent(in)::allpar,par
    real(dp),dimension(:),intent(out)::RVok
    real(dp),dimension(:,:),intent(out)::T0ok
    real(dp),dimension(:),intent(out)::resw
    real(dp),intent(out)::fitness

    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,i,lN
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp)::dt1,dt2
    logical::Hc
    !real(dp),parameter::resmax=huge(0._dp)
    !real(dp),parameter::resmax=1.e6_dp
    integer::cntRV
    real(dp),dimension(:),allocatable::RV_sim
    integer,dimension(:),allocatable::RV_stat
    real(dp),dimension(:,:),allocatable::gamma
    integer::cntT0,nTs
    real(dp),dimension(:,:),allocatable::T0_sim
    integer,dimension(:,:),allocatable::T0_stat
    logical::checkpar

    resw=zero
    fitness=zero
    nTs=maxval(nT0)

    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),&
        &w(NB),mA(NB),i(NB),lN(NB),clN(NB))

    ! from the allpar and par to keplerian elements: due to the LMdif approach/code
!     call par2kel(allpar,par,m,R,P,a,e,w,mA,i,lN)
    if(progtype.le.1)then
      call par2kel(allpar,par,m,R,P,a,e,w,mA,i,lN)
      checkpar = checkbounds(par)
    else
      call par2kel_fit(allpar,par,m,R,P,a,e,w,mA,i,lN,checkpar)
      if(.not.checkpar)then
        if(nRV.gt.0) RVok=zero
        if(nTs.gt.0) T0ok=zero
  !       resw=sqrt(resmax)/real(ndata,dp)
        resw=set_max_residuals(ndata)
        fitness=resmax
        return
      end if
      checkpar = checkbounds_fit(par)
    end if

    if(.not.checkpar)then
      if(nRV.gt.0) RVok=zero
      if(nTs.gt.0) T0ok=zero
!       resw=sqrt(resmax)/real(ndata,dp)
      resw=set_max_residuals(ndata)
      fitness=resmax
      return
    end if

    
    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))
    ! IT CREATES THE INITIAL STATE VECTOR FROM KEPLERIAN ORBITAL ELEMENS IN THE ORBITAL PLANE
    call initial_state(P,a,e,mA,ra0)
    ! IT ROTATES THE STATE VECTORS FROM ORBITAL REF. SYSTEM TO THE OBSERVER REF. SYSTEM
    call orb2obs(ra0,-lN,-i,-w,ra1)

    cntRV=0
    if(nRV.gt.0)then
      allocate(RV_stat(nRV),RV_sim(nRV))
      RV_stat=0
      RV_sim=0._dp
    end if

    cntT0=0
    if(nTs.gt.0)then
      allocate(T0_stat(nTs,NB),T0_sim(nTs,NB))
      T0_stat=0
      T0_sim=0._dp
    end if

    dt1=tstart-tepoch
    dt2=dt1+tint
    if(dt1.lt.0._dp)then
      call ode_a(m,R,ra1,dt1,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
      if(.not.Hc)then
        if(abs(dt1).le.tint)then
          dt2=dt1+tint
          call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
              &cntT0,T0_stat,T0_sim,Hc)
        end if
      end if
    else
      dt2=dt1+tint
      call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
    end if
    if(Hc)then
      if(nRV.gt.0) RVok = zero
      if(nTs.gt.0) T0ok = zero
!       resw=sqrt(resmax)/real(ndata,dp)
      resw=set_max_residuals(ndata)
      fitness=resmax
    else
      !call setval(RVobs,RV_sim,T0obs,T0_sim,resw)
!       call setval_3(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
      call setval_4(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
      call check_max_residuals(resw,ndata)
      if(nRV.gt.0)then
        call setRVok(RV_sim,gamma,RVok)
      end if
      if(nTs.gt.0) T0ok=T0_sim
      fitness=sum(resw*resw)
    end if

    if(allocated(RV_sim)) deallocate(RV_stat,RV_sim)
    if(allocated(gamma)) deallocate(gamma)
    if(allocated(T0_sim)) deallocate(T0_stat,T0_sim)
    if(allocated(m)) deallocate(m,R,P,a,e,w,mA,i,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_parok

  ! subroutine called by the L-M to fit the parameters with observations in input (RV and T0)
  subroutine ode_boot(allpar,nm,nn,par,RV_obs,T0_obs,resw,iflag)
    integer,intent(in)::nm,nn
    real(dp),dimension(:),intent(in)::allpar,par,RV_obs
    real(dp),dimension(:,:),intent(in)::T0_obs
    real(dp),dimension(:),intent(out)::resw
    integer,intent(inout)::iflag

    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,i,lN
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp)::dt1,dt2
    logical::Hc
    !real(dp),parameter::resmax=huge(0._dp)
    !real(dp),parameter::resmax=1.e6_dp
    integer::cntRV
    real(dp),dimension(:),allocatable::RV_sim
    integer,dimension(:),allocatable::RV_stat
    real(dp),dimension(:,:),allocatable::gamma
    integer::cntT0,nTs
    real(dp),dimension(:,:),allocatable::T0_sim
    integer,dimension(:,:),allocatable::T0_stat
    logical::checkpar

    resw=zero
    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),&
        &w(NB),mA(NB),i(NB),lN(NB),clN(NB))

    ! from the allpar and par to keplerian elements: due to the LMdif approach/code
!     call par2kel(allpar,par,m,R,P,a,e,w,mA,i,lN)
    if(progtype.le.1)then
      call par2kel(allpar,par,m,R,P,a,e,w,mA,i,lN)
      checkpar = checkbounds(par)
    else
      call par2kel_fit(allpar,par,m,R,P,a,e,w,mA,i,lN,checkpar)
      if(.not.checkpar)then
  !       resw=sqrt(resmax)/real(ndata,dp)
          resw=set_max_residuals(ndata)
        return
      end if
      checkpar = checkbounds_fit(par)
    end if
    if(.not.checkpar)then
!       resw=sqrt(resmax)/real(ndata,dp)
        resw=set_max_residuals(ndata)
      return
    end if
      
    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))
    ! IT CREATES THE INITIAL STATE VECTOR FROM KEPLERIAN ORBITAL ELEMENS IN THE ORBITAL PLANE
    call initial_state(P,a,e,mA,ra0)
    ! IT ROTATES THE STATE VECTORS FROM ORBITAL REF. SYSTEM TO THE OBSERVER REF. SYSTEM
    call orb2obs(ra0,-lN,-i,-w,ra1)

    cntRV=0
    if(nRV.gt.0)then
      allocate(RV_stat(nRV),RV_sim(nRV))
      RV_stat=0
      RV_sim=zero
    end if

    cntT0=0
    nTs=maxval(nT0)
    if(nTs.gt.0)then
      allocate(T0_stat(nTs,NB),T0_sim(nTs,NB))
      T0_stat=0
      T0_sim=zero
    end if

    dt1=tstart-tepoch
    if(dt1.lt.0._dp)then
      call ode_a(m,R,ra1,dt1,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
      if(.not.Hc)then
        if(abs(dt1).le.tint)then
          dt2=dt1+tint
          call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
              &cntT0,T0_stat,T0_sim,Hc)
        end if
      end if
    else
      dt2=dt1+tint
      call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
    end if
    if(Hc)then
!       resw=sqrt(resmax)/real(ndata,dp)
      resw=set_max_residuals(ndata)
    else
      !call setval(RV_obs,RV_sim,T0_obs,T0_sim,resw)
!       call setval_3(RV_obs,RV_sim,eRVobs,gamma,T0_obs,T0_sim,eT0obs,resw)
      call setval_4(RV_obs,RV_sim,eRVobs,gamma,T0_obs,T0_sim,eT0obs,resw)
      call check_max_residuals(resw,ndata)
    end if

    if(allocated(RV_sim)) deallocate(RV_stat,RV_sim)
    if(allocated(gamma)) deallocate(gamma)
    if(allocated(T0_sim)) deallocate(T0_stat,T0_sim)
    if(allocated(m)) deallocate(m,R,P,a,e,w,mA,i,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_boot
  ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
  ! subroutine to write file and to screen the data ... what it writes
  ! depends on the option marked with 1 in arg.in file
  subroutine ode_out(cpuid,isim,wrtid,allpar,par,resw)
    integer,intent(in)::cpuid,isim,wrtid
    real(dp),dimension(:),intent(in)::allpar,par
    real(dp),dimension(:),intent(out)::resw

    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,i,lN
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp)::dt1,dt2
    logical::Hc
    !real(dp),parameter::resmax=huge(0._dp)
    !real(dp),parameter::resmax=1.e6_dp
    integer::cntRV
    real(dp),dimension(:),allocatable::RV_sim
    integer,dimension(:),allocatable::RV_stat
    real(dp),dimension(:,:),allocatable::gamma
    integer::cntT0,nTs
    real(dp),dimension(:,:),allocatable::T0_sim
    integer,dimension(:,:),allocatable::T0_stat
    logical::checkpar

    real(dp)::chi2r_RV,chi2r_T0,chi2wr_RV,chi2wr_T0,fitness,w_chi2r
    real(dp),dimension(:),allocatable::resw_temp
    
    ! units and file names to store and write to files
    integer::uorb,ucon
    character(512)::florb,fmorb,flcon,fmcon,fmele
    integer,dimension(:),allocatable::uele,utra
    character(512),dimension(:),allocatable::flele,fltra

    write(*,*)
    write(*,'(a)')" EXECUTING SIMPLE INTEGRATION AND WRITING FINAL FILES"
    write(*,'(a,i3)')" LM wrtid = ",wrtid
    write(*,*)
    
    resw=zero

    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),w(NB),mA(NB),i(NB),lN(NB),clN(NB))

    ! from the allpar and par to keplerian elements: due to the LMdif approach/code
    if(progtype.le.1)then
      call par2kel(allpar,par,m,R,P,a,e,w,mA,i,lN)
      checkpar = checkbounds(par)
    else
      call par2kel_fit(allpar,par,m,R,P,a,e,w,mA,i,lN,checkpar)
      if(.not.checkpar)then
        write(*,*)' checkpar is FALSE ... BAD'
  !       resw=sqrt(resmax)/real(ndata,dp)
        resw=set_max_residuals(ndata)
        return
      end if
      checkpar = checkbounds_fit(par)
    end if

    if(.not.checkpar)then
      write(*,*)' checkpar is FALSE ... BAD'
!       resw=sqrt(resmax)/real(ndata,dp)
      resw=set_max_residuals(ndata)
      return
    end if

    write(*,*)
    write(*,*)" par"
    write(*,*)par
    write(*,*)"  m = ",m
    write(*,*)"  P = ",P
    write(*,*)"  a = ",a
    write(*,*)"  e = ",e
    write(*,*)"  w = ",w
    write(*,*)" mA = ",mA
    write(*,*)"  i = ",i
    write(*,*)" lN = ",lN
    write(*,*)
    
    ! write orbital elements into a file
    call outElements(isim,wrtid,m,R,P,a,e,w,mA,i,lN)

    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))
    ! IT CREATES THE INITIAL STATE VECTOR FROM KEPLERIAN ORBITAL ELEMENS IN THE ORBITAL PLANE
    call initial_state(P,a,e,mA,ra0)
    ! IT ROTATES THE STATE VECTORS FROM ORBITAL REF. SYSTEM TO THE OBSERVER REF. SYSTEM
    call orb2obs(ra0,-lN,-i,-w,ra1) !OK 3T 1T 3T

    cntRV=0
    if(nRV.gt.0)then
      allocate(RV_stat(nRV),RV_sim(nRV))
      RV_stat=0
      RV_sim=zero
    end if

    if((idtra.ge.1).and.(idtra.le.NB)) call set_file_tra(cpuid,isim,wrtid,utra,fltra)
    cntT0=0
    nTs=maxval(nT0)
    if(nTs.gt.0)then
      allocate(T0_stat(nTs,NB),T0_sim(nTs,NB))
      T0_stat=0
      T0_sim=zero
    end if

    fmorb=trim(adjustl(fmtorbit()))
    if(wrtorb.eq.1) call set_file_orb(cpuid,isim,wrtid,uorb,florb)
    fmcon=trim(adjustl(fmtconst()))
    if(wrtconst.eq.1) call set_file_con(cpuid,isim,wrtid,ucon,flcon)
    fmele=fmtele()
    if(wrtel.eq.1) call set_file_elem(cpuid,isim,wrtid,uele,flele)

    dt1=tstart-tepoch
    dt2=dt1+tint
    if(dt1.lt.0._dp)then
!       write(*,'(a,g25.14,a)')" dt1 = ",dt1,&
!           &"< 0 => backward integration"
      call ode_a_o(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
          &m,R,ra1,dt1,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
      if(.not.Hc)then
        if(abs(dt1).le.tint)then
!           write(*,'(a,g25.14,a)')" dt2 =  ",dt2,&
!               &" => forward integration"
          call ode_a_o(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
              &m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
              &cntT0,T0_stat,T0_sim,Hc)
        end if
      end if
    else
!       write(*,'(a,g25.14,a)')" dt2 = ",dt2,&
!           &" => only forward integration"
      call ode_a_o(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
          &m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
    end if

    if((idtra.ge.1).and.(idtra.le.NB)) call close_tra(utra,fltra)
    if(wrtorb.eq.1) close(uorb)
    if(wrtconst.eq.1) close(ucon)
    if(wrtel.eq.1) call close_elem(uele,flele)

    if(Hc)then
!       resw=sqrt(resmax)/real(ndata,dp)
      resw=set_max_residuals(ndata)
    else
      !call setval(RVobs,RV_sim,T0obs,T0_sim,resw)
      !call setval_2(RVobs,RV_sim,eRVobs,T0obs,T0_sim,eT0obs,resw)
!       call setval_3(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
      call setval_4(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
      call check_max_residuals(resw,ndata)
    end if

    ! Reduced Chi Squares for report/summary
    chi2r_RV=zero
    chi2r_T0=zero
    chi2wr_RV=zero
    chi2wr_T0=zero
    
    if(cntRV.gt.0)then
      write(*,*)""
      write(*,'(a,i4)')" RADIAL VELOCITIES found: ",cntRV
      write(*,*)""
      call write_RV(gamma,RV_sim,RV_stat)
      call write_RV(cpuid,isim,wrtid,gamma,RV_sim,RV_stat)
      allocate(resw_temp(nRV))
      call set_RV_resw(RVobs,RV_sim,eRVobs,gamma,resw_temp)
      chi2r_RV=sum(resw_temp*resw_temp)*inv_dof
      w_chi2r=real(ndata,dp)/real(nRV,dp)
      chi2wr_RV=chi2r_RV*w_chi2r
      deallocate(resw_temp)
    end if

    if(cntT0.gt.0)then
      write(*,*)
      write(*,'(a,i5)')" T0 SIM found ",cntT0
      write(*,*)
      call write_T0(T0_sim,T0_stat)
      call write_T0(cpuid,isim,wrtid,T0_sim,T0_stat)
      allocate(resw_temp(sum(nT0)))
      call set_T0_resw(T0obs,T0_sim,eT0obs,resw_temp)
!       resw_temp=(T0obs-T0_sim)/eT0obs
      chi2r_T0=sum(resw_temp*resw_temp)*inv_dof
      w_chi2r=real(ndata,dp)/real(sum(nT0),dp)
      chi2wr_T0=chi2r_T0*w_chi2r
      deallocate(resw_temp)
    end if

    fitness = sum(resw*resw)
    call write_fitness_summary(cpuid,isim,wrtid,chi2r_RV,chi2wr_RV,chi2r_T0,chi2wr_T0,fitness)

    if(allocated(RV_sim)) deallocate(RV_stat,RV_sim)
    if(allocated(gamma)) deallocate(gamma)
    if(allocated(T0_sim)) deallocate(T0_stat,T0_sim)
    if(allocated(m))   deallocate(m,R,P,a,e,w,mA,i,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_out
  ! ------------------------------------------------------------------ !

end module ode_run



