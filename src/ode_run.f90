module ode_run
  use constants
  use parameters
  use parameters_conversion
  use celestial_mechanics
  use rotations,only:orb2obs
  use eq_motion,only:eqmastro
  use numerical_integrator,only:int_rk_a
  use transits
  use radial_velocities
  use gls_module,only:check_periodogram,check_periodogram_scale,check_and_write_periodogram
  use output_files
  use statistics,only:mean
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

  ! as setval but it accepts different arguments, i.e., gamma of RV, but for multiple set of RV
  ! and it allow to fit T0res or OCres.
  subroutine setval_5(RV_obs,RV_sim,e_RVobs,gamma,epoT0_obs,T0_obs,eT0_obs,T0_sim,resw,ocfit)
    real(dp),dimension(:),intent(in)::RV_obs,RV_sim,e_RVobs
    real(dp),dimension(:,:),allocatable,intent(out)::gamma
    integer,dimension(:,:),intent(in)::epoT0_obs
    real(dp),dimension(:,:),intent(in)::T0_obs,eT0_obs,T0_sim
    real(dp),dimension(:),intent(out)::resw
    integer,intent(in)::ocfit
    
    real(dp),dimension(:),allocatable::val,val_T0,val_oc

    allocate(val(ndata))
    resw=zero
    val=zero
    if(nRV.gt.0) call set_RV_resw(RV_obs,RV_sim,e_RVobs,gamma,val(1:nRV))
    if(ocfit.eq.1)then
      allocate(val_oc(sum(nT0)))
      val_oc=zero
      call set_oc_resw(epoT0_obs,T0_obs,eT0_obs,T0_sim,val_oc)
      val(nRV+1:ndata)=val_oc
      deallocate(val_oc)
    else if(ocfit.eq.2)then
      allocate(val_T0(sum(nT0)),val_oc(sum(nT0)))
      val_T0=zero
      call set_T0_resw(T0_obs,T0_sim,eT0_obs,val_T0)
      val_oc=zero
      call set_oc_resw(epoT0_obs,T0_obs,eT0_obs,T0_sim,val_oc)
      val(nRV+1:ndata)=sqrt_half*sqrt(val_T0*val_T0+val_oc*val_oc)
      deallocate(val_T0,val_oc)
    else
      allocate(val_T0(sum(nT0)))
      val_T0=zero
      call set_T0_resw(T0_obs,T0_sim,eT0_obs,val_T0)
      val(nRV+1:ndata)=val_T0
      deallocate(val_T0)
    end if
    call set_fitness(val,resw)
    deallocate(val)
    
    return
  end subroutine setval_5

  ! set only RV to the weighted residuals
  subroutine set_RV_resw(RV_obs,RV_sim,e_RVobs,gamma,resw)
    use statistics,only:wmean
    real(dp),dimension(:),intent(in)::RV_obs,RV_sim,e_RVobs
    real(dp),dimension(:,:),allocatable,intent(out)::gamma
    real(dp),dimension(:),intent(inout)::resw
    real(dp),dimension(:),allocatable::dRV,wi
    integer::j,aRV,bRV

    if(.not.allocated(gamma)) allocate(gamma(nRVset,2))
    gamma=zero
!     resw=zero
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
    
    real(dp)::obs,sim,eobs,weight_fit
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
!         resw(a)=resw(a)+weight_fit*((obs-sim)/eobs)
        resw(a)=((obs-sim)/eobs)
      end do
    end do
    
    return
  end subroutine set_T0_resw
  
  ! compute linear ephemeris for both obs and sim
  ! use O-Cobs/sim to compute residuals
  subroutine set_oc_resw(epoT0_obs,T0_obs,eT0_obs,T0_sim,resw)
    integer,dimension(:,:),intent(in)::epoT0_obs
    real(dp),dimension(:,:),intent(in)::T0_obs,eT0_obs,T0_sim
    real(dp),dimension(:),intent(out)::resw
    
    integer,dimension(:),allocatable::x
    real(dp),dimension(:),allocatable::y,ey
    real(dp)::Tref_sim,eTref_sim,Pref_sim,ePref_sim
    
    real(dp)::obs,sim,eobs
    
    integer::j,j1,a
    
    ! NB,nT0 are common variables
    resw=zero
    a=0
    do j=2,NB
      if(nT0(j).gt.0)then
        
        allocate(x(nT0(j)),y(nT0(j)),ey(nT0(j)))
        x=epoT0_obs(1:nT0(j),j)
        y=T0_sim(1:nT0(j),j)
        ey=mean(eT0_obs(1:nT0(j),j))
        Pref_sim=zero
        ePref_sim=zero
        Tref_sim=zero
        eTref_sim=zero
        call linfit(x,y,ey,Pref_sim,ePref_sim,Tref_sim,eTref_sim)
        deallocate(x,y,ey)
        
        do j1=1,nT0(j)
          a=a+1
          obs=T0_obs(j1,j) - (Tephem(j)+ epoT0_obs(j1,j)*Pephem(j))
          sim=T0_sim(j1,j) - (Tref_sim + epoT0_obs(j1,j)*Pref_sim )
!           sim=T0_sim(j1,j) - (Tephem(j)+ epoT0_obs(j1,j)*Pephem(j))
          eobs=eT0_obs(j1,j)
!           resw(a)=resw(a)+weight_fit*((obs-sim)/eobs)
          resw(a)=(obs-sim)/eobs
        end do
        
      end if
      
    end do
    
    return
  end subroutine set_oc_resw

  
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

!     Hc=.true.
    
!     Hc=Hillcheck(m,rin)
    if(do_hill_check) Hc=mutual_Hill_check(m,rin)
!     if(Hc) return
    if(.not.Hc) return

    allocate(X(NB),Y(NB),Z(NB),cX(NB),cY(NB),cR(NB),rmean(NB))
    X=0
    Y=0
    Z=0
    do j=2,NB
      X(j)=1+(j-1)*6
      Y(j)=2+(j-1)*6
      Z(j)=3+(j-1)*6
    end do
    cX=one
    cY=one
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
      if(do_hill_check) Hc=mutual_Hill_check(m,r2)
!       if(Hc) return
      if(.not.Hc) return

      ! RV check
      if(cntRV.lt.nRV) call check_RV(m,r1,dr,itime,hok,cntRV,RV_stat,RV_sim)

      ! T0 check
      do j=2,NB
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))
        rmean(j)=half*( rsky(r1(X(j):Y(j))) + rsky(r2(X(j):Y(j))) )
      end do

      if((idtra.gt.0).and.(idtra.le.NB))then
        if(cntT0.lt.sum(nT0))then
          if(idtra.eq.1)then
            do j=2,NB
              if(rmean(j).le.cR(j))then
                if(clN(j).eq.0)then
                  
                  if((cX(j).le.zero).and.(r1(Z(j)).gt.zero))then
                    call check_T0(j,m,R,r1,r2,&
                      &itime,hok,T0_stat,T0_sim,Hc) ! WARNING: CHANGE Hc WITHIN check_T0
!                     if(Hc)return
                    if(.not.Hc) return
                  end if
                
                else
                  
                  if((cY(j).le.zero).and.(r1(Z(j)).gt.zero))then
                    call check_T0(j,m,R,r1,r2,&
                      &itime,hok,T0_stat,T0_sim,Hc) ! WARNING: CHANGE Hc WITHIN check_T0
!                     if(Hc)return
                    if(.not.Hc) return
                  end if
                
                end if
              end if
            end do
          else
            if(rmean(idtra).le.cR(idtra))then
              if(clN(idtra).eq.0)then
                
                if((cX(idtra).le.zero).and.(r1(Z(idtra)).gt.zero))then
                  call check_T0(idtra,m,R,r1,r2,&
                    &itime,hok,T0_stat,T0_sim,Hc)
!                   if(Hc)return
                  if(.not.Hc) return
                end if
                
              else
              
                if((cY(idtra).le.zero).and.(r1(Z(idtra)).gt.zero))then
                  call check_T0(idtra,m,R,r1,r2,&
                    &itime,hok,T0_stat,T0_sim,Hc)
!                   if(Hc)return
                  if(.not.Hc) return
                end if
                
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

!     Hc=.false.
!     Hc=.true.
    
!     Hc=Hillcheck(m,rin)
    if(do_hill_check) Hc=mutual_Hill_check(m,rin)
!     write(*,*)' mutual_Hill_check = ',Hc
!     if(Hc) return
    if(.not.Hc) return

    allocate(X(NB),Y(NB),Z(NB),cX(NB),cY(NB),cR(NB),rmean(NB))
    X=0
    Y=0
    Z=0
    do j=2,NB
      X(j)=1+(j-1)*6
      Y(j)=2+(j-1)*6
      Z(j)=3+(j-1)*6
    end do
    cX=one
    cY=one
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

!     write(*,'(a)',advance='NO')' mutual_Hill_check: '
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
      if(do_hill_check) Hc=mutual_Hill_check(m,r2)
!       if(Hc) write(*,'(a7,i6,a6,l)')'iter = ',j1,' Hc = ',Hc
!       if(Hc) return
      if(.not.Hc) return

      ! RV check
      if(cntRV.lt.nRV) call check_RV(m,r1,dr,itime,hok,cntRV,RV_stat,RV_sim)

      ! T0 check (to compare and all)
      do j=2,NB
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))
        rmean(j)=half*( rsky(r1(X(j):Y(j))) + rsky(r2(X(j):Y(j))) )
      end do
      if((idtra.gt.0).and.(idtra.le.NB))then
        !if(cntT0.lt.sum(nT0))then
        if(idtra.eq.1)then
          do j=2,NB
            if(rmean(j).le.cR(j))then
              if(clN(j).eq.0)then
              
                if((cX(j).le.zero).and.(r1(Z(j)).gt.zero))then
                  if(nT0(j).gt.0) call check_T0(j,m,R,r1,r2,&
                    &itime,hok,T0_stat,T0_sim,Hc)
!                   if(Hc)return
                  if(.not.Hc) return
                  jtra=jtra+1
                  call all_transits(jtra,j,m,R,r1,r2,&
                    &itime,hok,stat_tra,storetra)
                end if
              
              else
              
                if((cY(j).le.zero).and.(r1(Z(j)).gt.zero))then
                  if(nT0(j).gt.0)call check_T0(j,m,R,r1,r2,&
                    &itime,hok,T0_stat,T0_sim,Hc)
!                   if(Hc)return
                  if(.not.Hc) return
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
            
              if((cX(idtra).le.zero).and.&
                &(r1(Z(idtra)).gt.zero))then
                if(nT0(idtra).gt.0)call check_T0(idtra,m,R,r1,r2,&
                  &itime,hok,T0_stat,T0_sim,Hc)
!                 if(Hc)return
                if(.not.Hc) return
                jtra=jtra+1
                call all_transits(jtra,idtra,m,R,r1,r2,&
                  &itime,hok,stat_tra,storetra)
              end if
            
            else
            
              if((cY(idtra).le.zero).and.&
                &(r1(Z(idtra)).gt.zero))then
                if(nT0(idtra).gt.0)call check_T0(idtra,m,R,r1,r2,&
                  &itime,hok,T0_stat,T0_sim,Hc)
!                 if(Hc)return
                if(.not.Hc) return
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
      iftime=int(ftime+half)
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

    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
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
    logical::checkpar,gls_check
    
!     real(dp)::fit_scale,gls_scale

    Hc=.true.
    resw=zero
    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),w(NB),mA(NB),inc(NB),lN(NB),clN(NB))

!     fit_scale=one
!     gls_scale=one
    call convert_parameters(allpar,par,m,R,P,a,e,w,mA,inc,lN,checkpar)
    if(.not.checkpar)then
        resw=set_max_residuals(ndata)
      return
    end if

!     call convert_parameters_scale(allpar,par,m,R,P,a,e,w,mA,inc,lN,fit_scale)

    
    ! it is needed to define which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))
    ! IT CREATES THE INITIAL STATE VECTOR FROM KEPLERIAN ORBITAL ELEMENS IN THE ORBITAL PLANE
    call initial_state(P,a,e,mA,ra0)
    ! IT ROTATES THE STATE VECTORS FROM ORBITAL REF. SYSTEM TO THE OBSERVER REF. SYSTEM
    call orb2obs(ra0,-lN,-inc,-w,ra1)

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
      if(Hc)then
        if(abs(dt1).le.tint)then
          call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
            &cntT0,T0_stat,T0_sim,Hc)
        end if
      end if
    
    else
    
      call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
        &cntT0,T0_stat,T0_sim,Hc)
    
    end if

    gls_check=.true.
    if(.not.Hc) then
!       resw=sqrt(resmax)/real(ndata,dp)
        resw=set_max_residuals(ndata)
    else
      !call setval(RVobs,RV_sim,T0obs,T0_sim,resw)
      !call setval_2(RVobs,RV_sim,eRVobs,T0obs,T0_sim,eT0obs,resw)
      !call setval_3(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
!       call setval_4(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
      call setval_5(RVobs,RV_sim,eRVobs,gamma,epoT0obs,T0obs,eT0obs,T0_sim,resw,oc_fit)
      if(cntRV.gt.0) call check_periodogram(jdRV,RVobs-RV_sim,eRVobs,P,gls_check)
!       if(cntRV.gt.0) call check_periodogram_scale(jdRV,RVobs-RV_sim,eRVobs,P,gls_check,gls_scale)
      if(.not.gls_check)then
        resw=set_max_residuals(ndata)
      else
        call check_max_residuals(resw,ndata)
      end if
!       resw=resw*fit_scale*gls_scale
      call check_max_residuals(resw,ndata)
    end if

    if(allocated(RV_sim)) deallocate(RV_stat,RV_sim)
    if(allocated(gamma)) deallocate(gamma)
    if(allocated(T0_sim)) deallocate(T0_stat,T0_sim)
    if(allocated(m))   deallocate(m,R,P,a,e,w,mA,inc,lN,clN)
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

    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
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
    logical::checkpar,gls_check
    
!     real(dp)::fit_scale,gls_scale
    
    integer::ii
    
    Hc=.true.
    resw=zero
    fitness=zero
    nTs=maxval(nT0)

    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),&
        &w(NB),mA(NB),inc(NB),lN(NB),clN(NB))

!     fit_scale=one
!     gls_scale=one
    call convert_parameters(allpar,par,m,R,P,a,e,w,mA,inc,lN,checkpar)
    if(.not.checkpar)then
      if(nRV.gt.0) RVok=zero
      if(nTs.gt.0) T0ok=zero
      resw=set_max_residuals(ndata)
      fitness=resmax
      return
    end if
!     call convert_parameters_scale(allpar,par,m,R,P,a,e,w,mA,inc,lN,fit_scale)
    
    do ii=1,NB
      if(m(ii).lt.zero)then
        write(*,'(a,i3,a,10000(1x,f22.16))')' neg mass(ii=',ii,' ) : masses = ',m
        stop
      end if
    end do
    
    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))
    ! IT CREATES THE INITIAL STATE VECTOR FROM KEPLERIAN ORBITAL ELEMENS IN THE ORBITAL PLANE
    call initial_state(P,a,e,mA,ra0)
    ! IT ROTATES THE STATE VECTORS FROM ORBITAL REF. SYSTEM TO THE OBSERVER REF. SYSTEM
    call orb2obs(ra0,-lN,-inc,-w,ra1)

    cntRV=0
    if(nRV.gt.0)then
      allocate(RV_stat(nRV),RV_sim(nRV))
      RV_stat=0
      RV_sim=zero
    end if

    cntT0=0
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
      if(Hc)then
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
    if(.not.Hc) then
      if(nRV.gt.0) RVok = zero
      if(nTs.gt.0) T0ok = zero
      resw=set_max_residuals(ndata)
      fitness=resmax
    else
      !call setval(RVobs,RV_sim,T0obs,T0_sim,resw)
!       call setval_3(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
!       call setval_4(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
      call setval_5(RVobs,RV_sim,eRVobs,gamma,epoT0obs,T0obs,eT0obs,T0_sim,resw,oc_fit)
      if(cntRV.gt.0) call check_periodogram(jdRV,RVobs-RV_sim,eRVobs,P,gls_check)
!       if(cntRV.gt.0) call check_periodogram_scale(jdRV,RVobs-RV_sim,eRVobs,P,gls_check,gls_scale)
!       resw=resw*fit_scale*gls_scale
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
    if(allocated(m)) deallocate(m,R,P,a,e,w,mA,inc,lN,clN)
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

    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
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
    logical::checkpar,gls_check
    
!     real(dp)::fit_scale,gls_scale
    
    Hc=.true.
    resw=zero
    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),&
        &w(NB),mA(NB),inc(NB),lN(NB),clN(NB))

!     fit_scale=one
!     gls_scale=one
    call convert_parameters(allpar,par,m,R,P,a,e,w,mA,inc,lN,checkpar)
    if(.not.checkpar)then
!       resw=sqrt(resmax)/real(ndata,dp)
        resw=set_max_residuals(ndata)
      return
    end if
!     call convert_parameters_scale(allpar,par,m,R,P,a,e,w,mA,inc,lN,fit_scale)

      
    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))
    ! IT CREATES THE INITIAL STATE VECTOR FROM KEPLERIAN ORBITAL ELEMENS IN THE ORBITAL PLANE
    call initial_state(P,a,e,mA,ra0)
    ! IT ROTATES THE STATE VECTORS FROM ORBITAL REF. SYSTEM TO THE OBSERVER REF. SYSTEM
    call orb2obs(ra0,-lN,-inc,-w,ra1)

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
    if(dt1.lt.zero)then
      call ode_a(m,R,ra1,dt1,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
      if(Hc)then
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
    if(.not.Hc) then
!       resw=sqrt(resmax)/real(ndata,dp)
      resw=set_max_residuals(ndata)
    else
      !call setval(RV_obs,RV_sim,T0_obs,T0_sim,resw)
!       call setval_3(RV_obs,RV_sim,eRVobs,gamma,T0_obs,T0_sim,eT0obs,resw)
!       call setval_4(RV_obs,RV_sim,eRVobs,gamma,T0_obs,T0_sim,eT0obs,resw)
      call setval_5(RV_obs,RV_sim,eRVobs,gamma,epoT0obs,T0_obs,eT0obs,T0_sim,resw,oc_fit)
      if(cntRV.gt.0) call check_periodogram(jdRV,RVobs-RV_sim,eRVobs,P,gls_check)
!       if(cntRV.gt.0) call check_periodogram_scale(jdRV,RVobs-RV_sim,eRVobs,P,gls_check,gls_scale)
!       resw=resw*fit_scale*gls_scale
      call check_max_residuals(resw,ndata)
    end if

    if(allocated(RV_sim)) deallocate(RV_stat,RV_sim)
    if(allocated(gamma)) deallocate(gamma)
    if(allocated(T0_sim)) deallocate(T0_stat,T0_sim)
    if(allocated(m)) deallocate(m,R,P,a,e,w,mA,inc,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_boot
  ! ------------------------------------------------------------------ !

  
  subroutine ode_full_args(m,R,rin,time,step_in,clN,&
    &cntRV,tRV,RV_stat,RV_sim,id_transit_body,transit_flag,dur_check,&
    &cntT0,n_T0,T0_num,T0_stat,T0_sim,Hc)
    real(dp),dimension(:),intent(in)::m,R,rin
    real(dp),intent(in)::time,step_in
    integer,dimension(:),intent(in)::clN
    integer,intent(inout)::cntRV
    real(dp),dimension(:),intent(in)::tRV
    integer,dimension(:),intent(inout)::RV_stat
    real(dp),dimension(:),intent(inout)::RV_sim
    integer,intent(in)::id_transit_body
    logical,dimension(:),intent(in)::transit_flag
    integer,intent(in)::dur_check
    integer,intent(inout)::cntT0
    integer,dimension(:),intent(in)::n_T0
    integer,dimension(:,:),intent(in)::T0_num
    integer,dimension(:,:),intent(inout)::T0_stat
    real(dp),dimension(:,:),intent(inout)::T0_sim
    logical,intent(inout)::Hc

    integer::n_body,nb_dim
    real(dp),dimension(:),allocatable::dr,r1,r2,err
    integer,dimension(:),allocatable::X,Y,Z
    real(dp),dimension(:),allocatable::cX,cY,cR,rmean
    real(dp)::hw,hok,hnext,itime
    integer::j,j1 !,nj
    integer::n_rv

    
    if(do_hill_check) Hc=mutual_Hill_check(m,rin)
    if(.not.Hc) return
    
    n_body=size(m)
    allocate(X(n_body),Y(n_body),Z(n_body),cX(n_body),cY(n_body),&
      &cR(n_body),rmean(n_body))
    X=0
    Y=0
    Z=0
    do j=2,n_body
      X(j)=1+(j-1)*6
      Y(j)=2+(j-1)*6
      Z(j)=3+(j-1)*6
    end do
    cX=one
    cY=one
    cR=zero
    cR(2:n_body)=1.5_dp*(R(1)+R(2:n_body))*RsunAU ! RsunAU from constants module
    rmean=9.e3_dp

    nb_dim=6*n_body
    allocate(dr(nb_dim),r1(nb_dim),r2(nb_dim),err(nb_dim))
    hw=step_in
    if(time.lt.zero) hw=-hw
    itime=zero
    r1=rin
    r2=zero
    err=zero

    n_rv=size(tRV)
    
    j1=0
    integration: do
      j1=j1+1
      if(abs(itime+hw).gt.abs(time)) hw=time-itime
      call eqmastro(m,r1,dr)
      call int_rk_a(m,r1,dr,hw,hok,hnext,r2,err)
      
      if(do_hill_check) Hc=mutual_Hill_check(m,r2)
      if(.not.Hc) return

      ! RV check
      if(cntRV.lt.n_rv) call check_RV&
        &(m,r1,dr,itime,hok,cntRV,tRV,RV_stat,RV_sim)
      
      ! T0 check
      do j=2,n_body
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))
        rmean(j)=half*( rsky(r1(X(j):Y(j))) + rsky(r2(X(j):Y(j))) )
      end do
      
      if((id_transit_body.gt.0).and.(id_transit_body.le.n_body))then
        if(cntT0.lt.sum(n_T0))then
          if(id_transit_body.eq.1)then
            do j=2,n_body
              if(rmean(j).le.cR(j))then
                if(clN(j).eq.0)then
                  
                  if((cX(j).le.zero).and.(r1(Z(j)).gt.zero))then
                    call check_T0(j,m,R,r1,r2,&
                      &itime,hok,transit_flag,dur_check,n_T0,&
                      &T0_num,T0_stat,T0_sim,Hc)
                    if(.not.Hc) return
                  end if
                
                else
                  
                  if((cY(j).le.zero).and.(r1(Z(j)).gt.zero))then
                    call check_T0(j,m,R,r1,r2,&
                      &itime,hok,transit_flag,dur_check,n_T0,&
                      &T0_num,T0_stat,T0_sim,Hc)
                    if(.not.Hc) return
                  end if
                
                end if
              end if
            end do
          else
            if(rmean(id_transit_body).le.cR(id_transit_body))then
              if(clN(id_transit_body).eq.0)then
                
                if((cX(id_transit_body).le.zero).and.&
                  &(r1(Z(id_transit_body)).gt.zero))then
                  call check_T0(j,m,R,r1,r2,&
                    &itime,hok,transit_flag,dur_check,n_T0,&
                    &T0_num,T0_stat,T0_sim,Hc)
                  if(.not.Hc) return
                end if
                
              else
              
                if((cY(id_transit_body).le.zero).and.&
                  &(r1(Z(id_transit_body)).gt.zero))then
                  call check_T0(j,m,R,r1,r2,&
                    &itime,hok,transit_flag,dur_check,n_T0,&
                    &T0_num,T0_stat,T0_sim,Hc)
                  if(.not.Hc) return
                end if
                
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
  end subroutine ode_full_args
  
  
  subroutine orbits_to_data(t_start,t_epoch,step_in,t_int,&
    &m,R,P,ecc,argp,mA,inc,lN,&
    &tRV,RV_sim,&
    &id_transit_body,transit_flag,dur_check,n_T0,T0_num,T0_sim)
    real(dp),intent(in)::t_start,t_epoch,step_in,t_int
    real(dp),dimension(:),intent(in)::m,R,P,ecc,argp,mA,inc,lN
    real(dp),dimension(:),intent(in)::tRV
    real(dp),dimension(:),allocatable,intent(out)::RV_sim
    integer,intent(in)::id_transit_body
    logical,dimension(:),intent(in)::transit_flag
    integer,intent(in)::dur_check
    integer,dimension(:),intent(in)::n_T0
    integer,dimension(:,:),intent(in)::T0_num
    real(dp),dimension(:,:),allocatable,intent(out)::T0_sim
    
    integer::n_body,nb_dim
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp),dimension(:),allocatable::sma
    real(dp)::dt1,dt2
    logical::Hc
    integer::cntRV,n_rv
    integer,dimension(:),allocatable::RV_stat
    integer::cntT0,nTs
    integer,dimension(:,:),allocatable::T0_stat
    logical::checkpar,gls_check
    
    integer::ii
    
    Hc=.true.
    
    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)
    
    n_body=size(m)
    allocate(sma(n_body))
    sma=zero
    call semax_vec(m(1),m(2:n_body),P(2:n_body), sma(2:n_body))
    
    nb_dim=6*n_body
    if(.not.allocated(ra0)) allocate(ra0(nb_dim),ra1(nb_dim))
    ! IT CREATES THE INITIAL STATE VECTOR FROM KEPLERIAN ORBITAL ELEMENS IN THE ORBITAL PLANE
    call initial_state(P,sma,ecc,mA,ra0)
    ! IT ROTATES THE STATE VECTORS FROM ORBITAL REF. SYSTEM TO THE OBSERVER REF. SYSTEM
    call orb2obs(ra0,-lN,-inc,-argp,ra1)

    cntRV=0
    n_rv=size(tRV)
    if(n_rv.gt.0)then
      allocate(RV_stat(n_rv),RV_sim(n_rv))
      RV_stat=0
      RV_sim=zero
    end if

    nTs=maxval(n_T0)
    cntT0=0
    if(nTs.gt.0)then
      allocate(T0_stat(nTs,n_body),T0_sim(nTs,n_body))
      T0_stat=0
      T0_sim=zero
    end if

    dt1=t_start-t_epoch
    dt2=dt1+t_int
    if(dt1.lt.zero)then
!       call ode_a(m,R,ra1,dt1,clN,cntRV,RV_stat,RV_sim,&
!           &cntT0,T0_stat,T0_sim,Hc)
      call ode_full_args(m,R,ra1,dt1,step_in,clN,&
        &cntRV,tRV,RV_stat,RV_sim,id_transit_body,&
        &transit_flag,dur_check,&
        &cntT0,n_T0,T0_num,T0_stat,T0_sim,Hc)
      if(Hc)then
        if(abs(dt1).le.t_int)then
          dt2=dt1+t_int
!           call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
!               &cntT0,T0_stat,T0_sim,Hc)
          call ode_full_args(m,R,ra1,dt2,step_in,clN,&
            &cntRV,tRV,RV_stat,RV_sim,id_transit_body,&
            &transit_flag,dur_check,&
            &cntT0,n_T0,T0_num,T0_stat,T0_sim,Hc)
        end if
      end if
    else
      dt2=dt1+t_int
!       call ode_a(m,R,ra1,dt2,clN,cntRV,RV_stat,RV_sim,&
!           &cntT0,T0_stat,T0_sim,Hc)
      call ode_full_args(m,R,ra1,dt2,step_in,clN,&
        &cntRV,tRV,RV_stat,RV_sim,id_transit_body,&
        &transit_flag,dur_check,&
        &cntT0,n_T0,T0_num,T0_stat,T0_sim,Hc)
    end if
    
    if(.not.Hc) then
      if(n_rv.gt.0) RV_sim = zero
      if(nTs.gt.0) T0_sim = zero
    end if

    if(allocated(RV_stat)) deallocate(RV_stat)
    if(allocated(T0_stat)) deallocate(T0_stat)
    if(allocated(ra0)) deallocate(ra0,ra1)

    
    return
  end subroutine orbits_to_data
  
  
    ! ------------------------------------------------------------------ !
  ! subroutine to write file and to screen the data ... what it writes
  ! depends on the option isim and wrtid/lmon in arg.in file
  subroutine ode_out(cpuid,isim,wrtid,allpar,par,resw,fit_scale,gls_scale)
    integer,intent(in)::cpuid,isim,wrtid
    real(dp),dimension(:),intent(in)::allpar,par
    real(dp),dimension(:),intent(out)::resw
    
    real(dp),optional,intent(out)::fit_scale,gls_scale

    
    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp)::dt1,dt2
    logical::Hc
    !real(dp),parameter::resmax=huge(zero)
    !real(dp),parameter::resmax=1.e6_dp
    integer::cntRV
    real(dp),dimension(:),allocatable::RV_sim
    integer,dimension(:),allocatable::RV_stat
    real(dp),dimension(:,:),allocatable::gamma
    integer::cntT0,nTs
    real(dp),dimension(:,:),allocatable::T0_sim
    integer,dimension(:,:),allocatable::T0_stat
    logical::checkpar,gls_check

    real(dp)::chi2r_RV,chi2r_T0,chi2wr_RV,chi2wr_T0,chi2r_oc,fitness,w_chi2r
    real(dp),dimension(:),allocatable::resw_temp
    
    ! units and file names to store and write to files
    integer::uorb,ucon
    character(512)::florb,fmorb,flcon,fmcon,fmele
    integer,dimension(:),allocatable::uele,utra
    character(512),dimension(:),allocatable::flele,fltra

    integer::i_par
        
    write(*,'(a)')''
    write(*,'(a)')" EXECUTING SIMPLE INTEGRATION AND WRITING FINAL FILES"
!     write(*,'(a,i3)')" LM on[1]/off[0] = ",wrtid
    write(*,'(a)')''
    
    resw=zero
    Hc=.true.
    
    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),w(NB),mA(NB),inc(NB),lN(NB),clN(NB))

    if(present(fit_scale))then
      fit_scale=one
      gls_scale=one
      call convert_parameters_scale(allpar,par,m,R,P,a,e,w,mA,inc,lN,fit_scale)
    else
      call convert_parameters(allpar,par,m,R,P,a,e,w,mA,inc,lN,checkpar,.true.)
      if(.not.checkpar)then
        write(*,*)' checkpar after convert_parameters is FALSE ... BAD'
  !       resw=sqrt(resmax)/real(ndata,dp)
        resw=set_max_residuals(ndata)
        return
      end if
    end if

    
!     write(*,*)
!     write(*,'(a12,4(1x,a25))')'parid','value','min','max','allpar(id)'
!     do i_par=1,nfit
!       write(*,'(a12,10(1x,f25.15))')parid(i_par),par(i_par),minpar(i_par),maxpar(i_par),allpar(idall(i_par))
!     end do
!     write(*,*)
!     write(*,'(a12,1000(f25.15,1x))')"  m = ",m
!     write(*,'(a12,1000(f25.15,1x))')"  P = ",P
!     write(*,'(a12,1000(f25.15,1x))')"  a = ",a
!     write(*,'(a12,1000(f25.15,1x))')"  e = ",e
!     write(*,'(a12,1000(f25.15,1x))')"  w = ",w
!     write(*,'(a12,1000(f25.15,1x))')" mA = ",mA
!     write(*,'(a12,1000(f25.15,1x))')"  i = ",inc
!     write(*,'(a12,1000(f25.15,1x))')" lN = ",lN
!     write(*,*)
!     write(*,'(a)')' CHECK FIT/KEP.ELEM. BOUNDARIES'
    if(present(fit_scale)) write(*,'(a,es23.16)')' fit_scale = ',fit_scale
!     write(*,*)
    flush(6)
    
    ! write orbital elements into a file
    call outElements(isim,wrtid,m,R,P,a,e,w,mA,inc,lN)

    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))
    ! IT CREATES THE INITIAL STATE VECTOR FROM KEPLERIAN ORBITAL ELEMENS IN THE ORBITAL PLANE
    call initial_state(P,a,e,mA,ra0)
    ! IT ROTATES THE STATE VECTORS FROM ORBITAL REF. SYSTEM TO THE OBSERVER REF. SYSTEM
    call orb2obs(ra0,-lN,-inc,-w,ra1) !OK 3T 1T 3T

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
    if(dt1.lt.zero)then
!       write(*,'(a,g25.14,a)')" dt1 = ",dt1,&
!           &"< 0 => backward integration"
      call ode_a_o(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
          &m,R,ra1,dt1,clN,cntRV,RV_stat,RV_sim,&
          &cntT0,T0_stat,T0_sim,Hc)
      if(Hc)then
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

    write(*,*)
    if(.not.Hc) then
      write(*,'(a,l2,a)')' flag Hc = ',Hc,' : BAD ==> PROBLEM DURING INTEGRATION'
    else
      write(*,'(a,l2,a)')' flag Hc = ',Hc,' : OK  ==> NO PROBLEM DURING INTEGRATION'
    end if
    write(*,*)
    flush(6)
    
    if((idtra.ge.1).and.(idtra.le.NB)) call close_tra(utra,fltra)
    if(wrtorb.eq.1) close(uorb)
    if(wrtconst.eq.1) close(ucon)
    if(wrtel.eq.1) call close_elem(uele,flele)

    if(.not.Hc) then
!       resw=sqrt(resmax)/real(ndata,dp)
      resw=set_max_residuals(ndata)
    else
      !call setval(RVobs,RV_sim,T0obs,T0_sim,resw)
      !call setval_2(RVobs,RV_sim,eRVobs,T0obs,T0_sim,eT0obs,resw)
!       call setval_3(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
!       call setval_4(RVobs,RV_sim,eRVobs,gamma,T0obs,T0_sim,eT0obs,resw)
      call setval_5(RVobs,RV_sim,eRVobs,gamma,epoT0obs,T0obs,eT0obs,T0_sim,resw,oc_fit)
      if(present(fit_scale))then
        if(cntRV.gt.0) call check_periodogram_scale(jdRV,RVobs-RV_sim,eRVobs,P,gls_check,gls_scale)
        write(*,*)
        write(*,'(a)')' CHECK GLS RESIDUALS-RV'
        write(*,'(a,es23.16)')' gls_scale = ',gls_scale
        write(*,*)
        flush(6)
        resw=resw*fit_scale*gls_scale
      else
        if(cntRV.gt.0) call check_periodogram(jdRV,RVobs-RV_sim,eRVobs,P,gls_check)
      end if
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
      ! 2016-04-08: added gls check
      if(present(gls_scale))then
        call check_and_write_periodogram(cpuid,isim,wrtid,jdRV,RVobs-RV_sim,eRVobs,P,gls_check,gls_scale)
      else
        call check_and_write_periodogram(cpuid,isim,wrtid,jdRV,RVobs-RV_sim,eRVobs,P,gls_check)
      end if
      if(.not.gls_check)resw=set_max_residuals(ndata)
    else
    write(*,*)
      write(*,'(a)')' RADIAL VELOCITIES NOT FOUND'
      write(*,*)
    end if

    if(cntT0.gt.0)then
      write(*,*)
      write(*,'(a,i5)')" T0 SIM found ",cntT0
      write(*,*)
      call write_T0(T0_sim,T0_stat)
      call write_T0(cpuid,isim,wrtid,T0_sim,T0_stat)
      allocate(resw_temp(sum(nT0)))
      resw_temp=zero
      call set_T0_resw(T0obs,T0_sim,eT0obs,resw_temp)
!       resw_temp=(T0obs-T0_sim)/eT0obs
      chi2r_T0=sum(resw_temp*resw_temp)*inv_dof
      w_chi2r=real(ndata,dp)/real(sum(nT0),dp)
      chi2wr_T0=chi2r_T0*w_chi2r
!       deallocate(resw_temp)
!       allocate(resw_temp(sum(nT0)))
      resw_temp=zero
      call set_oc_resw(epoT0obs,T0obs,eT0obs,T0_sim,resw_temp)
      chi2r_oc=sum(resw_temp*resw_temp)*inv_dof
      deallocate(resw_temp)
    else
      write(*,*)
      write(*,'(a)')' TRANSIT TIMES NOT FOUND'
      write(*,*)
    end if

    fitness = sum(resw*resw)
    if(present(fit_scale))then
      call write_fitness_summary(cpuid,isim,wrtid,chi2r_RV,chi2wr_RV,chi2r_T0,chi2wr_T0,chi2r_oc,fitness,fit_scale,gls_scale)
    else
      call write_fitness_summary(cpuid,isim,wrtid,chi2r_RV,chi2wr_RV,chi2r_T0,chi2wr_T0,chi2r_oc,fitness)
    end if

    if(allocated(RV_sim)) deallocate(RV_stat,RV_sim)
    if(allocated(gamma)) deallocate(gamma)
    if(allocated(T0_sim)) deallocate(T0_stat,T0_sim)
    if(allocated(m))   deallocate(m,R,P,a,e,w,mA,inc,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_out
  ! ------------------------------------------------------------------ !

end module ode_run



