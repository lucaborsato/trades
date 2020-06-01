module ode_run
  use constants
  use custom_type
  use parameters
  use parameters_conversion
  use linear_ephem
  use celestial_mechanics
!   use rotations,only:orb2obs
  use eq_motion,only:eqmastro
  use numerical_integrator,only:int_rk_a,rkck_a
  use transits
  use radial_velocities
  use gls_module,only:check_periodogram,check_periodogram_scale,check_and_write_periodogram
  use output_files
  use sorting,only:indexx
  use statistics,only:mean
  implicit none

  contains

    ! set only RV to the weighted residuals
!   subroutine set_RV_resw(RV_obs,RV_sim,e_RVobs,gamma,resw)
  subroutine set_RV_resw(obsRV,simRV,resw)
    use statistics,only:wmean
!     real(dp),dimension(:),intent(in)::RV_obs,RV_sim,e_RVobs
!     real(dp),dimension(:,:),allocatable,intent(out)::gamma
    type(dataRV),intent(in)::obsRV
    type(dataRV),intent(inout)::simRV

    real(dp),dimension(:),intent(inout)::resw
    real(dp),dimension(:),allocatable::xRV,wi
    integer::j,aRV,bRV,nset,nRVs

    nset=obsRV%nRVset
    if(.not.allocated(simRV%gamma)) allocate(simRV%gamma(nset,2))
    simRV%gamma=zero
    if(obsRV%nRV.gt.0)then
      aRV=0
      bRV=0
      do j=1,nset
        nRVs=obsRV%nRVsingle(j)
        aRV=bRV+1
        bRV=bRV+nRVs
        allocate(xRV(nRVs),wi(nRVs))
        xRV=zero
        xRV=obsRV%RV(aRV:bRV)-(simRV%RV(aRV:bRV)+simRV%trend(aRV:bRV))
        wi=one/(obsRV%eRV(aRV:bRV)*obsRV%eRV(aRV:bRV))
        simRV%gamma(j,1)=wmean(xRV,wi)
        simRV%gamma(j,2)=sqrt(one/sum(wi))
        xRV=xRV-simRV%gamma(j,1)
        resw(aRV:bRV)=xRV/obsRV%eRV(aRV:bRV)
        deallocate(xRV,wi)
      end do
    end if

    return
  end subroutine set_RV_resw

  ! set only T0 to the weighted residuals
  subroutine set_T0_resw(obsT0,simT0,resw)
!     real(dp),dimension(:,:),intent(in)::T0_obs,T0_sim,eT0_obs
    type(dataT0),dimension(:),intent(in)::obsT0
    type(dataT0),dimension(:),intent(inout)::simT0
    real(dp),dimension(:),intent(out)::resw

    integer::j,j1,a,b,nT0

    ! NB is a common variable
    resw=zero

    a=0
    b=0
    do j=2,NB
      j1=j-1
      nT0=obsT0(j1)%nT0
      if(nT0.gt.0)then
        a=a+1
        b=b+nT0
        resw(a:b)=(obsT0(j1)%T0-simT0(j1)%T0)/obsT0(j1)%eT0
        a=b
      end if
    end do


    return
  end subroutine set_T0_resw

  ! set only dur to the weighted residuals
  subroutine set_dur_resw(obsT0,simT0,resw)
!     real(dp),dimension(:,:),intent(in)::T0_obs,T0_sim,eT0_obs
    type(dataT0),dimension(:),intent(in)::obsT0
    type(dataT0),dimension(:),intent(inout)::simT0
    real(dp),dimension(:),intent(out)::resw

    integer::j,j1,a,b,nd

    ! NB is a common variable
    resw=zero

    a=0
    b=0
    do j=2,NB
      j1=j-1
      nd=obsT0(j1)%nDur
      if(nd.gt.0)then
        a=a+1
        b=b+nd
        resw(a:b)=(obsT0(j1)%dur-simT0(j1)%dur)/obsT0(j1)%edur
        a=b
      end if
    end do

    return
  end subroutine set_dur_resw

  ! compute linear ephemeris for both obs and sim
  ! use O-Cobs/sim to compute residuals
!   subroutine set_oc_resw(epoT0_obs,T0_obs,eT0_obs,T0_sim,resw)
  subroutine set_oc_resw(obsT0,simT0,resw)
!     integer,dimension(:,:),intent(in)::epoT0_obs
!     real(dp),dimension(:,:),intent(in)::T0_obs,eT0_obs,T0_sim
    type(dataT0),dimension(:),intent(in)::obsT0
    type(dataT0),dimension(:),intent(inout)::simT0
    real(dp),dimension(:),intent(out)::resw

!     integer,dimension(:),allocatable::x
!     real(dp),dimension(:),allocatable::y,ey
!     real(dp)::Tref_sim,eTref_sim,Pref_sim,ePref_sim

!     real(dp)::obs,sim,eobs

    integer::j,nTx,a,b

    ! NB is a common variables
    resw=zero

    a=0
    b=0
    do j=2,NB
      nTx=obsT0(j-1)%nT0

      if(nTx.gt.0)then
        a=a+1
        b=b+nTx
        ! it computes the linear ephemeris from simulated data
        ! call set_ephem(simT0(j-1))
        call set_ephem_simT0(simT0(j-1))
        call compute_oc_one_planet(simT0(j-1))
        resw(a:b)=(obsT0(j-1)%oc-simT0(j-1)%oc)/obsT0(j-1)%eT0
        a=b
      end if
    end do

    return
  end subroutine set_oc_resw


  subroutine set_fitness(oDataIn,val,resw)
    type(dataObs),intent(in)::oDataIn
    real(dp),dimension(:),intent(in)::val
    real(dp),dimension(:),intent(out)::resw

    real(dp),dimension(:),allocatable::resa1,resa2,resb1,resb2
    integer::nRV,nTTs,nDurs

    nRV=oDataIn%obsRV%nRV
    nTTs=oDataIn%nTTs
    nDurs=oDataIn%nDurs

    ! set the values scaled by the k_chi2r and k_chi2wr
    ! They are set in such a way the sum of resw^2 is:
    resw=zero
    if(k_chi2wr.eq.zero)then
      ! Chi2r (Normal way)
      resw=val*sqrt(oDataIn%inv_dof) ! chi2r = chi2 / dof

    else if(k_chi2r.eq.zero)then
      ! Chi2wr (Only weighted chi2r)

      ! only RV, no TT
      if(nRV.ne.0.and.nTTs.eq.0)then
        resw=val*k_b

      ! only TT,no RV
      else if(nRV.eq.0.and.nTTs.ne.0)then
        resw(1:nTTs)=val(1:nTTS)*k_b(1)
        if(durcheck.eq.1) resw(1+nTTs:nTTS+nDurs)=val(1+nTTs:nTTS+nDurs)*k_b(2)

      else
        resw(1:nRV)=val(1:nRV)*k_b(1)
        resw(1+nRV:nRV+nTTs)=val(1+nRV:nRV+nTTs)*k_b(2)
        if(durcheck.eq.1) &
          &resw(1+nRV+nTTs:nRV+nTTs+nDurs)=val(1+nRV+nTTs:nRV+nTTs+nDurs)*k_b(3)

      end if

    else
      ! fitness = Chi2r*k_chi2r + Chi2wr*k_chi2wr

      ! NO T0: nRV!=0, nTTs==0
      if(nRV.ne.0.and.nTTs.eq.0)then
        allocate(resa1(nRV),resb1(nRV))
        resa1=val*k_a
        resb1=val*k_b
        resw=sqrt(resa1*resa1 + resb1*resb1)
        deallocate(resa1,resb1)

      ! NO RV: nRV==0, nT0!=0
      else if(nRV.eq.0.and.nTTs.ne.0)then
        allocate(resa2(nTTs),resb2(nTTs))
        resa2=val(1:nTTs)*k_a
        resb2=val(1:nTTs)*k_b(1)
        resw(1:nTTs)=sqrt(resa2*resa2 + resb2*resb2)
        if(durcheck.eq.1)then
          resa2=val(1+nTTs:nTTs+nDurs)*k_a
          resb2=val(1+nTTs:nTTs+nDurs)*k_b(2)
          resw(1+nTTs:nTTs+nDurs)=sqrt(resa2*resa2 + resb2*resb2)
        end if
        deallocate(resa2,resb2)

      else

        allocate(resa1(nRV),resb1(nRV),resa2(nTTs),resb2(nTTs))
        ! RV
        resa1=val(1:nRV)*k_a
        resb1=val(1:nRV)*k_b(1)
        resw(1:nRV)=sqrt(resa1*resa1 + resb1*resb1)
        ! TT
        resa2=val(1+nRV:nRV+nTTs)*k_a
        resb2=val(1+nRV:nRV+nTTs)*k_b(2)
        resw(nRV+1:nRV+nTTs)=sqrt(resa2*resa2 + resb2*resb2)
        ! Dur
        if(durcheck.eq.1)then
          resa2=zero
          resa2=val(1+nRV+nTTs:nRV+nTTs+nDurs)*k_a
          resb2=zero
          resb2=val(1+nRV+nTTs:nRV+nTTs+nDurs)*k_b(3)
          resw(1+nRV+nTTs:nRV+nTTs+nDurs)=sqrt(resa2*resa2 + resb2*resb2)
        end if
        deallocate(resa1,resa2,resb1,resb2)
      end if

    end if

    return
  end subroutine set_fitness
  ! ------------------------------------------------------------------ !


  ! setting properly the weighted residuals: RV and T0 or OC
  !   subroutine setval_1(RV_obs,RV_sim,e_RVobs,gamma,epoT0_obs,T0_obs,eT0_obs,T0_sim,resw,ocfit)
  subroutine setval_1(oDataIn,simRV,simT0,resw,ocfit)
    type(dataObs),intent(in)::oDataIn
!     type(dataRV),intent(in)::simRV
!     type(dataT0),dimension(:),intent(in)::simT0
    type(dataRV),intent(inout)::simRV
    type(dataT0),dimension(:),intent(inout)::simT0
    real(dp),dimension(:),intent(out)::resw
    integer,intent(in)::ocfit

    real(dp),dimension(:),allocatable::val,val_T0,val_dur,val_oc
    integer::nRV,nTTs,nDurs

    nRV=oDataIn%obsRV%nRV
    nTTs=oDataIn%nTTs
    nDurs=oDataIn%nDurs
    allocate(val(oDataIn%ndata))

    resw=zero
    val=zero

    if(nRV.gt.0) call set_RV_resw(oDataIn%obsRV,simRV,val(1:nRV))

    if(ocfit.eq.1)then
      ! fit only O-C

      allocate(val_oc(nTTs))
      val_oc=zero
      call set_oc_resw(oDataIn%obsT0,simT0,val_oc)
      val(nRV+1:nRV+nTTs)=val_oc
      deallocate(val_oc)

      if(durcheck.eq.1)then
        allocate(val_dur(nDurs))
        val_dur=zero
        call set_dur_resw(oDataIn%obsT0,simT0,val_dur)
        val(nRV+nTTs+1:nRV+nTTs+nDurs)=val_dur
        deallocate(val_dur)
      end if

    else if(ocfit.eq.2)then
      ! fit T0 and O-C and weight it half

      allocate(val_T0(nTTs),val_oc(nTTs))
      val_T0=zero
      call set_T0_resw(oDataIn%obsT0,simT0,val_T0)
      val_oc=zero
      call set_oc_resw(oDataIn%obsT0,simT0,val_oc)
      val(nRV+1:nRV+nTTs)=sqrt_half*sqrt(val_T0*val_T0+val_oc*val_oc)
      deallocate(val_T0,val_oc)

      if(durcheck.eq.1)then
        allocate(val_dur(nDurs))
        val_dur=zero
        call set_dur_resw(oDataIn%obsT0,simT0,val_dur)
        val(nRV+nTTs+1:nRV+nTTs+nDurs)=val_dur
        deallocate(val_dur)
      end if

    else
      ! fit only T0, default way

      allocate(val_T0(nTTs))
      val_T0=zero
      call set_T0_resw(oDataIn%obsT0,simT0,val_T0)
      val(nRV+1:nRV+nTTs)=val_T0
      deallocate(val_T0)

      if(durcheck.eq.1)then
        allocate(val_dur(nDurs))
        val_dur=zero
        call set_dur_resw(oDataIn%obsT0,simT0,val_dur)
        val(nRV+nTTs+1:nRV+nTTs+nDurs)=val_dur
        deallocate(val_dur)
      end if

    end if
    call set_fitness(obsData,val,resw)
!     resw=val
    deallocate(val)

    return
  end subroutine setval_1



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

    if(ndata.gt.0)then
      resw_max=maxval(resw)
      resw_max_possible=set_max_residuals(ndata)
      if(resw_max.gt.resw_max_possible) resw=resw_max_possible
    end if

    return
  end subroutine check_max_residuals


  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! performs the integration, integration time to be provided as argument
!   subroutine ode_a(m,R,rin,time,clN,cntRV,RV_stat,RV_sim,cntT0,T0_stat,T0_sim,Hc)
subroutine ode_a(m,R,rin,time,clN,simRV,simT0,Hc)
    real(dp),dimension(:),intent(in)::m,R,rin
    real(dp),intent(in)::time
    integer,dimension(:),intent(in)::clN
!     integer,intent(inout)::cntRV
!     integer,dimension(:),intent(inout)::RV_stat
!     real(dp),dimension(:),intent(inout)::RV_sim
!     integer,intent(inout)::cntT0
!     real(dp),dimension(:,:),intent(inout)::T0_sim
!     integer,dimension(:,:),intent(inout)::T0_stat
    type(dataRV),intent(inout)::simRV
    type(dataT0),dimension(:),intent(inout)::simT0
    logical,intent(out)::Hc

    real(dp),dimension(:),allocatable::dr,r1,r2,err
    integer,dimension(:),allocatable::X,Y,Z
    real(dp),dimension(:),allocatable::cX,cY,cR,rmean
    real(dp)::hw,hok,hnext,itime
    integer::j,j1 !,nj
    integer::nRV,nTTs,nDurs

    nRV=obsData%obsRV%nRV
    nTTs=obsData%nTTs
    nDurs=obsData%nDurs

    Hc = separation_mutual_Hill_check(m,R,rin,do_hill_check)
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

      Hc = separation_mutual_Hill_check(m,R,r2,do_hill_check)
      if(.not.Hc) return

      ! RV check
      if((nRV.gt.0).and.(simRV%nRV.lt.nRV)) call check_RV(m,r1,dr,itime,hok,simRV)

      ! T0 check
      do j=2,NB
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))
        rmean(j)=half*( rsky(r1(X(j):Y(j))) + rsky(r2(X(j):Y(j))) )
      end do

      if((idtra.gt.0).and.(idtra.le.NB))then ! if 1
        ! if((sum(simT0(:)%nT0).lt.nTTs)then ! if 2
          if(idtra.eq.1)then ! if 3
            do j=2,NB
              if(rmean(j).le.cR(j))then ! if 4
                if(clN(j).eq.0)then ! if 5

                  if((cX(j).le.zero).and.(r1(Z(j)).gt.zero))then ! if 6
                    call check_T0(j,m,R,r1,r2,itime,hok,simT0,Hc)
                    if(.not.Hc) return
                  end if ! end if 6

                else

                  if((cY(j).le.zero).and.(r1(Z(j)).gt.zero))then ! if 7
                    call check_T0(j,m,R,r1,r2,itime,hok,simT0,Hc)
                    if(.not.Hc) return
                  end if ! end if 7

                end if ! end if 5
              end if ! end if 4
            end do
          else
            if(rmean(idtra).le.cR(idtra))then ! if 8
              if(clN(idtra).eq.0)then ! if 9

                if((cX(idtra).le.zero).and.(r1(Z(idtra)).gt.zero))then ! if 10
                  call check_T0(idtra,m,R,r1,r2,itime,hok,simT0,Hc)
                  if(.not.Hc) return
                end if ! end if 10

              else

                if((cY(idtra).le.zero).and.(r1(Z(idtra)).gt.zero))then ! if 11
                  call check_T0(idtra,m,R,r1,r2,itime,hok,simT0,Hc)
                  if(.not.Hc) return
                end if ! end if 11

              end if ! end if 9
            end if ! end if 8
          end if ! end if 3
        ! end if ! end if 2
      end if ! end if 1

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
    &m,R,rin,time,clN,simRV,simT0,Hc)
    integer,intent(in)::uorb,ucon
    integer,dimension(:),intent(in)::uele,utra
    character(*),intent(in)::fmorb,fmcon,fmele
    real(dp),dimension(:),intent(in)::m,R,rin
    real(dp),intent(in)::time
    integer,dimension(:),intent(in)::clN

    type(dataRV),intent(inout)::simRV
    type(dataT0),dimension(:),intent(inout)::simT0

    logical,intent(inout)::Hc

    real(dp),dimension(:),allocatable::dr,r1,r2,err
    real(dp),dimension(:),allocatable::drdt_save,r_save
    integer,dimension(:),allocatable::X,Y,Z
    real(dp),dimension(:),allocatable::cX,cY,cR,rmean
    real(dp)::hw,hok,hnext,itime
    integer::j,j1,j2,jtra
    real(dp),dimension(:,:),allocatable::storeorb,storecon,storetra
    integer,dimension(:,:),allocatable::stat_tra
    integer::Norb
    real(dp)::itwrt,Etot,Eold,htot,hold

    real(dp)::step_save,step_write

    integer::nRV

    nRV=obsData%obsRV%nRV

    Hc = separation_mutual_Hill_check(m,R,rin,do_hill_check)
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
    rmean=9.e8_dp

    Norb=NBDIM+3

    allocate(dr(NBDIM),r1(NBDIM),r2(NBDIM),err(NBDIM))
    allocate(drdt_save(NBDIM),r_save(NBDIM))
    hw=step_0
    if(time.lt.zero) hw=-hw
    itime=zero
    r1=rin
    r2=zero
    err=zero

    j1=0
    j2=1

    step_write=wrttime
    if(time.lt.zero) step_write=-step_write
    itwrt=step_write

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

    integration: do
      j1=j1+1

      if(abs(itime+hw).gt.abs(time)) hw=time-itime
      call eqmastro(m,r1,dr) ! computes the eq. of motion
      call int_rk_a(m,r1,dr,hw,hok,hnext,r2,err) ! computes the next orbit step

      ! check if it has to compute or not 'something'
      if((wrtorb.eq.1).or.(wrtel.eq.1).or.(wrtconst.eq.1))then

        if(abs(itime+hok).ge.(abs(itwrt)))then

          saveloop: do

            j2=j2+1
!             computes the proper step -> itwrt will be updated in the loop
            step_save=itwrt-itime
!             computes state vector from r1 to the new time step (smaller then hok)
            drdt_save=dr
            r_save=r1
            call rkck_a(m,r1,drdt_save,step_save,r_save,err)

            ! check and store orbit or kep. elements
            if((wrtorb.eq.1).or.(wrtel.eq.1))then
              call store_orb(j2,itime+step_save,m,r_save,storeorb) ! save r_save!!
              if(j2.eq.DIMMAX)then ! write into file if reach max dimension
                if(wrtorb.eq.1) call write_file(DIMMAX,uorb,fmorb,storeorb)
                if(wrtel.eq.1) call write_elem(DIMMAX,uele,fmele,m,storeorb)
                storeorb=zero ! reset the storeorb variable
              end if
            end if

            ! check and store constants of motion
            if(wrtconst.eq.1)then
              call compute_con(m,r_save,Etot,htot) ! if selected it computes the Energy and Angular momentum
              call store_con(j2,itime+step_save,Etot,Eold,htot,hold,storecon)
              if(j2.eq.DIMMAX)then ! write into file if reach max dimension
                call write_file(DIMMAX,ucon,fmcon,storecon)
                storecon=zero
              end if
            end if

            if(j2.eq.DIMMAX) j2=0

            itwrt=itwrt+step_write
            if(abs(itwrt).gt.abs(itime+hok)) exit saveloop

          end do saveloop

        end if

      end if
      ! ===============================

      Hc = separation_mutual_Hill_check(m,R,r2,do_hill_check)
      if(.not.Hc)then
        ! write(*,'(a,l)')' INLOOP Hc: ',Hc
        return
      end if

      ! RV check
      if((nRV.gt.0).and.(simRV%nRV.lt.nRV)) call check_RV(m,r1,dr,itime,&
        &hok,simRV)

      ! T0 check (to compare and all)
      do j=2,NB
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))
        rmean(j)=half*( rsky(r1(X(j):Y(j))) + rsky(r2(X(j):Y(j))) )
      end do
      if((idtra.gt.0).and.(idtra.le.NB))then
        if(idtra.eq.1)then
          do j=2,NB
            if(rmean(j).le.cR(j))then
              if(clN(j).eq.0)then

                if((cX(j).le.zero).and.(r1(Z(j)).gt.zero))then
                  if(obsData%obsT0(j-1)%nT0.gt.0) &
                    &call check_T0(j,m,R,r1,r2,itime,hok,simT0,Hc)
                  if(.not.Hc) return
                  jtra=jtra+1
                  call all_transits(jtra,j,m,R,r1,r2,&
                    &itime,hok,stat_tra,storetra)
                end if

              else

                if((cY(j).le.zero).and.(r1(Z(j)).gt.zero))then
                  if(obsData%obsT0(j-1)%nT0.gt.0) &
                    &call check_T0(j,m,R,r1,r2,itime,hok,simT0,Hc)
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
                if(obsData%obsT0(idtra-1)%nT0.gt.0) &
                    &call check_T0(idtra,m,R,r1,r2,itime,hok,simT0,Hc)
                if(.not.Hc) return
                jtra=jtra+1
                call all_transits(jtra,idtra,m,R,r1,r2,&
                  &itime,hok,stat_tra,storetra)
              end if

            else

              if((cY(idtra).le.zero).and.&
                &(r1(Z(idtra)).gt.zero))then
                if(obsData%obsT0(idtra-1)%nT0.gt.0) &
                    &call check_T0(idtra,m,R,r1,r2,itime,hok,simT0,Hc)
                if(.not.Hc) return
                jtra=jtra+1
                call all_transits(jtra,idtra,m,R,r1,r2,&
                  &itime,hok,stat_tra,storetra)
              end if

            end if
          end if
        end if

        if(jtra.eq.DIMMAX)then
          call write_tra(jtra,utra,stat_tra,storetra)
          jtra=0
          stat_tra=0
          storetra=zero
        end if

      end if

      itime=itime+hok

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
    deallocate(dr,r1,r2,err,drdt_save,r_save)

    return
  end subroutine ode_a_o
  ! ------------------------------------------------------------------ !


    ! ------------------------------------------------------------------ !
  ! subroutine that integrates orbit, computes, stores and
  ! writes orbits, orbital elements, constants of motion, transit time and RV into files
!   ! it is called by ode_integrates, that check in which direction (in time) integrates
  ! it doesn't need obsRV.dat or NB#_observations.dat files, it computes only all the possible transits
  ! and the RVs are stored in the #_#_rotorbit.dat file.
  subroutine ode_a_orbit(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
      &m,R,rin,time,clN,Hc)
    integer,intent(in)::uorb,ucon
    integer,dimension(:),intent(in)::uele,utra
    character(*),intent(in)::fmorb,fmcon,fmele
    real(dp),dimension(:),intent(in)::m,R,rin
    real(dp),intent(in)::time
    integer,dimension(:),intent(in)::clN
    logical,intent(inout)::Hc

    real(dp),dimension(:),allocatable::dr,r1,r2,err
    real(dp),dimension(:),allocatable::drdt_save,r_save
    integer,dimension(:),allocatable::X,Y,Z
    real(dp),dimension(:),allocatable::cX,cY,cR,rmean
    real(dp)::hw,hok,hnext,itime
    integer::j,j1,j2,jtra
    real(dp),dimension(:,:),allocatable::storeorb,storecon,storetra
    integer,dimension(:,:),allocatable::stat_tra
    integer::Norb
    real(dp)::itwrt,Etot,Eold,htot,hold

    real(dp)::step_save,step_write

    Hc = separation_mutual_Hill_check(m,R,rin,do_hill_check)

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
    rmean=9.e8_dp

    Norb=NBDIM+3

    allocate(dr(NBDIM),r1(NBDIM),r2(NBDIM),err(NBDIM))
    allocate(drdt_save(NBDIM),r_save(NBDIM))
    hw=step_0
    if(time.lt.zero) hw=-hw
    itime=zero
    r1=rin
    r2=zero
    err=zero

    j1=0
    j2=1

    step_write=wrttime
    if(time.lt.zero) step_write=-step_write
    itwrt=step_write

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

    integration: do
      j1=j1+1

      if(abs(itime+hw).gt.abs(time)) hw=time-itime
      call eqmastro(m,r1,dr) ! computes the eq. of motion
      call int_rk_a(m,r1,dr,hw,hok,hnext,r2,err) ! computes the next orbit step

      ! ===============================
      ! NEW VERSION
      ! check if it has to compute or not 'something'
      if((wrtorb.eq.1).or.(wrtel.eq.1).or.(wrtconst.eq.1))then

        if(abs(itime+hok).ge.(abs(itwrt)))then

          saveloop: do

            j2=j2+1
!             computes the proper step -> itwrt will be updated in the loop
            step_save=itwrt-itime
!             computes state vector from r1 to the new time step (smaller then hok)
            drdt_save=dr
            r_save=r1
            call rkck_a(m,r1,drdt_save,step_save,r_save,err)

            ! check and store orbit or kep. elements
            if((wrtorb.eq.1).or.(wrtel.eq.1))then
              call store_orb(j2,itime+step_save,m,r_save,storeorb) ! save r_save!!
              if(j2.eq.DIMMAX)then ! write into file if reach max dimension
                if(wrtorb.eq.1) call write_file(DIMMAX,uorb,fmorb,storeorb)
                if(wrtel.eq.1) call write_elem(DIMMAX,uele,fmele,m,storeorb)
                storeorb=zero ! reset the storeorb variable
              end if
            end if

            ! check and store constants of motion
            if(wrtconst.eq.1)then
              call compute_con(m,r_save,Etot,htot) ! if selected it computes the Energy and Angular momentum
              call store_con(j2,itime+step_save,Etot,Eold,htot,hold,storecon)
              if(j2.eq.DIMMAX)then ! write into file if reach max dimension
                call write_file(DIMMAX,ucon,fmcon,storecon)
                storecon=zero
              end if
            end if

            if(j2.eq.DIMMAX) j2=0

            itwrt=itwrt+step_write
            if(abs(itwrt).gt.abs(itime+hok)) exit saveloop

          end do saveloop

        end if

      end if
      ! ===============================

      Hc = separation_mutual_Hill_check(m,R,r2,do_hill_check)

      ! T0 check (to compare and all)
      do j=2,NB
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))
        rmean(j)=half*( rsky(r1(X(j):Y(j))) + rsky(r2(X(j):Y(j))) )
      end do
      if((idtra.gt.0).and.(idtra.le.NB))then
        if(idtra.eq.1)then
          do j=2,NB
            if(rmean(j).le.cR(j))then
              if(clN(j).eq.0)then

                if((cX(j).le.zero).and.(r1(Z(j)).gt.zero))then
                  jtra=jtra+1
                  call all_transits(jtra,j,m,R,r1,r2,&
                    &itime,hok,stat_tra,storetra)
                end if

              else

                if((cY(j).le.zero).and.(r1(Z(j)).gt.zero))then
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
                jtra=jtra+1
                call all_transits(jtra,idtra,m,R,r1,r2,&
                  &itime,hok,stat_tra,storetra)
              end if

            else

              if((cY(idtra).le.zero).and.&
                &(r1(Z(idtra)).gt.zero))then
                jtra=jtra+1
                call all_transits(jtra,idtra,m,R,r1,r2,&
                  &itime,hok,stat_tra,storetra)
              end if

            end if
          end if
        end if

        if(jtra.eq.DIMMAX)then
          call write_tra(jtra,utra,stat_tra,storetra)
          jtra=0
          stat_tra=0
          storetra=zero
        end if

      end if

      itime=itime+hok

      if(abs(itime).ge.abs(time)) exit integration
      hw=hnext
      r1=r2

    end do integration

    if((idtra.ge.1).and.(idtra.le.NB)) call write_tra(jtra,utra,&
      &stat_tra,storetra)

    if((wrtorb.eq.1).or.(wrtel.eq.1))then
      if(j2.gt.0)then
        if(wrtorb.eq.1) call write_file(j2,uorb,fmorb,storeorb)
        if(wrtel.eq.1) call write_elem(j2,uele,fmele,m,storeorb)
      end if
      deallocate(storeorb)
    end if
    if(wrtconst.eq.1)then
      if(j2.gt.0)then
        call write_file(j2,ucon,fmcon,storecon)
      end if
      deallocate(storecon)
    end if

    deallocate(X,Y,Z,cX,cY,cR,rmean)
    deallocate(dr,r1,r2,err,drdt_save,r_save)

    return
  end subroutine ode_a_orbit
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

    type(dataRV)::simRV

    type(dataT0),dimension(:),allocatable::simT0

    logical::checkpar,gls_check

    integer::ndata,nRV,nTTs,nT0,ibd

    ndata=obsData%ndata
    nRV=obsData%obsRV%nRV
    nTTs=obsData%nTTs

    Hc=.true.
    resw=zero
    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),w(NB),mA(NB),inc(NB),lN(NB),clN(NB))

    call convert_parameters(allpar,par,m,R,P,a,e,w,mA,inc,lN,checkpar)
    if(.not.checkpar)then
      resw=set_max_residuals(ndata)
      return
    end if

    ! it is needed to define which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))

    ! NEW VERSION 2017-11-21
    call kepelements2statevector(m,a,e,mA,w,inc,lN,ra0)
    ra1=ra0

    if(nRV.gt.0)then
      call init_dataRV(nRV,simRV)
      simRV%nRV=0
    end if

    if(nTTs.gt.0)then
      allocate(simT0(NB-1))
      do ibd=1,NB-1
        nT0=obsData%obsT0(ibd)%nT0
        call init_dataT0(nT0,simT0(ibd),durcheck)
        simT0(ibd)%nT0=0
        simT0(ibd)%nDur=0
        !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
        simT0(ibd)%eT0=one/s24h ! 1s in day
        simT0(ibd)%edur=one/60.0_dp ! 1s in min
      end do
    end if

    dt1=tstart-tepoch
    dt2=dt1+tint
    if(dt1.lt.zero)then

      call ode_a(m,R,ra1,dt1,clN,simRV,simT0,Hc)
      if(Hc)then
        if(abs(dt1).le.tint)then
          call ode_a(m,R,ra1,dt2,clN,simRV,simT0,Hc)
        end if
      end if

    else

      call ode_a(m,R,ra1,dt2,clN,simRV,simT0,Hc)

    end if

    gls_check=.true.
    if(.not.Hc) then
        resw=set_max_residuals(ndata)
    else

      ! compute RV trend
      if(rv_trend_order.gt.0)then
        call addRVtrend(simRV%jd, par(nfit-rv_trend_order:nfit), simRV%trend)
      end if

      call setval_1(obsData,simRV,simT0,resw,oc_fit)

      if(simRV%nRV.gt.0) call check_periodogram(obsData%obsRV%jd,&
        &obsData%obsRV%RV-simRV%RV-simRV%trend,obsData%obsRV%eRV,&
        &P,gls_check)

      if(.not.gls_check)then
        resw=set_max_residuals(ndata)
      else
        call check_max_residuals(resw,ndata)
      end if
      call check_max_residuals(resw,ndata)

    end if

    call deallocate_dataRV(simRV)
    if(nTTs.gt.0)then
      do ibd=1,NB-1
        call deallocate_dataT0(simT0(ibd))
      end do
    end if

    if(allocated(m))   deallocate(m,R,P,a,e,w,mA,inc,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_lm
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! subroutines called by the bootstrap to calculate the new set of T_0 and RV from
  ! the fitted parameters
  subroutine ode_parok(allpar,par,simRV,simT0,resw,fitness)
    real(dp),dimension(:),intent(in)::allpar,par

    type(dataRV),intent(inout)::simRV
    type(dataT0),dimension(:),intent(inout)::simT0

    real(dp),dimension(:),intent(out)::resw
    real(dp),intent(out)::fitness

    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp)::dt1,dt2
    logical::Hc

    logical::checkpar,gls_check

    integer::ndata,nRV,nTTs,nT0

    integer::ii,ibd

    ndata=obsData%ndata
    nRV=obsData%obsRV%nRV
    nTTs=obsData%nTTs

    Hc=.true.
    resw=zero
    fitness=zero

    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),&
        &w(NB),mA(NB),inc(NB),lN(NB),clN(NB))

    call convert_parameters(allpar,par,m,R,P,a,e,w,mA,inc,lN,checkpar)

    if(nRV.gt.0)then
      call init_dataRV(nRV,simRV)
      simRV%nRV=0
    end if
!
    if(nTTs.gt.0)then
      do ibd=1,NB-1
        nT0=obsData%obsT0(ibd)%nT0
        call init_dataT0(nT0,simT0(ibd),durcheck)
        simT0(ibd)%nT0=0
        simT0(ibd)%nDur=0
        !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
        simT0(ibd)%eT0=1./s24h ! 1s in day
        simT0(ibd)%edur=1./60.0_dp ! 1s in min
      end do
    end if

    if(.not.checkpar)then
      if(nRV.gt.0) simRV%RV=zero

      if(nTTs.gt.0)then
        do ibd=1,NB-1
          simT0(ibd)%T0=zero
        end do
      end if

      resw=set_max_residuals(ndata)
      fitness=resmax
      return
    end if

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

    ! NEW VERSION 2017-11-21
    call kepelements2statevector(m,a,e,mA,w,inc,lN,ra0)
    ra1=ra0

    dt1=tstart-tepoch
    dt2=dt1+tint
    if(dt1.lt.zero)then

      call ode_a(m,R,ra1,dt1,clN,simRV,simT0,Hc)

      if(Hc)then
        if(abs(dt1).le.tint)then
          dt2=dt1+tint
          call ode_a(m,R,ra1,dt2,clN,simRV,simT0,Hc)
        end if
      end if

    else

      dt2=dt1+tint
      call ode_a(m,R,ra1,dt2,clN,simRV,simT0,Hc)

    end if
    if(.not.Hc) then
      if(nRV.gt.0) simRV%RV=zero

      if(nTTs.gt.0)then
        do ibd=1,NB-1
          simT0(ibd)%T0=zero
        end do
      end if

      resw=set_max_residuals(ndata)
      fitness=resmax

    else

      ! compute RV trend
      if(rv_trend_order.gt.0)then
        call addRVtrend(simRV%jd, par(nfit-rv_trend_order:nfit), simRV%trend)
      end if
      call setval_1(obsData,simRV,simT0,resw,oc_fit)

      if(simRV%nRV.gt.0) call check_periodogram(obsData%obsRV%jd,&
        &obsData%obsRV%RV-simRV%RV-simRV%trend,obsData%obsRV%eRV,&
        &P,gls_check)

      call check_max_residuals(resw,ndata)

      fitness=sum(resw*resw)
    end if

    if(allocated(m)) deallocate(m,R,P,a,e,w,mA,inc,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_parok

  ! subroutine called by the L-M to fit the parameters with observations in input (RV and T0)
  subroutine ode_boot(allpar,nm,nn,par,oDataIn,resw,iflag)
    integer,intent(in)::nm,nn
    real(dp),dimension(:),intent(in)::allpar,par
    type(dataObs),intent(in)::oDataIn

    real(dp),dimension(:),intent(out)::resw
    integer,intent(inout)::iflag

    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp)::dt1,dt2
    logical::Hc

    type(dataRV)::simRV
    type(dataT0),dimension(:),allocatable::simT0

    integer::ndata,nRV,nTTs,nT0,ibd

    logical::checkpar,gls_check

    ndata=oDataIn%ndata
    nRV=oDataIn%obsRV%nRV
    nTTs=oDataIn%nTTs

    Hc=.true.
    resw=zero
    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),&
        &w(NB),mA(NB),inc(NB),lN(NB),clN(NB))

    call convert_parameters(allpar,par,m,R,P,a,e,w,mA,inc,lN,checkpar)
    if(.not.checkpar)then
        resw=set_max_residuals(ndata)
      return
    end if

    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))

    ! NEW VERSION 2017-11-21
    call kepelements2statevector(m,a,e,mA,w,inc,lN,ra0)
    ra1=ra0

    if(nRV.gt.0)then
      call init_dataRV(nRV,simRV)
      simRV%nRV=0
    end if

    if(nTTs.gt.0)then
      allocate(simT0(NB-1))
      do ibd=1,NB-1
        nT0=obsData%obsT0(ibd)%nT0
        call init_dataT0(nT0,simT0(ibd),durcheck)
        simT0(ibd)%nT0=0
        simT0(ibd)%nDur=0
        !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
        simT0(ibd)%eT0=1./s24h ! 1s in day
        simT0(ibd)%edur=1./60.0_dp ! 1s in min
      end do
    end if

    dt1=tstart-tepoch
    if(dt1.lt.zero)then
      call ode_a(m,R,ra1,dt1,clN,simRV,simT0,Hc)
      if(Hc)then
        if(abs(dt1).le.tint)then
          dt2=dt1+tint
          call ode_a(m,R,ra1,dt2,clN,simRV,simT0,Hc)
        end if
      end if
    else
      dt2=dt1+tint
      call ode_a(m,R,ra1,dt2,clN,simRV,simT0,Hc)
    end if
    if(.not.Hc) then
      resw=set_max_residuals(ndata)
    else

      ! compute RV trend
      if(rv_trend_order.gt.0)then
        call addRVtrend(simRV%jd, par(nfit-rv_trend_order:nfit), simRV%trend)
      end if
      call setval_1(oDataIn,simRV,simT0,resw,oc_fit)

      if(simRV%nRV.gt.0) call check_periodogram(oDataIn%obsRV%jd,&
        &oDataIn%obsRV%RV-simRV%RV-simRV%trend,oDataIn%obsRV%eRV,&
        &P,gls_check)

      call check_max_residuals(resw,ndata)

    end if

    call deallocate_dataRV(simRV)
    if(nTTs.gt.0)then
      do ibd=1,NB-1
        call deallocate_dataT0(simT0(ibd))
      end do
    end if

    if(allocated(m)) deallocate(m,R,P,a,e,w,mA,inc,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_boot
  ! ------------------------------------------------------------------ !


    subroutine ode_full_args(m,R,rin,time,step_in,clN,simRV,&
      &id_transit_body,transit_flag,dur_check,simT0,Hc)

    real(dp),dimension(:),intent(in)::m,R,rin
    real(dp),intent(in)::time,step_in
    integer,dimension(:),intent(in)::clN

    type(dataRV),intent(inout)::simRV

    integer,intent(in)::id_transit_body
    logical,dimension(:),intent(in)::transit_flag
    integer,intent(in)::dur_check

    type(dataT0),dimension(:),intent(inout)::simT0

    logical,intent(inout)::Hc

    integer::n_body,nb_dim
    real(dp),dimension(:),allocatable::dr,r1,r2,err
    integer,dimension(:),allocatable::X,Y,Z
    real(dp),dimension(:),allocatable::cX,cY,cR,rmean
    real(dp)::hw,hok,hnext,itime
    integer::j,j1 !,nj

    integer::nRV,nTTs

    Hc = separation_mutual_Hill_check(m,R,rin,do_hill_check)
    if(.not.Hc) return

    nRV=obsData%obsRV%nRV
    nTTs=obsData%nTTs

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

    j1=0
    integration: do
      j1=j1+1
      if(abs(itime+hw).gt.abs(time)) hw=time-itime
      call eqmastro(m,r1,dr)
      call int_rk_a(m,r1,dr,hw,hok,hnext,r2,err)

      Hc = separation_mutual_Hill_check(m,R,r2,do_hill_check)
      if(.not.Hc) return

      ! RV check
      if(simRV%nRV.lt.nRV) call check_RV(m,r1,dr,itime,hok,simRV)

      ! T0 check
      do j=2,n_body
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))
        rmean(j)=half*( rsky(r1(X(j):Y(j))) + rsky(r2(X(j):Y(j))) )
      end do

      if((id_transit_body.gt.0).and.(id_transit_body.le.n_body))then
        if(sum(simT0(:)%nT0).lt.nTTs)then
          if(id_transit_body.eq.1)then
            do j=2,n_body
              if(rmean(j).le.cR(j))then
                if(clN(j).eq.0)then

                  if((cX(j).le.zero).and.(r1(Z(j)).gt.zero))then
                      call check_T0(j,m,R,r1,r2,itime,hok,&
                        &transit_flag,dur_check,obsData,simT0,Hc)
                    if(.not.Hc) return
                  end if

                else

                  if((cY(j).le.zero).and.(r1(Z(j)).gt.zero))then
                    call check_T0(j,m,R,r1,r2,itime,hok,&
                      &transit_flag,dur_check,obsData,simT0,Hc)
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
                  call check_T0(id_transit_body,m,R,r1,r2,itime,hok,&
                    &transit_flag,dur_check,obsData,simT0,Hc)
                  if(.not.Hc) return
                end if

              else

                if((cY(id_transit_body).le.zero).and.&
                  &(r1(Z(id_transit_body)).gt.zero))then
                  call check_T0(id_transit_body,m,R,r1,r2,itime,hok,&
                    &transit_flag,dur_check,obsData,simT0,Hc)
                  if(.not.Hc) return
                end if

              end if
            end if
          end if
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
    &id_transit_body,transit_flag,dur_check,T0_sim)

    real(dp),intent(in)::t_start,t_epoch,step_in,t_int
    real(dp),dimension(:),intent(in)::m,R,P,ecc,argp,mA,inc,lN

    real(dp),dimension(:),intent(in)::tRV
    real(dp),dimension(:),allocatable,intent(out)::RV_sim

    integer,intent(in)::id_transit_body
    logical,dimension(:),intent(in)::transit_flag
    integer,intent(in)::dur_check

    real(dp),dimension(:,:),allocatable,intent(out)::T0_sim


    type(dataRV)::simRV
    type(dataT0),dimension(:),allocatable::simT0

    integer::n_body,nb_dim
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp),dimension(:),allocatable::sma
    real(dp)::dt1,dt2
    logical::Hc

    integer::nRV,nTTs
    integer::ibd,nT0

    nRV=size(tRV)
    nTTs=obsData%nTTs

    Hc=.true.

    ! it is needed to define which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    n_body=size(m)
    allocate(sma(n_body))
    sma=zero
    call semax_vec(m(1),m(2:n_body),P(2:n_body), sma(2:n_body))

    nb_dim=6*n_body
    if(.not.allocated(ra0)) allocate(ra0(nb_dim),ra1(nb_dim))

    ! NEW VERSION 2017-11-21
    call kepelements2statevector(m,sma,ecc,mA,argp,inc,lN,ra0)
    ra1=ra0

    if(nRV.gt.0)then
      call init_dataRV(nRV,simRV)
      simRV%nRV=0
      allocate(RV_sim(nRV))
      RV_sim=zero
      simRV%jd=tRV
    end if

    if(nTTs.gt.0)then
      allocate(simT0(n_body-1))
      do ibd=1,n_body-1
        nT0=obsData%obsT0(ibd)%nT0
        call init_dataT0(nT0,simT0(ibd),dur_check)
        simT0(ibd)%nT0=0
        simT0(ibd)%nDur=0
        !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
        simT0(ibd)%eT0=1./s24h ! 1s in day
        simT0(ibd)%edur=1./60.0_dp ! 1s in min
      end do
      allocate(T0_sim(maxval(obsData%obsT0(:)%nT0),n_body))
      T0_sim=zero
    end if

    dt1=t_start-t_epoch
    dt2=dt1+t_int
    if(dt1.lt.zero)then
      call ode_full_args(m,R,ra1,dt1,step_in,clN,simRV,&
        &id_transit_body,transit_flag,dur_check,simT0,Hc)

      if(Hc)then
        if(abs(dt1).le.t_int)then
          dt2=dt1+t_int
          call ode_full_args(m,R,ra1,dt2,step_in,clN,simRV,&
            &id_transit_body,transit_flag,dur_check,simT0,Hc)

        end if
      end if
    else
      dt2=dt1+t_int
      call ode_full_args(m,R,ra1,dt2,step_in,clN,simRV,&
        &id_transit_body,transit_flag,dur_check,simT0,Hc)

    end if

    if(.not.Hc) then
        if(nRV.gt.0)then
          simRV%RV=zero
          simRV%RV_stat=0
          RV_sim=zero
        end if
        if(nTTs.gt.0)then
          do ibd=1,n_body-1
            if(simT0(ibd)%nT0.gt.0)then
              simT0(ibd)%T0=zero
              simT0(ibd)%T0_stat=0
              if(dur_check.eq.1)then
                simT0(ibd)%dur=zero
                simT0(ibd)%dur_stat=0
              end if
              T0_sim(:,ibd+1)=zero
            end if
          end do
        end if

    else

      if(nRV.gt.0) RV_sim=simRV%RV
      if(nTTs.gt.0)then
        do ibd=1,n_body-1
          nT0=simT0(ibd)%nT0
          if(nT0.gt.0) T0_sim(1:nT0,ibd+1)=simT0(ibd)%T0
        end do
      end if

    end if

    if(allocated(ra0)) deallocate(ra0,ra1)

    call deallocate_dataRV(simRV)
    if(nTTs.gt.0)then
      do ibd=1,n_body-1
        call deallocate_dataT0(simT0(ibd))
      end do
    end if
      
    return
  end subroutine orbits_to_data


  ! ------------------------------------------------------------------ !
  ! subroutine to write file and to screen the data ... what it writes
  ! depends on the option isim and wrtid/lmon in arg.in file
  subroutine ode_out(cpuid,isim,wrtid,allpar,par,resw,fit_scale,gls_scale,to_screen)
    integer,intent(in)::cpuid,isim,wrtid
    real(dp),dimension(:),intent(in)::allpar,par
    real(dp),dimension(:),intent(out)::resw
    logical,optional,intent(in)::to_screen
    real(dp),optional,intent(out)::fit_scale,gls_scale

    real(dp),dimension(:),allocatable::m,R,P,sma,ecc,w,mA,inc,lN
    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp)::dt1,dt2
    logical::Hc

    type(dataRV)::simRV
    type(dataT0),dimension(:),allocatable::simT0

    logical::checkpar,gls_check

    integer::ndata,nRV,nTTs,nDurs,ibd,nT0

    real(dp)::chi2r_RV,chi2r_T0,chi2r_dur,chi2r_oc
    real(dp)::chi2wr_RV,chi2wr_T0,fitness,w_chi2r
    real(dp),dimension(:),allocatable::resw_temp

    ! units and file names to store and write to files
    integer::uorb,ucon
    character(512)::florb,fmorb,flcon,fmcon,fmele
    integer,dimension(:),allocatable::uele,utra
    character(512),dimension(:),allocatable::flele,fltra

    write(*,'(a)')''
    write(*,'(a)')" EXECUTING SIMPLE INTEGRATION AND WRITING FINAL FILES"
    write(*,'(a)')''

    ndata=obsData%ndata
    nRV=obsData%obsRV%nRV
    nTTs=obsData%nTTs
    nDurs=obsData%nDurs

    resw=zero
    Hc=.true.
    checkpar=.true.

    write(*,'(a)')' fitting parameters'
    write(*,'(a)')trim(paridlist)
    write(*,'(1000(1x,es23.16))')par

    allocate(m(NB),R(NB),P(NB),sma(NB),ecc(NB),w(NB),mA(NB),inc(NB),lN(NB),clN(NB))

    if(present(fit_scale))then
      fit_scale=one
      gls_scale=one
      call convert_parameters_scale(allpar,par,m,R,P,sma,ecc,w,mA,inc,lN,fit_scale)
    else
      call convert_parameters(allpar,par,m,R,P,sma,ecc,w,mA,inc,lN,checkpar)
      if(.not.checkpar)then
        write(*,'(a)')' checkpar after convert_parameters is FALSE ... BAD'
      if(ndata > 0) resw=set_max_residuals(ndata)
      end if
    end if


    if(present(fit_scale)) write(*,'(a,es23.16)')' fit_scale = ',fit_scale
    flush(6)

    ! write orbital elements into a file
    call outElements(isim,wrtid,m,R,P,sma,ecc,w,mA,inc,lN)

    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))

    ! NEW VERSION 2017-11-21
    call kepelements2statevector(m,sma,ecc,mA,w,inc,lN,ra0)
    ra1=ra0

    if(nRV.gt.0)then
      call init_dataRV(nRV,simRV)
      simRV%nRV=0
    end if

    if((idtra.ge.1).and.(idtra.le.NB)) call set_file_tra(cpuid,isim,wrtid,utra,fltra)

    allocate(simT0(NB-1))
    if(nTTs.gt.0)then
        do ibd=1,NB-1
          nT0=obsData%obsT0(ibd)%nT0
          call init_dataT0(nT0,simT0(ibd),durcheck)
          simT0(ibd)%nT0=0
          simT0(ibd)%nDur=0
          !  let's set error on TT and Dur to 1 sec! It is needed to avoid divisione by zero
          simT0(ibd)%eT0=1./s24h ! 1s in day
          simT0(ibd)%edur=1./60.0_dp ! 1s in min
        end do
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
        call ode_a_o(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
          &m,R,ra1,dt1,clN,simRV,simT0,Hc)
      if(Hc)then
        if(abs(dt1).le.tint)then
          call ode_a_o(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
            &m,R,ra1,dt2,clN,simRV,simT0,Hc)
        end if
      end if
    else
      call ode_a_o(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
        &m,R,ra1,dt2,clN,simRV,simT0,Hc)
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

    if(ndata.gt.0)then

      if(.not.Hc) then
        resw=set_max_residuals(ndata)

        if(.not.allocated(simRV%gamma)) &
          &allocate(simRV%gamma(obsData%obsRV%nRVset,2))
        simRV%gamma=zero
      else
        ! compute RV trend
        if(rv_trend_order.gt.0)then
          call addRVtrend(simRV%jd, par(nfit-rv_trend_order:nfit), simRV%trend)
        end if
        call setval_1(obsData,simRV,simT0,resw,oc_fit)

        if(present(fit_scale))then
          if(simRV%nRV.gt.0) call check_periodogram_scale(obsData%obsRV%jd,&
            &obsData%obsRV%RV-simRV%RV-simRV%trend,obsData%obsRV%eRV,&
            P,gls_check,gls_scale)
          write(*,*)
          write(*,'(a)')' CHECK GLS RESIDUALS-RV'
          write(*,'(a,es23.16)')' gls_scale = ',gls_scale
          write(*,*)
          flush(6)
          resw=resw*fit_scale*gls_scale
        else
          if(simRV%nRV.gt.0) call check_periodogram(obsData%obsRV%jd,&
            &obsData%obsRV%RV-simRV%RV-simRV%trend,obsData%obsRV%eRV,&
            &P,gls_check)
        end if
        call check_max_residuals(resw,ndata)
      end if

      ! Reduced Chi Squares for report/summary
      chi2r_RV=zero
      chi2r_T0=zero
      chi2r_dur=zero
      chi2wr_RV=zero
      chi2wr_T0=zero

      if(simRV%nRV.gt.0)then

        write(*,*)""
        write(*,'(a,i4)')" RADIAL VELOCITIES found: ",simRV%nRV
        write(*,*)""

        if(present(to_screen)) call write_RV(simRV)
        call write_RV(cpuid,isim,wrtid,simRV)

        allocate(resw_temp(nRV))
        call set_RV_resw(obsData%obsRV,simRV,resw_temp)
        chi2r_RV=sum(resw_temp*resw_temp)*obsData%inv_dof
        w_chi2r=real(ndata,dp)/real(nRV,dp)
        chi2wr_RV=chi2r_RV*w_chi2r
        deallocate(resw_temp)
        write(*,'(a)')' DEBUG: set_RV_resw, chi2r_RV, deallocated resw_temp'
        ! 2016-04-08: added gls check
        if(present(gls_scale))then
          call check_and_write_periodogram(cpuid,isim,wrtid,&
            &obsData%obsRV%jd,obsData%obsRV%RV-simRV%RV-simRV%trend,&
            &obsData%obsRV%eRV,P,gls_check,gls_scale)
        else
          call check_and_write_periodogram(cpuid,isim,wrtid,&
            &obsData%obsRV%jd,obsData%obsRV%RV-simRV%RV-simRV%trend,&
            &obsData%obsRV%eRV,P,gls_check)
        end if
        if(.not.gls_check)resw=set_max_residuals(ndata)

        write(*,'(a)')' DEBUG: gls check'

      else

        write(*,*)
        write(*,'(a)')' RADIAL VELOCITIES NOT FOUND'
        write(*,*)

      end if

      if(allocated(simT0).and.sum(simT0(:)%nT0).gt.0)then
        write(*,*)
        write(*,'(a,i5)')" T0 SIM found ",sum(simT0(:)%nT0)
        write(*,*)
        if(present(to_screen)) call write_T0(simT0)
        call write_T0(cpuid,isim,wrtid,simT0)

        allocate(resw_temp(nTTs))
        resw_temp=zero
        call set_T0_resw(obsData%obsT0,simT0,resw_temp)
        chi2r_T0=sum(resw_temp*resw_temp)*obsData%inv_dof
        w_chi2r=real(ndata,dp)/real(nTTs,dp)
        chi2wr_T0=chi2r_T0*w_chi2r
        resw_temp=zero
        call set_oc_resw(obsData%obsT0,simT0,resw_temp)
        chi2r_oc=sum(resw_temp*resw_temp)*obsData%inv_dof
        ! duration
        if(durcheck.eq.1)then
          resw_temp=zero
          call set_dur_resw(obsData%obsT0,simT0,resw_temp)
          chi2r_dur=sum(resw_temp*resw_temp)*obsData%inv_dof
        end if
        deallocate(resw_temp)
      else
        chi2r_T0=zero
        chi2r_oc=zero
        chi2r_dur=zero
        write(*,*)
        write(*,'(a)')' TRANSIT TIMES NOT FOUND'
        write(*,*)
      end if

      fitness = sum(resw*resw)
      write(*,'(a)')' DEBUG: fitness'

      if(present(to_screen))then
        if(present(fit_scale))then
          call write_fitness_summary(cpuid,isim,wrtid,chi2r_RV,chi2wr_RV,chi2r_T0,&
            &chi2r_dur,chi2wr_T0,chi2r_oc,fitness,&
            &fit_scale,gls_scale,to_screen=to_screen)
        else
          call write_fitness_summary(cpuid,isim,wrtid,chi2r_RV,chi2wr_RV,chi2r_T0,&
            &chi2r_dur,chi2wr_T0,chi2r_oc,fitness,to_screen=to_screen)
        end if
      else
        if(present(fit_scale))then
          call write_fitness_summary(cpuid,isim,wrtid,chi2r_RV,chi2wr_RV,chi2r_T0,&
            &chi2r_dur,chi2wr_T0,chi2r_oc,fitness,fit_scale,gls_scale)
        else
          call write_fitness_summary(cpuid,isim,wrtid,chi2r_RV,chi2wr_RV,chi2r_T0,&
            &chi2r_dur,chi2wr_T0,chi2r_oc,fitness)
        end if
      end if

    else ! ndata <= 0

      write(*,'(a)')' ndata == 0: NO FITNESS SUMMARY (SCREEN AND FILES)'
      flush(6)

    end if  ! ndata

    write(*,'(a)')' DEBUG: deallocating simRV'
    call deallocate_dataRV(simRV)
    write(*,'(a)')' DEBUG: deallocating simT0'
    if(nTTs.gt.0)then
      do ibd=1,NB-1
        call deallocate_dataT0(simT0(ibd))
      end do
    end if

    write(*,'(a)')' DEBUG: ode_out deallocated simRV and simT0...now deallocating kep/ra'

    if(allocated(m))   deallocate(m,R,P,sma,ecc,w,mA,inc,lN,clN)
    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_out

  ! ------------------------------------------------------------------ !

    subroutine ode_integrates(cpuid,isim,wrtid,m,R,P,sma,ecc,w,mA,inc,lN)
    integer,intent(in)::cpuid,isim,wrtid
    real(dp),dimension(:),intent(in)::m,R,P,sma,ecc,w,mA,inc,lN

    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1
    real(dp)::dt1,dt2
    logical::Hc

    ! units and file names to store and write to files
    integer::uorb,ucon
    character(512)::florb,fmorb,flcon,fmcon,fmele
    integer,dimension(:),allocatable::uele,utra
    character(512),dimension(:),allocatable::flele,fltra

    write(*,'(a)')''
    write(*,'(a)')" EXECUTING SIMPLE INTEGRATION AND WRITING FINAL FILES"
    write(*,'(a)')''

!     resw=zero
    Hc=.true.

    ! write orbital elements into a file
    call outElements(isim,wrtid,m,R,P,sma,ecc,w,mA,inc,lN)

    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))

    ! NEW VERSION 2017-11-21
    call kepelements2statevector(m,sma,ecc,mA,w,inc,lN,ra0)
    ra1=ra0

    if((idtra.ge.1).and.(idtra.le.NB)) call set_file_tra(cpuid,isim,wrtid,utra,fltra)

    fmorb=trim(adjustl(fmtorbit()))
    if(wrtorb.eq.1) call set_file_orb(cpuid,isim,wrtid,uorb,florb)
    fmcon=trim(adjustl(fmtconst()))
    if(wrtconst.eq.1) call set_file_con(cpuid,isim,wrtid,ucon,flcon)
    fmele=fmtele()
    if(wrtel.eq.1) call set_file_elem(cpuid,isim,wrtid,uele,flele)

    write(*,'(a)')' INITIALISED ORBIT AND OUTPUT FILES'
    write(*,'(a)')' RUNNING INTEGRATION ...'
    flush(6)

    dt1=tstart-tepoch
    dt2=dt1+tint
    if(dt1.lt.zero)then
      call ode_a_orbit(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
          &m,R,ra1,dt1,clN,Hc)
      if(abs(dt1).le.tint)then
        call ode_a_orbit(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
            &m,R,ra1,dt2,clN,Hc)
      end if
    else
      call ode_a_orbit(uorb,ucon,uele,utra,fmorb,fmcon,fmele,&
          &m,R,ra1,dt2,clN,Hc)
    end if

    write(*,'(a)')' COMPLETED'
    if(.not.Hc) then
      write(*,'(a)')'WARNING'
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

    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_integrates
  ! ------------------------------------------------------------------ !

! ==============================================================================
! ==============================================================================
! ORBITAL PARAMETERS TO ALL TRANSIT TIMES OF ALL PLANETS AND RV MODEL
! NO T0 FIT/DATA, NO RV FIT/DATA, NO OUTPUT FILES
! ==============================================================================
! ==============================================================================

  ! ------------------------------------------------------------------ !
  subroutine ode_b_orbit(m,R,rin,time_int,wrt_time,clN,&
    &last_tra,ttra_full,dur_full,id_ttra_full,stats_ttra,&
    &last_rv,time_rv_nmax,rv_nmax,stats_rv,&
    &Hc)
    real(dp),dimension(:),intent(in)::m,R,rin
    real(dp),intent(in)::time_int,wrt_time
    integer,dimension(:),intent(in)::clN

    ! transit times variables to be updated!!
    integer,intent(inout)::last_tra ! if first call it is zero, otherwise it is the last transit position
    real(dp),dimension(:),intent(inout)::ttra_full,dur_full
    integer,dimension(:),intent(inout)::id_ttra_full
    logical,dimension(:),intent(inout)::stats_ttra

    ! radial velocities variables to be updated!!
    integer,intent(inout)::last_rv
    real(dp),dimension(:),intent(inout)::time_rv_nmax,rv_nmax
    logical,dimension(:),intent(inout)::stats_rv

    logical,intent(out)::Hc

    integer::ntra_full

    real(dp)::step_rv,rv_temp,ttra_temp,dur_tra_temp
    logical::check_ttra
    integer::nrv_max

    real(dp),dimension(:),allocatable::dr,r1,r2,err
    integer,dimension(:),allocatable::X,Y,Z
    real(dp),dimension(:),allocatable::cX,cY,cR,rmean
    real(dp)::hw,hok,hnext,iter_time,iter_write,step_write
    integer::j,step_num

    ntra_full=size(ttra_full)
    nrv_max=size(rv_nmax)

    Hc=.true.
!     if(do_hill_check) Hc=mutual_Hill_check(m,rin)
    Hc = separation_mutual_Hill_check(m,R,rin,do_hill_check)
    if(.not.Hc) return

    ! init the state vector
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
    rmean=9.e9_dp

    allocate(dr(NBDIM),r1(NBDIM),r2(NBDIM),err(NBDIM))

    ! set initial stepsize to the value in the arg.in file
    hw=step_0
    if(time_int.lt.zero) hw=-hw ! reverse if backward integration

    ! set the iteration time: it will be updated with each step
    iter_time=zero

    ! set the initial state vector
    r1=rin
    r2=zero
    err=zero

    step_num=0 ! step counter

    ! set the iteration write time: updated each wrt_time passed
    step_write=wrt_time
    if(time_int.lt.zero) step_write=-step_write
    iter_write=step_write

    integration: do
      step_num=step_num+1

      if(abs(iter_time+hw).gt.abs(time_int)) hw=time_int-iter_time ! if last step exceeds integration time create new stepsize
      call eqmastro(m,r1,dr) ! computes the eq. of motion
      call int_rk_a(m,r1,dr,hw,hok,hnext,r2,err) ! computes the next orbit step

      Hc = separation_mutual_Hill_check(m,R,r2,do_hill_check)
      if(.not.Hc) return

      ! check if it passes the iter_write and compute the rv at proper time and update last_rv
      if(abs(iter_time+hok).ge.(abs(iter_write)))then

        rvloop: do

          last_rv=last_rv+1
          if(last_rv.gt.nrv_max)then
                Hc=.false.
                return
          end if
!           computes the proper step
          step_rv=iter_write-iter_time
          call calcRV(m,r1,dr,step_rv,rv_temp) ! computes the rv
!           save time rv status
          time_rv_nmax(last_rv)=iter_write
          rv_nmax(last_rv)=rv_temp
          stats_rv(last_rv)=.true.
!           update next writing time
          iter_write=iter_write+step_write
          if(abs(iter_write).gt.abs(iter_time+hok)) exit rvloop

        end do rvloop

      end if

      ! transit times!!
      ! loop on planets -> checks for all planets
      do j=2,NB
        ! transit time check variables
        cX(j)=r1(X(j))*r2(X(j))
        cY(j)=r1(Y(j))*r2(Y(j))

        ! check the transit criteria for each body

        check_ttra=.false. ! initialise to .false.
        ttra_temp=-9.e10_dp     !               zero
        dur_tra_temp=-9.e10_dp  !               zero

        if(clN(j).eq.0)then ! condition to check X

          if((cX(j).le.zero).and.(r1(Z(j)).gt.zero))then
            call transit_time(j,m,R,r1,r2,iter_time,hok,ttra_temp,dur_tra_temp,check_ttra)
            if(check_ttra)then
              last_tra=last_tra+1
              if(last_tra.gt.ntra_full)then
                Hc=.false.
                return
              end if
              ttra_full(last_tra)=ttra_temp
              dur_full(last_tra)=dur_tra_temp
              id_ttra_full(last_tra)=j
              stats_ttra(last_tra)=.true.
            end if
          end if

        else ! condition to check Y

          if((cY(j).le.zero).and.(r1(Z(j)).gt.zero))then
            call transit_time(j,m,R,r1,r2,iter_time,hok,ttra_temp,dur_tra_temp,check_ttra)
            if(check_ttra)then
              last_tra=last_tra+1
              if(last_tra.gt.ntra_full)then
                Hc=.false.
                return
              end if
              ttra_full(last_tra)=ttra_temp
              dur_full(last_tra)=dur_tra_temp
              id_ttra_full(last_tra)=j
              stats_ttra(last_tra)=.true.
            end if
          end if

        end if ! end condition X,Y

!         end if ! end transit criteria

      end do

      ! update iteration time with the stepsize used (hok)
      iter_time=iter_time+hok

      ! check end of integration: time_int reached
      if(abs(iter_time).ge.abs(time_int)) exit integration
      ! update step and state vector
      hw=hnext
      r1=r2

    end do integration

    deallocate(X,Y,Z,cX,cY,cR,rmean)
    deallocate(dr,r1,r2,err)

    return
  end subroutine ode_b_orbit
  ! ------------------------------------------------------------------ !



  ! ------------------------------------------------------------------ !
  subroutine ode_all_ttra_rv(wrt_time,m,R,P,a,e,w,mA,inc,lN,&
    &ttra_full,dur_full,id_ttra_full,stats_ttra,&
    &time_rv_nmax,rv_nmax,stats_rv)
    real(dp),intent(in)::wrt_time
    real(dp),dimension(:),intent(in)::m,R,P,a,e,w,mA,inc,lN

    real(dp),dimension(:),intent(out)::ttra_full,dur_full
    integer,dimension(:),intent(out)::id_ttra_full
    logical,dimension(:),intent(out)::stats_ttra
    real(dp),dimension(:),intent(out)::time_rv_nmax,rv_nmax
    logical,dimension(:),intent(out)::stats_rv

    integer,dimension(:),allocatable::clN
    real(dp),dimension(:),allocatable::ra0,ra1

    integer::last_tra,last_rv

    real(dp)::dt1,dt2
    logical::Hc

    Hc=.true.

    ! it is needed to define the which is the alarm coordinate for the transit detection
    call lNset(lN,clN)

    NBDIM=6*NB
    if(.not.allocated(ra0)) allocate(ra0(NBDIM),ra1(NBDIM))

    ! NEW VERSION 2017-11-21
    call kepelements2statevector(m,a,e,mA,w,inc,lN,ra0)
    ra1=ra0

    ! prepare rv arrays
    time_rv_nmax=zero
    rv_nmax=zero
    stats_rv=.false.
    last_rv=0

    ! prepare ttra arrays
    ttra_full=-9.e10_dp
    dur_full=-9.e10_dp
    id_ttra_full=0
    stats_ttra=.false.
    last_tra=0

    dt1=tstart-tepoch
    dt2=dt1+tint

    if(dt1.lt.zero)then

      ! backward integration
      call ode_b_orbit(m,R,ra1,dt1,wrt_time,clN,&
        &last_tra,ttra_full,dur_full,id_ttra_full,stats_ttra,&
        &last_rv,time_rv_nmax,rv_nmax,stats_rv,Hc)
      if(.not.Hc)then
        if(allocated(ra0)) deallocate(ra0,ra1)
        return
      end if

      if(abs(dt1).le.tint)then
        ! forward integration
        call ode_b_orbit(m,R,ra1,dt2,wrt_time,clN,&
          &last_tra,ttra_full,dur_full,id_ttra_full,stats_ttra,&
          &last_rv,time_rv_nmax,rv_nmax,stats_rv,Hc)
        if(.not.Hc)then
          if(allocated(ra0)) deallocate(ra0,ra1)
          return
        end if
      end if

    else

      ! only forward integration
      call ode_b_orbit(m,R,ra1,dt2,wrt_time,clN,&
        &last_tra,ttra_full,dur_full,id_ttra_full,stats_ttra,&
        &last_rv,time_rv_nmax,rv_nmax,stats_rv,Hc)

    end if

    if(allocated(ra0)) deallocate(ra0,ra1)

    return
  end subroutine ode_all_ttra_rv
  ! ------------------------------------------------------------------ !

end module ode_run
