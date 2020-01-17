module bootstrap
  use constants,only:dp,sprec,zero,one,half
  use custom_type
  use parameters
  use Levenberg_Marquardt,only:lm_driver
  use gaussian
  use ode_run,only:ode_parok,ode_boot
  !$ use omp_lib
  implicit none

  contains

  ! bootstrap analysis
  ! given fitted parameters (parok) it computes the 'right' RV and T0 (call ode_parok)
  ! and then it runs nboot simulations (and fitting) for nboot series of ndata
  ! as a Gaussian N(mean,var) with mean the RVok and T0ok
  ! and var the square of the errors of RV and T0 observed
  subroutine strap_driver(isim,allpar,parok)
!     use init_trades,only:init_random_seed
    use random_trades
    integer,intent(in)::isim
    real(dp),dimension(:),intent(in)::allpar,parok

    real(dp),dimension(:),allocatable::Ms_boot,Rs_boot
    
!     real(dp),dimension(:),allocatable::RVok,RVboot
!     real(dp),dimension(:,:),allocatable::T0ok,T0boot
    type(dataRV)::okRV,bootRV
    type(dataT0),dimension(:),allocatable::okT0,bootT0
    type(dataObs)::bootData

    real(dp),dimension(:),allocatable::b_par,b_allpar
    real(dp),dimension(:,:),allocatable::storeboot

    real(dp),dimension(:),allocatable::resw,diag,sig
    integer,dimension(:),allocatable::iwa
    real(dp)::fitness,scale_errorbar
    integer::info,nRV,nTTs,nT0

    real(dp),dimension(:,:),allocatable::percentile
    integer,parameter::nperc=7
    integer::iboot,inb,cpu,chunk,ndata

    ndata=obsData%ndata
    allocate(resw(ndata))
    allocate(b_par(nfit),b_allpar(npar),diag(nfit),sig(nfit),iwa(nfit))
    
!     if(nRV.gt.0) allocate(RVok(nRV),RVboot(nRV))
    nRV=obsData%obsRV%nRV
!     call init_dataRV(nRV,okRV) ! it will be done within ode_parok
    call init_dataRV(nRV,bootRV)
    call init_dataRV(nRV,bootData%obsRV)
    
    nTTs=obsData%nTTs
!     if(nTTs.gt.0) allocate(T0ok(nTTs,NB),T0boot(nTTs,NB))
    if(nTTS.gt.0)then
      allocate(okT0(NB-1),bootT0(NB-1),bootData%obsT0(NB-1))
      do inb=1,NB-1
        nT0=obsData%obsT0(inb)%nT0
!         call init_dataT0(nT0,okT0(inb),durcheck) ! it will be done within ode_parok
        call init_dataT0(nT0,bootT0(inb),durcheck)
        call init_dataT0(nT0,bootData%obsT0(inb),durcheck)
      end do
    end if
    
    
    write(*,'(a)')''
    write(*,'(a,i6,a)')' BOOTSTRAP/MC ANALYSIS WITH ',nboot,' ITERATIONS'
    flush(6)
    
    !1. call subroutine to retrieve RVok and T0ok (ode_parok)
!     call ode_parok(allpar,parok,RVok,T0ok,resw,fitness)
    call ode_parok(allpar,parok,okRV,okT0,resw,fitness)

    call set_bootfile(isim)
    allocate(storeboot(nboot,nfit),percentile(nperc,nfit))

    scale_errorbar=one
    if(bootstrap_scaling) scale_errorbar=sqrt(fitness)
!     write(*,*)' Fitness value (Chi2r*k_chi2r + Chi2wr*k_chi2wr) = ', fitness, ' dof = ', dof
!     write(*,'(a,l,a,f20.12)')' Bootstrap scaling (',bootstrap_scaling,') = ',&
!       &scale_errorbar
!     write(*,*)
!     flush(6)
    
    call init_random_seed_clock(nboot)
    
    ! allocate and create nboot values of the M-R of the star with the provided errors
    allocate(Ms_boot(nboot),Rs_boot(nboot))
    call gaussdev(Ms_boot)
    Ms_boot = MR_star(1,1)+Ms_boot*MR_star(1,2)
    call gaussdev(Rs_boot)
    Rs_boot = MR_star(2,1)+Rs_boot*MR_star(2,2)
    call wrt_gaussian_star(isim,Ms_boot,Rs_boot)
    
    
    !2.  run nboot iterations

    chunk=1
!     !$ chunk=int(0.5_dp*(nboot/omp_get_num_threads()))

    !$omp parallel do schedule(DYNAMIC,chunk)&
    !$omp& default(private) NUM_THREADS(ncpu_in)&
    !$omp& shared(allpar, parok, ndata, nfit, obsData, storeboot,&
    !$omp& scale_errorbar, nboot, bootRV, bootDAta, bootT0, NB, nRV, &
    !$omp& Ms_boot, Rs_boot)
    
    bootstrap: do iboot=1,nboot

      cpu=1
      !$ cpu = omp_get_thread_num() + 1
!       write(*,'(2(a,i5),a)')" CPU ",cpu,&
!             &" BOOTSTRAP ITERATION NUMBER ",iboot," ... "
            
      !2.a generation of ndata (RVboot and T0boot) with Gaussian distribution N(mean,var)
      !    where mean = data (RVok or T0ok) and var = square error data (eRVobs, eT0obs)
      if(nRV.gt.0)then
!           RVboot=zero
!           call gaussdev(RVboot)
!           RVboot = RVok+RVboot*eRVobs*scale_errorbar
            call gaussdev(bootRV%RV)
            bootRV%RV=bootRV%RV*bootRV%eRV*scale_errorbar+okRV%RV
            bootData%obsRV=bootRV
      end if
      
!       if(nTTs.gt.0) T0boot=zero
      
      do inb=1,NB-1
!           nT0=nT0(inb)
!           if(nT0.gt.0)then
!             call gaussdev(T0boot(1:nT0,inb))
!             T0boot(1:nT0,inb)=T0ok(1:nT0,inb)+&
!                   &T0boot(1:nT0,inb)*eT0obs(1:nT0,inb)*&
!                   &scale_errorbar
!           end if
        nT0=okT0(inb)%nT0
        if(nT0.gt.0)then
          call gaussdev(bootT0(inb)%T0)
          bootT0(inb)%T0=bootT0(inb)%T0*bootT0(inb)%eT0*scale_errorbar+&
            &okT0(inb)%T0
        end if
      end do
      bootData%obsT0=bootT0
      
      !2.b lm_driver (_3) to retrieve temporary parameters fitting new RVboot and T0boot
      b_allpar = allpar
      ! use the Ms/Rs_boot star values
      b_allpar(1) = Ms_boot(iboot)
      b_allpar(2) = Rs_boot(iboot)
      
      b_par = parok
      
      resw=zero
      diag=zero
      sig=zero
      iwa=0
      info=0
      
      call lm_driver(ode_boot,b_allpar,ndata,nfit,b_par,&
            &bootData,resw,diag,sig,info,iwa)

      storeboot(iboot,:)=b_par
      write(*,'(2(a,i5))')" CPU ",cpu,&
            &" DONE bootstrap iteration ",iboot
      flush(6)

    end do bootstrap
    !$omp end parallel do

!     call wrt_storeboot(isim,parok,storeboot)
    call wrt_storeboot(isim,storeboot)
!     write(*,'(a)')" Written bootstrap file"

    !3.  sort all new parameters and compute the percentiles
!     call calc_percentile(parok,storeboot,percentile)
    call calc_percentile(storeboot,percentile)

    !4.  write 2 files, one with all fitted parameters, one with percentiles
    call wrt_percentile(isim,percentile)

    deallocate(storeboot)
    deallocate(Ms_boot,Rs_boot)
!     if(allocated(RVok)) deallocate(RVok,RVboot)
!     if(allocated(T0ok)) deallocate(T0ok,T0boot)
    call deallocate_dataRV(okRV)
    call deallocate_dataRV(bootRV)
    
    do inb=1,NB-1
      call deallocate_dataT0(bootT0(inb))
    end do

    call deallocate_dataObs(bootData)
    
    deallocate(resw)

    return
  end subroutine strap_driver

  ! it prepares files for bootstrap analysis
  subroutine set_bootfile(isim)
    use init_trades,only:get_unit
    use convert_type,only:string
    integer,intent(in)::isim
    integer::cpuid,uboot
    character(512)::flboot
    logical::fstat

    flboot=trim(path)//trim(adjustl(string(isim)))//"_bootstrap_sim.dat"
    flboot=trim(adjustl(flboot))
    inquire(file=trim(flboot),exist=fstat)
    cpuid=1
    !$ cpuid=omp_get_thread_num()+1
    uboot=get_unit(cpuid)
    if(fstat)then
      open(uboot,file=trim(flboot),status='replace')
    else
      open(uboot,file=trim(flboot),status='new')
    end if
    close(uboot)

    return
  end subroutine set_bootfile

  ! it writes into file the bootstrap analysis
!   subroutine wrt_storeboot(isim,parok,storeboot)
  subroutine wrt_storeboot(isim,storeboot)
    use init_trades,only:get_unit
    use convert_type,only:string
    integer,intent(in)::isim
!     real(dp),dimension(:),intent(in)::parok
    real(dp),dimension(:,:),intent(in)::storeboot
    real(dp),dimension(:),allocatable::storevec
!     real(dp),dimension(2)::val,valabs
    integer::cpuid,uboot,ifit,iboot
    character(512)::flboot,fmt

    allocate(storevec(nfit))
    storevec=zero
    fmt=adjustl("(i6,1x,100000(1x,"//trim(sprec)//"))")
    flboot=trim(path)//trim(adjustl(string(isim)))//"_bootstrap_sim.dat"
    flboot=trim(adjustl(flboot))
    cpuid=1
    !$ cpuid=omp_get_thread_num()+1
    uboot=get_unit(cpuid)
    open(uboot,file=trim(flboot))
    write(uboot,'(a,a)')"# nboot "//trim(paridlist)
    do iboot = 1,nboot
      storevec=storeboot(iboot,:)
!       do ifit=1,nfit
!           if((parid(ifit)(1:1).eq.'w').or.&
!               &(parid(ifit)(1:2).eq.'mA').or.&
!               &(parid(ifit)(1:2).eq.'lN'))then
!             val(1)=storevec(ifit)-parok(ifit)
!             val(2)=storevec(ifit)-360._dp-parok(ifit) 
!             valabs=abs(val)
!             storevec(ifit)=val(minloc(valabs,dim=1))+parok(ifit)
!           end if
!       end do
      write(uboot,trim(fmt))iboot,storevec
    end do
    close(uboot)
    deallocate(storevec)

    return
  end subroutine wrt_storeboot

  ! it computes the distribution of the parameters
  ! and it returns the percentiles
!   subroutine calc_percentile(parok,store,perc)
  subroutine calc_percentile(store,perc)
    use sorting,only:sort
!     real(dp),dimension(:),intent(in)::parok
    real(dp),dimension(:,:),intent(in)::store
    real(dp),dimension(:,:),intent(out)::perc

    real(dp),dimension(:),allocatable::vec
    real(dp),parameter::hundred=100._dp
    ! percentiles of a distribution, they are not the percentiles of the variance distribution 
    real(dp),parameter::p1=0.13_dp/hundred, p2=2.28_dp/hundred,&
                      &p3=15.87_dp/hundred, p4=half,& ! this the median
                      &p5=84.13_dp/hundred, p6=97.72_dp/hundred,&
                      &p7=99.87_dp/hundred
    integer,dimension(7)::intperc
!     real(dp),dimension(2)::val,valabs
    integer::ifit,iper,nperc!,iboot
    
    nperc=size(perc(:,1))

    intperc(1)=int(p1*nboot+half) !  0.13%
    if(intperc(1).eq.0) intperc(1)=1
    intperc(2)=int(p2*nboot+half) !  2.28%
    intperc(3)=int(p3*nboot+half) ! 15.87%
    intperc(4)=int(p4*nboot+half) ! 50.00%
    intperc(5)=int(p5*nboot+half) ! 84.13%
    intperc(6)=int(p6*nboot+half) ! 97.72%
    intperc(7)=int(p7*nboot+half) ! 99.87%
    if(intperc(1).eq.1) intperc(7)=nboot
    allocate(vec(nboot))
!     write(*,*)"ifit iper inteperc(iper) vec(intperc(iper)) perc(iper,ifit)"
    do ifit=1,nfit
      vec=store(:,ifit)
!       if((parid(ifit)(1:1).eq.'w').or.&
!             &(parid(ifit)(1:2).eq.'mA').or.&
!             &(parid(ifit)(1:2).eq.'lN'))then
!           do iboot=1,nboot
!             val(1)=vec(iboot)-parok(ifit)
!             val(2)=vec(iboot)-360._dp-parok(ifit)
!             valabs=abs(val)
!             vec(iboot)=val(minloc(valabs,dim=1))+parok(ifit)
!           end do
!       end if
      call sort(vec)
      do iper=1,nperc
          perc(iper,ifit)=vec(intperc(iper))
!           write(*,*)ifit,iper,intperc(iper),vec(intperc(iper)),perc(iper,ifit)
      end do
    end do
    deallocate(vec)

    return
  end subroutine calc_percentile

  ! it writes percentiles of the distribution
  ! it would be better to analyze the distribution
  ! of the parameters in the isim_boostrap.dat file
  ! if there are some angles or in case of bimodal distribution
  subroutine wrt_percentile(isim,perc)
    use init_trades,only:get_unit
    use convert_type,only:string
    integer,intent(in)::isim
    real(dp),dimension(:,:),intent(in)::perc
    integer::cpuid,uboot,ip
    character(512)::flboot,fmt

    fmt=adjustl("(10000(1x,"//trim(sprec)//"))")
    flboot=trim(path)//trim(adjustl(string(isim)))//"_perc_boot.dat"
    flboot=trim(adjustl(flboot))
    cpuid=1
    !$ cpuid=omp_get_thread_num()+1
    uboot=get_unit(cpuid)
    open(uboot,file=trim(flboot))
!     write(*,*)" writing ",trim(flboot)
    do ip=1,7
      write(uboot,trim(fmt))perc(ip,:)
!       write(*,'(i5,1x,10000(f25.14,1x))')ip,perc(ip,:)
    end do
    close(uboot)

    return
  end subroutine wrt_percentile

  subroutine wrt_gaussian_star(isim,Ms_boot,Rs_boot)
    use init_trades,only:get_unit
    use convert_type,only:string
    integer,intent(in)::isim
    real(dp),dimension(:),intent(in)::Ms_boot,Rs_boot
    integer::cpuid,uboot,ip
    character(512)::flboot,fmt

    fmt=adjustl("(10000(1x,"//trim(sprec)//"))")
    flboot=trim(path)//trim(adjustl(string(isim)))//"_gaussian_star_boot.dat"
    flboot=trim(adjustl(flboot))
    cpuid=1
    !$ cpuid=omp_get_thread_num()+1
    uboot=get_unit(cpuid)
    open(uboot,file=trim(flboot))
    
    write(uboot,'(4(a,es23.16),a)')'# Mstar=N(',MR_star(1,1),' , ',MR_star(1,2),') Rstar=N(',MR_star(2,1),' , ',MR_star(2,2),') '
!     write(*,*)" writing ",trim(flboot)
    do ip=1,nboot
      write(uboot,trim(fmt))Ms_boot(ip),Rs_boot(ip)
    end do
    close(uboot)

    return
  end subroutine wrt_gaussian_star


end module bootstrap

