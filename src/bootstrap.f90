module bootstrap
  use constants,only:dp,sprec,zero,one
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

    real(dp),dimension(:),allocatable::RVok,RVboot
    real(dp),dimension(:,:),allocatable::T0ok,T0boot

    real(dp),dimension(:),allocatable::b_par,b_allpar
    real(dp),dimension(:,:),allocatable::storeboot

    real(dp),dimension(:),allocatable::resw,diag,sig
    integer,dimension(:),allocatable::iwa
    real(dp)::fitness,scale_errorbar
    integer::info,nTs,iT0

    real(dp),dimension(:,:),allocatable::percentile
    integer,parameter::nperc=7
    integer::iboot,inb,cpu

    allocate(resw(ndata))
    allocate(b_par(nfit),b_allpar(npar),diag(nfit),sig(nfit),iwa(nfit))
    if(nRV.gt.0) allocate(RVok(nRV),RVboot(nRV))
    nTs=maxval(nT0)
    if(nTs.gt.0) allocate(T0ok(nTs,NB),T0boot(nTs,NB))

    !1. call subroutine to retrieve RVok and T0ok (ode_parok)
    call ode_parok(allpar,parok,RVok,T0ok,resw,fitness)

    call set_bootfile(isim)
    allocate(storeboot(nboot,nfit),percentile(nperc,nfit))

    scale_errorbar=one
    if(bootstrap_scaling) scale_errorbar=sqrt(fitness)
    write(*,*)' Fitness value (Chi2r*k_chi2r + Chi2wr*k_chi2wr) = ', fitness, ' dof = ', dof
    write(*,'(a,l,a,f20.12)')' Bootstrap scaling (',bootstrap_scaling,') = ',&
      &scale_errorbar
    write(*,*)
    flush(6)
    
    call init_random_seed_clock(nboot)
    !2.  run nboot iterations

    !$omp parallel do default(private)&
    !$omp& shared(allpar, parok, ndata, nfit, eT0obs, eRVobs, storeboot,&
    !$omp& scale_errorbar, nboot, RVok, T0ok, NB, nRV, nT0, nTs)
    bootstrap: do iboot=1,nboot
      cpu=1
      !$ cpu = omp_get_thread_num() + 1
      write(*,'(2(a,i5),a)')" CPU ",cpu,&
            &" BOOTSTRAP ITERATION NUMBER ",iboot," ... "
      !2.a generation of ndata (RVboot and T0boot) with Gaussian distribution N(mean,var)
      !    where mean = data (RVok or T0ok) and var = square error data (eRVobs, eT0obs)
      if(nRV.gt.0)then
          RVboot=zero
          call gaussdev(RVboot)
          RVboot = RVok+RVboot*eRVobs
      end if
      if(nTs.gt.0) T0boot=zero
      do inb=2,NB
          iT0=nT0(inb)
          if(iT0.gt.0)then
            call gaussdev(T0boot(1:iT0,inb))
            T0boot(1:iT0,inb)=T0ok(1:iT0,inb)+&
                  &T0boot(1:iT0,inb)*eT0obs(1:iT0,inb)*&
                  &scale_errorbar
          end if
      end do
      
      !2.b lm_driver (_3) to retrieve temporary parameters fitting new RVboot and T0boot
      b_allpar = allpar
      b_par = parok
      resw=zero
      diag=zero
      sig=zero
      iwa=0
      info=0
      call lm_driver(ode_boot,b_allpar,ndata,nfit,b_par,&
            &RVboot,T0boot,resw,diag,sig,info,iwa)

      storeboot(iboot,:)=b_par
      write(*,'(2(a,i5))')" CPU ",cpu,&
            &" DONE bootstrap iteration ",iboot
      flush(6)

    end do bootstrap
    !$omp end parallel do

    call wrt_storeboot(isim,parok,storeboot)
    write(*,'(a)')" Written bootstrap file"

    !3.  sort all new parameters and compute the percentiles
    call calc_percentile(parok,storeboot,percentile)

    !4.  write 2 files, one with all fitted parameters, one with percentiles
    call wrt_percentile(isim,percentile)

    deallocate(storeboot)
    if(allocated(RVok)) deallocate(RVok,RVboot)
    if(allocated(T0ok)) deallocate(T0ok,T0boot)
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
  subroutine wrt_storeboot(isim,parok,storeboot)
    use init_trades,only:get_unit
    use convert_type,only:string
    integer,intent(in)::isim
    real(dp),dimension(:),intent(in)::parok
    real(dp),dimension(:,:),intent(in)::storeboot
    real(dp),dimension(:),allocatable::storevec
    real(dp),dimension(2)::val,valabs
    integer::cpuid,uboot,iboot,ifit
    character(512)::flboot,fmt

    allocate(storevec(nfit))
    storevec=zero
    fmt=adjustl("(i6,1x,100000"//trim(sprec)//")")
    flboot=trim(path)//trim(adjustl(string(isim)))//"_bootstrap_sim.dat"
    flboot=trim(adjustl(flboot))
    cpuid=1
    !$ cpuid=omp_get_thread_num()+1
    uboot=get_unit(cpuid)
    open(uboot,file=trim(flboot))
    write(uboot,'(a,a)')"# nboot "//trim(paridlist)
    do iboot = 1,nboot
      storevec=storeboot(iboot,:)
      do ifit=1,nfit
          if((parid(ifit)(1:1).eq.'w').or.&
              &(parid(ifit)(1:2).eq.'mA').or.&
              &(parid(ifit)(1:2).eq.'lN'))then
            val(1)=storevec(ifit)-parok(ifit)
            val(2)=storevec(ifit)-360._dp-parok(ifit) 
            valabs=abs(val)
            storevec(ifit)=val(minloc(valabs,dim=1))+parok(ifit)
          end if
      end do
      write(uboot,trim(fmt))iboot,storevec
    end do
    close(uboot)
    deallocate(storevec)

    return
  end subroutine wrt_storeboot

  ! it computes the distribution of the parameters
  ! and it returns the percentiles
  subroutine calc_percentile(parok,store,perc)
    use sorting,only:sort
    real(dp),dimension(:),intent(in)::parok
    real(dp),dimension(:,:),intent(in)::store
    real(dp),dimension(:,:),intent(out)::perc

    real(dp),dimension(:),allocatable::vec
    ! percentiles of a distribution, they are not the percentiles of the variance distribution 
    real(dp),parameter::p1=0.13_dp/100._dp, p2=2.28_dp/100._dp,&
                      &p3=15.87_dp/100._dp, p4=0.5_dp,& ! this the median
                      &p5=84.13_dp/100._dp, p6=97.72_dp/100._dp,&
                      &p7=99.87_dp/100._dp
    integer,dimension(7)::intperc
    real(dp),dimension(2)::val,valabs
    integer::ifit,iper,iboot,nperc
    
    nperc=size(perc(:,1))

    intperc(1)=int(p1*nboot+0.5_dp) !  0.13%
    if(intperc(1).eq.0) intperc(1)=1
    intperc(2)=int(p2*nboot+0.5_dp) !  2.28%
    intperc(3)=int(p3*nboot+0.5_dp) ! 15.87%
    intperc(4)=int(p4*nboot+0.5_dp) ! 50.00%
    intperc(5)=int(p5*nboot+0.5_dp) ! 84.13%
    intperc(6)=int(p6*nboot+0.5_dp) ! 97.72%
    intperc(7)=int(p7*nboot+0.5_dp) ! 99.87%
    if(intperc(1).eq.1) intperc(7)=nboot
    allocate(vec(nboot))
    write(*,*)"ifit iper inteperc(iper) vec(intperc(iper)) perc(iper,ifit)"
    do ifit=1,nfit
      vec=store(:,ifit)
      if((parid(ifit)(1:1).eq.'w').or.&
            &(parid(ifit)(1:2).eq.'mA').or.&
            &(parid(ifit)(1:2).eq.'lN'))then
          do iboot=1,nboot
            val(1)=vec(iboot)-parok(ifit)
            val(2)=vec(iboot)-360._dp-parok(ifit)
            valabs=abs(val)
            vec(iboot)=val(minloc(valabs,dim=1))+parok(ifit)
          end do
      end if
      call sort(vec)
      do iper=1,nperc
          perc(iper,ifit)=vec(intperc(iper))
          write(*,*)ifit,iper,intperc(iper),vec(intperc(iper)),perc(iper,ifit)
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

    fmt=adjustl("(10000"//trim(sprec)//")")
    flboot=trim(path)//trim(adjustl(string(isim)))//"_perc_boot.dat"
    flboot=trim(adjustl(flboot))
    cpuid=1
    !$ cpuid=omp_get_thread_num()+1
    uboot=get_unit(cpuid)
    open(uboot,file=trim(flboot))
    write(*,*)" writing ",trim(flboot)
    do ip=1,7
      write(uboot,trim(fmt))perc(ip,:)
      write(*,'(i5,1x,10000(f25.14,1x))')ip,perc(ip,:)
    end do
    close(uboot)

    return
  end subroutine wrt_percentile


end module bootstrap

