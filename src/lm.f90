module Levenberg_Marquardt
  use constants
  ! Luca Borsato 2014
  ! Levenberg-Marquardt by MINPACK converted to Fortran90
  ! I left some original code commented.
  ! Some subroutines have been duplicated and modified for different purposes
  ! MINPACK routines which are used by both LMDIF & LMDER

  implicit none
  !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)

  !dpIVATE
  !PUBLIC::dp,lmdif1,lmdif,lmder1,lmder,enorm
  !PUBLIC::lmdif1,lmdif2,lmdif,lmder1,lmder,enorm

  interface lm_driver
    module procedure lm_driver_3,lm_driver_4s,lm_driver_4p
  end interface lm_driver

  interface lmdif
    module procedure lmdif_2,lmdif_3
  end interface lmdif

  interface fdjac2
    module procedure fdjac2_1,fdjac2_2,fdjac2_3
  end interface fdjac2

  interface get_best
    module procedure get_best_1,get_best_1_b,get_best_2
  end interface get_best

  private
  public::lm_driver,lmdif

  contains

  ! it updates the allpar vector from parameters fitted in par
!   subroutine allpar_update(par,allpar)
!     use parameters,only:npar,tofit
!     real(dp),dimension(:),intent(in)::par
!     real(dp),dimension(:)::allpar
!     integer::j1,j2
! 
!     j2=0
!     do j1=1,npar
!       if(tofit(j1).eq.1)then
!           j2=j2+1
!           allpar(j1)=par(j2)
!       end if
!     end do
! 
!     return
!   end subroutine allpar_update

!   ! it writes allpar vector in sequence...not used...
!   subroutine writeallpar(allpar)
!     use parameters,only:npar,tofit
!     real(dp),dimension(:),intent(in)::allpar
!     integer::j1
! 
!     do j1=1,npar
!       if(tofit(j1).eq.1)then
!           write(*,'(g25.14)',advance='no')allpar(j1)
!       end if
!     end do
!     write(*,*)""
! 
!     return
!   end subroutine writeallpar
  
  ! it calculates sigmas of the parameters from the covariance matrix
  subroutine get_covar(fjac,iwa,m,n,diag,sig,cntsig)
    real(dp),dimension(:,:),intent(in)::fjac
    integer,dimension(:),intent(in)::iwa
    integer,intent(in)::m,n
    real(dp),dimension(:),intent(out)::diag,sig
    integer,intent(out)::cntsig
    real(dp),dimension(:,:),allocatable::r
    integer,dimension(:),allocatable::ipvt
    real(dp),dimension(:),allocatable::wa2
    integer::j,ldr

    allocate(r(m,n))
    allocate(ipvt(n))
    allocate(wa2(n))

    ! calculates covariance matrix each run
    r=fjac
    ipvt=iwa
    ldr=m
    wa2=zero
    diag=zero
    sig=zero
    cntsig=0
    call covar(n,r,ldr,ipvt,TOLERANCE,wa2)
    do j=1,n
      diag(j)=r(j,j)
      sig(j)=sqrt(diag(j))
      if(sig(j).le.TOLERANCE)cntsig=cntsig+1
    end do

    deallocate(r,ipvt,wa2)

    return
  end subroutine get_covar

  ! subroutine to perform lm analysis but using a subroutine with more arguments
  ! (needed for bootstrap analysis)
  ! USED BY:
  ! BOOTSTRAP
  subroutine lm_driver_3(fcn,allpar,m,n,x,RV_obs,T0_obs,fvec,diag,sig,info,iwa)
    use parameters,only:maxfev,nprint,lmtols
    use parameters_conversion,only:param_adj
    use init_trades,only:get_unit
    use convert_type,only:string
    !$ use omp_lib
    integer,intent(IN)::m
    integer,intent(IN)::n
    real(dp),dimension(:),intent(in)::RV_obs
    real(dp),dimension(:,:),intent(in)::T0_obs
    real(dp),dimension(:),intent(inout)::allpar
    real(dp),dimension(:),intent(INOUT)::x
    real(dp),dimension(:),intent(OUT)::fvec,diag,sig
    integer,intent(OUT)::info
    integer,dimension(:),intent(OUT)::iwa

    real(dp),dimension(:),allocatable::wallpar
    !integer::it,itmax,cpuid,iflag
    integer::cpuid,iflag
    real(dp)::wftol,wxtol,wgtol,wepsfcn,fitness

    interface
      subroutine fcn(allpar,m,n,x,RV_obs,T0_obs,fvec,iflag)
        use constants
        implicit none
        integer,intent(IN)::m,n
        real(dp),dimension(:),intent(IN)::allpar,x,RV_obs
        real(dp),dimension(:,:),intent(in)::T0_obs
        real(dp),dimension(:),intent(OUT)::fvec
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface

    real(dp),parameter::factor = 100._dp

    ! added for the covar subroutine
    real(dp),dimension(:,:),allocatable::fjac
    real(dp),dimension(:),allocatable::wa
    integer::mode,cntsig,nfev

!     write(*,'(a)')' **************************'
!     write(*,'(a)')' *** CALLED lm_driver_3 ***'
!     write(*,'(a)')' **************************'
!     flush(6)
    
    info = 0
    iflag = 0 ! added by Luca in order to avoid negative value
    ! check the input parameters for errors.
    !if (n <= 0 .or. m < n .or. tol < zero) return
    if (n <= 0 .or. m < n) return

    !     call lmdif.
    allocate(fjac(m,n),wa(2*n))
    fjac=zero
    wa=zero
    allocate(wallpar(size(allpar)))
    wallpar=allpar
    cpuid = 1
    !$ cpuid = omp_get_thread_num() + 1

    ! working tolerances from the lmtols array with the rigth values 
    ! already determined in previous fitting
    wxtol = lmtols(cpuid,1)
    wftol = lmtols(cpuid,2)
    wgtol = lmtols(cpuid,3)
    wepsfcn = lmtols(cpuid,4)

!     write(*,'(a,i3,a)')" lm_driver from CPU ",cpuid,&
!         &" setted tolerances: "
!     write(*,'(a,g25.14)')   " ftol   = ",wftol
!     write(*,'(a,g25.14)')   " xtol   = ",wxtol
!     write(*,'(a,g25.14)')   " gtol   = ",wgtol
!     write(*,'(2(a,g25.14))')" epsfcn = ",wepsfcn,&
!         &" -> eps = sqrt(epsfcn) = ",sqrt(wepsfcn)
    mode = 1
!     write(*,'(a)')""
! 
!     write(*,'(a)')" LM iteration"
!     !write(*,'(2(a,g25.14))') "  setted ftol = ",wftol," xtol = ",wxtol
!     write(*,'(a,1000g25.14)')"  Input parameters:  ",x
    ! LM call
    call lmdif(fcn,wallpar,m,n,x,RV_obs,T0_obs,fvec,wftol,wxtol,wgtol,&
        &maxfev,wepsfcn,wa,&
        &mode,factor,nprint,info,nfev,fjac,iwa,wa(n+1:))

    cntsig=0 ! this counter is used to determine if a sigma in sig vector has zero value
    call get_covar(fjac,iwa,m,n,diag,sig,cntsig)
    !write(*,'(2g25.15,i3)')diag(1),sig(1),cntsig

!     call param_adj(x,sig) ! adjusting parameters, i.e., angles [0, 360] deg
    
!     write(*,'(a,1000(g25.14))')" Output parameters: ",x
!     write(*,'(a,1000(g25.14))')"  Sigma parameters: ",sig
!     fitness_x_dof = sum(fvec*fvec)
    
    fitness = sum(fvec*fvec)
!     write(*,'(a)')""
!     write(*,'(2(a,g25.14),a,i4)')" Fitness*dof = ",fitness*real(dof,dp)," Fitness = ",&
!         &fitness," dof = ",dof
    if(info == 8)then
      write(*,'(a)')"** info = 8 will be set to 4 **"
      info = 4
    end if
!     write(*,'(a,i3)')" LM info = ",info
!     call allpar_update(x,wallpar)
!     write(*,*)""
    deallocate(fjac,wa)
    deallocate(wallpar)

    return

  end subroutine lm_driver_3

  ! it initializes the all the possible values to use
  ! as epsfcn, the parameter that determines the initial Jacobian
  ! neps log steps from valMin to valMax = sqrt(valMin)
  function initEpsfcn(neps) result(lmepsfcn)
    real(dp),dimension(:),allocatable::lmepsfcn
    integer,intent(in)::neps
    real(dp)::valMin,lvalMin,valMax,lvalMax,dleps,leps
    integer::i

    if(prec.le.15)then
      valMin=epsilon(one)
    else
      valMin=epsilon(1.d0)
    end if
    lvalMin=log10(valMin)
    valMax=sqrt(valMin)
    lvalMax=log10(valMax)

    allocate(lmepsfcn(neps))
    if(neps.gt.1)then
      lmepsfcn=valMin
      dleps=(lvalMax-lvalMin)/(neps-1)
      do i=1,neps
          leps=lvalMin+((i-1)*dleps)
          lmepsfcn(i)=10._dp**(leps)
      end do
    else
      lmepsfcn(1)=epsilon(one)
    end if

    return
  end function initEpsfcn

  ! it determines the minimun of a vector,
  ! used to determine the min fitness
  function get_best_1(fitness) result(best)
    integer::best
    integer,dimension(1)::b1
    real(dp),dimension(:),intent(in)::fitness

    b1=minloc(fitness)
    best=b1(1)
    if(best.gt.size(fitness).or.best.le.0)best=1
    
    return
  end function get_best_1

  ! it determines the minimun of a vector,
  ! used to determine the min fitness
  function get_best_1_b(fitness,infos) result(best)
    use sorting,only:indexx
    integer::best
    real(dp),dimension(:),intent(in)::fitness
    integer,dimension(:),intent(in)::infos
    integer,dimension(:),allocatable::idx
    integer::ninfo,i
    
    best=1
    ninfo=size(infos)
    allocate(idx(ninfo))
    call indexx(fitness,idx)
    selinfo: do i=1,ninfo
      if((infos(idx(i)).gt.0).and.(infos(idx(i)).lt.4))then
          best=idx(i)
          exit selinfo
      end if
    end do selinfo
    deallocate(idx)
    if(best.gt.ninfo.or.best.le.0)best=1

    return
  end function get_best_1_b

  ! it determines the minimun of a section of an array,
  ! used to determine the min fitness
  function get_best_2(fitness) result(best)
    integer::best
    integer,dimension(1)::b1
    real(dp),dimension(:,:),intent(in)::fitness
    real(dp),dimension(:),allocatable::fitness1
    integer,dimension(2)::n1

    n1=shape(fitness)
    allocate(fitness1(n1(1)))
    b1=minloc(fitness1)
    best=b1(1)
    if(best.gt.n1(1).or.best.le.0)best=1
    deallocate(fitness1)

    return
  end function get_best_2

  
  ! SUBROUTINE lm_driver_4 WITH SEMI-AUTOMATIC SELECTION OF TOLERANCES
  ! varying epsfcn
  ! to use in GRID Search
  ! USED BY:
  ! LM (not BOOTSTRAP)
  subroutine lm_driver_4s(cpuid,isim,fcn,allpar,m,n,x,fvec,diag,sig,info,iwa)
    use parameters,only:path,maxfev,nprint,npar,ndata,dof,&
        &xtol,ftol,gtol,lmtols,progtype,paridlist,sig_paridlist
    use parameters_conversion,only:param_adj
    use init_trades,only:get_unit
    use convert_type,only:string
    !$ use omp_lib
    integer,intent(in)::cpuid,isim
    integer,intent(IN)::m
    integer,intent(IN)::n
    real(dp),dimension(:),intent(INOUT)::allpar,x
    real(dp),dimension(:),intent(OUT)::fvec,diag,sig
    integer,intent(OUT)::info
    integer,dimension(:),intent(OUT)::iwa

    integer::iflag

    real(dp),dimension(:),allocatable::inpar

    integer,parameter::neps=10 ! number of epsfcn values
    !integer,parameter::neps=2
    real(dp),dimension(:),allocatable::lmepsfcn
    integer::ieps

    real(dp)::fitness,fitness_x_dof
    integer::ulm,cpuid2
    character(512)::fllm
    character(80)::wrtfmt
    
    ! working variables
    real(dp),dimension(:),allocatable::wallpar,wpar,wfvec
    real(dp),dimension(:),allocatable::wa
    integer::mode,nfev
    real(dp),parameter::factor = 100._dp
    real(dp),dimension(:,:),allocatable::fjac

    real(dp),dimension(:),allocatable::wdiag,wsig
    integer::cntsig

    ! summary variables, store all the different outputs of the LM for each epsfcn value
    real(dp),dimension(:,:),allocatable::lmfitness,lmpar,lmfvec,&
        &lmdiag,lmsig
    integer,dimension(:,:),allocatable::lmiwa
    integer,dimension(:),allocatable::lmstat,infos
    
    integer::best

    interface
      subroutine fcn(allpar,m,n,x,fvec,iflag)
        use constants
        implicit none
        integer,intent(IN)::m,n
        real(dp),dimension(:),intent(IN)::allpar,x
        real(dp),dimension(:),intent(OUT)::fvec
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface

    write(*,'(a)')' ***************************'
    write(*,'(a)')' *** CALLED lm_driver_4s ***'
    write(*,'(a)')' ***************************'
    flush(6)


    info = 0
    iflag = 0 ! added by Luca, avoid negative value
    !     check the input parameters for errors.
    if (n <= 0 .or. m < n) return

    allocate(inpar(n))
    inpar=x
    
    ! initializes epsfcn values
    lmepsfcn=initEpsfcn(neps)
    
    ! makes a call to the function so it is possible to write the results as
    ! header of a summary file #_lmfit.log
    call fcn(allpar,m,n,x,fvec,iflag)
!     fitness_x_dof=sum(fvec*fvec)
!     fitness=fitness_x_dof/real(dof,dp)

    fitness=sum(fvec*fvec)
    fitness_x_dof=fitness*real(dof,dp)

    ulm = get_unit(cpuid)
    fllm = trim(path)//trim(adjustl(string(isim)))//"_lmfit.log"
    open(ulm,file=trim(fllm))
    write(ulm,'(a)')"# lmsim cntstat info "//trim(paridlist)//&
        &trim(sig_paridlist)//" Fitness*dof Fitness xtol ftol gtol epsfcn"
    wrtfmt = "(a,i4,i5,i3,1000("//trim(sprec)//",1x))"
    write(ulm,wrtfmt)"# ",0,0,0,x,(x-x),fitness_x_dof,fitness,one,one,one,one

    ! allocates and initializes variables to store all the results
    allocate(lmfitness(neps,2),lmpar(neps,n),lmfvec(neps,ndata))
    allocate(lmdiag(neps,n),lmsig(neps,n))
    allocate(lmiwa(neps,n))
    allocate(lmstat(neps),infos(neps))
    lmfitness=zero
    lmpar=zero
    lmfvec=zero
    lmdiag=zero
    lmsig=zero
    lmiwa=0
    lmstat=0
    infos=0

    ! allocate variable for LM computation
    allocate(wallpar(npar),wpar(n),wfvec(m))
    allocate(wa(2*n))
    allocate(fjac(m,n))
    allocate(wdiag(n),wsig(n))

    !$omp parallel do schedule(DYNAMIC,1)&
    !$omp& default(private) &
    !$omp& shared(cpuid,allpar,m,n,x,lmepsfcn,xtol,ftol,gtol,maxfev,&
    !$omp& nprint,lmfitness,dof,lmpar,lmfvec,lmdiag,lmsig,lmiwa,lmstat,&
    !$omp& infos)
    do ieps=1,neps
      cpuid2=cpuid
      !$ cpuid2 = omp_get_thread_num() + 1
      ! initialize every time to default values [these are private!]
      wallpar=allpar
      wpar=x
      wfvec=zero
      wa=zero
      mode=1
      info=0
      nfev=0
      fjac=zero
      iwa=0
      fitness_x_dof=zero
      fitness=zero

      ! call the LM
      call lmdif(fcn,wallpar,m,n,wpar,wfvec,ftol,xtol,gtol,maxfev,&
            &lmepsfcn(ieps),wa,mode,factor,nprint,info,nfev,&
            fjac,iwa,wa(n+1:))

      if(info == 8)then
          write(*,'(a)')"** info = 8 will be set to 4 **"
          info = 4
      end if
!       write(*,'(2(a,i3))')" subCPU ",cpuid2," of CPU ",cpuid,&
!             &" END LM. INFO = ",info
!       write(*,*)"" 

      ! calculates sigmas and adjust parameters
      wdiag=zero
      wsig=zero
      cntsig=0
      call get_covar(fjac,iwa,m,n,wdiag,wsig,cntsig)
!       call param_adj(wpar,wsig)

      ! save LM outputs
!       fitness_x_dof=sum(wfvec*wfvec)
!       fitness=fitness_x_dof/real(dof,dp)

      fitness=sum(wfvec*wfvec)
      fitness_x_dof=fitness*real(dof,dp)

      lmfitness(ieps,1)=fitness_x_dof
      lmfitness(ieps,2)=fitness
      lmstat(ieps)=cntsig
      infos(ieps)=info
      lmpar(ieps,:)=wpar
      lmfvec(ieps,:)=wfvec
      lmdiag(ieps,:)=wdiag
      lmsig(ieps,:)=wsig
      lmiwa(ieps,:)=iwa
      
      flush(6)
    end do
    !$omp end parallel do

    ! write final otuputs into a file
    wrtfmt=adjustl("(i4,i5,i3,10000("//trim(sprec)//",1x))")
    do ieps=1,neps
      write(ulm,trim(wrtfmt))&
            &ieps,lmstat(ieps),infos(ieps),lmpar(ieps,:),lmsig(ieps,:),&
            &lmfitness(ieps,:),xtol,ftol,gtol,lmepsfcn(ieps)
    end do
    close(ulm)

    ! select best simulation, it is the output of the subroutine
    !best=get_best(lmfitness(:,2))
    best=get_best(lmfitness(:,2),infos)
    x=lmpar(best,:)
    diag=lmdiag(best,:)
    sig=lmsig(best,:)
    fvec=lmfvec(best,:)
    iwa=lmiwa(best,:)
    info=infos(best)
    if(progtype.eq.1)then
      ! case GRID search
      lmtols(cpuid,4)=lmepsfcn(best)
    else
      ! case LM, PIKAIA, and PSO
      lmtols(:,4)=lmepsfcn(best)
    end if
!     call allpar_update(x,allpar)

    
    deallocate(inpar)
    if(allocated(lmepsfcn)) deallocate(lmepsfcn)
    deallocate(wallpar,wpar,wfvec)
    deallocate(wa)
    deallocate(fjac)
    deallocate(wdiag,wsig)
    deallocate(lmfitness,lmpar,lmfvec,lmdiag,lmsig)
    deallocate(lmiwa)
    deallocate(lmstat,infos)

    return
    !     last card of subroutine lmdif1.
  end subroutine lm_driver_4s

  ! SUBROUTINE lm_driver_4 WITH SEMI-AUTOMATIC SELECTION OF TOLERANCES
  ! varying epsfcn
  ! parallel version, to use with only LM, after PIKAIA or PSO
  subroutine lm_driver_4p(isim,fcn,allpar,m,n,x,fvec,diag,sig,info,iwa)
    use parameters,only:path,maxfev,nprint,npar,ndata,dof,&
        &xtol,ftol,gtol,lmtols,progtype,paridlist,sig_paridlist
    use parameters_conversion,only:param_adj
    use init_trades,only:get_unit
    use convert_type,only:string
    use ode_run,only:ode_lm
    !$ use omp_lib
    integer,intent(in)::isim
    integer,intent(IN)::m
    integer,intent(IN)::n
    real(dp),dimension(:),intent(INOUT)::allpar,x
    real(dp),dimension(:),intent(OUT)::fvec,diag,sig
    integer,intent(OUT)::info
    integer,dimension(:),intent(OUT)::iwa

    integer::iflag

    real(dp),dimension(:),allocatable::inpar

    integer,parameter::neps=10
    !integer,parameter::neps=2
    real(dp),dimension(:),allocatable::lmepsfcn
    integer::ieps

    real(dp)::fitness_x_dof,fitness
    integer::ulm,cpuid
    character(512)::fllm
    character(80)::wrtfmt

    real(dp),dimension(:),allocatable::wallpar,wpar,wfvec
    real(dp),dimension(:),allocatable::wa,qtf
    integer::mode,nfev
    real(dp),parameter::factor = 100._dp
    real(dp),dimension(:,:),allocatable::fjac

    real(dp),dimension(:),allocatable::wdiag,wsig
    integer::cntsig

    real(dp),dimension(:,:),allocatable::lmfitness,lmpar,lmfvec,&
        &lmdiag,lmsig
    integer,dimension(:,:),allocatable::lmiwa
    integer,dimension(:),allocatable::lmstat,infos
    
    integer::best

    interface
      subroutine fcn(allpar,m,n,x,fvec,iflag)
        use constants
        implicit none
        integer,intent(IN)::m,n
        real(dp),dimension(:),intent(IN)::allpar,x
        real(dp),dimension(:),intent(OUT)::fvec
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface

    write(*,'(a)')' ***************************'
    write(*,'(a)')' *** CALLED lm_driver_4p ***'
    write(*,'(a)')' ***************************'
    flush(6)

    info = 0
    iflag = 0 ! added by Luca, avoid negative value
    !     check the input parameters for errors.
    if (n <= 0 .or. m < n) return

    allocate(inpar(n))
    inpar=x
    
!     write(*,'(a)')" READY TO INITIALIZE EACH EPSFCN VALUE"
    lmepsfcn=initEpsfcn(neps)
    
    fitness_x_dof=zero
    fitness=zero
    fvec=zero
!     write(*,'(a)')" READY TO RUN FCN THE FIRST TIME"
    call fcn(allpar,m,n,x,fvec,iflag)
!     fitness_x_dof=sum(fvec*fvec)
!     fitness=fitness_x_dof/real(dof,dp)

    fitness=sum(fvec*fvec)
    fitness_x_dof=fitness*real(dof,dp)

    cpuid = 1
    ulm = get_unit(cpuid)
    fllm = trim(path)//trim(adjustl(string(isim)))//"_lmfit.log"
    open(ulm,file=trim(fllm))
    write(ulm,'(a)')"# lmsim cntstat info "//trim(paridlist)//&
        &trim(sig_paridlist)//" Fitness*dof Fitness xtol ftol gtol epsfcn"
    wrtfmt = adjustl("(a,i4,i5,i3,1000("//trim(sprec)//",1x))")
    write(ulm,trim(wrtfmt))&
        &"# ",0,0,0,x,(x-x),fitness_x_dof,fitness,one,one,one,one
    flush(6)
    
    ! allocate variables to store all the results
    allocate(lmfitness(neps,2),lmpar(neps,n),lmfvec(neps,ndata))
    allocate(lmdiag(neps,n),lmsig(neps,n))
    allocate(lmiwa(neps,n))
    allocate(lmstat(neps),infos(neps))
    lmfitness=zero
    lmpar=zero
    lmfvec=zero
    lmdiag=zero
    lmsig=zero
    lmiwa=0
    lmstat=0
    infos=0

    ! allocate variable for LM computation
    allocate(wallpar(npar),wpar(n),wfvec(m))
    allocate(wa(n),qtf(n))
    allocate(fjac(m,n))
    allocate(wdiag(n),wsig(n))

!     write(*,'(a)')" READY TO RUN LM FOR EACH EPSFCN VALUE"
!     write(*,'(a)')""

    !$omp parallel do schedule(DYNAMIC,1) &
    !$omp& default(private) &
    !$omp& shared(allpar,m,n,dof,x,lmepsfcn,xtol,ftol,gtol,maxfev,&
    !$omp& nprint,lmfitness,lmpar,lmfvec,lmdiag,lmsig,lmiwa,lmstat,&
    !$omp& infos) &
    !$omp& private(wallpar,wpar,wfvec,wa,qtf)
    
    do ieps=1,neps
      !$ cpuid = omp_get_thread_num() + 1

      ! initialize every time to default values
      wallpar=allpar
      wpar=x
      wfvec=zero
      wa=zero
      mode=1
      info=0
      nfev=0
      fjac=zero
      iwa=0
      qtf=zero
      fitness_x_dof=zero
      fitness=zero
      
      call lmdif(fcn,wallpar,m,n,wpar,wfvec,ftol,xtol,gtol,maxfev,&
            &lmepsfcn(ieps),wa,mode,factor,nprint,info,nfev,&
            fjac,iwa,qtf)

      if(info == 8)then
          write(*,'(a)')"** info = 8 will be set to 4 **"
          info = 4
      end if
!       write(*,'(2(a,i3))')" -- CPU ",cpuid," END LM. INFO = ",info
!       write(*,*)"" 

      wdiag=zero
      wsig=zero
      cntsig=0
      !write(*,'(a)')" GET_COVAR "
      call get_covar(fjac,iwa,m,n,wdiag,wsig,cntsig)
      !write(*,'(a)')" DONE "
      
!       call param_adj(wpar,wsig)
      !write(*,'(a)')" PAR ADJUSTED "

      ! save all the LM outputs
!       fitness_x_dof=sum(wfvec*wfvec)
!       fitness=fitness_x_dof/real(dof,dp)
      fitness=sum(wfvec*wfvec)
      fitness_x_dof=fitness*real(dof,dp)
      
      lmfitness(ieps,1)=fitness_x_dof
      lmfitness(ieps,2)=fitness
      lmstat(ieps)=cntsig
      infos(ieps)=info
      lmpar(ieps,:)=wpar
      lmfvec(ieps,:)=wfvec
      lmdiag(ieps,:)=wdiag
      lmsig(ieps,:)=wsig
      lmiwa(ieps,:)=iwa

    end do
    !$omp end parallel do

    ! write outputs of all the LM into a file
    wrtfmt=adjustl("(i4,i5,i3,10000("//trim(sprec)//",1x))")
    do ieps=1,neps
      write(ulm,trim(wrtfmt))ieps,lmstat(ieps),infos(ieps),&
            &lmpar(ieps,:),lmsig(ieps,:),&
            &lmfitness(ieps,:),xtol,ftol,gtol,lmepsfcn(ieps)
    end do
    close(ulm)

    ! select best simulation, output of lmdif
    !best=get_best(lmfitness(:,2))
    best=get_best(lmfitness(:,2),infos)
    if(best.eq.0)best=1
    x=lmpar(best,:)
    diag=lmdiag(best,:)
    sig=lmsig(best,:)
    fvec=lmfvec(best,:)
    iwa=lmiwa(best,:)
    info=infos(best)
    ! it should not be used in the grid search, but a check is better than a bug
    if(progtype.eq.1)then
      ! case GRID search
      !$ cpuid = omp_get_thread_num() + 1
      lmtols(cpuid,4)=lmepsfcn(best)
    else
      ! case LM, PIKAIA, and PSO
      lmtols(:,4)=lmepsfcn(best)
    end if

!     call allpar_update(x,allpar)

    deallocate(inpar)
    if(allocated(lmepsfcn)) deallocate(lmepsfcn)
    deallocate(wallpar,wpar,wfvec)
    deallocate(wa,qtf)
    deallocate(fjac)
    deallocate(wdiag,wsig)
    deallocate(lmfitness,lmpar,lmfvec,lmdiag,lmsig)
    deallocate(lmiwa)
    deallocate(lmstat,infos)

    return
    !     last card of subroutine lmdif1.
  end subroutine lm_driver_4p

  ! this is the true lmdif, but renamed and modified at my will :-)
  subroutine lmdif_2(fcn,allpar,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn, &
      diag,mode,factor,nprint,info,nfev,fjac,ipvt,qtf)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:59

    integer,intent(IN)::m,n
    real(dp),dimension(:),intent(INOUT)::allpar,x
    real(dp),dimension(:),intent(OUT)::fvec
    real(dp),intent(IN)::ftol,xtol
    real(dp),intent(INOUT)::gtol
    integer,intent(INOUT)::maxfev
    real(dp),intent(INOUT)::epsfcn
    real(dp),dimension(:),intent(OUT)::diag
    integer,intent(IN)::mode
    real(dp),intent(IN)::factor
    integer,intent(IN)::nprint
    integer,intent(OUT)::info
    integer,intent(OUT)::nfev
    real(dp),dimension(:,:),intent(OUT)::fjac   ! fjac(ldfjac,n)
    integer,dimension(:),intent(OUT)::ipvt
    real(dp),dimension(:),intent(OUT)::qtf

    ! EXTERNAL fcn

    interface
      subroutine fcn(allpar,m,n,x,fvec,iflag)
        use constants
        implicit none
        integer,intent(IN)::m,n
        real(dp),dimension(:),intent(IN)::allpar,x
        real(dp),dimension(:),intent(OUT)::fvec
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface

    !     **********

    !     subroutine lmdif

    !     the purpose of lmdif is to minimize the sum of the squares of
    !     m nonlinear functions in n variables by a modification of
    !     the levenberg-marquardt algorithm. The user must provide a
    !     subroutine which calculates the functions. The jacobian is
    !     then calculated by a forward-difference approximation.

    !     the subroutine statement is

    !       subroutine lmdif(fcn,allpar,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
    !                        diag,mode,factor,nprint,info,nfev,fjac,
    !                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)

    !     where

    !       fcn is the name of the user-supplied subroutine which
    !         calculates the functions. fcn must be declared
    !         in an external statement in the user calling
    !         program,and should be written as follows.

    !         subroutine fcn(m,n,x,fvec,iflag)
    !         integer m,n,iflag
    !         REAL (dp) x(:),fvec(m)
    !         ----------
    !         calculate the functions at x and
    !         return this vector in fvec.
    !         ----------
    !         return
    !         end

    !         the value of iflag should not be changed by fcn unless
    !         the user wants to terminate execution of lmdif.
    !         in this case set iflag to a negative integer.

    !       m is a positive integer input variable set to the number
    !         of functions.

    !       n is a positive integer input variable set to the number
    !         of variables. n must not exceed m.

    !       x is an array of length n. on input x must contain
    !         an initial estimate of the solution vector. on output x
    !         contains the final estimate of the solution vector.

    !       fvec is an output array of length m which contains
    !         the functions evaluated at the output x.

    !       ftol is a nonnegative input variable. termination
    !         occurs when both the actual and predicted relative
    !         reductions in the sum of squares are at most ftol.
    !         therefore,ftol measures the relative error desired
    !         in the sum of squares.

    !       xtol is a nonnegative input variable. termination
    !         occurs when the relative error between two consecutive
    !         iterates is at most xtol. therefore,xtol measures the
    !         relative error desired in the approximate solution.

    !       gtol is a nonnegative input variable. termination
    !         occurs when the cosine of the angle between fvec and
    !         any column of the jacobian is at most gtol in absolute
    !         value. therefore,gtol measures the orthogonality
    !         desired between the function vector and the columns
    !         of the jacobian.

    !       maxfev is a positive integer input variable. termination
    !         occurs when the number of calls to fcn is at least
    !         maxfev by the end of an iteration.

    !       epsfcn is an input variable used in determining a suitable
    !         step length for the forward-difference approximation. this
    !         approximation assumes that the relative errors in the
    !         functions are of the order of epsfcn. if epsfcn is less
    !         than the machine precision,it is assumed that the relative
    !         errors in the functions are of the order of the machine
    !         precision.

    !       diag is an array of length n. if mode = 1 (see
    !         below),diag is internally set. if mode = 2,diag
    !         must contain positive entries that serve as
    !         multiplicative scale factors for the variables.

    !       mode is an integer input variable. if mode = 1,the
    !         variables will be scaled internally. if mode = 2,
    !         the scaling is specified by the input diag. other
    !         values of mode are equivalent to mode = 1.

    !       factor is a positive input variable used in determining the
    !         initial step bound. this bound is set to the product of
    !         factor and the euclidean norm of diag*x if nonzero,or else
    !         to factor itself. in most cases factor should lie in the
    !         interval (.1,100.). 100. is a generally recommended value.

    !       nprint is an integer input variable that enables controlled
    !         printing of iterates if it is positive. in this case,
    !         fcn is called with iflag = 0 at the beginning of the first
    !         iteration and every nprint iterations thereafter and
    !         immediately prior to return,with x and fvec available
    !         for printing. if nprint is not positive,no special calls
    !         of fcn with iflag = 0 are made.

    !       info is an integer output variable.  If the user has terminated
    !         execution,info is set to the (negative) value of iflag.
    !         See description of fcn.  Otherwise,info is set as follows.

    !         info = 0  improper input parameters.

    !         info = 1  both actual and predicted relative reductions
    !                   in the sum of squares are at most ftol.

    !         info = 2  relative error between two consecutive iterates
    !                   is at most xtol.

    !         info = 3  conditions for info = 1 and info = 2 both hold.

    !         info = 4  the cosine of the angle between fvec and any column of
    !                   the Jacobian is at most gtol in absolute value.

    !         info = 5  number of calls to fcn has reached or exceeded maxfev.

    !         info = 6  ftol is too small. no further reduction in
    !                   the sum of squares is possible.

    !         info = 7  xtol is too small. no further improvement in
    !                   the approximate solution x is possible.

    !         info = 8  gtol is too small. fvec is orthogonal to the
    !                   columns of the jacobian to machine precision.

    !       nfev is an integer output variable set to the number of calls to fcn.

    !       fjac is an output m by n array. the upper n by n submatrix
    !         of fjac contains an upper triangular matrix r with
    !         diagonal elements of nonincreasing magnitude such that

    !                t     t           t
    !               p *(jac *jac)*p = r *r,

    !         where p is a permutation matrix and jac is the final calculated
    !         Jacobian.  Column j of p is column ipvt(j) (see below) of the
    !         identity matrix. the lower trapezoidal part of fjac contains
    !         information generated during the computation of r.

    !       ldfjac is a positive integer input variable not less than m
    !         which specifies the leading dimension of the array fjac.

    !       ipvt is an integer output array of length n. ipvt
    !         defines a permutation matrix p such that jac*p = q*r,
    !         where jac is the final calculated jacobian,q is
    !         orthogonal (not stored),and r is upper triangular
    !         with diagonal elements of nonincreasing magnitude.
    !         column j of p is column ipvt(j) of the identity matrix.

    !       qtf is an output array of length n which contains
    !         the first n elements of the vector (q transpose)*fvec.

    !       wa1,wa2,and wa3 are work arrays of length n.
    
    !       wa4 is a work array of length m.
    
    !     subprograms called
    
    !       user-supplied ...... fcn
    
    !       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
    
    !       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
    
    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more
    
    !     **********
    integer::i,iflag,iter,j,l
    real(dp)::actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm, &
        par,pnorm,prered,ratio,sum2,temp,temp1,temp2,xnorm
    !real(dp),dimension(n)::wa1,wa2,wa3
    !real(dp),dimension(m)::wa4
    real(dp),dimension(:),allocatable::wa1,wa2,wa3
    real(dp),dimension(:),allocatable::wa4
    real(dp),parameter::p1 = 0.1_dp,p5 = 0.5_dp, &
        p25 = 0.25_dp,p75 = 0.75_dp,p0001 = 0.0001_dp

    epsmch=epsilon(epsmch)

    info=0
    iflag=0
    nfev=0

    allocate(wa1(n),wa2(n),wa3(n))
    allocate(wa4(m))

!     write(*,'(a)')""
!     write(*,'(a)')" -- IN lmdif_2 --"
!     write(*,'(a,g25.15)')" ftol = ",ftol
!     write(*,'(a,g25.15)')" xtol = ",xtol
!     write(*,'(a,g25.15)')" gtol = ",gtol
!     write(*,'(a,g25.15)')" epsfcn = ",epsfcn
!     write(*,'(a,i10)')" maxfev = ",maxfev
!     write(*,'(a)')""

    if(n<=0)then
      go to 300
    else if(m<n)then
      go to 300
      !else if(ldfjac<m)then
      !go to 300
    else if(ftol<zero)then
      go to 300
    else if(xtol<zero)then
      go to 300
    else if(gtol<zero)then
      go to 300
    else if(maxfev<=0)then
      go to 300
    else if(factor<=zero)then
      go to 300
    endif

    if(mode==2)then
      do j=1,n
          if(diag(j)<=zero) go to 300
      end do
    endif
    !
    !  Evaluate the function at the starting point and calculate its norm.
    !
    iflag=1
    call fcn (allpar,m,n,x,fvec,iflag)
    nfev=1

    if (iflag<0) go to 300

    !fnorm=enorm(m,fvec)
    fnorm=enorm3(m,fvec)

    !
    !  Initialize Levenberg-Marquardt parameter and iteration counter.
    !
    par=zero
    iter=1
    !
    !  Beginning of the outer loop.
    !
30  continue
    !
    !  Calculate the jacobian matrix.
    !
    iflag=2
    !call fdjac2(fcn,allpar,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn)
    call fdjac2(fcn,allpar,m,n,x,fvec,fjac,iflag,epsfcn)
    nfev = nfev + n


    if ( iflag < 0 ) go to 300
    !
    !  If requested, call FCN to enable printing of iterates.
    !
    if ( 0 < nprint ) then
      iflag = 0
      if (mod(iter-1,nprint)==0) call fcn(allpar,m,n,x,fvec,iflag)
      if ( iflag < 0 ) go to 300
    end if
    !
    !  Compute the QR factorization of the jacobian.
    !
    !call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2)
    call qrfac(m,n,fjac,.true.,ipvt,wa1,wa2) ! Luca

    !
    !  On the first iteration and if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
    if ( iter == 1 ) then

      if ( mode /= 2 ) then
          diag(1:n) = wa2(1:n)
          do j = 1, n
            if ( wa2(j) == zero ) diag(j) = one
          end do
      end if
      !
      !  On the first iteration, calculate the norm of the scaled X
      !  and initialize the step bound DELTA.
      !
      wa3(1:n) = diag(1:n) * x(1:n)
      !xnorm = enorm ( n, wa3 )
      xnorm = enorm3( n, wa3 )
      delta = factor * xnorm
      if ( delta == zero ) delta = factor

    end if

    !
    !  Form Q' * FVEC and store the first N components in QTF.
    !
    wa4(1:m) = fvec(1:m)

    do j = 1, n

      if ( fjac(j,j) /= zero ) then
          sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
          temp = - sum2 / fjac(j,j)
          wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
      end if

      fjac(j,j) = wa1(j)
      qtf(j) = wa4(j)

    end do


    !
    !  Compute the norm of the scaled gradient.
    !
    gnorm = zero

    if ( fnorm /= zero ) then

      do j = 1, n

          l = ipvt(j)

          if ( wa2(l) /= zero ) then
            sum2 = zero
            do i = 1, j
                sum2 = sum2 + fjac(i,j) * ( qtf(i) / fnorm )
            end do
            gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
          end if

      end do

    end if

    !
    !  Test for convergence of the gradient norm.
    !
    if ( gnorm <= gtol ) then
      info = 4
      go to 300
    end if
    !
    !  Rescale if necessary.
    !
    if ( mode /= 2 ) then
      do j = 1, n
          diag(j) = max ( diag(j), wa2(j) )
      end do
    end if

    !
    !  Beginning of the inner loop.
    !
200 continue
    !
    !  Determine the Levenberg-Marquardt parameter.
    !
    !call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
    call lmpar(n,fjac,ipvt,diag,qtf,delta,par,wa1,wa2)

    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
    wa1(1:n) = -wa1(1:n)
    wa2(1:n) = x(1:n) + wa1(1:n)
    wa3(1:n) = diag(1:n) * wa1(1:n)

    !pnorm = enorm ( n, wa3 )
    pnorm = enorm3( n, wa3 )

    !
    !  On the first iteration, adjust the initial step bound.
    !
    if ( iter == 1 ) delta = min ( delta, pnorm )
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
    iflag = 1
    call fcn(allpar,m,n,wa2,wa4,iflag)
    nfev = nfev + 1

    if ( iflag < 0 ) go to 300
    !fnorm1 = enorm ( m, wa4 )
    fnorm1 = enorm3( m, wa4 )
    
    !
    !  Compute the scaled actual reduction.
    !
    if ( p1 * fnorm1 < fnorm ) then
      actred = one - ( fnorm1 / fnorm )**2
    else
      actred = -one
    end if

    !
    !  Compute the scaled predicted reduction and the scaled directional derivative.
    !
    do j = 1, n
      wa3(j) = zero
      l = ipvt(j)
      temp = wa1(l)
      wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
    end do

    !temp1 = enorm ( n, wa3 ) / fnorm
    temp1 = enorm3( n, wa3 ) / fnorm
    temp2 = ( sqrt ( par ) * pnorm ) / fnorm
    prered = temp1**2 + temp2**2 / p5
    dirder = - ( temp1**2 + temp2**2 )
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
    ratio = zero
    if ( prered /= zero ) ratio = actred / prered

    !
    !  Update the step bound.
    !
    if ( ratio <= p25 ) then

      if ( actred >= zero ) temp = p5
      if ( actred < zero ) temp = p5 * dirder / ( dirder + p5 * actred )
      if ( p1 * fnorm1 >= fnorm .or. temp < p1 ) temp = p1
      delta = temp * min ( delta, pnorm / p1  )
      par = par / temp

    else

      if ( par == zero .or. ratio >= p75 ) then
          delta = 2._dp * pnorm
          par = p5 * par
      end if

    end if

    !
    !  Test for successful iteration.
    !

    !
    !  Successful iteration. update X, FVEC, and their norms.
    !
    !if ( 0.0001D+00 <= ratio ) then
    if ( p0001 <= ratio ) then
      x(1:n) = wa2(1:n)
      wa2(1:n) = diag(1:n) * x(1:n)
      fvec(1:m) = wa4(1:m)
      !xnorm = enorm ( n, wa2 )
      xnorm = enorm3( n, wa2 )
      fnorm = fnorm1
      iter = iter + 1
    end if
    !
    !  Tests for convergence.
    !
    if ( abs ( actred) <= ftol .and. prered <= ftol &
        .and. p5 * ratio <= one ) info = 1

    if ( delta <= xtol * xnorm ) info = 2

    if ( abs ( actred) <= ftol .and. prered <= ftol &
        .and. p5 * ratio <= one .and. info == 2 ) info = 3

    if ( info /= 0 ) go to 300
    !
    !  Tests for termination and stringent tolerances.
    !
    if ( maxfev <= nfev ) info = 5

    if ( abs ( actred) <= epsmch .and. prered <= epsmch &
        .and. p5 * ratio <= one ) info = 6

    if ( delta <= epsmch * xnorm ) info = 7

    if ( gnorm <= epsmch ) info = 8

    if ( info /= 0 ) go to 300
    !
    !  End of the inner loop.  Repeat if iteration unsuccessful.
    !
    if ( ratio < p0001 ) go to 200
    !
    !  End of the outer loop.
    !
    go to 30

300 continue
    
    ! add by Luca to check value of fvec and fitness_x_dof, fitness
    ! ! call fcn(allpar,m,n,x,fvec,iflag)
    !call callodeout(allpar,x,fvec)
    !write(*,'(a)')""
    !write(*,'(a)')" ================IN lmdif_2======================="
    !write(*,'(a,i6,a)')" called fcn ",nfev," times"
    !write(*,'(a)')" parameters:"
    !write(*,'(1000g25.15)')x(1:n)
    !write(*,'(2(a,g25.15))')" fitness_x_dof = ",sum(fvec*fvec),&
    !     " fitness = ",sum(fvec*fvec)/real((m-n),dp)
    !write(*,'(a)')" ================IN lmdif_2======================="
    !write(*,'(a)')""

    !
    !  Termination, either normal or user imposed.
    !
    if ( iflag < 0 ) info = iflag

    deallocate(wa1,wa2,wa3)
    deallocate(wa4)
    iflag = 0

    if ( 0 < nprint ) call fcn(allpar,m,n,x,fvec,iflag)

    return
  end subroutine lmdif_2

  ! as lmdif_2 but with more arguments needed by the fcn to use
  subroutine lmdif_3(fcn,allpar,m,n,x,RV_obs,T0_obs,fvec,&
      &ftol,xtol,gtol,maxfev,epsfcn,&
      &diag,mode,factor,nprint,info,nfev,fjac,ipvt,qtf)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:59

    ! N.B. Arguments LDFJAC,WA1,WA2,WA3 & WA4 have been removed.

    integer,intent(IN)::m,n
    real(dp),dimension(:),intent(INOUT)::allpar,x
    real(dp),dimension(:),intent(in)::RV_obs
    real(dp),dimension(:,:),intent(in)::T0_obs
    real(dp),dimension(:),intent(OUT)::fvec
    real(dp),intent(IN)::ftol,xtol
    real(dp),intent(INOUT)::gtol
    integer,intent(INOUT)::maxfev
    real(dp),intent(INOUT)::epsfcn
    real(dp),dimension(:),intent(OUT)::diag
    integer,intent(IN)::mode
    real(dp),intent(IN)::factor
    integer,intent(IN)::nprint
    integer,intent(OUT)::info
    integer,intent(OUT)::nfev
    real(dp),dimension(:,:),intent(OUT)::fjac   ! fjac(ldfjac,n)
    integer,dimension(:),intent(OUT)::ipvt
    real(dp),dimension(:),intent(OUT)::qtf

    ! EXTERNAL fcn

    interface
      subroutine fcn(allpar,m,n,x,RV_obs,T0_obs,fvec,iflag)
        use constants
        implicit none
        integer,intent(IN)::m,n
        real(dp),dimension(:),intent(IN)::allpar,x,RV_obs
        real(dp),dimension(:,:),intent(in)::T0_obs
        real(dp),dimension(:),intent(OUT)::fvec
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface

    !     **********
    integer::i,iflag,iter,j,l
    real(dp)::actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm, &
        par,pnorm,prered,ratio,sum2,temp,temp1,temp2,xnorm
    !real(dp),dimension(n)::wa1,wa2,wa3
    !real(dp),dimension(m)::wa4
    real(dp),dimension(:),allocatable::wa1,wa2,wa3
    real(dp),dimension(:),allocatable::wa4
    real(dp),parameter::p1 = 0.1_dp,p5 = 0.5_dp, &
        p25 = 0.25_dp,p75 = 0.75_dp,p0001 = 0.0001_dp

    !     epsmch is the machine precision.

    epsmch = epsilon(zero)

    info = 0
    iflag = 0
    nfev = 0

    allocate(wa1(n),wa2(n),wa3(n))
    allocate(wa4(m))

!     write(*,'(a)')""
!     write(*,'(a)')" -- IN lmdif_3 --"
!     write(*,'(a,g25.15)')" ftol = ",ftol
!     write(*,'(a,g25.15)')" xtol = ",xtol
!     write(*,'(a,g25.15)')" gtol = ",gtol
!     write(*,'(a,g25.15)')" epsfcn = ",epsfcn
!     write(*,'(a)')""

    !     check the input parameters for errors.
    !write(*,*)" 10 "

    if(n<=0)then
      go to 300
    else if(m<n)then
      go to 300
      !else if(ldfjac<m)then
      !go to 300
    else if(ftol<zero)then
      go to 300
    else if(xtol<zero)then
      go to 300
    else if(gtol<zero)then
      go to 300
    else if(maxfev<=0)then
      go to 300
    else if(factor<=zero)then
      go to 300
    endif

    if(mode==2)then
      do j=1,n
          if(diag(j)<=zero) go to 300
      end do
    endif
    !
    !  Evaluate the function at the starting point and calculate its norm.
    !
    iflag=1
    call fcn(allpar,m,n,x,RV_obs,T0_obs,fvec,iflag)
    nfev = 1

    if (iflag < 0) go to 300

    !fnorm=enorm(m,fvec)
    fnorm=enorm3(m,fvec)
    !
    !  Initialize Levenberg-Marquardt parameter and iteration counter.
    !
    par=zero
    iter=1
    !
    !  Beginning of the outer loop.
    !
30  continue
    !
    !  Calculate the jacobian matrix.
    !
    iflag=2
    call fdjac2(fcn,allpar,m,n,x,RV_obs,T0_obs,fvec,fjac,iflag,epsfcn)
    nfev = nfev + n
    if (iflag < 0) go to 300
    !
    !        if requested, call fcn to enable printing of iterates.
    !
    if ( 0 < nprint ) then
      iflag = 0
      if (mod(iter-1,nprint) == 0) &
            &call fcn(allpar,m,n,x,RV_obs,T0_obs,fvec,iflag)
      if ( iflag < 0 ) go to 300
    end if
    !
    !  Compute the QR factorization of the jacobian.
    !
    !call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2)
    call qrfac(m,n,fjac,.true.,ipvt,wa1,wa2) ! Luca
    !
    !  On the first iteration and if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
    if ( iter == 1 ) then

      if ( mode /= 2 ) then
          diag(1:n) = wa2(1:n)
          do j = 1, n
            if ( wa2(j) == zero ) diag(j) = one
          end do
      end if
      !
      !  On the first iteration, calculate the norm of the scaled X
      !  and initialize the step bound DELTA.
      !
      wa3(1:n) = diag(1:n) * x(1:n)
      !xnorm = enorm ( n, wa3 )
      xnorm = enorm3( n, wa3 )
      delta = factor * xnorm
      if ( delta == zero ) delta = factor
    end if
    !
    !  Form Q' * FVEC and store the first N components in QTF.
    !
    wa4(1:m) = fvec(1:m)

    do j = 1, n

      if ( fjac(j,j) /= zero ) then
          sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
          temp = - sum2 / fjac(j,j)
          wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
      end if

      fjac(j,j) = wa1(j)
      qtf(j) = wa4(j)

    end do

    !
    !  Compute the norm of the scaled gradient.
    !
    gnorm = zero

    if ( fnorm /= zero ) then

      do j = 1, n

          l = ipvt(j)

          if ( wa2(l) /= zero ) then
            sum2 = zero
            do i = 1, j
                sum2 = sum2 + fjac(i,j) * ( qtf(i) / fnorm )
            end do
            gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
          end if

      end do

    end if
    !
    !  Test for convergence of the gradient norm.
    !
    if ( gnorm <= gtol ) then
      info = 4
      go to 300
    end if
    !
    !  Rescale if necessary.
    !
    if ( mode /= 2 ) then
      do j = 1, n
          diag(j) = max ( diag(j), wa2(j) )
      end do
    end if
    !
    !  Beginning of the inner loop.
    !
200 continue
    !
    !  Determine the Levenberg-Marquardt parameter.
    !
    call lmpar(n,fjac,ipvt,diag,qtf,delta,par,wa1,wa2)
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
    wa1(1:n) = -wa1(1:n)
    wa2(1:n) = x(1:n) + wa1(1:n)
    wa3(1:n) = diag(1:n) * wa1(1:n)

    !pnorm = enorm ( n, wa3 )
    pnorm = enorm3( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
    if ( iter == 1 ) delta = min ( delta, pnorm )
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
    iflag = 1
    call fcn(allpar,m,n,wa2,RV_obs,T0_obs,wa4,iflag)
    nfev = nfev + 1
    if ( iflag < 0 ) go to 300
    !fnorm1 = enorm ( m, wa4 )
    fnorm1 = enorm3( m, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
    if ( p1 * fnorm1 < fnorm ) then
      actred = one - ( fnorm1 / fnorm )**2
    else
      actred = -one
    end if
    !
    !  Compute the scaled predicted reduction and the scaled directional derivative.
    !
    do j = 1, n
      wa3(j) = zero
      l = ipvt(j)
      temp = wa1(l)
      wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
    end do

    !temp1 = enorm ( n, wa3 ) / fnorm
    temp1 = enorm3( n, wa3 ) / fnorm
    temp2 = ( sqrt ( par ) * pnorm ) / fnorm
    prered = temp1**2 + temp2**2 / p5
    dirder = - ( temp1**2 + temp2**2 )
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
    ratio = zero
    if ( prered /= zero ) ratio = actred / prered
    !
    !  Update the step bound.
    !
    if ( ratio <= p25 ) then

      if ( actred >= zero ) temp = p5
      if ( actred < zero ) temp = p5 * dirder / ( dirder + p5 * actred )
      if ( p1 * fnorm1 >= fnorm .or. temp < p1 ) temp = p1
      delta = temp * min ( delta, pnorm / p1  )
      par = par / temp

    else

      if ( par == zero .or. ratio >= p75 ) then
          delta = 2._dp * pnorm
          par = p5 * par
      end if

    end if
    !
    !  Test for successful iteration.
    !

    !
    !  Successful iteration. update X, FVEC, and their norms.
    !
    !if ( 0.0001D+00 <= ratio ) then
    if ( p0001 <= ratio ) then
      x(1:n) = wa2(1:n)
      wa2(1:n) = diag(1:n) * x(1:n)
      fvec(1:m) = wa4(1:m)
      !xnorm = enorm ( n, wa2 )
      xnorm = enorm3( n, wa2 )
      fnorm = fnorm1
      iter = iter + 1
    end if
    !
    !  Tests for convergence.
    !
    if ( abs ( actred) <= ftol .and. prered <= ftol &
        .and. p5 * ratio <= one ) info = 1

    if ( delta <= xtol * xnorm ) info = 2

    if ( abs ( actred) <= ftol .and. prered <= ftol &
        .and. p5 * ratio <= one .and. info == 2 ) info = 3

    if ( info /= 0 ) go to 300
    !
    !  Tests for termination and stringent tolerances.
    !
    if ( maxfev <= nfev ) info = 5

    if ( abs ( actred) <= epsmch .and. prered <= epsmch &
        .and. p5 * ratio <= one ) info = 6

    if ( delta <= epsmch * xnorm ) info = 7

    if ( gnorm <= epsmch ) info = 8

    if ( info /= 0 ) go to 300
    !
    !  End of the inner loop.  Repeat if iteration unsuccessful.
    !
    if ( ratio < p0001 ) go to 200
    !
    !  End of the outer loop.
    !
    go to 30

300 continue
    !
    !  Termination, either normal or user imposed.
    !
    if ( iflag < 0 ) info = iflag

    iflag = 0
    ! Luca
    deallocate(wa1,wa2,wa3)
    deallocate(wa4)

    if (nprint > 0) call fcn(allpar,m,n,x,RV_obs,T0_obs,fvec,iflag)

    return
  end subroutine lmdif_3


  subroutine lmder1(fcn,m,n,x,fvec,fjac,tol,info,ipvt)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:54

    ! N.B. Arguments LDFJAC,WA & LWA have been removed.

    integer,intent(IN)::m
    integer,intent(IN)::n
    real(dp),intent(INOUT)::x(:)
    real(dp),intent(INOUT)::fvec(:)
    real(dp),intent(INOUT)::fjac(:,:)    ! fjac(ldfjac,n)
    real(dp),intent(IN)::tol
    integer,intent(OUT)::info
    integer,intent(INOUT)::ipvt(:)


    ! EXTERNAL fcn

    interface
      subroutine fcn(m,n,x,fvec,fjac,iflag)
        use constants
        implicit none
        !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)
        integer,intent(IN)::m,n
        real(dp),intent(IN)::x(:)
        real(dp),intent(OUT)::fvec(:)
        real(dp),intent(OUT)::fjac(:,:)
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface

    !     **********

    !     subroutine lmder1

    !     The purpose of lmder1 is to minimize the sum of the squares of
    !     m nonlinear functions in n variables by a modification of the
    !     levenberg-marquardt algorithm.  This is done by using the more
    !     general least-squares solver lmder.  The user must provide a
    !     subroutine which calculates the functions and the jacobian.

    !     the subroutine statement is

    !       subroutine lmder1(fcn,m,n,x,fvec,fjac,tol,info,ipvt)

    !     where

    !       fcn is the name of the user-supplied subroutine which
    !         calculates the functions and the jacobian.  fcn must
    !         be declared in an interface statement in the user
    !         calling program,and should be written as follows.

    !         subroutine fcn(m,n,x,fvec,fjac,iflag)
    !         integer::m,n,ldfjac,iflag
    !         REAL (dp)::x(:),fvec(:),fjac(:,:)
    !         ----------
    !         if iflag = 1 calculate the functions at x and
    !         return this vector in fvec. do not alter fjac.
    !         if iflag = 2 calculate the jacobian at x and
    !         return this matrix in fjac. do not alter fvec.
    !         ----------
    !         return
    !         end

    !         the value of iflag should not be changed by fcn unless
    !         the user wants to terminate execution of lmder1.
    !         in this case set iflag to a negative integer.

    !       m is a positive integer input variable set to the number of functions.

    !       n is a positive integer input variable set to the number
    !         of variables.  n must not exceed m.

    !       x is an array of length n. on input x must contain
    !         an initial estimate of the solution vector. on output x
    !         contains the final estimate of the solution vector.

    !       fvec is an output array of length m which contains
    !         the functions evaluated at the output x.

    !       fjac is an output m by n array. the upper n by n submatrix
    !         of fjac contains an upper triangular matrix r with
    !         diagonal elements of nonincreasing magnitude such that

    !                t     t           t
    !               p *(jac *jac)*p = r *r,

    !         where p is a permutation matrix and jac is the final calculated
    !         Jacobian.  Column j of p is column ipvt(j) (see below) of the
    !         identity matrix.  The lower trapezoidal part of fjac contains
    !         information generated during the computation of r.

    !       ldfjac is a positive integer input variable not less than m
    !         which specifies the leading dimension of the array fjac.

    !       tol is a nonnegative input variable. termination occurs
    !         when the algorithm estimates either that the relative
    !         error in the sum of squares is at most tol or that
    !         the relative error between x and the solution is at most tol.

    !       info is an integer output variable.  If the user has terminated
    !         execution,info is set to the (negative) value of iflag.
    !         See description of fcn.  Otherwise,info is set as follows.

    !         info = 0  improper input parameters.

    !         info = 1  algorithm estimates that the relative error
    !                   in the sum of squares is at most tol.

    !         info = 2  algorithm estimates that the relative error
    !                   between x and the solution is at most tol.

    !         info = 3  conditions for info = 1 and info = 2 both hold.

    !         info = 4  fvec is orthogonal to the columns of the
    !                   jacobian to machine precision.

    !         info = 5  number of calls to fcn with iflag = 1 has reached 100*(n+1).

    !         info = 6  tol is too small.  No further reduction in
    !                   the sum of squares is possible.

    !         info = 7  tol is too small.  No further improvement in
    !                   the approximate solution x is possible.

    !       ipvt is an integer output array of length n. ipvt
    !         defines a permutation matrix p such that jac*p = q*r,
    !         where jac is the final calculated jacobian,q is
    !         orthogonal (not stored),and r is upper triangular
    !         with diagonal elements of nonincreasing magnitude.
    !         column j of p is column ipvt(j) of the identity matrix.

    !       wa is a work array of length lwa.

    !       lwa is a positive integer input variable not less than 5*n+m.

    !     subprograms called

    !       user-supplied ...... fcn

    !       minpack-supplied ... lmder

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    integer::maxfev,mode,nfev,njev,nprint
    real(dp)::ftol,gtol,xtol,wa(2*n)
    real(dp),parameter::factor = 100._dp,zero = 0.0_dp

    info = 0

    !     check the input parameters for errors.

    if ( n <= 0 .or. m < n .or. tol < zero ) GO TO 10

    !     call lmder.

    maxfev = 100*(n + 1)
    ftol = tol
    xtol = tol
    gtol = zero
    mode = 1
    nprint = 0
    call lmder(fcn,m,n,x,fvec,fjac,ftol,xtol,gtol,maxfev, &
        wa,mode,factor,nprint,info,nfev,njev,ipvt,wa(n+1:) )
    if (info == 8) info = 4

10  return

    !     last card of subroutine lmder1.

  end subroutine lmder1



  subroutine lmder(fcn,m,n,x,fvec,fjac,ftol,xtol,gtol,maxfev,&
      diag,mode,factor,nprint,info,nfev,njev,ipvt,qtf)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:50

    ! N.B. Arguments LDFJAC,WA1,WA2,WA3 & WA4 have been removed.

    integer,intent(IN)::m
    integer,intent(IN)::n
    real(dp),intent(INOUT)::x(:)
    real(dp),intent(OUT)::fvec(m)
    real(dp),intent(OUT)::fjac(:,:)    ! fjac(ldfjac,n)
    real(dp),intent(IN)::ftol
    real(dp),intent(IN)::xtol
    real(dp),intent(INOUT)::gtol
    integer,intent(INOUT)::maxfev
    real(dp),intent(OUT)::diag(:)
    integer,intent(IN)::mode
    real(dp),intent(IN)::factor
    integer,intent(IN)::nprint
    integer,intent(OUT)::info
    integer,intent(OUT)::nfev
    integer,intent(OUT)::njev
    integer,intent(OUT)::ipvt(:)
    real(dp),intent(OUT)::qtf(:)

    interface
      subroutine fcn(m,n,x,fvec,fjac,iflag)
        use constants
        implicit none
        !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)
        integer,intent(IN)::m,n
        real(dp),intent(IN)::x(:)
        real(dp),intent(OUT)::fvec(:)
        real(dp),intent(OUT)::fjac(:,:)
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface


    !     **********

    !     subroutine lmder

    !     the purpose of lmder is to minimize the sum of the squares of
    !     m nonlinear functions in n variables by a modification of
    !     the levenberg-marquardt algorithm. the user must provide a
    !     subroutine which calculates the functions and the jacobian.

    !     the subroutine statement is

    !       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
    !                        maxfev,diag,mode,factor,nprint,info,nfev,
    !                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)

    !     where

    !       fcn is the name of the user-supplied subroutine which
    !         calculates the functions and the jacobian. fcn must
    !         be declared in an external statement in the user
    !         calling program,and should be written as follows.

    !         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
    !         integer m,n,ldfjac,iflag
    !         REAL (dp) x(:),fvec(m),fjac(ldfjac,n)
    !         ----------
    !         if iflag = 1 calculate the functions at x and
    !         return this vector in fvec. do not alter fjac.
    !         if iflag = 2 calculate the jacobian at x and
    !         return this matrix in fjac. do not alter fvec.
    !         ----------
    !         return
    !         end

    !         the value of iflag should not be changed by fcn unless
    !         the user wants to terminate execution of lmder.
    !         in this case set iflag to a negative integer.

    !       m is a positive integer input variable set to the number
    !         of functions.

    !       n is a positive integer input variable set to the number
    !         of variables. n must not exceed m.

    !       x is an array of length n. on input x must contain
    !         an initial estimate of the solution vector. on output x
    !         contains the final estimate of the solution vector.

    !       fvec is an output array of length m which contains
    !         the functions evaluated at the output x.

    !       fjac is an output m by n array. the upper n by n submatrix
    !         of fjac contains an upper triangular matrix r with
    !         diagonal elements of nonincreasing magnitude such that

    !                t     t           t
    !               p *(jac *jac)*p = r *r

    !         where p is a permutation matrix and jac is the final calculated
    !         jacobian.  Column j of p is column ipvt(j) (see below) of the
    !         identity matrix.  The lower trapezoidal part of fjac contains
    !         information generated during the computation of r.

    !       ldfjac is a positive integer input variable not less than m
    !         which specifies the leading dimension of the array fjac.

    !       ftol is a nonnegative input variable.  Termination occurs when both
    !         the actual and predicted relative reductions in the sum of squares
    !         are at most ftol.   Therefore,ftol measures the relative error
    !         desired in the sum of squares.

    !       xtol is a nonnegative input variable. termination
    !         occurs when the relative error between two consecutive
    !         iterates is at most xtol. therefore,xtol measures the
    !         relative error desired in the approximate solution.

    !       gtol is a nonnegative input variable.  Termination occurs when the
    !         cosine of the angle between fvec and any column of the jacobian is
    !         at most gtol in absolute value.  Therefore,gtol measures the
    !         orthogonality desired between the function vector and the columns
    !         of the jacobian.

    !       maxfev is a positive integer input variable.  Termination occurs when
    !         the number of calls to fcn with iflag = 1 has reached maxfev.

    !       diag is an array of length n.  If mode = 1 (see below),diag is
    !         internally set.  If mode = 2,diag must contain positive entries
    !         that serve as multiplicative scale factors for the variables.

    !       mode is an integer input variable.  if mode = 1,the
    !         variables will be scaled internally.  if mode = 2,
    !         the scaling is specified by the input diag.  other
    !         values of mode are equivalent to mode = 1.

    !       factor is a positive input variable used in determining the
    !         initial step bound. this bound is set to the product of
    !         factor and the euclidean norm of diag*x if nonzero,or else
    !         to factor itself. in most cases factor should lie in the
    !         interval (.1,100.).100. is a generally recommended value.

    !       nprint is an integer input variable that enables controlled printing
    !         of iterates if it is positive.  In this case,fcn is called with
    !         iflag = 0 at the beginning of the first iteration and every nprint
    !         iterations thereafter and immediately prior to return,with x,fvec,
    !         and fjac available for printing.  fvec and fjac should not be
    !         altered.  If nprint is not positive,no special calls of fcn with
    !         iflag = 0 are made.

    !       info is an integer output variable.  If the user has terminated
    !         execution,info is set to the (negative) value of iflag.
    !         See description of fcn.  Otherwise,info is set as follows.

    !         info = 0  improper input parameters.

    !         info = 1  both actual and predicted relative reductions
    !                   in the sum of squares are at most ftol.

    !         info = 2  relative error between two consecutive iterates
    !                   is at most xtol.

    !         info = 3  conditions for info = 1 and info = 2 both hold.

    !         info = 4  the cosine of the angle between fvec and any column of
    !                   the jacobian is at most gtol in absolute value.

    !         info = 5  number of calls to fcn with iflag = 1 has reached maxfev.

    !         info = 6  ftol is too small.  No further reduction in
    !                   the sum of squares is possible.

    !         info = 7  xtol is too small.  No further improvement in
    !                   the approximate solution x is possible.

    !         info = 8  gtol is too small.  fvec is orthogonal to the
    !                   columns of the jacobian to machine precision.

    !       nfev is an integer output variable set to the number of
    !         calls to fcn with iflag = 1.

    !       njev is an integer output variable set to the number of
    !         calls to fcn with iflag = 2.

    !       ipvt is an integer output array of length n.  ipvt
    !         defines a permutation matrix p such that jac*p = q*r,
    !         where jac is the final calculated jacobian,q is
    !         orthogonal (not stored),and r is upper triangular
    !         with diagonal elements of nonincreasing magnitude.
    !         column j of p is column ipvt(j) of the identity matrix.

    !       qtf is an output array of length n which contains
    !         the first n elements of the vector (q transpose)*fvec.

    !       wa1,wa2,and wa3 are work arrays of length n.

    !       wa4 is a work array of length m.

    !     subprograms called

    !       user-supplied ...... fcn

    !       minpack-supplied ... dpmpar,enorm,lmpar,qrfac

    !       fortran-supplied ... ABS,MAX,MIN,SQRT,mod

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    integer::i,iflag,iter,j,l
    real(dp)::actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm, &
        par,pnorm,prered,ratio,sum,temp,temp1,temp2,xnorm
    real(dp)::wa1(n),wa2(n),wa3(n),wa4(m)
    real(dp),parameter::one = 1.0_dp,p1 = 0.1_dp,p5 = 0.5_dp, &
        p25 = 0.25_dp,p75 = 0.75_dp,p0001 = 0.0001_dp,&
        zero = 0.0_dp

    !     epsmch is the machine precision.

    epsmch = epsilon(zero)

    info = 0
    iflag = 0
    nfev = 0
    njev = 0

    !     check the input parameters for errors.

    if (n <= 0 .or. m < n .or. ftol < zero .or. xtol < zero .or. gtol < zero  &
        .or. maxfev <= 0 .or. factor <= zero) GO TO 300
    if (mode /= 2) GO TO 20
    do  j = 1,n
      if (diag(j) <= zero) GO TO 300
    end do

    !     evaluate the function at the starting point and calculate its norm.

20  iflag = 1
    call fcn(m,n,x,fvec,fjac,iflag)
    nfev = 1
    if (iflag < 0) GO TO 300
    !fnorm = enorm(m,fvec)
    fnorm = enorm3(m,fvec)

    !     initialize levenberg-marquardt parameter and iteration counter.

    par = zero
    iter = 1

    !     beginning of the outer loop.

    !        calculate the jacobian matrix.

30  iflag = 2
    call fcn(m,n,x,fvec,fjac,iflag)
    njev = njev + 1
    if (iflag < 0) GO TO 300

    !        if requested,call fcn to enable printing of iterates.

    if (nprint <= 0) GO TO 40
    iflag = 0
    if (mod(iter-1,nprint) == 0) call fcn(m,n,x,fvec,fjac,iflag)
    if (iflag < 0) GO TO 300

    !        compute the qr factorization of the jacobian.

40  call qrfac(m,n,fjac,.true.,ipvt,wa1,wa2)

    !        on the first iteration and if mode is 1,scale according
    !        to the norms of the columns of the initial jacobian.

    if (iter /= 1) GO TO 80
    if (mode == 2) GO TO 60
    do  j = 1,n
      diag(j) = wa2(j)
      if (wa2(j) == zero) diag(j) = one
    end do

    !        on the first iteration,calculate the norm of the scaled x
    !        and initialize the step bound delta.

60  wa3(1:n) = diag(1:n)*x(1:n)
    !xnorm = enorm(n,wa3)
    xnorm = enorm3(n,wa3)
    delta = factor*xnorm
    if (delta == zero) delta = factor

    !        form (q transpose)*fvec and store the first n components in qtf.

80  wa4(1:m) = fvec(1:m)
    do  j = 1,n
      if (fjac(j,j) == zero) GO TO 120
      sum = dot_product( fjac(j:m,j),wa4(j:m) )
      temp = -sum/fjac(j,j)
      do  i = j,m
          wa4(i) = wa4(i) + fjac(i,j)*temp
      end do
120    fjac(j,j) = wa1(j)
      qtf(j) = wa4(j)
    end do

    !        compute the norm of the scaled gradient.

    gnorm = zero
    if (fnorm == zero) GO TO 170
    do  j = 1,n
      l = ipvt(j)
      if (wa2(l) == zero) cycle
      sum = zero
      do  i = 1,j
          sum = sum + fjac(i,j)*(qtf(i)/fnorm)
      end do
      gnorm = max(gnorm,abs(sum/wa2(l)))
    end do

    !        test for convergence of the gradient norm.

170 if (gnorm <= gtol) info = 4
    if (info /= 0) GO TO 300

    !        rescale if necessary.

    if (mode == 2) GO TO 200
    do  j = 1,n
      diag(j) = max(diag(j),wa2(j))
    end do

    !        beginning of the inner loop.

    !           determine the levenberg-marquardt parameter.

200 call lmpar(n,fjac,ipvt,diag,qtf,delta,par,wa1,wa2)

    !           store the direction p and x + p. calculate the norm of p.

    do  j = 1,n
      wa1(j) = -wa1(j)
      wa2(j) = x(j) + wa1(j)
      wa3(j) = diag(j)*wa1(j)
    end do
    !pnorm = enorm(n,wa3)
    pnorm = enorm3(n,wa3)

    !           on the first iteration,adjust the initial step bound.

    if (iter == 1) delta = min(delta,pnorm)

    !           evaluate the function at x + p and calculate its norm.

    iflag = 1
    call fcn(m,n,wa2,wa4,fjac,iflag)
    nfev = nfev + 1
    if (iflag < 0) GO TO 300
    !fnorm1 = enorm(m,wa4)
    fnorm1 = enorm3(m,wa4)

    !           compute the scaled actual reduction.

    actred = -one
    if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

    !           compute the scaled predicted reduction and
    !           the scaled directional derivative.

    do  j = 1,n
      wa3(j) = zero
      l = ipvt(j)
      temp = wa1(l)
      do  i = 1,j
          wa3(i) = wa3(i) + fjac(i,j)*temp
      end do
    end do
    !temp1 = enorm(n,wa3)/fnorm
    temp1 = enorm3(n,wa3)/fnorm
    temp2 = (sqrt(par)*pnorm)/fnorm
    prered = temp1**2 + temp2**2/p5
    dirder = -(temp1**2 + temp2**2)

    !           compute the ratio of the actual to the predicted reduction.

    ratio = zero
    if (prered /= zero) ratio = actred/prered

    !           update the step bound.

    if (ratio <= p25) then
      if (actred >= zero) temp = p5
      if (actred < zero) temp = p5*dirder/(dirder + p5*actred)
      if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
      delta = temp*min(delta,pnorm/p1)
      par = par/temp
    else
      if (par /= zero .and. ratio < p75) GO TO 260
      delta = pnorm/p5
      par = p5*par
    end if

    !           test for successful iteration.

260 if (ratio < p0001) GO TO 290

    !           successful iteration. update x,fvec,and their norms.

    do  j = 1,n
      x(j) = wa2(j)
      wa2(j) = diag(j)*x(j)
    end do
    fvec(1:m) = wa4(1:m)
    !xnorm = enorm(n,wa2)
    xnorm = enorm3(n,wa2)
    fnorm = fnorm1
    iter = iter + 1

    !           tests for convergence.

290 if (abs(actred) <= ftol .and. prered <= ftol .and. p5*ratio <= one) info = 1
    if (delta <= xtol*xnorm) info = 2
    if (abs(actred) <= ftol .and. prered <= ftol  &
        .and. p5*ratio <= one .and. info == 2) info = 3
    if (info /= 0) GO TO 300

    !           tests for termination and stringent tolerances.

    if (nfev >= maxfev) info = 5
    if (abs(actred) <= epsmch .and. prered <= epsmch  &
        .and. p5*ratio <= one) info = 6
    if (delta <= epsmch*xnorm) info = 7
    if (gnorm <= epsmch) info = 8
    if (info /= 0) GO TO 300

    !           end of the inner loop. repeat if iteration unsuccessful.

    if (ratio < p0001) GO TO 200

    !        end of the outer loop.

    GO TO 30

    !     termination,either normal or user imposed.

300 if (iflag < 0) info = iflag
    iflag = 0
    if (nprint > 0) call fcn(m,n,x,fvec,fjac,iflag)
    return

    !     last card of subroutine lmder.

  end subroutine lmder


  subroutine lmpar(n,r,ipvt,diag,qtb,delta,par,x,sdiag)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:46:12

    ! N.B. Arguments LDR,WA1 & WA2 have been removed.

    integer,intent(IN)::n
    real(dp),dimension(:,:),intent(INOUT)::r
    integer,dimension(:),intent(IN)::ipvt
    real(dp),dimension(:),intent(IN)::diag,qtb
    real(dp),intent(IN)::delta
    real(dp),intent(OUT)::par
    real(dp),dimension(:),intent(OUT)::x,sdiag


    !     **********

    !     subroutine lmpar

    !     given an m by n matrix a,an n by n nonsingular diagonal
    !     matrix d,an m-vector b,and a positive number delta,
    !     the problem is to determine a value for the parameter
    !     par such that if x solves the system

    !           a*x = b ,    sqrt(par)*d*x = 0 ,

    !     in the least squares sense,and dxnorm is the euclidean
    !     norm of d*x,then either par is zero and

    !           (dxnorm-delta) <= 0.1*delta ,

    !     or par is positive and

    !           abs(dxnorm-delta) <= 0.1*delta .

    !     this subroutine completes the solution of the problem
    !     if it is provided with the necessary information from the
    !     qr factorization,with column pivoting,of a. that is,if
    !     a*p = q*r,where p is a permutation matrix,q has orthogonal
    !     columns,and r is an upper triangular matrix with diagonal
    !     elements of nonincreasing magnitude,then lmpar expects
    !     the full upper triangle of r,the permutation matrix p,
    !     and the first n components of (q transpose)*b. on output
    !     lmpar also provides an upper triangular matrix s such that

    !            t   t                   t
    !           p *(a *a + par*d*d)*p = s *s .

    !     s is employed within lmpar and may be of separate interest.

    !     only a few iterations are generally needed for convergence
    !     of the algorithm. if,however,the limit of 10 iterations
    !     is reached,then the output par will contain the best
    !     value obtained so far.

    !     the subroutine statement is

    !       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2)

    !     where

    !       n is a positive integer input variable set to the order of r.

    !       r is an n by n array. on input the full upper triangle
    !         must contain the full upper triangle of the matrix r.
    !         on output the full upper triangle is unaltered,and the
    !         strict lower triangle contains the strict upper triangle
    !         (transposed) of the upper triangular matrix s.

    !       ldr is a positive integer input variable not less than n
    !         which specifies the leading dimension of the array r.

    !       ipvt is an integer input array of length n which defines the
    !         permutation matrix p such that a*p = q*r. column j of p
    !         is column ipvt(j) of the identity matrix.

    !       diag is an input array of length n which must contain the
    !         diagonal elements of the matrix d.

    !       qtb is an input array of length n which must contain the first
    !         n elements of the vector (q transpose)*b.

    !       delta is a positive input variable which specifies an upper
    !         bound on the euclidean norm of d*x.

    !       par is a nonnegative variable. on input par contains an
    !         initial estimate of the levenberg-marquardt parameter.
    !         on output par contains the final estimate.

    !       x is an output array of length n which contains the least
    !         squares solution of the system a*x = b,sqrt(par)*d*x = 0,
    !         for the output par.

    !       sdiag is an output array of length n which contains the
    !         diagonal elements of the upper triangular matrix s.

    !       wa1 and wa2 are work arrays of length n.

    !     subprograms called

    !       minpack-supplied ... dpmpar,enorm,qrsolv

    !       fortran-supplied ... ABS,MAX,MIN,SQRT

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    integer,dimension(2)::nldr
    integer::iter,j,k,l,nsing,ldr
    real(dp)::dxnorm,dwarf,fp,gnorm,parc,parl,paru,sum2,temp
    real(dp),dimension(n)::wa1,wa2
    real(dp),parameter::p1 = 0.1_dp,p001 = 0.001_dp

    nldr=shape(r)
    ldr=nldr(1)
    !
    !  DWARF is the smallest positive magnitude.
    !
    dwarf = tiny ( dwarf )
    !
    !  Compute and store in X the Gauss-Newton direction.
    !
    !  If the jacobian is rank-deficient, obtain a least squares solution.
    !
    nsing = n

    do j = 1, n
      wa1(j) = qtb(j)
      if ( r(j,j) == zero .and. nsing == n ) then
          nsing = j - 1
      end if
      if ( nsing < n ) then
          wa1(j) = zero
      end if
    end do

    do k = 1, nsing
      j = nsing - k + 1
      wa1(j) = wa1(j) / r(j,j)
      temp = wa1(j)
      wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
    end do

    do j = 1, n
      l = ipvt(j)
      x(l) = wa1(j)
    end do
    !
    !  Initialize the iteration counter.
    !  Evaluate the function at the origin, and test
    !  for acceptance of the Gauss-Newton direction.
    !
    iter = 0
    wa2(1:n) = diag(1:n) * x(1:n)
    !dxnorm = enorm ( n, wa2 )
    dxnorm = enorm3( n, wa2 )
    fp = dxnorm - delta

    if ( fp <= p1 * delta ) then
      if ( iter == 0 ) then
          par = zero
      end if
      return
    end if
    !
    !  If the jacobian is not rank deficient, the Newton
    !  step provides a lower bound, PARL, for the zero of
    !  the function.
    !
    !  Otherwise set this bound to zero.
    !
    parl = zero

    if ( n <= nsing ) then

      do j = 1, n
          l = ipvt(j)
          wa1(j) = diag(l) * ( wa2(l) / dxnorm )
      end do

      do j = 1, n
          sum2 = dot_product ( wa1(1:j-1), r(1:j-1,j) )
          wa1(j) = ( wa1(j) - sum2 ) / r(j,j)
      end do

      !temp = enorm ( n, wa1 )
      temp = enorm3( n, wa1 )
      parl = ( ( fp / delta ) / temp ) / temp

    end if
    !
    !  Calculate an upper bound, PARU, for the zero of the function.
    !
    do j = 1, n
      sum2 = dot_product ( qtb(1:j), r(1:j,j) )
      l = ipvt(j)
      wa1(j) = sum2 / diag(l)
    end do

    !gnorm = enorm ( n, wa1 )
    gnorm = enorm3( n, wa1 )
    paru = gnorm / delta

    if ( paru == zero ) then
      paru = dwarf / min ( delta, p1 )
    end if
    !
    !  If the input PAR lies outside of the interval (PARL, PARU),
    !  set PAR to the closer endpoint.
    !
    par = max ( par, parl )
    par = min ( par, paru )
    if ( par == zero ) then
      par = gnorm / dxnorm
    end if
    !
    !  Beginning of an iteration.
    !
    do

      iter = iter + 1
      !
      !  Evaluate the function at the current value of PAR.
      !
      if ( par == zero ) then
          par = max ( dwarf, p001 * paru )
      end if

      wa1(1:n) = sqrt ( par ) * diag(1:n)

      !call qrsolv ( n, r, ldr, ipvt, wa1, qtb, x, sdiag )
      call qrsolv ( n, r, ipvt, wa1, qtb, x, sdiag )

      wa2(1:n) = diag(1:n) * x(1:n)
      !dxnorm = enorm ( n, wa2 )
      dxnorm = enorm3( n, wa2 )
      temp = fp
      fp = dxnorm - delta
      !
      !  If the function is small enough, accept the current value of PAR.
      !
      if ( abs ( fp ) <= p1 * delta ) then
          exit
      end if
      !
      !  Test for the exceptional cases where PARL
      !  is zero or the number of iterations has reached 10.
      !
      if ( parl == zero .and. fp <= temp .and. temp < zero ) then
          exit
      else if ( iter == 10 ) then
          exit
      end if
      !
      !  Compute the Newton correction.
      !
      do j = 1, n
          l = ipvt(j)
          wa1(j) = diag(l) * ( wa2(l) / dxnorm )
      end do

      do j = 1, n
          wa1(j) = wa1(j) / sdiag(j)
          temp = wa1(j)
          wa1(j+1:n) = wa1(j+1:n) - r(j+1:n,j) * temp
      end do

      !temp = enorm ( n, wa1 )
      temp = enorm3( n, wa1 )
      parc = ( ( fp / delta ) / temp ) / temp
      !
      !  Depending on the sign of the function, update PARL or PARU.
      !
      if ( zero < fp ) then
          parl = max ( parl, par )
      else if ( fp < zero ) then
          paru = min ( paru, par )
      end if
      !
      !  Compute an improved estimate for PAR.
      !
      par = max ( parl, par + parc )
      !
      !  End of an iteration.
      !
    end do
    !
    !  Termination.
    !
    if ( iter == 0 ) then
      par = zero
    end if

    return

  end subroutine lmpar



  subroutine qrfac(m,n,a,pivot,ipvt,rdiag,acnorm)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:46:17

    ! N.B. Arguments LDA,LIPVT & WA has been removed.

    integer,intent(IN)::m
    integer,intent(IN)::n
    real(dp),intent(INOUT)::a(:,:)
    logical,intent(IN)::pivot
    !integer,intent(OUT)::ipvt(:)
    integer,dimension(:),intent(out)::ipvt
    real(dp),dimension(:),intent(OUT)::rdiag,acnorm

    !     **********

    !     subroutine qrfac

    !     this subroutine uses householder transformations with column
    !     pivoting (optional) to compute a qr factorization of the
    !     m by n matrix a. that is,qrfac determines an orthogonal
    !     matrix q,a permutation matrix p,and an upper trapezoidal
    !     matrix r with diagonal elements of nonincreasing magnitude,
    !     such that a*p = q*r. the householder transformation for
    !     column k,k = 1,2,...,min(m,n),is of the form

    !                           t
    !           i - (1/u(k))*u*u

    !     where u has zeros in the first k-1 positions. the form of
    !     this transformation and the method of pivoting first
    !     appeared in the corresponding linpack subroutine.

    !     the subroutine statement is

    !       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)

    !     where

    !       m is a positive integer input variable set to the number of rows of a.

    !       n is a positive integer input variable set to the number
    !         of columns of a.

    !       a is an m by n array. on input a contains the matrix for
    !         which the qr factorization is to be computed.  on output
    !         the strict upper trapezoidal part of a contains the strict
    !         upper trapezoidal part of r,and the lower trapezoidal
    !         part of a contains a factored form of q (the non-trivial
    !         elements of the u vectors described above).

    !       lda is a positive integer input variable not less than m
    !         which specifies the leading dimension of the array a.

    !       pivot is a logical input variable.  if pivot is set true,
    !         then column pivoting is enforced.  if pivot is set false,
    !         then no column pivoting is done.

    !       ipvt is an integer output array of length lipvt.  ipvt
    !         defines the permutation matrix p such that a*p = q*r.
    !         column j of p is column ipvt(j) of the identity matrix.
    !         if pivot is false,ipvt is not referenced.

    !       lipvt is a positive integer input variable. if pivot is false,
    !         then lipvt may be as small as 1. if pivot is true,then
    !         lipvt must be at least n.

    !       rdiag is an output array of length n which contains the
    !         diagonal elements of r.

    !       acnorm is an output array of length n which contains the
    !         norms of the corresponding columns of the input matrix a.
    !         If this information is not needed,then acnorm can coincide
    !         with rdiag.

    !       wa is a work array of length n. if pivot is false,then wa
    !         can coincide with rdiag.

    !     subprograms called

    !       minpack-supplied ... dpmpar,enorm

    !       fortran-supplied ... MAX,SQRT,MIN

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    integer::j,k,kmax,minmn,i4_temp
    real(dp)::ajnorm,epsmch,temp
    real(dp),dimension(:),allocatable::r8_temp
    real(dp),dimension(n)::wa
    real(dp),parameter::p05 = 0.05_dp

    !     epsmch is the machine precision.

    epsmch = epsilon(zero)
    !
    !  Compute the initial column norms and initialize several arrays.
    !
    do j = 1, n
      !acnorm(j) = enorm ( m, a(1:m,j) )
      acnorm(j) = enorm3( m, a(1:m,j) )
    end do

    rdiag(1:n) = acnorm(1:n)
    wa(1:n) = acnorm(1:n)

    if ( pivot ) then
      do j = 1, n
          ipvt(j) = j
      end do
    end if

    !
    !  Reduce A to R with Householder transformations.
    !
    minmn = min ( m, n )

    allocate(r8_temp(m))
    do j = 1, minmn
      !
      !  Bring the column of largest norm into the pivot position.
      !
      if ( pivot ) then

          kmax = j

          do k = j, n
            if ( rdiag(k) > rdiag(kmax) ) then
                kmax = k
            end if
          end do

          if ( kmax /= j ) then

            r8_temp(1:m) = a(1:m,j)
            a(1:m,j)     = a(1:m,kmax)
            a(1:m,kmax)  = r8_temp(1:m)

            rdiag(kmax) = rdiag(j)
            wa(kmax) = wa(j)

            i4_temp    = ipvt(j)
            ipvt(j)    = ipvt(kmax)
            ipvt(kmax) = i4_temp

          end if

      end if
      !
      !  Compute the Householder transformation to reduce the
      !  J-th column of A to a multiple of the J-th unit vector.
      !
      !ajnorm = enorm ( m-j+1, a(j,j) ) ! original ... a(j,j) error Rank mismatch
      !ajnorm = enorm ( m-j+1, a(j:,j) )
      ajnorm = enorm3( m-j+1, a(j:,j) )

      if ( ajnorm /= zero ) then

          if ( a(j,j) < zero ) ajnorm = -ajnorm

          a(j:m,j) = a(j:m,j) / ajnorm
          a(j,j) = a(j,j) + one
          !
          !  Apply the transformation to the remaining columns and update the norms.
          !
          do k = j+1, n

            temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)

            a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

            if ( pivot .and. rdiag(k) /= zero ) then

                temp = a(j,k) / rdiag(k)
                rdiag(k) = rdiag(k) * sqrt ( max ( zero, one-temp**2 ) )

                if ( p05 * ( rdiag(k) / wa(k) )**2 <= epsmch ) then
                  !rdiag(k) = enorm ( m-j, a(j+1:,k) )
                  rdiag(k) = enorm3( m-j, a(j+1:,k) )
                  wa(k) = rdiag(k)
                end if

            end if

          end do

      end if

      rdiag(j) = -ajnorm

    end do
    deallocate(r8_temp)

    return
  end subroutine qrfac



  subroutine qrsolv(n,r,ipvt,diag,qtb,x,sdiag)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:46:21

    ! N.B. Arguments LDR & WA has been removed.

    integer,intent(IN)::n
    real(dp),intent(INOUT)::r(:,:)
    integer,intent(IN)::ipvt(:)
    real(dp),intent(IN)::diag(:)
    real(dp),intent(IN)::qtb(:)
    real(dp),intent(OUT)::x(:)
    real(dp),intent(OUT)::sdiag(:)


    !     **********

    !     subroutine qrsolv

    !     given an m by n matrix a,an n by n diagonal matrix d,
    !     and an m-vector b,the problem is to determine an x which
    !     solves the system

    !           a*x = b ,    d*x = 0 ,

    !     in the least squares sense.

    !     this subroutine completes the solution of the problem
    !     if it is provided with the necessary information from the
    !     qr factorization,with column pivoting,of a. that is,if
    !     a*p = q*r,where p is a permutation matrix,q has orthogonal
    !     columns,and r is an upper triangular matrix with diagonal
    !     elements of nonincreasing magnitude,then qrsolv expects
    !     the full upper triangle of r,the permutation matrix p,
    !     and the first n components of (q transpose)*b. the system
    !     a*x = b,d*x = 0,is then equivalent to

    !                  t       t
    !           r*z = q *b , p *d*p*z = 0 ,

    !     where x = p*z. if this system does not have full rank,
    !     then a least squares solution is obtained. on output qrsolv
    !     also provides an upper triangular matrix s such that

    !            t   t               t
    !           p *(a *a + d*d)*p = s *s .

    !     s is computed within qrsolv and may be of separate interest.

    !     the subroutine statement is

    !       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)

    !     where

    !       n is a positive integer input variable set to the order of r.

    !       r is an n by n array. on input the full upper triangle
    !         must contain the full upper triangle of the matrix r.
    !         on output the full upper triangle is unaltered,and the
    !         strict lower triangle contains the strict upper triangle
    !         (transposed) of the upper triangular matrix s.

    !       ldr is a positive integer input variable not less than n
    !         which specifies the leading dimension of the array r.

    !       ipvt is an integer input array of length n which defines the
    !         permutation matrix p such that a*p = q*r. column j of p
    !         is column ipvt(j) of the identity matrix.

    !       diag is an input array of length n which must contain the
    !         diagonal elements of the matrix d.

    !       qtb is an input array of length n which must contain the first
    !         n elements of the vector (q transpose)*b.

    !       x is an output array of length n which contains the least
    !         squares solution of the system a*x = b,d*x = 0.

    !       sdiag is an output array of length n which contains the
    !         diagonal elements of the upper triangular matrix s.

    !       wa is a work array of length n.

    !     subprograms called

    !       fortran-supplied ... ABS,SQRT

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    integer::i,j,jp1,k,kp1,l,nsing
    real(dp)::COS,cotan,qtbpj,SIN,sum,TAN,temp,wa(n)
    real(dp),parameter::p5 = 0.5_dp,p25 = 0.25_dp,zero = 0.0_dp

    !     copy r and (q transpose)*b to preserve input and initialize s.
    !     in particular,save the diagonal elements of r in x.

    do  j = 1,n
      do  i = j,n
          r(i,j) = r(j,i)
      end do
      x(j) = r(j,j)
      wa(j) = qtb(j)
    end do

    !     eliminate the diagonal matrix d using a givens rotation.

    do  j = 1,n

      !        prepare the row of d to be eliminated,locating the
      !        diagonal element using p from the qr factorization.

      l = ipvt(j)
      if (diag(l) == zero) cycle
      sdiag(j:n) = zero
      sdiag(j) = diag(l)

      !        the transformations to eliminate the row of d
      !        modify only a single element of (q transpose)*b
      !        beyond the first n,which is initially zero.

      qtbpj = zero
      do  k = j,n

          !           determine a givens rotation which eliminates the
          !           appropriate element in the current row of d.

          if (sdiag(k) == zero) cycle
          if (abs(r(k,k)) < abs(sdiag(k))) then
            cotan = r(k,k)/sdiag(k)
            SIN = p5/sqrt(p25 + p25*cotan**2)
            COS = SIN*cotan
          else
            TAN = sdiag(k)/r(k,k)
            COS = p5/sqrt(p25 + p25*TAN**2)
            SIN = COS*TAN
          end if

          !           compute the modified diagonal element of r and
          !           the modified element of ((q transpose)*b,0).

          r(k,k) = COS*r(k,k) + SIN*sdiag(k)
          temp = COS*wa(k) + SIN*qtbpj
          qtbpj = -SIN*wa(k) + COS*qtbpj
          wa(k) = temp

          !           accumulate the tranformation in the row of s.

          kp1 = k + 1
          do  i = kp1,n
            temp = COS*r(i,k) + SIN*sdiag(i)
            sdiag(i) = -SIN*r(i,k) + COS*sdiag(i)
            r(i,k) = temp
          end do
      end do

      !        store the diagonal element of s and restore
      !        the corresponding diagonal element of r.

      sdiag(j) = r(j,j)
      r(j,j) = x(j)
    end do

    !     solve the triangular system for z. if the system is
    !     singular,then obtain a least squares solution.

    nsing = n
    do  j = 1,n
      if (sdiag(j) == zero .and. nsing == n) nsing = j - 1
      if (nsing < n) wa(j) = zero
    end do

    do  k = 1,nsing
      j = nsing - k + 1
      sum = zero
      jp1 = j + 1
      do  i = jp1,nsing
          sum = sum + r(i,j)*wa(i)
      end do
      wa(j) = (wa(j) - sum)/sdiag(j)
    end do

    !     permute the components of z back to components of x.

    do  j = 1,n
      l = ipvt(j)
      x(l) = wa(j)
    end do
    return

    !     last card of subroutine qrsolv.

  end subroutine qrsolv


  function enorm(n,x) result(fn_val)

    ! substitute the original enorm ==> enorm2
    integer,intent(in)::n
    real(dp),dimension(:),intent(in)::x
    real(dp)::fn_val

    fn_val = sqrt( sum( x(1:n)**2 ) )
    
    return
  end function enorm

  function enorm2(n,x) result(fn_val)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:34

    integer,intent(IN)::n
    real(dp),intent(IN)::x(:)
    real(dp)::fn_val


    !     **********

    !     function enorm

    !     given an n-vector x,this function calculates the euclidean norm of x.

    !     the euclidean norm is computed by accumulating the sum of squares in
    !     three different sums.  The sums of squares for the small and large
    !     components are scaled so that no overflows occur.  Non-destructive
    !     underflows are permitted.  Underflows and overflows do not occur in the
    !     computation of the unscaled sum of squares for the intermediate
    !     components.  The definitions of small,intermediate and large components
    !     depend on two constants,rdwarf and rgiant.  The main restrictions on
    !     these constants are that rdwarf**2 not underflow and rgiant**2 not
    !     overflow.  The constants given here are suitable for every known computer.

    !     the function statement is

    !       REAL (dp) function enorm(n,x)

    !     where

    !       n is a positive integer input variable.

    !       x is an input array of length n.

    !     subprograms called

    !       fortran-supplied ... ABS,SQRT

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    integer::i
    real(dp)::agiant,floatn,s1,s2,s3,xabs,x1max,x3max
    real(dp),parameter::rdwarf = 3.834E-20_dp, &
        &rgiant = 1.304E+19_dp

    s1 = zero
    s2 = zero
    s3 = zero
    x1max = zero
    x3max = zero
    floatn = real(n,dp)
    agiant = rgiant/floatn
    do  i = 1,n
      xabs = abs(x(i))
      if (xabs > rdwarf .and. xabs < agiant) GO TO 70
      if (xabs <= rdwarf) GO TO 30

      !              sum for large components.

      if (xabs <= x1max) GO TO 10
      s1 = one + s1*(x1max/xabs)**2
      x1max = xabs
      GO TO 20

10     s1 = s1 + (xabs/x1max)**2

20     GO TO 60

      !              sum for small components.

30     if (xabs <= x3max) GO TO 40
      s3 = one + s3*(x3max/xabs)**2
      x3max = xabs
      GO TO 60

40     if (xabs /= zero) s3 = s3 + (xabs/x3max)**2

60     cycle

      !           sum for intermediate components.

70     s2 = s2 + xabs**2
    end do

    !     calculation of norm.

    if (s1 == zero) GO TO 100
    fn_val = x1max*sqrt(s1 + (s2/x1max)/x1max)
    GO TO 120

100 if (s2 == zero) GO TO 110
    if (s2 >= x3max) fn_val = sqrt(s2*(one + (x3max/s2)*(x3max*s3)))
    if (s2 < x3max) fn_val = sqrt(x3max*((s2/x3max) + (x3max*s3)))
    GO TO 120

110 fn_val = x3max*sqrt(s3)

120 return

    !     last card of function enorm.

  end function enorm2

  function enorm3 ( n, x )
    
    !*****************************************************************************80
    !
    !! ENORM2 computes the Euclidean norm of a vector.
    !
    !  Discussion:
    !
    !    This routine was named ENORM.  It has been renamed "ENORM2",
    !    and a simplified routine has been substituted.
    !
    !    The Euclidean norm is computed by accumulating the sum of
    !    squares in three different sums.  The sums of squares for the
    !    small and large components are scaled so that no overflows
    !    occur.  Non-destructive underflows are permitted.  Underflows
    !    and overflows do not occur in the computation of the unscaled
    !    sum of squares for the intermediate components.
    !
    !    The definitions of small, intermediate and large components
    !    depend on two constants, RDWARF and RGIANT.  The main
    !    restrictions on these constants are that RDWARF**2 not
    !    underflow and RGIANT**2 not overflow.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, is the length of the vector.
    !
    !    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
    !
    !    Output, real ( kind = 8 ) ENORM2, the Euclidean norm of the vector.
    !
    integer::n

    real(dp)::agiant
    real(dp)::enorm3
    integer::i
    real(dp)::rdwarf
    real(dp)::rgiant
    real(dp)::s1
    real(dp)::s2
    real(dp)::s3
    real(dp)::x(n)
    real(dp)::xabs
    real(dp)::x1max
    real(dp)::x3max

    rdwarf = sqrt ( tiny ( rdwarf ) )
    rgiant = sqrt ( huge ( rgiant ) )

    s1 = zero
    s2 = zero
    s3 = zero
    x1max = zero
    x3max = zero
    agiant = rgiant / real ( n, dp)

    do i = 1, n

      xabs = abs ( x(i) )
      if ( xabs <= rdwarf ) then
          if ( x3max < xabs ) then
            s3 = one + s3 * ( x3max / xabs )**2
            x3max = xabs
          else if ( xabs /= zero ) then
            s3 = s3 + ( xabs / x3max )**2
          end if
      else if ( agiant <= xabs ) then
          if ( x1max < xabs ) then
            s1 = one + s1 * ( x1max / xabs )**2
            x1max = xabs
          else
            s1 = s1 + ( xabs / x1max )**2
          end if
      else
          s2 = s2 + xabs**2
      end if

    end do
    !
    !  Calculation of norm.
    !
    if ( s1 /= zero ) then
      enorm3 = x1max * sqrt ( s1 + ( s2 / x1max ) / x1max )
    else if ( s2 /= zero ) then
      if ( x3max <= s2 ) then
          enorm3 = sqrt ( s2 * ( one + ( x3max / s2 ) * ( x3max * s3 ) ) )
      else
          enorm3 = sqrt ( x3max * ( ( s2 / x3max ) + ( x3max * s3 ) ) )
      end if
    else
      enorm3 = x3max * sqrt ( s3 )
    end if

    return
  end function enorm3
  ! inserted enorm3 f90 version of 77 enorm3

  subroutine fdjac2_1(fcn,m,n,x,fvec,fjac,iflag,epsfcn)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:44

    ! N.B. Arguments LDFJAC & WA have been removed.

    integer,intent(IN)::m
    integer,intent(IN)::n
    real(dp),dimension(:),intent(INOUT)::x
    real(dp),dimension(:),intent(IN)::fvec
    real(dp),dimension(:,:),intent(OUT)::fjac    ! fjac(ldfjac,n)
    integer,intent(inout)::iflag
    real(dp),intent(IN)::epsfcn

    interface
      subroutine fcn(m,n,x,fvec,iflag)
        use constants
        implicit none
        !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)
        integer,intent(IN)::m,n
        real(dp),intent(IN)::x(:)
        real(dp),intent(OUT)::fvec(:)
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface
    
    !*****************************************************************************80
    !
    !! FDJAC2 estimates an M by N jacobian matrix using forward differences.
    !
    !  Discussion:
    !
    !    This subroutine computes a forward-difference approximation
    !    to the M by N jacobian matrix associated with a specified
    !    problem of M functions in N variables.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Parameters:
    !
    !    Input, external FCN, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !
    !      subroutine fcn ( m, n, x, fvec, iflag )
    !      integer::n
    !      real fvec(m)
    !      integer::iflag
    !      real x(n)
    !
    !    The value of IFLAG should not be changed by FCN unless
    !    the user wants to terminate execution of the routine.
    !    In this case set IFLAG to a negative integer.
    !
    !    Input, integer::M, is the number of functions.
    !
    !    Input, integer::N, is the number of variables.  
    !    N must not exceed M.
    !
    !    Input, real(dp)::X(N), the point where the jacobian is evaluated.
    !
    !    Input, real(dp)::FVEC(M), the functions evaluated at X.
    !
    !    Output, real(dp)::FJAC(LDFJAC,N), the M by N approximate
    !    jacobian matrix.
    !
    !    Input, integer::LDFJAC, the leading dimension of FJAC, 
    !    which must not be less than M.
    !
    !    Output, integer::IFLAG, is an error flag returned by FCN.  
    !    If FCN returns a nonzero value of IFLAG, then this routine returns 
    !    immediately to the calling program, with the value of IFLAG.
    !
    !    Input, real(dp)::EPSFCN, is used in determining a suitable
    !    step length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    
    integer::j
    real(dp)::eps,epsmch,h,temp
    real(dp),dimension(m)::wa(m)

    epsmch = epsilon ( epsmch )

    eps = sqrt ( max ( epsfcn, epsmch ) )

    do j = 1, n

      temp = x(j)
      h = eps * abs ( temp )
      if ( h == zero ) then
          h = eps
      end if

      x(j) = temp + h
      call fcn ( m, n, x, wa, iflag )

      if ( iflag < 0 ) exit

      x(j) = temp
      fjac(1:m,j) = ( wa(1:m) - fvec(1:m) ) / h

    end do

    return
  end subroutine fdjac2_1

  subroutine fdjac2_2(fcn,allpar,m,n,x,fvec,fjac,iflag,epsfcn)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:44

    ! N.B. Arguments LDFJAC & WA have been removed.

    integer,intent(IN)::m,n
    real(dp),dimension(:),intent(INOUT)::allpar,x
    real(dp),dimension(:),intent(IN)::fvec
    real(dp),dimension(:,:),intent(OUT)::fjac    ! fjac(ldfjac,n)
    integer,intent(inout)::iflag
    real(dp),intent(IN)::epsfcn

    interface
      subroutine fcn(allpar,m,n,x,fvec,iflag)
        use constants
        implicit none
        integer,intent(IN)::m,n
        real(dp),dimension(:),intent(IN)::allpar,x
        real(dp),dimension(:),intent(OUT)::fvec
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface

    integer::j
    real(dp)::eps,epsmch,h,temp
    real(dp),dimension(m)::wa

    epsmch = epsilon ( epsmch )
    eps = sqrt ( max ( epsfcn, epsmch ) )

    do j = 1, n

      temp = x(j)
      h = eps * abs ( temp )
      if ( h == zero ) h = eps

      x(j) = temp + h
      call fcn(allpar,m,n,x,wa,iflag)

      if ( iflag < 0 ) exit

      x(j) = temp
      fjac(1:m,j) = ( wa(1:m) - fvec(1:m) ) / h

    end do

    return
  end subroutine fdjac2_2

  subroutine fdjac2_3(fcn,allpar,m,n,x,RV_obs,T0_obs,fvec,fjac,iflag,epsfcn)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:44

    ! N.B. Arguments LDFJAC & WA have been removed.

    integer,intent(IN)::m,n
    real(dp),dimension(:),intent(inout)::allpar,x
    real(dp),dimension(:),intent(in)::RV_obs
    real(dp),dimension(:,:),intent(in)::T0_obs
    real(dp),dimension(:),intent(IN)::fvec
    real(dp),dimension(:,:),intent(OUT)::fjac    ! fjac(ldfjac,n)
    integer,intent(inout)::iflag
    real(dp),intent(IN)::epsfcn

    interface
      subroutine fcn(allpar,m,n,x,RV_obs,T0_obs,fvec,iflag)
        use constants
        implicit none
        integer,intent(IN)::m,n
        real(dp),dimension(:),intent(IN)::allpar,x,RV_obs
        real(dp),dimension(:,:),intent(in)::T0_obs
        real(dp),dimension(:),intent(OUT)::fvec
        integer,intent(inout)::iflag
      end subroutine fcn
    end interface

    integer::j
    real(dp)::eps,epsmch,h,temp
    real(dp),dimension(m)::wa

    epsmch = epsilon(zero)

    eps = sqrt(max(epsfcn,epsmch))

    do  j = 1,n

      temp = x(j)
      h = eps*abs(temp)
      if (h == zero) h = eps

      x(j) = temp + h
      call fcn(allpar,m,n,x,RV_obs,T0_obs,wa,iflag)

      if (iflag < 0) exit

      x(j) = temp
      fjac(1:m,j) = ( wa(1:m) - fvec(1:m) ) / h

    end do

    return
  end subroutine fdjac2_3

  subroutine covar(n,r,ldr,ipvt,tol,wa)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2013-01-12  Time: 18:09:08

    integer,intent(in)::n
    real(dp),dimension(:,:),intent(out)::r !(ldr,n)
    integer,intent(inout)::ldr
    integer,dimension(:),intent(in)::ipvt
    real(dp),intent(in)::tol
    real(dp),dimension(:),intent(out)::wa !(n)


    !     **********

    !     subroutine covar

    !     given an m by n matrix a,the problem is to determine
    !     the covariance matrix corresponding to a,defined as

    !              t
    !     inverse(a *a) .

    !     this subroutine completes the solution of the problem
    !     if it is provided with the necessary information from the
    !     qr factorization,with column pivoting,of a. that is,if
    !     a*p = q*r,where p is a permutation matrix,q has orthogonal
    !     columns,and r is an upper triangular matrix with diagonal
    !     elements of nonincreasing magnitude,then covar expects
    !     the full upper triangle of r and the permutation matrix p.
    !     the covariance matrix is then computed as

    !                t     t
    !     p*inverse(r *r)*p  .

    !     if a is nearly rank deficient,it may be desirable to compute
    !     the covariance matrix corresponding to the linearly independent
    !     columns of a. to define the numerical rank of a,covar uses
    !     the tolerance tol. if l is the largest integer such that

    !     abs(r(l,l)) .gt. tol*abs(r(1,1)) ,

    !     then covar computes the covariance matrix corresponding to
    !     the first l columns of r. for k greater than l,column
    !     and row ipvt(k) of the covariance matrix are set to zero.

    !     the subroutine statement is

    !     subroutine covar(n,r,ldr,ipvt,tol,wa)

    !     where

    !     n is a positive integer input variable set to the order of r.

    !     r is an n by n array. on input the full upper triangle must
    !     contain the full upper triangle of the matrix r. on output
    !     r contains the square symmetric covariance matrix.

    !     ldr is a positive integer input variable not less than n
    !     which specifies the leading dimension of the array r.

    !     ipvt is an integer input array of length n which defines the
    !     permutation matrix p such that a*p = q*r. column j of p
    !     is column ipvt(j) of the identity matrix.

    !     tol is a nonnegative input variable used to define the
    !     numerical rank of a in the manner described above.

    !     wa is a work array of length n.

    !     subprograms called

    !     fortran-supplied ... dabs

    !     argonne national laboratory. minpack project. august 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    integer::i,ii,j,jj,k,km1,l
    logical::sing
    real(dp)::temp,tolr
    real(dp),parameter::one=1._dp,zero=0._dp

    !     form the inverse of r in the full upper triangle of r.
    !write(*,'(a)')" covar 1"
    tolr = tol*abs(r(1,1))
    l = 0
    do  k = 1,n
      if (abs(r(k,k)) <= tolr) exit
      r(k,k) = one/r(k,k)
      km1 = k - 1
      if (km1 < 1) GO TO 30
      do  j = 1,km1
          temp = r(k,k)*r(j,k)
          r(j,k) = zero
          do  i = 1,j
            r(i,k) = r(i,k) - temp*r(i,j)
          end do
      end do
30     continue
      l = k
    end do
    ! 50 CONTINUE
    !write(*,'(a)')" covar 2"
    !     form the full upper triangle of the inverse of (r transpose)*r
    !     in the full upper triangle of r.

    if (l < 1) GO TO 110
    do  k = 1,l
      km1 = k - 1
      if (km1 < 1) GO TO 80
      do  j = 1,km1
          temp = r(j,k)
          do  i = 1,j
            r(i,j) = r(i,j) + temp*r(i,k)
          end do
      end do
80     continue
      temp = r(k,k)
      do  i = 1,k
          r(i,k) = temp*r(i,k)
      end do
    end do
110 continue
    !write(*,'(a)')" covar 3"
    !     form the full lower triangle of the covariance matrix
    !     in the strict lower triangle of r and in wa.

    do j = 1,n
      !write(*,'(a,i4,1x)',advance='no')" j = ",j
      jj = ipvt(j)
      !write(*,'(a,i4,1x)',advance='no')" jj = ipvt(j) = ",jj
      sing = j > l
      !write(*,'(a,L2,1x)',advance='no')" sing = ",sing
      do  i = 1,j
          if (sing) r(i,j) = zero
          ii = ipvt(i)
          if (ii > jj) r(ii,jj) = r(i,j)
          if (ii < jj) r(jj,ii) = r(i,j)
      end do
      !write(*,'(a,g25.14)',advance='no')" r(j,j) = ",r(j,j)
      wa(jj) = r(j,j)
    end do
    !write(*,'(a)')" covar 4"
    !     symmetrize the covariance matrix in r.
    do  j = 1,n
      do  i = 1,j
          r(i,j) = r(j,i)
      end do
      !write(*,'(a,i4,a,g25.14)')"j=",j," wa(j)=",wa(j)
      r(j,j) = wa(j)
    end do
    !write(*,'(a)')" covar 5"
    return

    !     last card of subroutine covar.

  end subroutine covar


end module Levenberg_Marquardt
