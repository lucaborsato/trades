module Genetic_Algorithm
  use constants,only:dp,sprec,TOLERANCE,zero,one
  use parameters,only:path,wrtAll,seed_pik,ctrl,dof,inv_dof,ndata,nfit,npar,parid,minpar,maxpar
  use parameters_conversion
  use init_trades,only:get_unit
  use random_trades
  use convert_type,only:string
  use fitness_module
  implicit none

  !     Common block to make iseed visible to rninit (and to save
  !     it between calls)
  ! COMMON /rnseed/ iseed
  integer, save  :: iseed

  contains

  ! initializes the variables to store best simulations
  subroutine bestpik_init(ng,idbest,bestpik)
    use parameters,only:nfit
    integer,intent(in)::ng
    real(dp),dimension(:,:),allocatable,intent(out)::bestpik
    integer,dimension(:),allocatable,intent(out)::idbest
    integer::nstore

    nstore=nfit+3
    if(.not.allocated(bestpik)) allocate(bestpik(nstore,ng),idbest(ng))
    bestpik=zero

    return
  end subroutine bestpik_init

  ! calculates the % of the simulations
  function getStat(i,n) result(out)
    integer::out
    integer,intent(in)::i,n
    real(dp)::fi,fn
    
    fi=real(i,dp)
    fn=real(n,dp)
    out=int((fi/fn)*100._dp)

    return
  end function getStat

    ! function to be use with PIKAIA genetic algorithm
!   function fpik(n,allpar,par) result(fn_val)
!     use ode_run,only:ode_lm
!     real(dp)::fn_val
!     integer,intent(in)::n
!     real(dp),dimension(:),intent(in)::allpar,par
!     real(dp),dimension(:),allocatable::wpar,wall
!     real(dp),dimension(:),allocatable::resw,resw2
! !     real(dp)::chi2
!     integer::iflag
! 
!     iflag=1
!     allocate(wpar(n),wall(npar))
!     wall=allpar
!     call norm2par(par,wpar,wall)
!     allocate(resw(ndata),resw2(ndata))
!     resw=0._dp
!     call ode_lm(wall,ndata,n,wpar,resw,iflag)
!     resw2=resw*resw
! !     chi2=sum(resw2)
! !     fn_val=1._dp/chi2
! !     fn_val=real(dof,dp)/sum(resw2)
!     fn_val=one/sum(resw2)
!     deallocate(wpar,wall,resw,resw2)
! 
!     return
!   end function fpik

  ! function to compute the fitness
  ! comment: now it uses the function in fitness_module module
!   function fpik(n,all_par,fitting_parameters_in) result(inv_fitness)
!     use ode_run,only:ode_lm
!     use derived_parameters_mod
!     real(dp)::inv_fitness
!     integer,intent(in)::n
!     real(dp),intent(in),dimension(:)::all_par
!     real(dp),intent(in),dimension(n)::fitting_parameters_in
!     real(dp),dimension(n)::fitting_parameters
!     real(dp)::fitness
!     logical::check
!     real(dp),dimension(:),allocatable::run_all_par
!     real(dp),dimension(:),allocatable::resw
!     integer::i,iflag
!     logical::check_status
!     
!     iflag=1
!     check=.true.
!     check_status=.true.
!     
!     ! needed only by PIKAIA: conversion from [0-1] to [minpar-maxpar] boundaries
!     allocate(run_all_par(npar))
!     run_all_par=all_par
!     call norm2par(fitting_parameters_in,fitting_parameters,run_all_par)
!     
!     checkloop: do i=1,nfit
!       if(fitting_parameters(i).lt.minpar(i))then
!         check=.false.
!         exit checkloop
!       else if(fitting_parameters(i).gt.maxpar(i))then
!         check=.false.
!         exit checkloop
!       end if
!     end do checkloop
! 
!     if(check)then
!     
!       if(check_derived) check_status=check_derived_parameters(fitting_parameters)
!       if(fix_derived) call fix_derived_parameters(fitting_parameters,run_all_par,check_status)
! 
!       if(check_status)then
!         allocate(resw(ndata))
!         resw=zero
!         call ode_lm(run_all_par,ndata,nfit,fitting_parameters,resw,iflag)
!         fitness=sum(resw*resw)
!         ! resw t.c. sum(resw^2) = fitness = Chi2r*K_chi2r + Chi2wr*K_chi2wr
!         deallocate(resw)
!         if (fitness.ge.resmax)then
! !           check=.false.
!           fitness=resmax
!         end if
!         
!       else ! check_status
!         fitness=resmax
!       end if
!     else
!       fitness=resmax
!     end if
! 
!     deallocate(run_all_par)
!     inv_fitness=one/fitness
!     
!     return
!   end function fpik

  ! function called by PIKAIA, based on fitness function in fitness_module
  function fpik(n,all_par,fitting_parameters_in) result(inv_fitness)
    real(dp)::inv_fitness
    integer,intent(in)::n
    real(dp),intent(in),dimension(:)::all_par
    real(dp),intent(in),dimension(n)::fitting_parameters_in

    real(dp),dimension(n)::fitting_parameters
    real(dp)::fitness
    real(dp),dimension(:),allocatable::run_all_par
    
    ! needed only by PIKAIA: conversion from [0-1] to [minpar-maxpar] boundaries
    allocate(run_all_par(npar))
    run_all_par=all_par
!     call norm2par(fitting_parameters_in,fitting_parameters,run_all_par)
    call norm2par(fitting_parameters_in,fitting_parameters)
    
    fitness=bound_fitness_function(run_all_par,fitting_parameters)
    inv_fitness=one/fitness
    
    deallocate(run_all_par)
    
    return
  end function fpik


  ! 2016-07-20 Luca Borsato own method to initialize population
  subroutine init_good_particles(nfit,npop,all_parameters,oldph)
    integer,intent(in)::nfit,npop            ! number of particle
    real(dp),dimension(:),intent(in)::all_parameters
    real(dp),dimension(:,:)::oldph ! particles
    real(dp),dimension(:),allocatable::xtemp,ptemp
    integer::iloop
    logical::check
    
    allocate(xtemp(nfit),ptemp(nfit))
    iloop=0
    initloop: do
        
        xtemp=zero
        ptemp=zero
        call random_number(xtemp)
        call norm2par(xtemp,ptemp)
        check=.true.
        check=check_only_boundaries(all_parameters,ptemp)
        if(check)then
!           write(*,*)' par (True) = ',ptemp
          iloop=iloop+1
          oldph(1:nfit,iloop)=xtemp(1:nfit)
          if(iloop.eq.npop)exit initloop
        end if
        
    end do initloop
    deallocate(xtemp,ptemp)
!     flush(6)
    
  end subroutine init_good_particles

  
  ! ------------------------------------------------------------------ !
  
  ! driver that initialize variabel, factors and calls the pikaia subroutine.
  subroutine ga_driver(iGlobal,ff,n,xall,x,f)
    !$ use omp_lib
    integer,intent(in)::iGlobal,n
    real(dp),dimension(:)::xall,x
    real(dp),intent(inout)::f

    integer::cpuid,STATUS

!     interface
!       function ff(n,xall,x) result(fn_val)
!         use constants,only:dp
!         implicit none
!         integer,intent(in)::n
!         real(dp),dimension(:),intent(in)::xall,x
!         real(dp)::fn_val
!       end function ff
!     end interface
    interface
      function ff(n,xall,x) result(fn_val)
        use constants,only:dp
        implicit none
        real(dp)::fn_val
        integer,intent(in)::n
        real(dp),intent(in),dimension(:)::xall
        real(dp),intent(in),dimension(n)::x
      end function ff
    end interface

    cpuid=1
    !$ cpuid = omp_get_thread_num() + 1

    !     First, initialize the random-number generator
    !seed = 123456
    write(*,*)
!     write(*,'(a,i8)')" Init random number with seed: ",seed_pik+iGlobal
! !     call rninit(seed_pik+iGlobal)
!     call init_random_seed(int(ctrl(1)),seed_pik+iGlobal)
!     flush(6)
    
    ! set initial fitting parameters to mean value of the boundaries
!     x = 0.5_dp*(minpar+maxpar)
    x=minpar
    
    call pikaia(ff,n,ctrl,xall,x,f,STATUS,iGlobal)
    !     Print the results
    write(*,'(a)') " in ga_driver "
    write(*,'(a,i3)') ' status: ', STATUS
    write(*,'(a,1000(g25.14))') '      x: ', x
    write(*,'(a,g25.14)') '      f: ', f
    write(*,20) ctrl
      
20  format('    ctrl: ', 6g11.6 / t11, 6g11.6)
    flush(6)

    return
  end subroutine ga_driver

  ! This is the PIKAIA subroutine
  ! modified by Luca Borsato 2014: added openMP directives and save/write simulations
  subroutine pikaia(ff,n,ctrl,xall,x,f,STATUS,iGlobal)
    !$ use omp_lib

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2001-07-09  Time: 15:54:13

    !====================================================================
    !  Optimization (maximization) of user-supplied "fitness" function ff
    !  over n-dimensional parameter space  x  using a basic genetic algorithm
    !  method.

    !  Paul Charbonneau & Barry Knapp
    !  High Altitude Observatory
    !  National Center for Atmospheric Research
    !  Boulder CO 80307-3000
    !  USA
    !  <paulchar@hao.ucar.edu>
    !  <knapp@hao.ucar.edu>

    !  Web site:
    !  http://www.hao.ucar.edu/public/research/si/pikaia/pikaia.html

    !  Version 1.0   [ 1995 December 01 ]

    !  Genetic algorithms are heuristic search techniques that incorporate in a
    !  computational setting, the biological notion of evolution by means of
    !  natural selection.  This subroutine implements the three basic operations
    !  of selection, crossover, and mutation, operating on "genotypes" encoded as
    !  strings.

    !  References:

    !     Charbonneau, Paul.  "Genetic Algorithms in Astronomy and Astrophysics."
    !        Astrophysical J. (Supplement), vol 101, in press (December 1995).

    !     Goldberg, David E.  Genetic Algorithms in Search, Optimization,
    !        & Machine Learning.  Addison-Wesley, 1989.

    !     Davis, Lawrence, ed.  Handbook of Genetic Algorithms.
    !        Van Nostrand Reinhold, 1991.
    !====================================================================
    !  USES: ff, urand, setctl, report, rnkpop, select, encode, decode,
    !        cross, mutate, genrep, stdrep, newpop, adjmut

    integer,intent(in)::n
    real(dp),intent(in out)::ctrl(12)
    real(dp),dimension(:),intent(in)::xall
    real(dp),dimension(:),intent(out)::x
    real(dp),intent(out)::f
    integer,intent(out)::STATUS
    integer,intent(in)::iGlobal

!     interface
!       function ff(n,xall,x) result(fn_val)
!         use constants,only:dp
!         implicit none
!         integer,intent(in)::n
!         real(dp),dimension(:),intent(in)::xall,x
!         real(dp)::fn_val
!       end function ff
!     end interface
    interface
      function ff(n,xall,x) result(fn_val)
        use constants,only:dp
        implicit none
        real(dp)::fn_val
        integer,intent(in)::n
        real(dp),intent(in),dimension(:)::xall
        real(dp),intent(in),dimension(n)::x
      end function ff
    end interface

    ! EXTERNAL ff

    !  Input:
    !  o Integer  n  is the parameter space dimension, i.e., the number
    !    of adjustable parameters.

    !  o Function  ff  is a user-supplied scalar function of n variables, which
    !    must have the calling sequence f = ff(n,xall,x), where x is a real parameter
    !    array of length n.  This function must be written so as to bound all
    !    parameters to the interval [0,1]; that is, the user must determine
    !    a priori bounds for the parameter space, and ff must use these bounds
    !    to perform the appropriate scalings to recover true parameter values in
    !    the a priori ranges.

    !    By convention, ff should return higher values for more optimal
    !    parameter values (i.e., individuals which are more "fit").
    !    For example, in fitting a function through data points, ff
    !    could return the inverse of chi**2.

    !    In most cases initialization code will have to be written
    !    (either in a driver or in a separate subroutine) which loads
    !    in data values and communicates with ff via one or more labeled
    !    common blocks.  An example exercise driver and fitness function
    !    are provided in the accompanying file, xpkaia.f.


    !  Input/Output:


    !  o Array  ctrl  is an array of control flags and parameters, to
    !    control the genetic behavior of the algorithm, and also printed
    !    output.  A default value will be used for any control variable
    !    which is supplied with a value less than zero.  On exit, ctrl
    !    contains the actual values used as control variables.  The
    !    elements of ctrl and their defaults are:

    !       ctrl( 1) - number of individuals in a population (default
    !                  is 100)
    !       ctrl( 2) - number of generations over which solution is
    !                  to evolve (default is 500)
    !       ctrl( 3) - number of significant digits (i.e., number of
    !                  genes) retained in chromosomal encoding (default
    !                  is 6)  (Note: This number is limited by the
    !                  machine floating point precision.  Most 32-bit
    !                  floating point representations have only 6 full
    !                  digits of precision.  To achieve greater preci-
    !                  sion this routine could be converted to double
    !                  precision, but note that this would also require
    !                  a double precision random number generator, which
    !                  likely would not have more than 9 digits of
    !                  precision if it used 4-byte integers internally.)
    !       ctrl( 4) - crossover probability; must be  <= 1.0 (default is 0.85)
    !       ctrl( 5) - mutation mode; 1/2=steady/variable (default is 2)
    !       ctrl( 6) - initial mutation rate; should be small (default is 0.005)
    !                  (Note: the mutation rate is the probability that any one
    !                  gene locus will mutate in any one generation.)
    !       ctrl( 7) - minimum mutation rate; must be >= 0.0 (default is 0.0005)
    !       ctrl( 8) - maximum mutation rate; must be <= 1.0 (default is 0.25)
    !       ctrl( 9) - relative fitness differential; range from 0
    !                  (none) to 1 (maximum).  (default is 1.)
    !       ctrl(10) - reproduction plan; 1/2/3=Full generational
    !                  replacement/Steady-state-replace-random/Steady-
    !                  state-replace-worst (default is 3)
    !       ctrl(11) - elitism flag; 0/1=off/on (default is 0)
    !                  (Applies only to reproduction plans 1 and 2)
    !       ctrl(12) - printed output 0/1/2=None/Minimal/Verbose (default is 0)


    ! Output:


    !  o Array  x(1:n)  is the "fittest" (optimal) solution found,
    !     i.e., the solution which maximizes fitness function ff

    !  o Scalar  f  is the value of the fitness function at x

    !  o Integer  status  is an indicator of the success or failure
    !     of the call to pikaia (0=success; non-zero=failure)


    ! Constants

    !INTEGER, PARAMETER :: nmax = 32, pmax = 128, dmax = 6
!     integer, parameter :: nmax = 50, pmax = 512, dmax = 12
    integer, parameter :: nmax = 110, pmax = 1024, dmax = 12

    !  o NMAX is the maximum number of adjustable parameters (n <= NMAX)

    !  o PMAX is the maximum population (ctrl(1) <= PMAX)

    !  o DMAX is the maximum number of Genes (digits) per Chromosome
    !        segement (parameter) (ctrl(3) <= DMAX)


    !     Local variables
    integer :: np, nd, ngen, imut, irep, ielite, ivrb, k, ip, ig, ip1,  &
        ip2, new, newtot
    real(dp) :: pcross, pmut, pmutmn, pmutmx, fdif

    real(dp) :: ph(nmax,2), oldph(nmax,pmax), newph(nmax,pmax)

    integer :: gn1(nmax*dmax), gn2(nmax*dmax)
    integer :: ifit(pmax), jfit(pmax)
    real(dp) :: fitns(pmax)

    !     User-supplied uniform random number generator
    ! real(dp) :: urand
    ! EXTERNAL urand

    real(dp),dimension(:,:),allocatable::bestpik
    integer,dimension(:),allocatable::idbest
    ! my counters for write into files
    integer::iwrt,wrtgen,igc,uall,ubest

    ! Function urand should not take any arguments.  If the user wishes to be able
    ! to initialize urand, so that the same sequence of random numbers can be
    ! repeated, this capability could be implemented with a separate subroutine,
    ! and called from the user's driver program.  An example urand function
    ! (and initialization subroutine) which uses the function ran0 (the "minimal
    ! standard" random number generator of Park and Miller [Comm. ACM 31, 1192-
    ! 1201, Oct 1988; Comm. ACM 36 No. 7, 105-110, July 1993]) is provided.

    ! ======================
    ! insert by Luca Borsato
    write(*,'(a)')" START PIKAIA SUBROUTINE "

    !     Set control variables from input and defaults
    call setctl(ctrl, n, np, ngen, nd, pcross, pmutmn, pmutmx, pmut, imut, fdif, &
        irep, ielite, ivrb, STATUS)
    ! ======================
    ! insert by Luca Borsato
    write(*,'(a)')" DONE setclt in PIKAIA"

    if (STATUS /= 0) then
      write (*, '(a)') ' Control vector (ctrl) argument(s) invalid'
      return
    end if

    !     Make sure locally-dimensioned arrays are big enough
    if (n > nmax .or. np > pmax .or. nd > dmax) then
      write (*, *) ' Number of parameters, population, or genes too large'
      write(*,'(3(a,i5))')" n = ",n," np = ",np," nd = ",nd
      STATUS = -1
      stop
      return
    end if

    write(*,'(a,i8)',advance='no')" Init random number with seed: ",seed_pik+iGlobal
!     call rninit(seed_pik+iGlobal)
    call init_random_seed_input(np*n,seed_pik+iGlobal)
    write(*,'(a)')" ... done. "
    flush(6)

    !     Compute initial (random but bounded) phenotypes
    write(*,'(a)',advance='no')" Computing initial random phenotypes..."
    call init_good_particles(n,np,xall,oldph)
    
    !$omp parallel do
    do  ip = 1, np
!       do  k = 1, n
! !           oldph(k,ip) = urand()
!         call random_number(oldph(k,ip))
!       end do

!       call random_number(oldph(:,ip))
      fitns(ip) = ff(n,xall,oldph(1:n,ip))
    end do
    !$omp end parallel do

    ! ======================
    ! insert by Luca Borsato
    write(*,'(a)')" done. "
    flush(6)

    !     Rank initial population by fitness order
    call rnkpop(np,fitns,ifit(1:np),jfit(1:np))
    write(*,'(a)')" Ranked initial population by fitness order"
    flush(6)
  
    iwrt=ivrb
    call bestpik_init(ngen,idbest,bestpik)
    wrtgen=10
    igc=1
    !call check_file_pik(uall)
    call check_file_pik(iGlobal,uall,ubest)

    !     Main Generation Loop
!     do  ig = 1, ngen
    ig=0
    gen_do: do
      ig=ig+1
      !        Main Population Loop
      newtot = 0
      do  ip = 1, np / 2
          !           1. pick two parents
          call select(np,jfit(1:np),fdif,ip1)
30        call select(np,jfit(1:np),fdif,ip2)
          if (ip1 == ip2) GO TO 30
          !           2. encode parent phenotypes
          call encode(n,nd,oldph(1:n,ip1),gn1)
          call encode(n,nd,oldph(1:n,ip2),gn2)
          !           3. breed
          call cross(n,nd,pcross,gn1,gn2)
          call mutate(n,nd,pmut,gn1)
          call mutate(n,nd,pmut,gn2)
          !           4. decode offspring genotypes
          call decode(n,nd,gn1,ph(1:n,1))
          call decode(n,nd,gn2,ph(1:n,2))
          !           5. insert into population
          if (irep == 1) then
            call genrep(nmax,n,np,ip,ph,newph)
          else
            call stdrep(ff,nmax,n,np,irep,ielite,xall,&
                  &ph,oldph,fitns,ifit,jfit,new)
            newtot = newtot + new
          end if
          !        End of Main Population Loop
      end do
      
      !        if running full generational replacement: swap populations
      if (irep == 1) call newpop(ff,ielite,nmax,n,np,xall,oldph,newph,ifit,&
            &jfit,fitns,newtot)
      !        adjust mutation rate?
      if (imut == 2)then
          call adjmut(np,fitns,ifit,pmutmn,pmutmx,pmut)
      end if
      !        print generation report to standard output?
      !if (ivrb > 0) call report(ivrb,nmax,n,np,nd,oldph,fitns,ifit,pmut,ig,newtot)
      !     End of Main Generation Loop
      
      ! save all data
      if(uall.gt.0) call write_allpik(uall,ig,ifit(1:np),oldph,fitns)
      ! stores the fittest for each generation
      !call bestpik_store(ig,np,ifit,oldph,fitns,idbest,bestpik)
      call bestpik_store(ubest,ig,np,ifit,oldph,fitns,idbest,bestpik)
!       write(*,'(a)')'All and Best written to file'
!       flush(6)
      
      if(getStat(ig,ngen).ge.wrtgen)then
          write(*,'(a,i4,a)')" Completed iterations: ",getStat(ig,ngen)," of 100%"
          write(*,'(a,i6,2(a,g25.14))')" DONE PIKAIA GENERATION ",ig,&
            & " fitness_x_dof = ",bestpik(n+2,ig),&
            & " fitness = ",bestpik(n+3,ig)
          write(*,'(a)')" Parameters "
          write(*,'(1000g20.8)')bestpik(1:n,ig)
          write(*,'(a)')""
          flush(6)
          igc = ig + 1
          wrtgen=wrtgen+10
      end if

!       if(iwrt.eq.ig)then
!           write(*,'(a,i6)')" DONE PIKAIA GENERATION ",ig
!           write(*,'(a,i6,2(a,g25.14))')" DONE PIKAIA GENERATION ",ig,&
!             & " fitness_x_dof = ",bestpik(n+2,ig),&
!             & " fitness = ",bestpik(n+3,ig)
!           write(*,'(a)')" Parameters "
!           write(*,'(1000g20.8)')bestpik(1:n,ig)
!           write(*,'(a)')""
!           flush(6)
!           iwrt=iwrt+ivrb
!       end if
!       if(ig.eq.1)stop

!     end do
      if(ig.ge.ngen)then
        if(bestpik(n+3,ig).ge.1000._dp)then
          write(*,'(a)')''
          write(*,'(a,i6,a)')' Too high fitness > 1000, continue for other ',ngen,' iterations'
          flush(6)
          ig=0
          wrtgen=10
          igc=1
          if(uall.gt.0) close(uall)
          close(ubest)
          call check_file_pik(iGlobal,uall,ubest) ! reset file
        else
          exit gen_do
        end if
      else if(bestpik(n+3,ig).le.one)then
        exit gen_do
      end if
    
    end do gen_do

    if(uall.gt.0) close(uall)
    !     Return best phenotype and its fitness
    do  k = 1, n
      x(k) = oldph(k,ifit(np))
    end do
    f = fitns(ifit(np))
    !call write_pik(ngen,storepik)
    close(ubest)

!     write(*,'(a,i6)')" DONE LAST PIKAIA GENERATION ",ngen
!     write(*,'(a)')" BEST " 
!     write(*,'(1(a,g20.8))')" fitness_x_dof = ",bestpik(n+2,ngen)
!     write(*,'(1(a,g20.8))')"       fitness = ",bestpik(n+3,ngen)
!     write(*,'(a)')" Parameters "
!     write(*,'(1000g20.8)')bestpik(1:n,ngen)
!     write(*,'(a)')""
    write(*,'(a,i6)')" DONE LAST PIKAIA GENERATION ",ig
    write(*,'(a)')" BEST " 
    write(*,'(a,i6,2(a,g25.14))')" DONE PIKAIA GENERATION ",ig,&
            & " fitness_x_dof = ",bestpik(n+2,ig),&
            & " fitness = ",bestpik(n+3,ig)
    write(*,'(a)')" Parameters "
    write(*,'(1000g20.8)')bestpik(1:n,ig)
    write(*,'(a)')""
    if(allocated(bestpik)) deallocate(idbest,bestpik)

    return
  end subroutine pikaia

  !********************************************************************

  subroutine setctl(ctrl,n,np,ngen,nd,pcross,pmutmn,pmutmx,pmut,  &
      imut,fdif,irep,ielite,ivrb,STATUS)
    !===================================================================
    !     Set control variables and flags from input and defaults
    !===================================================================

    !     Input
    !     Input/Output
    real(dp), intent(in out)  :: ctrl(12)
    integer, intent(in)   :: n

    !     Output
    integer, intent(out)  :: np
    integer, intent(out)  :: ngen
    integer, intent(out)  :: nd
    real(dp), intent(out)     :: pcross
    real(dp), intent(out)     :: pmutmn
    real(dp), intent(out)     :: pmutmx
    real(dp), intent(out)     :: pmut
    integer, intent(out)  :: imut
    real(dp), intent(out)     :: fdif
    integer, intent(out)  :: irep
    integer, intent(out)  :: ielite
    integer, intent(out)  :: ivrb
    integer, intent(out)  :: STATUS


    !     Local
    integer :: i
    real(dp), save  :: dfault(12) = (/ 100._dp, 500._dp, 5._dp, .85_dp, 2._dp, .005_dp, .0005_dp, .25_dp,  &
        1._dp, 1._dp, 1._dp, 0._dp /)

    do  i = 1, 12
      if (ctrl(i) < 0._dp) ctrl(i) = dfault(i)
    end do

    !  np = ctrl(1)
    !  ngen = ctrl(2)
    !  nd = ctrl(3)
    !  pcross = ctrl(4)
    !  imut = ctrl(5)
    !  pmut = ctrl(6)
    !  pmutmn = ctrl(7)
    !  pmutmx = ctrl(8)
    !  fdif = ctrl(9)
    !  irep = ctrl(10)
    !  ielite = ctrl(11)
    !  ivrb = ctrl(12)

    np = int(ctrl(1))
    ngen = int(ctrl(2))
    nd = int(ctrl(3))
    pcross = ctrl(4)
    imut = int(ctrl(5))
    pmut = ctrl(6)
    pmutmn = ctrl(7)
    pmutmx = ctrl(8)
    fdif = ctrl(9)
    irep = int(ctrl(10))
    ielite = int(ctrl(11))
    ivrb = int(ctrl(12))

    STATUS = 0

    !     Print a header
    if (ivrb > 0) then
      write (*,5000) ngen, np, n, nd, pcross, pmut, pmutmn, pmutmx, fdif
      if (imut == 1) write (*,5100) 'Constant'
      if (imut == 2) write (*,5100) 'Variable'
      if (irep == 1) write (*,5200) 'Full generational replacement'
      if (irep == 2) write (*,5200) 'Steady-state-replace-random'
      if (irep == 3) write (*,5200) 'Steady-state-replace-worst'
    end if

    !     Check some control values
    if (imut /= 1 .and. imut /= 2) then
      write (*,5300)
      STATUS = 5
    end if

    if (fdif > 1.) then
      write (*,5400)
      STATUS = 9
    end if

    if (irep /= 1 .and. irep /= 2 .and. irep /= 3) then
      write (*,5500)
      STATUS = 10
    end if

    if (pcross > 1.0 .or. pcross < 0.) then
      write (*,5600)
      STATUS = 4
    end if

    if (ielite /= 0 .and. ielite /= 1) then
      write (*,5700)
      STATUS = 11
    end if

    if (irep == 1 .and. imut == 1 .and. pmut > 0.5 .and. ielite == 0) then
      write (*,5800)
    end if

    if (irep == 1 .and. imut == 2 .and. pmutmx > 0.5 .and. ielite == 0) then
      write (*,5900)
    end if

    if (fdif < 0.33 .and. irep /= 3) then
      write (*,6000)
    end if

    if (mod(np,2) > 0) then
      np = np - 1
      write (*,6100) np
    end if

    return
5000 format (/' ', 60('*') /   &
        ' *', t16, 'PIKAIA Genetic Algorithm Report ', t60, '*' / &
        ' ', 60('*') //  &
        '   Number of Generations evolving: ', i4 /  &
        '       Individuals per generation: ', i4 /  &
        '    Number of Chromosome segments: ', i4 /  &
        '    Length of Chromosome segments: ', i4 /  &
        '            Crossover probability: ', f9.4 /   &
        '            Initial mutation rate: ', f9.4 /   &
        '            Minimum mutation rate: ', f9.4 /   &
        '            Maximum mutation rate: ', f9.4 /   &
        '    Relative fitness differential: ', f9.4)
5100 format ('                    Mutation Mode: '/ a)
5200 format ('                Reproduction Plan: '/ a)
5300 format (' ERROR: illegal value for imut (ctrl(5))')
5400 format (' ERROR: illegal value for fdif (ctrl(9))')
5500 format (' ERROR: illegal value for irep (ctrl(10))')
5600 format (' ERROR: illegal value for pcross (ctrl(4))')
5700 format (' ERROR: illegal value for ielite (ctrl(11))')
5800 format (' WARNinG: dangerously high value for pmut (ctrl(6));' /  &
        ' (Should enforce elitism with ctrl(11)=1.)')
5900 format (' WARNinG: dangerously high value for pmutmx (ctrl(8));' /  &
        ' (Should enforce elitism with ctrl(11)=1.)')
6000 format (' WARNinG: dangerously low value of fdif (ctrl(9))')
6100 format (' WARNinG: decreasing population size (ctrl(1)) to np='/ i4 )
  end subroutine setctl

  !********************************************************************

  subroutine report(ivrb, ndim, n, np, nd, oldph, fitns, ifit, pmut, ig, nnew)

    !     Write generation report to standard output

    !     Input:
    integer,  intent(in)  :: ivrb
    integer,  intent(in)  :: ndim
    integer,  intent(in)  :: n
    integer,  intent(in)  :: np
    integer,  intent(in)  :: nd
    real(dp), intent(in)      :: oldph(ndim, np)
    real(dp), intent(in)      :: fitns(np)
    integer, intent(in)   :: ifit(np)
    real(dp), intent(in)      :: pmut
    integer, intent(in)   :: ig
    integer, intent(in)   :: nnew

    !     Output: none

    !     Local
    real(dp), save  :: bestft = 0.0_dp, pmutpv = 0.0_dp
    integer  :: ndpwr, k
    logical  :: rpt

    rpt = .false.

    if (pmut /= pmutpv) then
      pmutpv = pmut
      rpt = .true.
    end if

    if (fitns(ifit(np)) /= bestft) then
      bestft = fitns(ifit(np))
      rpt = .true.
    end if

    if (rpt .or. ivrb >= 2) then

      !        Power of 10 to make integer genotypes for display
      ndpwr = nint(10.**nd)

      write (*, '(/i6, i6, f10.6, 4f12.6)') ig, nnew, pmut,  &
            fitns(ifit(np)), fitns(ifit(np-1)), fitns(ifit(np/2))
      do  k = 1, n
          write (*, '(22x, 3i16)') nint(ndpwr*oldph(k, ifit(np))),  &
              nint(ndpwr*oldph(k, ifit(np-1))), nint(ndpwr*oldph(k, ifit(np/2)))
      end do

    end if
    return
  end subroutine report

  !**********************************************************************
  !                         GENETICS MODULE
  !**********************************************************************

  !     ENCODE:    encodes phenotype into genotype
  !                called by: PIKAIA

  !     DECODE:    decodes genotype into phenotype
  !                called by: PIKAIA

  !     CROSS:     Breeds two offspring from two parents
  !                called by: PIKAIA

  !     MUTATE:    Introduces random mutation in a genotype
  !                called by: PIKAIA

  !     ADJMUT:    Implements variable mutation rate
  !                called by: PIKAIA

  !**********************************************************************

  subroutine encode(n, nd, ph, gn)
    !======================================================================
    !     encode phenotype parameters into integer genotype
    !     ph(k) are x, y coordinates [ 0 < x, y < 1 ]
    !======================================================================


    integer, intent(in)   :: n
    integer, intent(in)   :: nd
    real(dp), intent(in out)  :: ph(n)
    integer, intent(out)  :: gn(n*nd)

    !     Inputs:



    !     Output:


    !     Local:
    integer :: ip, i, j, ii
    real(dp) :: z

    z = 10._dp ** nd
    ii = 0
    do  i = 1, n
      ip = int(ph(i)*z)
      do  j = nd, 1, -1
          gn(ii+j) = mod(ip, 10)
          ip = ip / 10
      end do
      ii = ii + nd
    end do

    return
  end subroutine encode

  !**********************************************************************

  subroutine decode(n, nd, gn, ph)
    !======================================================================
    !     decode genotype into phenotype parameters
    !     ph(k) are x, y coordinates [ 0 < x, y < 1 ]
    !======================================================================


    integer, intent(in)  :: n
    integer, intent(in)  :: nd
    integer, intent(in)  :: gn(n*nd)
    real(dp), intent(out)    :: ph(n)

    !     Inputs:


    !     Output:


    !     Local:
    integer :: ip, i, j, ii
    real(dp) :: z

    z = 10._dp ** (-nd)
    ii = 0
    do  i = 1, n
      ip = 0
      do  j = 1, nd
          ip = 10 * ip + gn(ii+j)
      end do
      ph(i) = ip * z
      ii = ii + nd
    end do

    return
  end subroutine decode

  !**********************************************************************

  subroutine cross(n, nd, pcross, gn1, gn2)
    !======================================================================
    !     breeds two parent chromosomes into two offspring chromosomes
    !     breeding occurs through crossover starting at position ispl
    !======================================================================
    !     USES: urand

    !     Inputs:
    integer, intent(in)      :: n
    integer, intent(in)      :: nd
    real(dp), intent(in)         :: pcross

    !     Input/Output:
    integer, intent(in out)  :: gn1(n*nd)
    integer, intent(in out)  :: gn2(n*nd)

    !     Local:
    integer :: i, ispl, t

    !     Function
    ! real(dp) :: urand
    ! EXTERNAL urand


    !     Use crossover probability to decide whether a crossover occurs
    if (urand() < pcross) then

      !        Compute crossover point
      ispl = int(urand()*n*nd) + 1

      !        Swap genes at ispl and above
      do  i = ispl, n * nd
          t = gn2(i)
          gn2(i) = gn1(i)
          gn1(i) = t
      end do
    end if

    return
  end subroutine cross


  !**********************************************************************

  subroutine mutate(n, nd, pmut, gn)
    !======================================================================
    !     Mutations occur at rate pmut at all gene loci
    !======================================================================
    !     USES: urand

    !     Input:
    integer, intent(in)      :: n
    integer, intent(in)      :: nd
    real(dp), intent(in)         :: pmut

    !     Input/Output:
    integer, intent(in out)  :: gn(n*nd)


    !     Local:
    integer :: i

    !     Function:
    ! real(dp) :: urand
    ! EXTERNAL urand

    !     Subject each locus to mutation at the rate pmut
    do  i = 1, n * nd
      if (urand() < pmut) then
          gn(i) = int(urand()*10._dp)
      end if
    end do

    return
  end subroutine mutate

  !**********************************************************************

  subroutine adjmut(np, fitns, ifit, pmutmn, pmutmx, pmut)
    !======================================================================
    !     dynamical adjustment of mutation rate; criterion is relative
    !     difference in absolute fitnesses of best and median individuals
    !======================================================================

    !     Input:
    integer, intent(in)   :: np
    real(dp), intent(in)      :: fitns(:)
    integer, intent(in)   :: ifit(:)
    real(dp), intent(in)      :: pmutmn
    real(dp), intent(in)      :: pmutmx

    !     Input/Output:
    real(dp), intent(in out)  :: pmut

    !     Local:
    real(dp)  :: rdif
    real(dp), parameter  :: rdiflo = 0.05_dp, rdifhi = 0.25_dp, delta = 1.5_dp

    rdif = abs(fitns(ifit(np)) - fitns(ifit(np/2))) / (fitns(ifit(np)) +  &
        fitns(ifit(np/2)))
    if (rdif <= rdiflo) then
      pmut = min(pmutmx, pmut*delta)
    else if (rdif >= rdifhi) then
      pmut = max(pmutmn, pmut/delta)
    end if

    return
  end subroutine adjmut

  !**********************************************************************
  !                       REdpODUCTION MODULE
  !**********************************************************************

  !  SELECT:   Parent selection by roulette wheel algorithm
  !            called by: PIKAIA

  !  RNKPOP:   Ranks initial population
  !            called by: PIKAIA, NEWPOP

  !  GENREP:   Inserts offspring into population, for full
  !            generational replacement
  !            called by: PIKAIA

  !  STDREP:   Inserts offspring into population, for steady-state
  !            reproduction
  !            called by: PIKAIA
  !            calls:     FF

  !  NEWPOP:   Replaces old generation with new generation
  !            called by: PIKAIA
  !            calls:     FF, RNKPOP

  !**********************************************************************

  subroutine select(np, jfit, fdif, idad)
    !======================================================================
    !     Selects a parent from the population, using roulette wheel
    !     algorithm with the relative fitnesses of the phenotypes as
    !     the "hit" probabilities [see Davis 1991, chap. 1].
    !======================================================================
    !     USES: urand

    !     Input:
    integer, intent(in)   :: np
    integer, intent(in)   :: jfit(np)
    real(dp), intent(in)      :: fdif

    !     Output:
    integer, intent(out)  :: idad

    !     Local:
    integer :: np1, i
    real(dp) :: dice, rtfit

    !     Function:
    ! real(dp) :: urand
    ! EXTERNAL urand


    np1 = np + 1
    dice = urand() * np * np1
    rtfit = 0._dp
    do  i = 1, np
      rtfit = rtfit + np1 + fdif * (np1-2*jfit(i))
      if (rtfit >= dice) then
          idad = i
          GO TO 20
      end if
    end do
    !     Assert: loop will never exit by falling through

20  return
  end subroutine select

  !**********************************************************************

  subroutine rnkpop(n, arrin, indx, rank)
    !======================================================================
    !     Calls external sort routine to produce key index and rank order
    !     of input array arrin (which is not altered).
    !======================================================================
    !     USES: rqsort

    !     Input
    integer, intent(in)   :: n
    real(dp), intent(in)      :: arrin(:)

    !     Output
    integer, intent(out)  :: indx(:)
    integer, intent(out)  :: rank(:)


    !     Local
    integer :: i

    !     External sort subroutine
    ! EXTERNAL rqsort


    !     Compute the key index
    call rqsort(n, arrin, indx)

    !     ...and the rank order
    do  i = 1, n
      rank(indx(i)) = n - i + 1
    end do
    return
  end subroutine rnkpop

  !***********************************************************************

  subroutine genrep(ndim, n, np, ip, ph, newph)
    !=======================================================================
    !     full generational replacement: accumulate offspring into new
    !     population array
    !=======================================================================

    !     Input:
    integer, intent(in)  :: ndim
    integer, intent(in)  :: n
    integer, intent(in)  :: np
    integer, intent(in)  :: ip
    real(dp), intent(in)     :: ph(ndim, 2)

    !     Output:
    real(dp), intent(out)    :: newph(ndim, np)


    !     Local:
    integer :: i1, i2, k


    !     Insert one offspring pair into new population
    i1 = 2 * ip - 1
    i2 = i1 + 1
    do  k = 1, n
      newph(k, i1) = ph(k, 1)
      newph(k, i2) = ph(k, 2)
    end do

    return
  end subroutine genrep

  !**********************************************************************

  subroutine stdrep(ff, ndim, n, np, irep, ielite, xall, ph, oldph, fitns, ifit, jfit, nnew)
    !======================================================================
    !     steady-state reproduction: insert offspring pair into population
    !     only if they are fit enough (replace-random if irep=2 or
    !     replace-worst if irep=3).
    !======================================================================
    !     USES: ff, urand

    !     Input:
    integer, intent(in)      :: ndim
    integer, intent(in)      :: n
    integer, intent(in)      :: np
    integer, intent(in)      :: irep
    integer, intent(in)      :: ielite
    real(dp),dimension(:),intent(in)::xall
    real(dp), intent(in)         :: ph(ndim, 2)

    !     Input/Output:
    real(dp), intent(in out)     :: oldph(ndim, np)
    real(dp), intent(in out)     :: fitns(np)
    integer, intent(in out)  :: ifit(np)
    integer, intent(in out)  :: jfit(np)

    !     Output:
    integer, intent(out)     :: nnew

!     interface
!       function ff(n,xall,x) result(fn_val)
!         use constants,only:dp
!         implicit none
!         integer,intent(in)::n
!         real(dp),dimension(:),intent(in)::xall,x
!         real(dp)::fn_val
!       end function ff
!     end interface
    interface
      function ff(n,xall,x) result(fn_val)
        use constants,only:dp
        implicit none
        real(dp)::fn_val
        integer,intent(in)::n
        real(dp),intent(in),dimension(:)::xall
        real(dp),intent(in),dimension(n)::x
      end function ff
    end interface

    
    ! EXTERNAL ff

    !     Local:
    integer :: i, j, k, i1, if1
    real(dp) :: fit

    !     External function
    ! real(dp) :: urand
    ! EXTERNAL urand


    nnew = 0
    loop70:  do  j = 1, 2

      !        1. compute offspring fitness (with caller's fitness function)
      fit = ff(n,xall,ph(:, j))

      !        2. if fit enough, insert in population
      do  i = np, 1, -1
          if (fit > fitns(ifit(i))) then

            !              make sure the phenotype is not already in the population
            if (i < np) then
                do  k = 1, n
                  if (oldph(k, ifit(i+1)) /= ph(k, j)) GO TO 20
                end do
                cycle loop70
            end if

            !              offspring is fit enough for insertion, and is unique

            !              (i) insert phenotype at appropriate place in population
20           if (irep == 3) then
                i1 = 1
            else if (ielite == 0 .or. i == np) then
                i1 = int(urand()*np) + 1
            else
                i1 = int(urand()*(np-1)) + 1
            end if
            if1 = ifit(i1)
            fitns(if1) = fit
            do  k = 1, n
                oldph(k, if1) = ph(k, j)
            end do

            !              (ii) shift and update ranking arrays
            if (i < i1) then

                !                 shift up
                jfit(if1) = np - i
                do  k = i1 - 1, i + 1, -1
                  jfit(ifit(k)) = jfit(ifit(k)) - 1
                  ifit(k+1) = ifit(k)
                end do
                ifit(i+1) = if1
            else

                !                 shift down
                jfit(if1) = np - i + 1
                do  k = i1 + 1, i
                  jfit(ifit(k)) = jfit(ifit(k)) + 1
                  ifit(k-1) = ifit(k)
                end do
                ifit(i) = if1
            end if
            nnew = nnew + 1
            cycle loop70
          end if
      end do

    end do loop70

    return
  end subroutine stdrep

  !**********************************************************************

  subroutine newpop(ff, ielite, ndim, n, np, xall, oldph, newph, ifit, jfit, fitns, nnew)
    !======================================================================
    !     replaces old population by new; recomputes fitnesses & ranks
    !======================================================================
    !     USES: ff, rnkpop

    !     Input:
    integer, intent(in)   :: ielite
    integer, intent(in)   :: ndim
    integer, intent(in)   :: n
    integer, intent(in)   :: np
    real(dp),dimension(:),intent(in)::xall

    !     Input/Output:
    real(dp), intent(in out)  :: oldph(ndim, np)
    real(dp), intent(in out)  :: newph(ndim, np)

    !     Output:
    integer, intent(out)  :: ifit(np)
    integer, intent(out)  :: jfit(np)
    real(dp), intent(out)     :: fitns(np)
    integer, intent(out)  :: nnew

!     interface
!       function ff(n,xall,x) result(fn_val)
!         use constants,only:dp
!         implicit none
!         integer,intent(in)::n
!         real(dp),dimension(:),intent(in)::xall,x
!         real(dp)::fn_val
!       end function ff
!     end interface
    interface
      function ff(n,xall,x) result(fn_val)
        use constants,only:dp
        implicit none
        real(dp)::fn_val
        integer,intent(in)::n
        real(dp),intent(in),dimension(:)::xall
        real(dp),intent(in),dimension(n)::x
      end function ff
    end interface

    
    ! EXTERNAL ff

    !     Local:
    integer :: i, k

    nnew = np

    !     if using elitism, introduce in new population fittest of old
    !     population (if greater than fitness of the individual it is
    !     to replace)
    if (ielite == 1 .and. ff(n,xall,newph(1:n,1)) < fitns(ifit(np))) then
      do  k = 1, n
          newph(k, 1) = oldph(k, ifit(np))
      end do
      nnew = nnew - 1
    end if

    !     replace population
    !$omp parallel do
    do  i = 1, np
      do  k = 1, n
          oldph(k, i) = newph(k, i)
      end do

      !        get fitness using caller's fitness function
      fitns(i) = ff(n,xall,oldph(1:n,i))
    end do
    !$omp end parallel do


    !     compute new population fitness rank order
    call rnkpop(np, fitns, ifit(1:np), jfit(1:np))

    return
  end subroutine newpop

  !*********************************************************************

  function urand() result(fn_val)
    !=====================================================================
    !  Return the next pseudo-random deviate from a sequence which is
    !  uniformly distributed in the interval [0, 1]

    !  Uses the function ran0, the "minimal standard" random number
    !  generator of Park and Miller (Comm. ACM 31, 1192-1201, Oct 1988;
    !  Comm. ACM 36 No. 7, 105-110, July 1993).
    !=====================================================================

    !     Input - none

    !     Output
    real(dp)  :: fn_val

    !     Local
    ! inTEGER :: iseed
    ! real(dp) :: ran0
    ! EXTERNAL ran0

    !     Common block to make iseed visible to rninit (and to save
    !     it between calls)
    ! COMMON /rnseed/ iseed

!     fn_val = ran0()
      call random_number(fn_val)
    return
  end function urand

  !*********************************************************************

  subroutine rninit(seed)
    !=====================================================================
    !     Initialize random number generator urand with given seed
    !=====================================================================

    !     Input
    integer, intent(in)  :: seed

    !     Output - none

    !     Local
    ! inTEGER  :: iseed

    !     Common block to communicate with urand
    ! COMMON /rnseed/ iseed

    !     Set the seed value
    iseed = seed
    if (iseed <= 0) iseed = 123456
    return
  end subroutine rninit

  !*********************************************************************

  function ran0() result(fn_val)
    !=====================================================================
    !  "Minimal standard" pseudo-random number generator of Park and Miller.
    !  Returns a uniform random deviate r s.t. 0 < r < 1.0.
    !  Set seed to any non-zero integer value to initialize a sequence, then do
    !  not change seed between calls for successive deviates in the sequence.

    !  References:
    !     Park, S. and Miller, K., "Random Number Generators: Good Ones
    !        are Hard to Find", Comm. ACM 31, 1192-1201 (Oct. 1988)
    !     Park, S. and Miller, K., in "Remarks on Choosing and Implementing
    !        Random Number Generators", Comm. ACM 36 No. 7, 105-110 (July 1993)
    !=====================================================================
    ! *** Declaration section ***

    !     Output:
    real(dp)  :: fn_val

    !     Constants:

    integer, parameter  :: a = 48271, m = 2147483647, q = 44488, r = 3399

    real(dp), parameter :: scale = 1._dp/m, eps = 1.2E-7_dp, rnmx = 1._dp - eps

    !     Local:
    integer  :: j

    ! *** Executable section ***

    j = iseed / q
    iseed = a * (iseed - j*q) - r * j
    if (iseed < 0) iseed = iseed + m
    fn_val = min(iseed*scale, rnmx)

    return
  end function ran0

  !**********************************************************************

  subroutine rqsort(n, a, p)
    !======================================================================
    !  Return integer array p which indexes array a in increasing order.
    !  Array a is not disturbed.  The Quicksort algorithm is used.

    !  B. G. Knapp, 86/12/23

    !  Reference: N. Wirth, Algorithms and Data Structures/
    !  Prentice-Hall, 1986
    !======================================================================

    integer, intent(in)   :: n
    real(dp), intent(in)      :: a(:)
    integer, intent(out)  :: p(:)

    !     Constants

    integer, parameter :: lgn = 32, q = 11
    !        (LGN = log base 2 of maximum n;
    !         Q = smallest subfile to use quicksort on)

    !     Local:
    real(dp) :: x
    integer :: stackl(lgn), stackr(lgn), s, t, l, m, r, i, j

    !     Initialize the stack
    stackl(1) = 1
    stackr(1) = n
    s = 1

    !     Initialize the pointer array
    do  i = 1, n
      p(i) = i
    end do

20  if (s > 0) then
      l = stackl(s)
      r = stackr(s)
      s = s - 1

30     if ((r-l) < q) then

          !           Use straight insertion
          do  i = l + 1, r
            t = p(i)
            x = a(t)
            do  j = i - 1, l, -1
                if (a(p(j)) <= x) GO TO 50
                p(j+1) = p(j)
            end do
            j = l - 1
50           p(j+1) = t
          end do
      else

          !           Use quicksort, with pivot as median of a(l), a(m), a(r)
          m = (l+r) / 2
          t = p(m)
          if (a(t) < a(p(l))) then
            p(m) = p(l)
            p(l) = t
            t = p(m)
          end if
          if (a(t) > a(p(r))) then
            p(m) = p(r)
            p(r) = t
            t = p(m)
            if (a(t) < a(p(l))) then
                p(m) = p(l)
                p(l) = t
                t = p(m)
            end if
          end if

          !           Partition
          x = a(t)
          i = l + 1
          j = r - 1
70        if (i <= j) then
80           if (a(p(i)) < x) then
                i = i + 1
                GO TO 80
            end if
90           if (x < a(p(j))) then
                j = j - 1
                GO TO 90
            end if
            if (i <= j) then
                t = p(i)
                p(i) = p(j)
                p(j) = t
                i = i + 1
                j = j - 1
            end if
            GO TO 70
          end if

          !           Stack the larger subfile
          s = s + 1
          if (j-l > r-i) then
            stackl(s) = l
            stackr(s) = j
            l = i
          else
            stackl(s) = i
            stackr(s) = r
            r = j
          end if
          GO TO 30
      end if
      GO TO 20
    end if
    return
  end subroutine rqsort

  ! ---
  ! PIKAIA STORE AND WRIITING FILES
  ! ---
  
    ! check files for all and best simulations, creates new ones
  subroutine check_file_pik(iGlobal,uall,ubest)
    integer,intent(in)::iGlobal
    integer,intent(out)::uall,ubest
    character(512)::flpik,fmtpik,remove
    logical::exstat
    integer::ip

    ! check best simulation file
    flpik=trim(path)//trim(adjustl(string(iGlobal)))//"_bestpik_"//trim(adjustl(string(nfit)))//"par.out"
    flpik=trim(adjustl(flpik))
    inquire(file = trim(flpik), exist = exstat)
    if(exstat)then
      remove = "rm "//trim(flpik)
      remove = trim(adjustl(remove))
      call system(remove)
    end if
    ubest=get_unit(1)
    open(ubest,file=trim(flpik))
    fmtpik="# gen pop "
    do ip=1,nfit
      fmtpik=trim(adjustl(fmtpik))//&
            &" "//parid(ip)
    end do
!     fmtpik=trim(adjustl(fmtpik))//" invChi2r Chi2 Chi2r"
    fmtpik=trim(adjustl(fmtpik))//" inv_fitness fitness_x_dof fitness"
    write(ubest,'(a)')trim(adjustl(fmtpik))
    !close(ubest)

    ! create file to write all iterations all particles
    flpik=trim(path)//trim(adjustl(string(iGlobal)))//"_allpik_"//trim(adjustl(string(nfit)))//"par.out"
    inquire(file=trim(flpik), exist=exstat)
    if(exstat)then
      remove = "rm "//trim(flpik)
      remove = trim(adjustl(remove))
      call system(remove)
    end if
    if(wrtAll.eq.0)then
      uall=-1
    else
      uall=get_unit(1)
      open(uall,file=trim(adjustl(flpik)))
      write(uall,'(a)')trim(adjustl(fmtpik))
    end if

    return
  end subroutine check_file_pik

  ! saves the best simulations
  subroutine bestpik_store(ubest,igen,ipop,ifit,oldph,fitns,idbest,bestpik)
    integer,intent(in)::ubest,igen,ipop
    integer,dimension(:),intent(in)::ifit
    real(dp),dimension(:,:),intent(in)::oldph
    real(dp),dimension(:),intent(in)::fitns
    integer,dimension(:),intent(out)::idbest
    real(dp),dimension(:,:),intent(out)::bestpik
    integer::inp
    real(dp),dimension(:),allocatable::xpik,parpik
!     real(dp)::iChi2r,Chi2,Chi2r
    real(dp)::inv_fitness,fitness_x_dof,fitness
    character(512)::fmtw

    inp=ifit(ipop)
    allocate(xpik(nfit),parpik(nfit))
    xpik=oldph(:nfit,inp)
    parpik=zero
    call norm2par(xpik,parpik) ! convert parameters from [0, 1] to physical values
    bestpik(1:nfit,igen)=parpik

!         ! proper chi2 and chi2r to write
!     if(fitns(inp).ge.huge(zero))then
!       iChi2r=huge(zero)
!       Chi2=real(dof,dp)/iChi2r
!       Chi2r=Chi2
!     else if(fitns(inp).le.TOLERANCE)then
!       iChi2r=fitns(inp)
!       Chi2=huge(zero)
!       Chi2r=Chi2*inv_dof
!     else 
!       iChi2r=fitns(inp)
!       Chi2=real(dof,dp)/iChi2r
!       Chi2r=Chi2*inv_dof
!     end if
!     bestpik(nfit+1,igen)=iChi2r
!     bestpik(nfit+2,igen)=Chi2
!     bestpik(nfit+3,igen)=Chi2r
!     idbest(igen)=inp
    
    
    inv_fitness=fitns(inp)
    fitness=one/inv_fitness
    fitness_x_dof=real(dof,dp)*fitness
    bestpik(nfit+1,igen)=inv_fitness
    bestpik(nfit+2,igen)=fitness_x_dof
    bestpik(nfit+3,igen)=fitness
    idbest(igen)=inp
    
    ! added 2014-07-14: write best each iteration
    fmtw=adjustl("(i6,1x,i6,1000("//trim(sprec)//"))")
    write(ubest,trim(fmtw))igen,idbest(igen),bestpik(:,igen)
    flush(ubest)
    deallocate(xpik,parpik)

    return
  end subroutine bestpik_store

  ! write all the simulations into file
  subroutine write_allpik(uall,ig,ip,pp,fitns)
    integer,intent(in)::uall,ig
    real(dp),dimension(:,:),intent(in)::pp ! particles
    real(dp),dimension(:),intent(in)::fitns ! particles fitness
    integer,dimension(:),intent(in)::ip ! index of particles
    integer::np
    integer::j
    character(128)::wfmt
    real(dp)::inv_fitness,fitness_x_dof,fitness
!     real(dp)::iChi2r,Chi2,Chi2r,tChi2
    real(dp),dimension(:),allocatable::xpar

    np=size(ip)
    wfmt="(i6,1x,i4,1x,1000("//trim(sprec)//",1x))"

    allocate(xpar(nfit))

!     do j=1,np
!       iChi2r=fitns(j)
!       if(iChi2r.ge.huge(zero))then
!           tChi2=huge(zero)
!           Chi2=tChi2
!           Chi2r=Chi2
!       else if(iChi2r.le.TOLERANCE)then
!           tChi2=iChi2r
!           Chi2=huge(zero)
!           Chi2r=Chi2*inv_dof
!       else
!           tChi2=iChi2r
!           Chi2=real(dof,dp)/iChi2r
!           Chi2r=Chi2*inv_dof
!       end if
    do j=1,np
      inv_fitness=fitns(j)
      fitness=one/inv_fitness
      fitness_x_dof=real(dof,dp)*fitness
      xpar=zero
      call norm2par(pp(1:nfit,j),xpar)
!       write(uall,trim(wfmt))ig,ip(j),xpar,iChi2r,Chi2,Chi2r
      write(uall,trim(wfmt))ig,ip(j),xpar,inv_fitness,fitness_x_dof,fitness
      flush(uall)
    end do

    deallocate(xpar)

    return
  end subroutine write_allpik

end module Genetic_Algorithm
