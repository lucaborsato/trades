!
! -----------------------------------
! Module of POS optimization routines
! -----------------------------------
!    module util_random is requied
!    module util_sort   is requied
!    module util_qmc    is requied
!
!----------------------------------------------------------------------
module opti_pso
    use constants, only: dp, sprec, zero, half, one, two, three, TOLERANCE, TOL_dp
!   use parameters,only:path,wrtAll,seed_pso,ndata,nfit,dof,inv_dof,npar,parid,np_pso,nit_pso,wrt_pso,minpar,maxpar,population,population_fitness,pso_best_evolution,ncpu_in
!   use parameters,only:path,seed_pso,np_pso,nit_pso,wrt_pso,wrtAll,nGlobal,dof,minpar,maxpar,population,population_fitness,pso_best_evolution,inertia_in,self_in,swarm_in,randsearch_in,vmax_in,vrand_in
    use custom_type
    use parameters
    use parameters_conversion
    use init_trades, only: get_unit
    use random_trades
    use convert_type, only: string
    use fitness_module
!$  use omp_lib
    implicit none
    private

    !============= PUBLIC ROUTINES =============
    public :: pso_driver   ! driver by Luca Borsato
    public :: do_pso       ! simple routine with standard control parameters
    public :: do_psof      ! simple routine with configuration file
    public :: pso          ! primitive routine
    public :: pso_cbdummy  ! dummy callback function

    public :: evaluate_pso ! function that has to be call by other modules

    !============= PUBLIC PARAMETERS =============
    !----- control flags -----
    integer, public, parameter :: pso_maximize = 1  ! maximize object function
    integer, public, parameter :: pso_minimize = 0  ! minimize object function
    integer, public, parameter :: f_quiet = 1     ! restrain messages
    integer, public, parameter :: f_initp = 2     ! use initial parameters
    integer, public, parameter :: f_lrand = 4     ! local random search algorithm
    integer, public, parameter :: f_vlimit = 8     ! velocity limit algorithm
    integer, public, parameter :: f_vrand = 16    ! velocity perturbation algorithm
    integer, public, parameter :: f_gcenter = 32    ! evaluate center of gravity
    integer, public, parameter :: f_gjump = 64    ! move to center of gravity
    integer, public, parameter :: f_qmcinit = 128   ! initialize by Quasi-Monte Carlo

    !============= PRIVATE VARIABLES =============
    integer, parameter     :: max_nparam = 110      ! maximum number of parameter
    real(dp)               :: scale_a(max_nparam)   ! parameter scaling factor A
    real(dp)               :: scale_b(max_nparam)   ! parameter scaling factor B
    real(dp), parameter    :: constrained = -1.E37_dp   ! evaluation value of constrained particle
    real(dp), parameter    :: void = -1.E38_dp         ! evaluation value of undefined particle

    type type_p
        real(dp) :: x(max_nparam)          ! current position
        real(dp) :: v(max_nparam)          ! current velocity
        real(dp) :: ev                     ! current evaluation value
        real(dp) :: evbest                 ! best evaluation value of the particle
        real(dp) :: xbest(max_nparam)      ! best position
        integer :: constrained         ! number of constrained condition
    end type type_p

contains

    ! another function called by the PSO module, and now it uses the fitness function
    ! from the fitness_module module
    function evaluate_pso(allpar, par, pso_fitness) result(ir)
        integer::ir
        real(dp), dimension(:), intent(in)::allpar, par
        real(dp), intent(out)::pso_fitness
        ! real(dp)::fitness
        real(dp)::chi_square, reduced_chi_square, lnLikelihood, ln_const, bic

        call base_fitness_function(allpar, par,&
            &chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)
        ! pso_fitness = one/reduced_chi_square
        pso_fitness = lnLikelihood
        ir = 0

        return
    end function evaluate_pso

    subroutine pso_driver(iGlobal, evfunc, n_fit, allpar, minpar_in, maxpar_in, xpar, pso_fitness)
        integer, intent(in)::iGlobal, n_fit
        real(dp), dimension(:), intent(in)::allpar, minpar_in, maxpar_in
        real(dp), intent(out)::pso_fitness
!     real(dp),dimension(:)::xpar
        real(dp), dimension(:), intent(out)::xpar
        real(dp), dimension(:), allocatable::ipar

!     integer::ii

        interface
            !-------- evaluation function -------
            integer function evfunc(allpar_in, xpar_in, fitness)  ! return number of constrained condition
                use constants, only: dp
                real(dp), dimension(:), intent(in)::allpar_in, xpar_in    ! parameter set
                real(dp), intent(out) :: fitness       ! evaluation value
            end function evfunc
        end interface

        if (.not. allocated(ipar)) allocate (ipar(n_fit))
        if (.not. allocated(population)) allocate (population(n_fit, np_pso, nit_pso))
        if (.not. allocated(population_fitness)) allocate (population_fitness(np_pso, nit_pso))
        xpar = zero
        ipar = minpar_in
        population = zero
        population_fitness = zero

        call do_pso(n_fit, minpar_in, maxpar_in, ipar, np_pso, nit_pso, wrt_pso, 1,&
            &evfunc, allpar, xpar, pso_fitness, iGlobal)

        deallocate (ipar)

        return
    end subroutine pso_driver

    ! initialize bestpso
    subroutine bestpso_init(nit, idbest, bestpso)
!     use parameters,only:nfit
        use convert_type, only: string
        integer, intent(in)::nit
        real(dp), dimension(:, :), allocatable, intent(out)::bestpso
        integer, dimension(:), allocatable, intent(out)::idbest
        integer::nstore, n_fitness

        n_fitness = 2
        nstore = nfit + n_fitness
        if (.not. allocated(bestpso)) allocate (bestpso(nstore, nit), idbest(nit))
        bestpso = zero

        return
    end subroutine bestpso_init

    ! to determine the % of simulations
    function getStat(i, n) result(out)
        integer::out
        integer, intent(in)::i, n
        real(dp)::fi, fn

        fi = real(i, dp)
        fn = real(n, dp)
        out = int((fi/fn)*100._dp)

        return
    end function getStat

    ! --------------------------------------------------------------------------

    !======================================================================
    ! PSO Simplified Wrapper Routine (with standard control parameters)
    !======================================================================
    subroutine do_pso(n, xmin, xmax, xinit, np, it, ir, maximize, evfunc, xall, x, r, iGlobal)
        integer, intent(in):: n             ! number of parameter
        real(dp), intent(in):: xmin(1:n)     ! lower bounds
        real(dp), intent(in):: xmax(1:n)     ! upper bounds
        real(dp), intent(in):: xinit(1:n)    ! initial parameters
        integer, intent(in):: np            ! population size (number of particle)
        integer, intent(in):: it            ! maximum iteration count (<=0:unlimitted)
        integer, intent(in):: ir            ! interval of message printing (<=0:no mesage)
        integer, intent(in):: maximize      ! maximize-minimize flag (>0:max, <=0:min)
        real(dp), dimension(:), intent(in)::xall
        real(dp), intent(out)::x(1:n)        ! optimum parameter set
        real(dp), intent(out)::r             ! optimum evaluation value
        integer, intent(in)::iGlobal
        interface
            !-------- evaluation function -------
            integer function evfunc(xall_in, x_in, fitness)  ! return number of constrained condition
                use constants, only: dp
                real(dp), dimension(:), intent(in)::xall_in, x_in    ! parameter set
                real(dp), intent(out) :: fitness       ! evaluation value
            end function evfunc
        end interface

        !====== Local variables =======
        real(dp):: w, c1, c2, c3, vmax, vrand
        integer :: ie, flag

        !====== Setup control parameters =======
!     w = 0.9_dp            ! inertia parameter (0.9 is recommended) ! original
! !     w = 1.2_dp            ! inertia parameter (0.9 is recommended)
!     c1 = two           ! self intention parameter (2.0 is recommended)
!     c2 = two           ! swarm intention parameter (2.0 is recommended)
!     c3 = 1.e-5_dp         ! random search parameter (very small value is recommended) ! original value
!     !c3 = 1.e-4_dp         ! random search parameter (very small value is recommended)
!     vmax = half         ! limit of vector length (0.5-1.0 is recommended)
! !     vrand = 0.15_dp       ! velocity perturbation parameter (0.0-0.1 is recommended)
!     vrand = 0.07_dp
        ! 2017-05-18 add control parameters read from pso.opt file
        w = inertia_in
        c1 = self_in
        c2 = swarm_in
        c3 = randsearch_in
        vmax = vmax_in
        vrand = vrand_in

        ie = 0             ! maximum evaluation count (<=0 means unlimitted)

        !====== Setup control flags =======
        flag = 0           ! control flag
    !! flag = flag + f_quiet   ! restrain messages
        flag = flag + f_initp   ! use initial parameters
        flag = flag + f_lrand   ! local random search algorithm
        flag = flag + f_vlimit  ! velocity limit algorithm
        flag = flag + f_vrand   ! velocity perturbation algorithm
        flag = flag + f_gcenter ! evaluate center of gravity
    !! flag = flag + f_gjump   ! move to center of gravity
    !! flag = flag + f_qmcinit ! initialize by Quasi-Monte Carlo

        !====== Optimize =======
        call pso(n, xmin, xmax, xinit, np, w, c1, c2, c3, vmax, vrand,&
            &it, ie, ir, 0, flag, maximize, evfunc, pso_cbdummy, xall, x, r, iGlobal)

    end subroutine do_pso

    !======================================================================
    ! PSO Simplified Wrapper Routine (with configuration file)
    !======================================================================
    subroutine do_psof(n, xmin, xmax, xinit, maximize, evfunc, pfname, xall, x, r, iGlobal)
        integer, intent(in):: n             ! number of parameter
        real(dp), intent(in):: xmin(1:n)     ! lower bounds
        real(dp), intent(in):: xmax(1:n)     ! upper bounds
        real(dp), intent(in):: xinit(1:n)    ! initial parameters
        integer, intent(in):: maximize      ! maximize-minimize flag (>0:max, <=0:min)
        character(len=*)   :: pfname        ! configuration file name
        real(dp), dimension(:), intent(in)::xall
        real(dp), intent(out)::x(1:n)        ! optimum parameter set
        real(dp), intent(out)::r             ! optimum evaluation value
        integer, intent(in)::iGlobal
        interface
            !-------- evaluation function -------
            integer function evfunc(xall_in, x_in, fitness)  ! return number of constrained condition
                use constants, only: dp
                real(dp), dimension(:), intent(in)::xall_in, x_in    ! parameter set
                real(dp), intent(out) :: fitness       ! evaluation value
            end function evfunc
        end interface

        integer :: i, nn
        integer :: np, it, ie, ir, flag
        real(dp)    :: w, c1, c2, c3, vmax, vrand
        integer :: flag_quiet, flag_initp, flag_vlimit, flag_lrand, flag_vrand, &
            &            flag_gcenter, flag_gjump, flag_qmcinit
        namelist /opti_pso/ np, w, c1, c2, c3, vmax, vrand, it, ie, ir, &
            &            flag_quiet, flag_initp, flag_vlimit, flag_lrand, flag_vrand, &
            &            flag_gcenter, flag_gjump, flag_qmcinit

        open (1, file=pfname, status='old')
        read (1, opti_pso)
        close (1)

        flag = flag_quiet
        flag = flag + flag_initp*2**1
        flag = flag + flag_vlimit*2**2
        flag = flag + flag_lrand*2**3
        flag = flag + flag_vrand*2**4
        flag = flag + flag_gcenter*2**5
        flag = flag + flag_gjump*2**6
        flag = flag + flag_qmcinit*2**7

        nn = 0
        do i = 1, n
            if (xmin(i) /= xmax(i)) then
                nn = nn + 1
            end if
        end do
        c3 = c3*sqrt(real(nn, dp))
        vmax = vmax*sqrt(real(nn, dp))

        call pso(n, xmin, xmax, xinit, np, w, c1, c2, c3, vmax, vrand,&
            &it, ie, ir, 0, flag, maximize, evfunc, pso_cbdummy, xall, x, r, iGlobal)

    end subroutine do_psof

    !======================================================================
    ! PSO (primitive routine)
    !======================================================================
    subroutine pso(n, xmin, xmax, xinit, np, w, c1, c2, c3, vmax, vrand, &
        &it, ie, ir, ic, flag, maximize, evfunc, cbfunc,&
        &xall, xopti, ropti, iGlobal)
        !use parameters,only:nit_pso
        !use util_random, only : init_grand
        integer, intent(in):: n             ! number of parameter
        real(dp), intent(in):: xmin(1:n)     ! lower bounds
        real(dp), intent(in):: xmax(1:n)     ! upper bounds
        real(dp), intent(in):: xinit(1:n)    ! initial parameters
        integer, intent(in):: np            ! population size (number of particle)
        real(dp), intent(in):: w             ! inertia parameter
        real(dp), intent(in):: c1            ! self intention parameter
        real(dp), intent(in):: c2            ! swarm intention parameter
        real(dp), intent(in):: c3            ! random search parameter (f_lrand)
        real(dp), intent(in):: vmax          ! limit of vector (f_vlimit)
        real(dp), intent(in):: vrand         ! velocity parameter (f_vrand)
        integer, intent(in):: it            ! maximum iteration count (<=0:unlimitted)
        integer, intent(in):: ie            ! maximum evaluation count (<=0:unlimitted)
        integer, intent(in):: ir            ! interval of message printing (<=0:nomessage)
        integer, intent(in):: ic            ! interval of callback function (<=0:no call)
        integer, intent(in):: flag          ! control flag
        integer, intent(in):: maximize      ! maximize-minimize flag (>0:max, <=0:min)
        real(dp), dimension(:), intent(in)::xall
        real(dp), intent(out)::xopti(1:n)    ! optimum parameter set
        real(dp), intent(out)::ropti         ! optimum evaluation value
        integer, intent(in)::iGlobal
        interface
            !-------- evaluation function -------
            integer function evfunc(xall_in, x_in, fitness)  ! return number of constrained condition
                use constants, only: dp
                real(dp), dimension(:), intent(in)::xall_in, x_in    ! parameter set
                real(dp), intent(out) :: fitness       ! evaluation value
            end function evfunc
            !-------- callback function -------
            integer function cbfunc(it_in, iev, ev, x) ! return continuation flag(>=0:cont.,<0:quit)
                use constants, only: dp
                integer, intent(in) :: it_in    ! current iteration counter
                integer, intent(in) :: iev   ! current evaluation counter
                real(dp), intent(in) :: ev    ! current evaluation value
                real(dp), intent(in) :: x(:)  ! current parameter value
            end function cbfunc
        end interface

        !===== Local variables =====
        type(type_p) :: p(1:np)          ! particles
        integer      :: ip(1:np)         ! index of particles
        real(dp)     :: gebest           ! swarm best evaluation value
        real(dp)     :: gxbest(1:n)      ! swarm best position
        type(type_p) :: g                ! gravity center of good particles
        integer      :: evcount          ! counter of evaluation
        integer      :: sign             ! sign of evaluation value
        integer      :: best, worst      ! index of best/worst particle
        logical      :: exit_loop        ! exit flag
        integer      :: i, j             ! loop counter

        !-------------------------------------------
        ! variables for Luca Borsato's stuff
        !integer::ipso,nrows,nlast
!     real(dp),dimension(:,:),allocatable::bestpso
        integer, dimension(:), allocatable::idbest
        integer::wrt_it, uall, ubest !it_c
        !-------------------------------------------

        write (*, '(a,i4)') ' RUNNING PSO WITH NCPU = ', ncpu_in

        !===== Print PSO parameters =====
        if (flag_off(flag, f_quiet)) then
            call print_param(n, np, w, c1, c2, c3, vmax, it, ie, flag)
!      write(*,'(a)')'PSO --- cycle / best / mean / constrained(%) / evaluation ---'
        end if

        !===== Check arguments =====
        call check_param(n, xmin, xmax, xinit, np, w, c1, c2, c3, vmax, it, flag)

        !===== Initialize scaling function =====
        call init_scaling(n, xmin, xmax)

    !!....!===== Initialize random number generator =====
!     write(*,'(a)')" Initializing seeds in PSO"
        call init_random_seed_input(np*n, seed_pso + iGlobal - 1)
!     write(*,*)

        !===== Initialize sign of evaluation =====
        if (maximize > 0) then
            sign = 1
        else
            sign = -1
        end if

        !===== Initialize particles =====
        evcount = 0
!     call init_particles(n, xinit, np, flag, p)
!     call init_good_particles(n,np,xall,p)
        call init_population_zero(n, np, p)
        write (*, '(a)') ' PSO: PARTICLES INITITIALISED - show first and last particle positions'
        write (*, '(a)') ' -----'
        write (*, '(a)') ' -----'
        write (*, '(a)') ' first [0,1)'
        write (*, '(1000(1x,es23.16))') p(1)%x(1:n)
        write (*, '(a)') ' first [min,max)'
        write (*, '(1000(1x,es23.16))') unscaling(n, p(1)%x(1:n))
        ! write (*, '(a, 1x, es23.16)') ' PSO fitness',p(1)%ev
        write (*, '(a)') ' last [0,1)'
        write (*, '(1000(1x,es23.16))') p(np)%x(1:n)
        write (*, '(a)') ' last [min,max)'
        write (*, '(1000(1x,es23.16))') unscaling(n, p(np)%x(1:n))
        ! write (*, '(a, 1x, es23.16)') ' PSO fitness',p(np)%ev
        write (*, '(a)') ' -----'
        write (*, '(a)') ' -----'
        gebest = void

        ! added by Luca Borsato, stuff to save in write to a file best particles
        call check_file_pso(iGlobal, uall, ubest)
        call bestpso_init(it, idbest, pso_best_evolution)
        wrt_it = 10
!     it_c=1

        write (*, '(a)') ' PSO: START MAIN LOOP'
        !***********************************************************
        !*********************** Main Loop *************************
        i = 0
        do

            !===== Evaluate all particles =====
            !$omp parallel do NUM_THREADS(ncpu_in) schedule(DYNAMIC,1)&
            !$omp& shared(n,np,xall,p)
            do j = 1, np
!           write(*,*)' ---- '
!           write(*,*)''
                p(j)%ev = evaluate(n, xall, p(j)%x, sign, evfunc, evcount)
                if (abs(p(j)%ev) > p(j)%evbest) then
                    p(j)%evbest = p(j)%ev
                    p(j)%xbest(1:n) = p(j)%x(1:n)
                end if
!           write(*,*)''
!           write(*,*)' ---- '
            end do
            !$omp end parallel do

            !===== Sort particles by evaluation value =====
            call psort(np, p, ip)

            !===== Update swarm best position =====
            best = ip(1)
            if (p(best)%ev > gebest) then
                gxbest(1:n) = p(best)%x(1:n)
                gebest = p(best)%ev
            end if

            !===== Calculate gravity center of good particles =====
            if (flag_on(flag, f_gcenter) .or. flag_on(flag, f_gjump)) then
                g%x(1:n) = gcenter(n, 3, p, ip)
                g%ev = evaluate(n, xall, g%x, sign, evfunc, evcount)
            end if

            !===== Update swarm best position by gravity center =====
            if (flag_on(flag, f_gcenter)) then
                if (g%ev > gebest) then
                    gxbest(1:n) = g%x(1:n)
                    gebest = g%ev
                end if
            end if

            !===== Move worst particle to gravity center =====
            if (flag_on(flag, f_gjump)) then
                worst = ip(np)
                p(worst)%v(1:n) = g%x(1:n) - p(worst)%x(1:n)
                p(worst)%x(1:n) = g%x(1:n)
                p(worst)%ev = g%ev
                if (p(worst)%ev > p(worst)%evbest) then
                    p(worst)%xbest(1:n) = p(worst)%x(1:n)
                    p(worst)%evbest = p(worst)%ev
                end if
            end if

            ! added by Luca...the original is after exit_loop
            i = i + 1

            !===== Check exit conditions =====
            if ((it > 0 .and. i >= it) .or. (ie > 0 .and. evcount >= ie)) then
                exit_loop = .true.
            else
                exit_loop = .false.
            end if

            !===== Call callback function =====
            if ((ic > 0) .and. (mod(i, ic) == 0 .or. exit_loop)) then
                if (cbfunc(i, evcount, gebest, unscaling(n, gxbest(1:n))) < 0) then
                    exit_loop = .true.
                end if
            end if

            !===== Print intermediate conditions =====
            if (flag_off(flag, f_quiet) .and. (ir > 0) &
                  &                              .and. (mod(i, ir) == 0 .or. exit_loop)) then
                write (*, '(a)') ""
                print *, 'PSO ---', i, gebest*sign, evmean(np, p)*sign, &
                  & infeasible(np, p), evcount
                write (*, '(10(a,es23.16))')&
                  &" Best Swarm fitness = ", gebest
                  write (*, '(a)') " With parameters: "
                write (*, '(10000(es23.16,1x))') unscaling(n, gxbest)
                write (*, '(a)') " "
            end if

            ! save all data
!       write(*,*)'iter=',i,' saving to iter = iter = ',i
            if (uall .gt. 0) call write_allpso(uall, i, ip, p)
            call save_population(i, ip, p) ! it modifies population(:,:,i), population_fitness(:,i)
            ! save and store best by Luca Borsato
            !call bestpso_store(i,best,gxbest,gebest,idbest,pso_best_evolution)
            call bestpso_store(ubest, i, best, gxbest, gebest, idbest, pso_best_evolution)
            if (getStat(i, it) .ge. wrt_it) then
                write (*, '(a,i4,a)') " Completed iterations: ", getStat(i, it), " of 100%"
                !call write_bestpso(it_c,i,idbest,bestpso)
                write (*, '(a)') ' Best position'
                write (*, '(a,i4)') '             index best = ', best
                write (*, '(a)') ' Best fitness'
                write (*, '(a,es23.16)') '     fitness(best) = ', gebest
                write (*, *)
!           it_c=i
                wrt_it = wrt_it + 10
            end if

            !===== Exit loop =====
            if (exit_loop) exit

            ! 2017-02-16 -- PARALLEL --
            !===== Move particles =====
            !$omp parallel do NUM_THREADS(ncpu_in) schedule(DYNAMIC,1)&
            !$omp& shared(np,n,gxbest,w,c1,c2,c3,vmax,vrand,flag,p)
            do j = 1, np
                call p_move(n, gxbest, w, c1, c2, c3, vmax, vrand, flag, p(j))
            end do
            !$omp end parallel do
            ! 2017-02-16 -- PARALLEL --

        end do
        !*********************** Main Loop *************************

        !***********************************************************
        !if(it_c.lt.it) call write_bestpso(it_c,it,idbest,bestpso)
        if (allocated(idbest)) deallocate (idbest)
        close (ubest)
        if (uall .gt. 0) close (uall)

        !===== Unscaling parameters =====
        xopti = unscaling(n, gxbest)
        ropti = gebest*sign

        !===== Print exit message =====
        if (flag_off(flag, f_quiet)) then
            print *, 'PSO --- finish ---'
            print *
        end if

    end subroutine pso

    !----------------------------------------------------------------------
    ! Move Particle
    !----------------------------------------------------------------------
    subroutine p_move(n, gxbest, w, c1, c2, c3, vmax, vrand, flag, p)
    !!use util_random, only : grand, vngrand
        integer, intent(in):: n             ! number of parameter
        real(dp), intent(in):: gxbest(1:n)   ! swarm best position
        real(dp), intent(in):: w             ! inertia parameter
        real(dp), intent(in):: c1            ! self intention parameter
        real(dp), intent(in):: c2            ! swarm intention parameter
        real(dp), intent(in):: c3            ! random search parameter
        real(dp), intent(in):: vmax          ! limit of vector
        real(dp), intent(in):: vrand         ! velocity perturbation parameter
        integer, intent(in)    :: flag          ! control flag
        type(type_p), intent(inout) :: p             ! current particle
        real(dp)    :: r1, r2, r3
        real(dp), dimension(n)::x, v, v0, v1, v2, v3
        real(dp)    :: vv

        !===== Prep. random number =====
    !!r1 = grand()
    !!r2 = grand()
        r1 = random_scalar()
        r2 = random_scalar()
        r3 = one
        v3 = zero

        !---- Velocity random expansion (optional) ----
        if (flag_on(flag, f_vrand)) then
      !!r3 = (one - vrand) + grand() * vrand
            r3 = (one - vrand) + random_scalar()*vrand
        end if

        !---- Local random search (optional) ----
        if (flag_on(flag, f_lrand)) then
      !!v3(1:n) = c3 * vngrand(n)
            call random_number(v3)
            v3 = c3*v3
        end if

        !===== Calculate velocity vector =====
        v0(1:n) = w*p%v(1:n)
        v1(1:n) = c1*r1*(p%xbest(1:n) - p%x(1:n))
        v2(1:n) = c2*r2*(gxbest(1:n) - p%x(1:n))
        v(1:n) = v0(1:n) + v1(1:n) + v2(1:n)

        !===== Limit velocity length (optional) =====
        if (flag_on(flag, f_vlimit)) then
            vv = vabs(n, v)
            if (vv > vmax) then
                v(1:n) = v(1:n)*vmax/vv
            end if
        end if

        !===== Update velocity vector =====
        v(1:n) = r3*v(1:n) + v3(1:n)

        !===== perturbation =====
        vv = vabs(n, v)
        if (vv == zero) then
      !!v(1:n) = c3 * vngrand(n)
            call random_number(v)
            v = c3*v
        end if

        !===== Update position =====
        x(1:n) = p%x(1:n) + v(1:n)

        !===== Reflection =====
        call reflect(n, x, v)

        !===== Update particle =====
        p%x(1:n) = x(1:n)
        p%v(1:n) = v(1:n)

    end subroutine p_move

    !----------------------------------------------------------------------
    ! Reflection
    !----------------------------------------------------------------------
    subroutine reflect(n, x, v)
        integer, intent(in)    :: n
        real(dp), intent(inout) :: x(1:n), v(1:n)
        integer :: i
        logical :: retry
        retry = .true.
        do while (retry)
            retry = .false.
            do i = 1, n
                if (x(i) < zero) then
                    x(i) = -one*x(i)
                    !x(i) = zero
                    v(i) = -one*v(i)
                    retry = .true.
                else if (x(i) > one) then
                    x(i) = two - one*x(i)
                    !x(i) = one
                    v(i) = -one*v(i)
                    retry = .true.
                end if
            end do
        end do
    end subroutine reflect

    !----------------------------------------------------------------------
    ! Calculate Vector Length
    !----------------------------------------------------------------------
    function vabs(n, v) result(r)
        integer, intent(in) :: n
        real(dp), intent(in) :: v(:)
        real(dp)                :: r
        integer :: i
        r = zero
        do i = 1, n
            r = r + v(i)**2
        end do
        r = sqrt(r)
    end function vabs

    !----------------------------------------------------------------------
    ! Check Arguments
    !----------------------------------------------------------------------
    subroutine check_param(n, xmin, xmax, xinit, np, w, c1, c2, c3, vmax, it, flag)
        integer, intent(in)  :: n, np, it, flag
        real(dp), intent(in)  :: xmin(1:n), xmax(1:n), xinit(1:n)
        real(dp), intent(in)  :: w, c1, c2, c3, vmax
        integer :: i, err

        err = 0

        if (n > max_nparam) then
            print *, 'pso: error: n > max_nparam', n
            err = err + 1
        end if
        if (np < 1) then
            print *, 'pso: error: np < 1', np
            err = err + 1
        end if
        if (w < zero) then
            print *, 'pso: error: w < 0', w
            err = err + 1
        end if
        if (c1 < zero) then
            print *, 'pso: error: c1 < 0', c1
            err = err + 1
        end if
        if (c2 < zero) then
            print *, 'pso: error: c2 < 0', c2
            err = err + 1
        end if
        if (c3 < zero) then
            print *, 'pso: error: c3 < 0', c3
            err = err + 1
        end if
        if (vmax <= zero) then
            print *, 'pso: error: Rmax <= 0', vmax
            err = err + 1
        end if
        if (it < 0) then
            print *, 'pso: error: it < 0', it
            err = err + 1
        end if

        do i = 1, n
            if (xmax(i) < xmin(i)) then
                print *, 'pso: error: xmax < xmin', i
                err = err + 1
            end if
            if (flag_on(flag, f_initp)) then
                if (xinit(i) < xmin(i)) then
                    print *, 'pso: error: xinit < xmin', i, xinit(i), xmin(i)
                    err = err + 1
                end if
                if (xinit(i) > xmax(i)) then
                    print *, 'pso: error: xinit > xmax', i, xinit(i), xmax(i)
                    err = err + 1
                end if
            end if
        end do

        if (err > 0) then
            stop
        end if
    end subroutine check_param

    !----------------------------------------------------------------------
    ! Print PSO Parameters
    !----------------------------------------------------------------------
    subroutine print_param(n, np, w, c1, c2, c3, vmax, it, ie, flag)
        integer, intent(in) :: n, np, it, ie, flag
        real(dp), intent(in) :: w, c1, c2, c3, vmax
        write (*, *)
        write (*, *) '****** PSO parameters ******'
        write (*, *) 'n    :', n
        write (*, *) 'np   :', np
        write (*, *) 'w    :', w
        write (*, *) 'c1   :', c1
        write (*, *) 'c2   :', c2
        write (*, *) 'c3   :', c3
        write (*, *) 'Vmax :', vmax
        write (*, *) 'it   :', it
        write (*, *) 'ie   :', ie
        write (*, *) 'flag_quiet   :', chk_flag(flag, f_quiet)
        write (*, *) 'flag_initp   :', chk_flag(flag, f_initp)
        write (*, *) 'flag_vlimit  :', chk_flag(flag, f_vlimit)
        write (*, *) 'flag_lrand   :', chk_flag(flag, f_lrand)
        write (*, *) 'flag_vrand   :', chk_flag(flag, f_vrand)
        write (*, *) 'flag_gcenter :', chk_flag(flag, f_gcenter)
        write (*, *) 'flag_gjump   :', chk_flag(flag, f_gjump)
        write (*, *) 'flag_qmcinit :', chk_flag(flag, f_qmcinit)
        write (*, *)
    end subroutine print_param

    !----------------------------------------------------------------------
    ! Initialize Particles
    !----------------------------------------------------------------------
    subroutine init_particles(n, xinit, np, flag, p)
        !use util_random, only : vgrand
        use util_qmc, only: hammersley
        integer, intent(in)    :: n             ! number of parameter
        real(dp), intent(in)    :: xinit(1:n)    ! initial parameters
        integer, intent(in)    :: np            ! number of particle
        integer, intent(in)    :: flag          ! control flag
        type(type_p), intent(out)   :: p(1:np)       ! particles
        integer :: i, j
        integer :: na, ip(1:n)
        real(dp)    :: x(1:np*n)

        if (flag_on(flag, f_qmcinit)) then
!       write(*,'(a)')' flag on --> hammersley'
            na = 0                           ! number of effective dimension
            do j = 1, n
                if (scale_a(j) > zero) then      ! not fixed parameter
                    na = na + 1
                    ip(j) = na
                end if
            end do
            call hammersley(np, na, x)
            do i = 1, np
                do j = 1, n
                    if (scale_a(j) > zero) then
                        p(i)%x(j) = x((i - 1)*na + ip(j))
                    else
                        p(i)%x(j) = zero
                    end if
                end do
            end do
        else
!       write(*,'(a)')' flag off --> random_number'
            do i = 1, np
                !p(i)%x(1:n) = vgrand(n)
!           call random_number(p(i)%x(1:n))
                call random_number(p(i)%x(:))
            end do
        end if

        if (flag_on(flag, f_initp)) then
            p(1)%x(1:n) = scaling(n, xinit)
            p(1)%v(1:n) = zero
            p(1)%evbest = void
            p(1)%xbest = p(1)%x
        end if

        do i = 1, np
            p(i)%v(1:n) = zero
            p(i)%evbest = void
            p(i)%xbest(1:n) = p(i)%x(1:n)
        end do

    end subroutine init_particles

    ! 2016-07-20 Luca Borsato own method to initialize population
    subroutine init_good_particles(n_fit, n_pop, all_parameters, p)
        integer, intent(in)::n_fit, n_pop            ! number of particle
        real(dp), dimension(:), intent(in)::all_parameters
        type(type_p), dimension(:), intent(out)::p ! particles
        real(dp), dimension(:), allocatable::xtemp, ptemp
        integer::iloop
        logical::check

!     write(*,*)' size of p = ',size(p)
!     flush(6)
        allocate (xtemp(n_fit), ptemp(n_fit))
        iloop = 0
        initloop: do

            xtemp = zero
            call random_number(xtemp)
            ptemp = unscaling(n_fit, xtemp)
            check = .true.
            check = check_only_boundaries(all_parameters, ptemp)
            if (check) then
!           write(*,*)' par (True) = ',unscaling(n_fit,xtemp)
                iloop = iloop + 1
                p(iloop)%x(1:n_fit) = xtemp(1:n_fit)
                p(iloop)%v(1:n_fit) = zero
                p(iloop)%evbest = void
                p(iloop)%xbest(1:n_fit) = p(iloop)%x(1:n_fit)
                if (iloop .eq. n_pop) exit initloop
            end if

        end do initloop
        deallocate (xtemp)

        return
    end subroutine init_good_particles

    !2017-02-15 AIMED TO SIMPLICITY
    subroutine init_population_zero(n_fit, n_pop, p)
        integer, intent(in)::n_fit, n_pop            ! number of particle
        type(type_p), dimension(:), intent(out)::p ! particles
        real(dp), dimension(:), allocatable::xtemp ! position in [0,1)
        integer::iloop

        allocate (xtemp(n_fit))
        initloop: do iloop = 1, n_pop
            xtemp = zero
            call random_number(xtemp)
            p(iloop)%x(1:n_fit) = xtemp(1:n_fit)
            p(iloop)%v(1:n_fit) = zero
            p(iloop)%evbest = void
            p(iloop)%xbest(1:n_fit) = p(iloop)%x(1:n_fit)
        end do initloop
        deallocate (xtemp)

        return
    end subroutine init_population_zero

    !----------------------------------------------------------------------
    ! Call Evaluation Function
    !----------------------------------------------------------------------
    function evaluate(n, xall, x, sign, evfunc, evcount) result(r)
        integer, intent(in)    :: n
        real(dp), dimension(:), intent(in):: xall, x
        integer, intent(in)    :: sign
        integer, intent(inout) :: evcount
        real(dp)                   :: r
        interface
            integer function evfunc(xall_in, x_in, fitness)
                use constants, only: dp
                real(dp), dimension(:), intent(in)::xall_in, x_in
                real(dp), intent(out) :: fitness
            end function evfunc
        end interface
        if (evfunc(xall, unscaling(n, x), r) == 0) then
            r = r*sign
        else
            r = constrained
        end if
        evcount = evcount + 1
    end function evaluate

    !----------------------------------------------------------------------
    ! Sort Particles by Evaluation
    !----------------------------------------------------------------------
    subroutine psort(np, p, ip)
        use util_sort, only: qsort
        integer, intent(in)  :: np
        type(type_p), intent(in)  :: p(1:np)
        integer, intent(out) :: ip(1:np)
        real(dp)    :: x(1:np)
        x(1:np) = p(1:np)%ev
        call qsort(x, ip, -1)
    end subroutine psort

    !----------------------------------------------------------------------
    ! Calculate Gravity Center of Particles
    !----------------------------------------------------------------------
    function gcenter(n, np, p, ip) result(x)
        integer, intent(in) :: n           ! number of parameter
        integer, intent(in) :: np          ! number of particle
        type(type_p), intent(in) :: p(:)        ! particles
        integer, intent(in) :: ip(:)       ! index of particles
        real(dp)::x(1:n)
        integer::i
        x(1:n) = zero
        do i = 1, np
            x(1:n) = x(1:n) + p(ip(i))%x(1:n)
        end do
        x(1:n) = x(1:n)/real(np, dp)
    end function gcenter

    !----------------------------------------------------------------------
    ! Calculate Mean Evaluation Value
    !----------------------------------------------------------------------
    function evmean(np, p) result(r)
        integer, intent(in) :: np
        type(type_p), intent(in) :: p(1:np)
        real(dp)                     :: x, r
        integer :: i, j
        r = zero
        j = 0
        do i = 1, np
            x = p(i)%ev
            if (x > constrained) then
                r = r + x
                j = j + 1
            end if
        end do
        if (j > 0) then
            r = r/real(j, dp)
        else
            r = constrained
        end if
    end function evmean

    !----------------------------------------------------------------------
    ! Calculate Ratio of Restricted Particles
    !----------------------------------------------------------------------
    function infeasible(np, p) result(r)
        integer, intent(in) :: np
        type(type_p), intent(in) :: p(1:np)
        integer                  :: r
        integer :: i, j
        j = 0
        do i = 1, np
            if (p(i)%ev <= constrained) then
                j = j + 1
            end if
        end do
        r = int(real(j, dp)/real(np, dp)*100.0_dp)
    end function infeasible

    !----------------------------------------------------------------------
    ! Initialize Parameter Scaling Functions
    !----------------------------------------------------------------------
    subroutine init_scaling(n, xmin, xmax)
        integer, intent(in) :: n
        real(dp), intent(in) :: xmin(1:n), xmax(1:n)
!     scale_a(1:n) = xmax(1:n) - xmin(1:n)
        scale_a(1:n) = abs(xmax(1:n) - xmin(1:n))
        scale_b(1:n) = xmin(1:n)
    end subroutine init_scaling

    !----------------------------------------------------------------------
    ! Parameter Scaling (normalize)
    !----------------------------------------------------------------------
    function scaling(n, x) result(px)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(1:n)
        real(dp)                :: px(1:n)
        px(1:n) = zero
        where (scale_a(1:n) /= zero) px(1:n) = (x(1:n) - scale_b(1:n))/scale_a(1:n)
    end function scaling

    !----------------------------------------------------------------------
    ! Parameter Unscaling (denormalize)
    !----------------------------------------------------------------------
    function unscaling(n, px) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: px(1:n)
        real(dp)                :: x(1:n)
        x(1:n) = px(1:n)*scale_a(1:n) + scale_b(1:n)
    end function unscaling

    !----------------------------------------------------------------------
    ! Check Flag
    !----------------------------------------------------------------------
    function chk_flag(flag, mask) result(r)
        integer, intent(in) :: flag, mask
        integer             :: r
        if (mask <= 0) then
            print *, 'error: invalid mask for chk_flag'
            stop
        end if
        r = mod(flag/mask, 2)
    end function chk_flag

    !----------------------------------------------------------------------
    ! Check Flag is ON
    !----------------------------------------------------------------------
    function flag_on(flag, mask) result(r)
        integer, intent(in) :: flag, mask
        logical             :: r
        if (chk_flag(flag, mask) > 0) then
            r = .true.
        else
            r = .false.
        end if
    end function flag_on

    !----------------------------------------------------------------------
    ! Check Flag is OFF
    !----------------------------------------------------------------------
    function flag_off(flag, mask) result(r)
        integer, intent(in) :: flag, mask
        logical             :: r
        r = .not. flag_on(flag, mask)
    end function flag_off

    !----------------------------------------------------------------------
    ! Dummy (sample) callback function
    !----------------------------------------------------------------------
    function pso_cbdummy(it, iev, ev, x) result(r)
        integer, intent(in) :: it    ! iteration counter
        integer, intent(in) :: iev   ! current evaluation counter
        real(dp), intent(in) :: ev    ! current evaluation value
        real(dp), intent(in) :: x(:)  ! current parameter value
        integer             :: r     ! continuation flag (>=0: cont., <0: quit)
        real(dp), parameter :: eps = 1.e-5_dp
        print *, 'PSO ---', it, ev, iev, x(:)
        if (ev < eps) then
            r = -1
        else
            r = 1
        end if
    end function pso_cbdummy

    ! ---
    ! PSO STORE AND WRITING FILES
    ! ---

    ! check if the pso files where store simulations exists and creates new ones
    subroutine check_file_pso(iGlobal, uall, ubest)
        integer, intent(in)::iGlobal
        integer, intent(out)::uall, ubest
        character(512)::flpso, fmtpso, remove
        logical::exstat
        integer::ip

        flpso = trim(path)//trim(adjustl(string(iGlobal)))//"_bestpso_"//trim(adjustl(string(nfit)))//"par.out"
        flpso = trim(adjustl(flpso))
        inquire (file=trim(flpso), exist=exstat)
        if (exstat) then
            remove = "rm "//trim(flpso)
            remove = trim(adjustl(remove))
            call system(remove)
        end if
        ubest = get_unit(1)
        open (ubest, file=trim(flpso))
        fmtpso = "# iter part "
        do ip = 1, nfit
            fmtpso = trim(adjustl(fmtpso))//&
                  &" "//parid(ip)
        end do
        fmtpso = trim(adjustl(fmtpso))//" pso_fitness bic"
        write (ubest, '(a)') trim(adjustl(fmtpso))
        !close(ubest)

        ! create file to write all iterations all particles
        flpso = trim(path)//trim(adjustl(string(iGlobal)))//"_allpso_"//trim(adjustl(string(nfit)))//"par.out"
        inquire (file=trim(flpso), exist=exstat)
        if (exstat) then
            remove = "rm "//trim(flpso)
            remove = trim(adjustl(remove))
            call system(remove)
        end if
        if (wrtAll .eq. 0) then
            uall = -1
        else
            uall = get_unit(1)
            open (uall, file=trim(adjustl(flpso)))
            write (uall, '(a)') trim(adjustl(fmtpso))
        end if

        return
    end subroutine check_file_pso

    subroutine bestpso_store(ubest, ipos, ipx, x_sc, pso_fitness, idbest, bestpso)
        integer, intent(in)::ubest, ipos, ipx
        real(dp), dimension(:), intent(in)::x_sc
        real(dp), intent(in)::pso_fitness
        integer, dimension(:), intent(out)::idbest
        real(dp), dimension(:, :)::bestpso
        real(dp), dimension(:), allocatable::x_unsc
        real(dp)::bic
        character(512)::fmtw

        allocate (x_unsc(nfit))
        x_unsc = unscaling(nfit, x_sc(1:nfit)) ! scale parameters to physical value
        bestpso(1:nfit, ipos) = x_unsc(1:nfit)
        deallocate (x_unsc)

        bic = -two*pso_fitness + bic_const ! bic_const global variable

        bestpso(nfit + 1, ipos) = pso_fitness
        bestpso(nfit + 2, ipos) = bic
        idbest(ipos) = ipx
        ! added 2014-07-14: write best each iteration
        fmtw = adjustl("(i6,1x,i6,1x,10000("//sprec//",1x))")
        write (ubest, trim(fmtw)) ipos, idbest(ipos), bestpso(:, ipos)
        flush (ubest)

        return
    end subroutine bestpso_store

    ! write best simulations into file
    subroutine write_bestpso(it_i, it_e, idbest, bestpso)
!$      use omp_lib
        integer, intent(in)::it_i, it_e
        integer, dimension(:), intent(in)::idbest
        real(dp), dimension(:, :), intent(in)::bestpso
        integer::cpuid, upso, ii
        character(512)::flpso, fmtw

        cpuid = 1
!$      cpuid = omp_get_thread_num() + 1
        flpso = trim(path)//"bestpso_"//trim(adjustl(string(nfit)))//"par.out"
        flpso = trim(adjustl(flpso))
        fmtw = adjustl("(i6,1x,i6,10000("//sprec//"))")
        upso = get_unit(cpuid)
        open (upso, file=trim(flpso), access='append')
        do ii = it_i, it_e
            write (upso, trim(fmtw)) ii, idbest(ii), bestpso(:, ii)
        end do
        close (upso)

        return
    end subroutine write_bestpso

    subroutine write_allpso(uall, it, ip, pp)
        integer, intent(in)::uall, it
        type(type_p), dimension(:), intent(in)::pp ! particles
        integer, dimension(:), intent(in)::ip ! index of particles
        integer::np
        integer::j
        character(128)::wfmt
        real(dp)::pso_fitness, bic
        real(dp), dimension(:), allocatable::xpar

        np = size(ip)
        wfmt = "(i6,1x,i4,1x,1000("//trim(sprec)//",1x))"

        allocate (xpar(nfit))

        do j = 1, np
            pso_fitness = pp(j)%ev
            ! fitness = one/pso_fitness
            ! fitness_x_dof = fitness*real(obsData%dof, dp)
            bic = -two*pso_fitness + bic_const ! bic_const global variable
            xpar = zero
            xpar = unscaling(nfit, pp(j)%x(:))
            write (uall, trim(wfmt)) it, ip(j), xpar, pso_fitness, bic
            flush (uall)
        end do

        deallocate (xpar)

        return
    end subroutine write_allpso

    subroutine save_population(it, ip, pp)
        integer, intent(in)::it
        integer, dimension(:), intent(in)::ip ! index of particles
        type(type_p), dimension(:), intent(in)::pp ! particles
        real(dp)::fitness
        integer::np
        integer::j
        np = size(ip)

        do j = 1, np
            fitness = pp(j)%ev ! get pso_fitness
            population_fitness(j, it) = fitness
            population(:, j, it) = unscaling(nfit, pp(j)%x(:))
        end do

        return
    end subroutine save_population

end module opti_pso
