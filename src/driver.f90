! collection of most of the drivers needed to run the different algorithms and summary
module driver
    use constants
    use custom_type
    use parameters
    use parameters_conversion
    use timing, only: timer
    use init_trades, only: initu, init_param
    use ode_run, only: ode_output, ode_lm
    use output_files, only: write_parameters, write_grid_summary
    use Levenberg_Marquardt
    use bootstrap, only: strap_driver
    use fitness_module
    use grid_search
    use Genetic_Algorithm, only: ga_driver, fpik
    use opti_pso, only: pso_driver, evaluate_pso

    implicit none

contains

    ! ============================================================================

    ! given the set of parameters to fit it integrates the orbits,
    ! computes the RVs and the T0s, and writes everything to screen and into files
    ! in this subroutine there will be not sigma/errors on the parameters
    subroutine write_summary_nosigma(cpuid, sim_id, lm_flag, all_parameters, fit_parameters,&
      &chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic)
        ! Input
        integer, intent(in)::cpuid, sim_id, lm_flag
        real(dp), dimension(:), intent(in)::all_parameters, fit_parameters
        ! Output
        real(dp), intent(out)::chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic

        ! Local variables
        real(dp), dimension(:), allocatable::resw
        integer::nobs ! == obsData%ndata == ndata

        nobs = obsData%ndata
        allocate (resw(nobs))
        resw = zero
        lnprior = zero

        ! integrates and write RVs and T0s
        call ode_output(cpuid, sim_id, lm_flag, all_parameters, fit_parameters, resw)
        if (nobs .gt. 0) call write_parameters(cpuid, sim_id, lm_flag, fit_parameters, resw)

        if (nobs .gt. 0) then
            chi_square = sum(resw*resw)
            if (chi_square .ge. resmax) then
                chi_square = resmax
            end if
            call set_fitness_values(fit_parameters, chi_square,&
              &reduced_chi_square, lnLikelihood, ln_const, bic)
            if (n_priors .gt. 0) then
                call ln_priors(all_parameters, fit_parameters,&
                  &priors_names, priors_values, lnprior)
            end if
        else
            chi_square = -one
            reduced_chi_square = -one
            lnLikelihood = one
            ln_const = zero
        end if
        if (allocated(resw)) deallocate (resw)

        write (*, *) "FITNESS SUMMARY:"
        write (*, *) "chi_square         = ", chi_square
        write (*, *) "reduced_chi_square = ", reduced_chi_square
        write (*, *) "lnLikelihood       = ", lnLikelihood
        write (*, *) "lnprior            = ", lnprior
        write (*, *) "ln_const           = ", ln_const
        write (*, *) "bic                = ", bic
        write (*, *)

        flush (6)

        return
    end subroutine write_summary_nosigma

    ! ============================================================================

    ! lmon == 1
    ! given the set of parameters to fit it integrates the orbits, computes the RVs and the T0s, and writes everything to screen and into files
    ! in this subroutine there will be sigma/errors on the parameters
    subroutine write_summary_sigma(cpuid, sim_id, lm_flag, all_parameters, fit_parameters, sigma_parameters,&
      &chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic)
        ! Input
        integer, intent(in)::cpuid, sim_id, lm_flag
        real(dp), dimension(:), intent(in)::all_parameters
        real(dp), dimension(:), intent(in)::fit_parameters, sigma_parameters
        ! Output
        real(dp), intent(out)::chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic

        ! Local variables
        real(dp), dimension(:), allocatable::resw

        integer::nobs ! == obsData%ndata == ndata
        nobs = obsData%ndata

        allocate (resw(nobs))
        resw = zero
        lnprior = zero
        ! integrates and write RVs and T0s
        call ode_output(cpuid, sim_id, lm_flag, all_parameters, fit_parameters, resw)
        if (nobs .gt. 0) call write_parameters(cpuid, sim_id, lm_flag, fit_parameters, resw)

        if (nobs .gt. 0) then
            chi_square = sum(resw*resw)
            if (chi_square .ge. resmax) then
                chi_square = resmax
            end if
            call set_fitness_values(fit_parameters, chi_square,&
              &reduced_chi_square, lnLikelihood, ln_const, bic)
            if (n_priors .gt. 0) then
                call ln_priors(all_parameters, fit_parameters,&
                  &priors_names, priors_values, lnprior)
            end if
        else
            chi_square = one
            reduced_chi_square = one
            lnLikelihood = one
            ln_const = one
        end if
        if (allocated(resw)) deallocate (resw)

        write (*, *) "FITNESS SUMMARY:"
        write (*, *) "chi_square         = ", chi_square
        write (*, *) "reduced_chi_square = ", reduced_chi_square
        write (*, *) "lnLikelihood       = ", lnLikelihood
        write (*, *) "lnprior            = ", lnprior
        write (*, *) "ln_const           = ", ln_const
        write (*, *) "bic                = ", bic
        write (*, *) ""

        flush (6)

        return
    end subroutine write_summary_sigma

    ! ============================================================================

    ! given the set of parameters it fits with the LM algorithm and then
    ! it integrates the orbits with the best-fit configuration
    ! and write the summary files
    subroutine run_and_write_lm(cpuid, sim_id, all_parameters, fit_parameters)
        use parameters_conversion, only: check_only_boundaries
        ! Input
        integer, intent(in)::cpuid, sim_id
        ! Input/Output
        real(dp), dimension(:), intent(inout)::all_parameters, fit_parameters

        ! Local variables
        real(dp), dimension(:), allocatable::resw
        real(dp), dimension(:), allocatable::covariance_parameters
        real(dp), dimension(:), allocatable::sigma_parameters
        integer::info
        integer, dimension(:), allocatable::iwa
        real(dp)::chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic

        integer::j1
        ! logical::check_lm

        integer::nobs ! == obsData%ndata == ndata
        nobs = obsData%ndata

        allocate (resw(nobs), iwa(nfit), covariance_parameters(nfit),&
          &sigma_parameters(nfit))

        write (*, '(a)') ''
        write (*, '(a)') ' STARTING LM FIT'
        flush (6)

        call lm_driver(cpuid, sim_id, ode_lm, all_parameters, nobs, nfit,&
          &fit_parameters, resw, covariance_parameters, sigma_parameters,&
          &info, iwa) ! IT CALLS L-M

        write (*, '(a,i4,a)') ' ENDED LM FIT (info = ', info, ' )'
        write (*, '(a)') ''

        write (*, '(a,es23.16)') ' NAME FITTED PARAMETERS  ( min | max ) fitness = sum(resw*resw) = ', sum(resw*resw)
        do j1 = 1, nfit
            write (*, '(a10,1x,f25.15,2(a,f25.15),a)') parid(j1), fit_parameters(j1), ' ( ', minpar(j1), ' | ', maxpar(j1), ' )'
        end do
        write (*, '(a)') ''
        write (*, '(a)') 'WRITE SUMMARY'
        write (*, '(a)') ''
        flush (6)

        call write_summary_sigma(cpuid, sim_id, 1, all_parameters, fit_parameters, sigma_parameters,&
            &chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic)

        deallocate (resw, iwa, covariance_parameters, sigma_parameters)

        flush (6)

        return
    end subroutine run_and_write_lm

    ! ============================================================================

    ! this is a simple call of the driver of the bootstrap in the 'boostrap' module.
    subroutine run_bootstrap(sim_id, all_parameters, fit_parameters)
        integer::sim_id
        real(dp), dimension(:)::all_parameters, fit_parameters

        call strap_driver(sim_id, all_parameters, fit_parameters)

        flush (6)

        return
    end subroutine run_bootstrap

    ! ============================================================================

    ! init and run grid search
    subroutine run_grid(all_parameters)
!$      use omp_lib
        real(dp), dimension(:), intent(in)::all_parameters

        integer::cpuid
        real(dp), dimension(:, :), allocatable::perturber_grid
        real(dp), dimension(:, :), allocatable::fitness_grid
        integer::n_grid, sim_id
!     integer::i

        real(dp), dimension(:), allocatable::cpu_all_parameters
        real(dp), dimension(:), allocatable::cpu_fit_parameters
        real(dp)::chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic

        real(dp), dimension(:, :), allocatable::original_grid_summary
        real(dp), dimension(:, :), allocatable::fitted_grid_summary ! used only if lmon==1

        cpuid = 1 ! set initially to 1

        ! read input file <-> idpert
        call read_parameters_grid(cpuid, perturber_parameters_grid)
        ! perturber_parameters_grid in 'parameters' module

        ! set properly the fields of the perturber_parameters_grid variable
        call set_parameters_grid(perturber_parameters_grid, n_grid)

!     ! print to screen to debug
!     do i=1,10
!       write(*,*)'name = ',perturber_parameters_grid(i)%name
!       write(*,*)'input values = ',perturber_parameters_grid(i)%input_values
!       write(*,*)'step type = ',perturber_parameters_grid(i)%step_type
!       write(*,*)'number of steps = ',perturber_parameters_grid(i)%n_steps
!       write(*,*)'step size = ',perturber_parameters_grid(i)%step_grid
!     end do

        ! create/build the full grid, with all the combination of the parameters of perturber body
        call build_grid(MR_star(1, 1), perturber_parameters_grid,&
          &perturber_grid, fitness_grid, n_grid)

        ! print to screen to debug
!     write(*,*)
!     write(*,*)'--------------'
!     write(*,*)'perturber_grid'
!     write(*,*)'--------------'
!     do i=1,n_grid
!       write(*,'(a,I4,a,10(F26.15))')'combination',i,' = ',perturber_grid(i,:)
!     end do

        write (*, *)
        write (*, *) 'n_grid = ', n_grid

        ! needed to update the max allowed value of the semi-major axis
        amax = 5.0_dp*maxval(perturber_grid(:, 4))
!     write(*,'(2(a,f20.12))')' amin = ',amin,' amax = ',amax

        ! ------------------------------------------
        ! here run the grid search loop with openmp!

        write (*, *)
        write (*, '(a)') ' RUN GRID SEARCH'
        write (*, *)
        flush (6)

        allocate (original_grid_summary(n_grid, npar+2))
        if (lmon .eq. 1) allocate (fitted_grid_summary(n_grid, npar+2))

        ! first initialise the openMP stuff

        !$omp parallel NUM_THREADS(ncpu_in) default(shared) &
        !$omp private(cpuid,ncpu)

!$      cpuid = omp_get_thread_num()+1
!$      ncpu = omp_get_num_threads()
        !$omp master
!$      call initu(nfiles, ncpu)
        !$omp end master
        !$omp barrier

        !
        !+++++++
        ! sketch
        !+++++++
        !
        !1. read perturber_grid(i_cpu) -> udpate all_parameters -> cpu_all_parameters
        !2. update fit_parameters with cpu_all_parameters
        !3. calculates the fitness without fit parameters -> for each cpu write file with cpu_all_parameters and fitness (n_files = n_cpu)
        !4. then writes summary for each combination
        !5. if lmon = 1 --> lm fit
        !6. new fit_parameters and update cpu_all_parameters -> for each cpu write file with cpu_all_parameters and fitness
        !7. if nboot > 0 --> bootstrap
        !8. write full summary of full parameters and fitness
        !+++++++
        !

        ! openMP private variables
        !$omp do private(sim_id,cpu_all_parameters,cpu_fit_parameters,&
        !$omp& chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic) &
        !$omp& schedule(dynamic,1)

        ! start the loop on the simulation grid
        grid_loop: do sim_id = 1, n_grid

            ! 0. allocate private variables and initialise them
            allocate (cpu_all_parameters(npar), cpu_fit_parameters(nfit))
            cpu_all_parameters = all_parameters

            ! 1. select proper set of parameters to use: perturber_grid to cpu_all_parameters
            call perturber_grid2parameters(sim_id, idpert, perturber_grid, cpu_all_parameters)
            ! 2. update from cpu_all_parameters to cpu_fit_parameters
            call init_param(cpu_all_parameters, cpu_fit_parameters)

            ! 3. integrates and calculates the fitness -> updated fitness_grid
            ! fitness=base_fitness_function(cpu_all_parameters,cpu_fit_parameters)
            ! fitness_x_dof=fitness*real(obsData%dof,dp)
            ! fitness_grid(sim_id,1)=fitness_x_dof
            ! fitness_grid(sim_id,2)=fitness
            call base_fitness_function(cpu_all_parameters, cpu_fit_parameters,&
              &chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic)
            fitness_grid(sim_id, 1) = chi_square
            fitness_grid(sim_id, 2) = reduced_chi_square

!       write(*,'(a,i6,a,i3,a,1000(F16.7))')' sim_id = ',sim_id,' (cpu = ',cpuid,') all_parameters = ',cpu_all_parameters
!       write(*,'(a,i6,a,i3,a,F16.5,a,1000(F16.7))')' sim_id = ',sim_id,' (cpu = ',cpuid,') fitness = ',fitness,' <-> fit_parameters = ',cpu_fit_parameters
!       flush(6)

            ! 4. save cpu_all_parameters and fitness to file, write summary
            call write_summary_nosigma(cpuid, sim_id, 0, cpu_all_parameters, cpu_fit_parameters,&
              &chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic)
            flush (6)

            original_grid_summary(sim_id, 1:npar) = cpu_all_parameters
            original_grid_summary(sim_id, npar+1:npar+2) = fitness_grid(sim_id, :)
            write (*, *)
            write (*, '(a,i4)') ' FINISHED NOFIT SIM NUMBER ', sim_id
            write (*, *)
            flush (6)

            ! 5. if LM yes
            if (lmon .eq. 1) then
                ! run lm and write summary
                call run_and_write_lm(cpuid, sim_id, cpu_all_parameters, cpu_fit_parameters)
                write (*, *)
                write (*, '(a,i4)') ' FINISHED FIT SIM NUMBER ', sim_id
                write (*, *)
                flush (6)
                ! 6. updates from cpu_fit_parameters to cpu_all_parameters and save them
                call update_parameters_fit2all(cpu_fit_parameters, cpu_all_parameters)
                ! TODO save cpu_all_parameters
                ! fitness=base_fitness_function(cpu_all_parameters,cpu_fit_parameters)
                ! fitness_x_dof=fitness*real(obsData%dof,dp)
                call base_fitness_function(cpu_all_parameters, cpu_fit_parameters,&
                  &chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic)
                fitted_grid_summary(sim_id, 1:npar) = cpu_all_parameters
                fitness_grid(sim_id, 1) = chi_square
                fitness_grid(sim_id, 2) = reduced_chi_square
                ! fitted_grid_summary(sim_id,npar+1) = fitness_x_dof
                ! fitted_grid_summary(sim_id,npar+2) = fitness
                flush (6)
                write (*, *)
                write (*, '(a,i4,a,es23.16)') ' REDUCED CHI SQUARE FOR SIM NUMBER ', sim_id, ': ', reduced_chi_square
                write (*, *)

            end if

            ! 7. if BOOTSTRAP yes
            if (nboot > 0) then
                call run_bootstrap(sim_id, cpu_all_parameters, cpu_fit_parameters)
                write (*, '(a,i3,a,i7)') 'CPU:', cpuid, ' END BOOTSTRAP FOR SIM #', sim_id
                flush (6)
            end if

            deallocate (cpu_all_parameters, cpu_fit_parameters)

        end do grid_loop
        ! end grid

        !$omp end do
        !$omp barrier

        !$omp end parallel
        ! ------------------------------------------

        call write_grid_summary(1, n_grid, 0, original_grid_summary)
        if (lmon .eq. 1) call write_grid_summary(1, n_grid, 1, fitted_grid_summary)

        write (*, *)
        write (*, '(a)') ' END GRID SEARCH'
        write (*, *)
        flush (6)

        if (allocated(perturber_grid)) deallocate (perturber_grid, fitness_grid)
        if (allocated(original_grid_summary)) deallocate (original_grid_summary)
        if (allocated(fitted_grid_summary)) deallocate (fitted_grid_summary)

        return
    end subroutine run_grid

    ! ============================================================================

    ! run PIKAIA/GENETIC algorithm
    subroutine run_pikaia(sim_id, all_parameters, fit_parameters)
        integer, intent(in)::sim_id
        real(dp), dimension(:), intent(in)::all_parameters
        real(dp), dimension(:)::fit_parameters ! inout, already allocated

        real(dp), dimension(:), allocatable::uniform_parameters
        real(dp)::pik_fitness

        pik_fitness = zero
        allocate (uniform_parameters(nfit))
        call ga_driver(sim_id, fpik, nfit, all_parameters, uniform_parameters,&
          &pik_fitness) ! GA DRIVER
!     call norm2par(uniform_parameters,fit_parameters,all_parameters) ! from [0-1] parameter values to physical values
        call norm2par(uniform_parameters, fit_parameters) ! from [0-1] parameter values to physical values
        ! norm2par in module 'parameters'
        deallocate (uniform_parameters)

        flush (6)

        return
    end subroutine run_pikaia

    ! ============================================================================

    ! run PARTICLE SWARM OPTIMIZATION (PSO) algorithm
    subroutine run_pso(sim_id, all_parameters, fit_parameters)
        integer, intent(in)::sim_id
        real(dp), dimension(:), intent(in)::all_parameters
        real(dp), dimension(:)::fit_parameters ! inout, already allocated

        real(dp)::pso_fitness
        pso_fitness = zero
        call pso_driver(sim_id, evaluate_pso, nfit, all_parameters,&
          &minpar, maxpar, fit_parameters, pso_fitness)
        ! minpar,maxpar in module 'parameters'

        flush (6)

        return
    end subroutine run_pso

    ! ============================================================================

end module driver

