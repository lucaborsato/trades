! Module to define the fitness function that should be call by all the algorithm.
! It is needed for consistency
!
module fitness_module
    use constants
    use parameters
    use parameters_conversion
    use utils
    use ode_run, only: ode_lm
    use derived_parameters_mod
    implicit none

contains

    subroutine base_fitness_function(run_all_parameters, fit_parameters,&
        &chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)
        ! Input
        real(dp), intent(in), dimension(:)::run_all_parameters
        real(dp), intent(in), dimension(:)::fit_parameters
        ! Output
        real(dp), intent(out)::chi_square, reduced_chi_square, lnLikelihood, ln_const, bic

        real(dp), dimension(:), allocatable::resw
        integer::iflag

        integer::nobs ! == obsData%ndata == ndata
        integer::nset !, ns, ne
        nobs = obsData%ndata
        nset = obsData%obsRV%nRVset

        allocate (resw(nobs))
        resw = zero
        call ode_lm(run_all_parameters, nobs, nfit, fit_parameters, resw, iflag)

        chi_square = sum(resw*resw)
        deallocate (resw)
        call set_fitness_values(fit_parameters, chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)

        return
    end subroutine base_fitness_function

    subroutine bound_fitness_function(all_parameters, fit_parameters,&
        &chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)
        ! Input
        real(dp), intent(in), dimension(:)::all_parameters
        real(dp), intent(in), dimension(:)::fit_parameters
        ! Output
        real(dp), intent(out)::chi_square, reduced_chi_square, lnLikelihood, ln_const, bic
        ! Local variable
        logical::check
        real(dp), dimension(:), allocatable::run_all_parameters
        integer::iflag
        logical::check_status

        iflag = 1
        check = .true.
        check_status = .true.

        check = check_only_boundaries(all_parameters, fit_parameters)

        if (check) then
            allocate (run_all_parameters(npar))
            run_all_parameters = all_parameters
            if (check_derived) check_status = check_derived_parameters(fit_parameters)
            if (fix_derived) call fix_derived_parameters(fit_parameters, run_all_parameters, check_status)
            if (check_status) then
                call base_fitness_function(run_all_parameters, fit_parameters,&
                  &chi_square, reduced_chi_square, lnLikelihood, ln_const, bic)
            else ! check_status
                chi_square = resmax
                call set_fitness_values(fit_parameters, chi_square,&
                  &reduced_chi_square, lnLikelihood, ln_const, bic)
            end if
            deallocate (run_all_parameters)

        else
            chi_square = resmax
            call set_fitness_values(fit_parameters, chi_square,&
              &reduced_chi_square, lnLikelihood, ln_const, bic)
        end if

        return
    end subroutine bound_fitness_function

end module fitness_module
