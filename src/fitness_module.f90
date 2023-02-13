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

! ===============================================================================================

    subroutine asymmetric_gaussian(fit, val, nsigma, psigma, out)
        real(dp), intent(in)::fit, val, nsigma, psigma
        real(dp), intent(out)::out

        real(dp)::delta

        out = zero
        delta = fit-val
        if (delta < 0.0) then
            out = -half*((delta*delta)/(nsigma*nsigma))
        else
            out = -half*((delta*delta)/(psigma*psigma))
        end if

        return
    end subroutine asymmetric_gaussian

! ===============================================================================================

    subroutine ln_priors(allpar, fit_parameters, names, values, lnp)
        real(dp), dimension(:), intent(in)::allpar, fit_parameters
        character(10), dimension(:), intent(in)::names
        real(dp), dimension(:, :), intent(in)::values
        real(dp), intent(out)::lnp

        real(dp)::xpen, xpar
        integer::npriors, iprior
        integer::ipar
        integer::cnt_prior

        real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc
        real(dp), dimension(:), allocatable::argp, meanA, inc, longN
        logical::checkpar ! not used, only because needed in convert_parameters

        npriors = size(names)
        lnp = zero
        xpen = zero
        cnt_prior = 0

        allocate (mass(NB), radius(NB), period(NB), sma(NB), ecc(NB))
        allocate (argp(NB), meanA(NB), inc(NB), longN(NB))

        call convert_parameters(allpar, fit_parameters,&
            &mass, radius, period, sma, ecc, argp, meanA, inc, longN,&
            &checkpar)

        priors_loop: do iprior = 1, npriors
            ! FIRST CHECK IF PRIORS NAMES IN FITTING PARAMETER ID (PARID)
            do ipar = 1, nfit
                xpen = zero
                if (names(iprior) .eq. parid(ipar)) then
                    call asymmetric_gaussian(fit_parameters(ipar),&
                        &values(iprior, 1), values(iprior, 2), values(iprior, 3),&
                        &xpen)
                    lnp = lnp+xpen
                    cnt_prior = cnt_prior+1
                end if
            end do
            if (cnt_prior .ge. npriors) then
                exit priors_loop
            end if
            ! THEN CHECK IF PRIORS NAMES IN PHYSICAL PARAMETERS
            do ipar = 1, NB
                ! mass
                if (names(iprior) .eq. mass_id(ipar)) then
                    call asymmetric_gaussian(mass(ipar)*Msear,&
                        &values(iprior, 1), values(iprior, 2), values(iprior, 3),&
                        &xpen)
                    lnp = lnp+xpen
                    cnt_prior = cnt_prior+1
                end if
                ! sma
                if (names(iprior) .eq. sma_id(ipar)) then
                    call asymmetric_gaussian(sma(ipar),&
                    &values(iprior, 1), values(iprior, 2), values(iprior, 3),&
                    &xpen)
                    lnp = lnp+xpen
                    cnt_prior = cnt_prior+1
                end if
                ! ecc
                if (names(iprior) .eq. ecc_id(ipar)) then
                    call asymmetric_gaussian(ecc(ipar),&
                    &values(iprior, 1), values(iprior, 2), values(iprior, 3),&
                    &xpen)
                    lnp = lnp+xpen
                    cnt_prior = cnt_prior+1
                end if
                ! argp
                if (names(iprior) .eq. argp_id(ipar)) then
                    call asymmetric_gaussian(argp(ipar),&
                    &values(iprior, 1), values(iprior, 2), values(iprior, 3),&
                    &xpen)
                    lnp = lnp+xpen
                    cnt_prior = cnt_prior+1
                end if
                ! meana
                if (names(iprior) .eq. meana_id(ipar)) then
                    call asymmetric_gaussian(meana(ipar),&
                    &values(iprior, 1), values(iprior, 2), values(iprior, 3),&
                    &xpen)
                    lnp = lnp+xpen
                    cnt_prior = cnt_prior+1
                end if
                ! inc
                if (names(iprior) .eq. inc_id(ipar)) then
                    call asymmetric_gaussian(inc(ipar),&
                    &values(iprior, 1), values(iprior, 2), values(iprior, 3),&
                    &xpen)
                    lnp = lnp+xpen
                    cnt_prior = cnt_prior+1
                end if
                ! longn
                if (names(iprior) .eq. longn_id(ipar)) then
                    call asymmetric_gaussian(longN(ipar),&
                    &values(iprior, 1), values(iprior, 2), values(iprior, 3),&
                    &xpen)
                    lnp = lnp+xpen
                    cnt_prior = cnt_prior+1
                end if

                if (cnt_prior .ge. npriors) then
                    exit priors_loop
                end if
            end do

            ! CHECK IF JITTER IN NORMAL UNIT (FITTING LOG2JITTER)
            if (nRVset .gt. 0) then
                do ipar = 1, nRVset
                    xpar = two**fit_parameters(nkel+ipar)
                    call asymmetric_gaussian(xpar,&
                    &values(iprior, 1), values(iprior, 2), values(iprior, 3),&
                    &xpen)
                    lnp = lnp+xpen
                    cnt_prior = cnt_prior+1
                end do
            end if

        end do priors_loop

        return
    end subroutine ln_priors

    subroutine base_fitness_function(run_all_parameters, fit_parameters,&
        &chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic)
        ! Input
        real(dp), intent(in), dimension(:)::run_all_parameters
        real(dp), intent(in), dimension(:)::fit_parameters
        ! Output
        real(dp), intent(out)::chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic

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
        lnprior = zero
        if (n_priors .gt. 0) then
            call ln_priors(run_all_parameters, fit_parameters,&
                &priors_names, priors_values, lnprior)
        end if

        return
    end subroutine base_fitness_function

    subroutine bound_fitness_function(all_parameters, fit_parameters,&
        &chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic)
        ! Input
        real(dp), intent(in), dimension(:)::all_parameters
        real(dp), intent(in), dimension(:)::fit_parameters
        ! Output
        real(dp), intent(out)::chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic
        ! Local variable
        logical::check
        real(dp), dimension(:), allocatable::run_all_parameters
        integer::iflag
        logical::check_status

        iflag = 1
        check = .true.
        check_status = .true.

        check = check_only_boundaries(all_parameters, fit_parameters)
        ! write(*,*)" check_only_boundaries : ",check
        lnprior = zero

        ! write(*,*)"=== "
        if (check) then
            ! write(*,*)"=== in bound_fitness_function (check is True)"
            allocate (run_all_parameters(npar))
            run_all_parameters = all_parameters
            if (check_derived) check_status = check_derived_parameters(fit_parameters)
            if (fix_derived) call fix_derived_parameters(fit_parameters, run_all_parameters, check_status)
            if (check_status) then
                call base_fitness_function(run_all_parameters, fit_parameters,&
                    &chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic)
                ! write(*,*)"=== in bound_fitness_function (check_status is True)"
                ! write(*,*)"=== chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic"
                ! write(*,*)chi_square, reduced_chi_square, lnLikelihood, lnprior, ln_const, bic
            else ! check_status
                ! write(*,*)"=== in bound_fitness_function (check_status is False)"
                ! write(*,*)"=== set statistics to bad values"
                chi_square = resmax
                call set_fitness_values(fit_parameters, chi_square,&
                  &reduced_chi_square, lnLikelihood, ln_const, bic)
            end if
            deallocate (run_all_parameters)
        else
            ! write(*,*)"=== in bound_fitness_function (check is False)"
            ! write(*,*)"=== set statistics to bad values"
            chi_square = resmax
            call set_fitness_values(fit_parameters, chi_square,&
              &reduced_chi_square, lnLikelihood, ln_const, bic)
        end if
        ! write(*,*)"=== "
        ! flush(6)

        return
    end subroutine bound_fitness_function

end module fitness_module
