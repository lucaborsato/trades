! Module to define the fitness function that should be call by all the algorithm.
! It is needed for consistency
!
module fitness_module
  use constants
  use parameters
  use parameters_conversion
  use ode_run,only:ode_lm
  use derived_parameters_mod
  implicit none

contains

  function base_fitness_function(run_all_parameters,fit_parameters) result(fitness)
    real(dp)::fitness
    real(dp),intent(in),dimension(:)::run_all_parameters
    real(dp),intent(in),dimension(:)::fit_parameters
    real(dp),dimension(:),allocatable::resw
    integer::iflag
  
    allocate(resw(ndata))
    resw=zero
    call ode_lm(run_all_parameters,ndata,nfit,fit_parameters,resw,iflag)
    fitness=sum(resw*resw)
    ! resw t.c. sum(resw^2) = fitness = Chi2r*K_chi2r + Chi2wr*K_chi2wr
    deallocate(resw)
    if(fitness.ge.resmax)then
      fitness=resmax
    end if
  
    return
  end function base_fitness_function

  function check_fit_boundaries(fit_parameters) result(check)
    logical::check
    real(dp),intent(in),dimension(:)::fit_parameters
    integer::i
    check=.true.
  
    checkloop: do i=1,nfit
      if(fit_parameters(i).lt.minpar(i))then
        check=.false.
        exit checkloop
      else if(fit_parameters(i).gt.maxpar(i))then
        check=.false.
        exit checkloop
      end if
    end do checkloop
  
    return
  end function check_fit_boundaries
  
  function bound_fitness_function(all_parameters,fit_parameters) result(fitness)
    real(dp)::fitness
    real(dp),intent(in),dimension(:)::all_parameters
    real(dp),intent(in),dimension(:)::fit_parameters
    logical::check
    real(dp),dimension(:),allocatable::run_all_parameters
!     real(dp),dimension(:),allocatable::resw
    integer::iflag
    logical::check_status
  
    iflag=1
    check=.true.
    check_status=.true.
    fitness=zero
    
    check=check_fit_boundaries(fit_parameters)

!     write(*,'(i2,2(L2,1x),es23.16)')0,check_status,check,fitness
    
    if(check)then

      allocate(run_all_parameters(npar))
      run_all_parameters=all_parameters
      if(check_derived) check_status=check_derived_parameters(fit_parameters)
      if(fix_derived) call fix_derived_parameters(fit_parameters,run_all_parameters,check_status)

      if(check_status)then
        fitness=base_fitness_function(run_all_parameters,fit_parameters)
      else ! check_status
        fitness=resmax
      end if
      deallocate(run_all_parameters)
      
    else
      fitness=resmax
    end if
  
!     write(*,'(i2,2(L2,1x),es23.16)')1,check_status,check,fitness
!     flush(6)
  
    return
  end function bound_fitness_function

 
end module fitness_module
