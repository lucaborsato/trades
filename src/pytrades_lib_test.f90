module f90trades
    use constants
    ! use parameters
    

    implicit none

    !! === WARNING === !!
    !! beware cannot be exposed non-parameter variables!!
    !! === WARNING === !!

    ! exposing variables in parameters module to trades_lib
    !f2py integer,parameter::dp
    !f2py real(dp), parameter::one

    !f2py real(dp), parameter::deg2rad
    !f2py real(dp), parameter::rad2deg


contains

    subroutine change_ncpu(new_ncpu)
        integer, intent(in)::new_ncpu

        ncpu = new_ncpu

        return
    end subroutine change_ncpu

    ! ============================================================================

    ! SET ARGS AND PARAMETERS WITH SINGLE SUBROUTINES TO BE CALLED BY PYTHON IF NEEDED

    ! --- subroutine useful to modify the working path of TRADES from python
    ! subroutine path_change(new_path)
    !     character(512), intent(in)::new_path

    !     path = trim(adjustl(new_path))

    !     return
    ! end subroutine path_change

    ! subroutine get_path(path_out)
    !     character(512), intent(out)::path_out

    !     path_out = trim(adjustl(path))

    !     return
    ! end subroutine get_path
!     ! ---

!     ! --- set number of bodies
!     subroutine set_n_bodies(n_body)
!         integer, intent(in)::n_body

!         NB = n_body
!         NBDIM = NB*6

!         return
!     end subroutine set_n_bodies

!     subroutine get_n_bodies(n_body)
!         integer, intent(out)::n_body

!         n_body = NB

!         return
!     end subroutine get_n_bodies
!     ! ---

!     ! --- set epcoh/reference time
!     subroutine set_epoch_time(t_epoch)
!         real(dp), intent(in)::t_epoch

!         tepoch = t_epoch

!         return
!     end subroutine set_epoch_time

!     ! --- set starting time
!     subroutine set_starting_time(t_starting)
!         real(dp), intent(in)::t_starting

!         tstart = t_starting

!         return
!     end subroutine set_starting_time

!     ! --- set integration time
!     subroutine set_integration_time(t_integration)
!         real(dp), intent(in)::t_integration

!         tint = t_integration

!         return
!     end subroutine set_integration_time
!     ! ============================================================================

    ! ! ============================================================================
    ! ! ============================================================================
    ! subroutine initialize_trades(path_in, sub_folder, n_threads_in)
    !     !f2py real(dp),dimension(:),allocatable::eRVobs
    !     !f2py real(dp),dimension(:,:),allocatable::eT0obs
         
    !     character*(*), intent(in)::path_in, sub_folder
    !     integer, intent(in),optional::n_threads_in

    !     ! Local
    !     real(dp), dimension(:), allocatable::mass, radius, period, sma, ecc, argp, meanA, inc, longN
    !     integer, dimension(:), allocatable::nset

    !     integer:: ipar

    !     write (*, *) "INITIALISING TRADES ..."

    !     call initu(nfiles, 1)

    !     ! IT READS THE COMMAND ARGUMENT THAT DEFINE THE PATH OF THE FOLDER WITH THE FILES
    !     path_0 = trim(adjustl(path_in))
    !     path = trim(adjustl(path_in))

    !     ! IT DEFINES THE STRING TO WRITE THE REAL WITH RIGHT DECIMAL: PRECISION
    !     sprec = 'es23.16'

    !     ! IT READS THE ARGUMENTS OF INTEGRATION AND STORE IN COMMON MODULE PARAMETERS.
    !     ! THE VARIBLES WILL NOT BE MODIFIED FURTHERMORE.
    !     call read_arg(1)
    !     write(*,'(a, f16.5)')"t_epoch = ",tepoch
    !     write(*,'(a, f16.5)')"t_start = ",tstart
    !     write(*,'(a, f16.5)')"t_int   = ",tint
    !     flush(6)

    !     n_bodies = NB ! needed to be used by python wrapper ... to check if I can avoid it
    !     if (allocated(e_bounds)) deallocate (e_bounds)
    !     allocate (e_bounds(2, NB))
    !     e_bounds(1, :) = TOL_dp
    !     e_bounds(2, :) = 1.0_dp-TOL_dp

    !     ! IT READS THE FILES AND THE NAMES OF THE BODIES AND DETERMINES THE PARAMETERS TO BE FITTED
    !     call read_list(1)
    !     flush(6)

    !     ! IT READS RV DATA
    !     nRV = 0
    !     call read_RVobs(1)
    !     flush(6)
    !     nRV = obsData%obsRV%nRV

    !     ! IT READS T0 DATA
    !     if (.not. allocated(obsData%obsT0)) allocate (obsData%obsT0(NB-1))
    !     call read_T0obs(1)
    !     if (idtra .ne. 0) then
    !         write (*, '(a,1000(i5,1x))') " T0 DATA: nT0 = ",&
    !             &obsData%obsT0(:)%nT0
    !         if (durcheck .eq. 1) write (*, '(a,1000(i5,1x))') " DUR DATA: nT0 = ",&
    !             &obsData%obsT0(:)%nT0
    !         write (*, '(a,1000(l2,1x))') ' DO TRANSIT FLAG: ', do_transit
    !     end if
    !     flush(6)

    !     ! IT SETS THE LINEAR EPHEMERIS FROM T0 DATA
    !     if (obsData%nTTs .gt. 0) then
    !         call set_ephem()
    !         call compute_oc(obsData%obsT0)
    !     end if
    !     flush(6)

    !     ! IT DETERMINES THE NDATA
    !     nTTs = obsData%nTTs
    !     nDurs = obsData%nDurs

    !     obsData%ndata = nRV+nTTs+nDurs
    !     obsData%nfree = 0 ! gamma now fitted

    !     ndata = obsData%ndata
    !     nfree = obsData%nfree

    !     write (*, '(a,a)') " NUMBER OF ORBITAL PARAMETERS TO FIT: nfit = ", trim(adjustl(string(nfit)))
    !     nfit = nfit+rv_trend_order+nRVset+nRVset ! 2 x nRVset: gamma + jitter for each RV dataset
    !     if (rv_trend_order .gt. 0) then
    !         write (*, '(a,i2)') " RV trend of order ", rv_trend_order
    !     end if
    !     if (nRVset .gt. 0) then
    !         write (*, '(a,i2)') " number RV dataset ", nRVset
    !     end if
    !     write (*, '(a,a)') " NUMBER OF PARAMETERS TO FIT: nfit = ", trim(adjustl(string(nfit)))
    !     obsData%dof = (obsData%ndata-nfit-obsData%nfree)
    !     dof = obsData%dof
    !     flush(6)

    !     if (dof .le. 0) then
    !         write (*, '(a,a)') ' FOUND dof <= 0 SO IT IS FORCED TO 1 IN CASE',&
    !         &' THE USER WANT TO SIMULATE/INTEGRATE AND NOT CHECK THE FIT.'
    !         obsData%dof = 1
    !         dof = obsData%dof
    !     end if
    !     flush(6)

    !     obsData%inv_dof = one/real(obsData%dof, dp)
    !     inv_dof = obsData%inv_dof

    !     if (obsData%ndata .gt. 0) then
    !         bic_const = real(nfit, dp)*log(real(obsData%ndata, dp))
    !     else
    !         bic_const = zero
    !     end if

    !     ! IT DEFINES THE ID OF THE PARAMETERS TO BE FITTED
    !     call idpar()
    !     flush(6)

    !     ! IT READS BOUNDARIES OF THE KEPLERIAN ELEMENTS
    !     call read_fullpar(1, mass, radius, period, sma, ecc, argp, meanA, inc, longN, system_parameters)
    !     flush(6)

    !     ! IT DETERMINES THE LN_ERR_CONST TO COMPUTE LOGLIKELIHOOD
    !     ln_err_const = get_lnec(obsData)

    !     ! IT SETS THE LIST OF THE PARAMETERS TO FIT
    !     call set_parid_list()
    !     ! IT SETS FITNESS PARAMETERS
    !     if (nRV .ne. 0 .and. nTTs .ne. 0) then
    !         if (durcheck .eq. 0) then
    !             allocate (nset(2))
    !             nset(1) = nRV
    !             nset(2) = nTTs
    !         else
    !             allocate (nset(3))
    !             nset(1) = nRV
    !             nset(2) = nTTs
    !             nset(3) = nDurs
    !         end if
    !     else if (nRV .ne. 0 .and. nTTs .eq. 0) then
    !         allocate (nset(1))
    !         nset(1) = nRV
    !     else if (nRV .eq. 0 .and. nTTs .ne. 0) then
    !         if (durcheck .eq. 0) then
    !             allocate (nset(1))
    !             nset(1) = nTTs
    !         else
    !             allocate (nset(2))
    !             nset(1) = nTTs
    !             nset(2) = nDurs
    !         end if
    !     else
    !         allocate (nset(1))
    !         nset(1) = 1
    !     end if
    !     deallocate (nset)
    !     flush(6)
        

    !     write (*, '(a)') ''
    !     write (*, '(a)') 'Initial-input Keplerian orbital elements: val'
    !     write (*, '(a, 1000(1x,es23.16))') "mass     [Msun] = ", mass
    !     write (*, '(a, 1000(1x,es23.16))') "radius   [Rsun] = ", radius
    !     write (*, '(a, 1000(1x,es23.16))') "period   [days] = ", period
    !     write (*, '(a, 1000(1x,es23.16))') "sma      [au]   = ", sma
    !     write (*, '(a, 1000(1x,es23.16))') "ecc             = ", ecc
    !     write (*, '(a, 1000(1x,es23.16))') "argp     [deg]  = ", argp
    !     write (*, '(a, 1000(1x,es23.16))') "meana    [deg]  = ", meana
    !     write (*, '(a, 1000(1x,es23.16))') "inc      [deg]  = ", inc
    !     write (*, '(a, 1000(1x,es23.16))') "longn    [deg]  = ", longn
    !     write (*, '(a)') ''
    !     flush(6)

    !     write (*, '(a23,3(1x,a23))') 'Full System Parameters', "value", '( par_min ,', 'par_max )'
    !     do ipar = 1, npar
    !         write (*, '(a23,1x,es23.16,2(a,es23.16),a)') all_names_list(ipar), system_parameters(ipar),&
    !             &' ( ', par_min(ipar), ' , ', par_max(ipar), ' )'
    !     end do
    !     flush(6)

    !     ! read priors within fortran!
    !     call read_priors(1)
        
    !     str_len = len(parid(1))

    !     ! check if there are derived parameters to compute and to check
    !     call init_derived_parameters(1, path)

    !     ! deallocated variables not needed anymore
    !     if (allocated(mass)) deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)

    !     path = trim(adjustl(path_in))//trim(adjustl(sub_folder))

    !     return
    ! end subroutine initialize_trades
    ! ! ============================================================================

    ! subroutine get_system_parameters_and_minmax(npar, sys_params, sys_par_min, sys_par_max)
    !     ! Input
    !     integer,intent(in)::npar

    !     ! Output
    !     real(dp),dimension(npar),intent(out)::sys_params, sys_par_min, sys_par_max

    !     sys_params = system_parameters
    !     sys_par_min = par_min
    !     sys_par_max = par_max
    
    !     return
    ! end subroutine get_system_parameters_and_minmax

    ! ! ============================================================================

    ! subroutine get_default_fitting_parameters_and_minmax(fit_params, fit_par_minmax, n_fit)
    !     ! Input
    !     integer,intent(in)::n_fit

    !     ! Output
    !     real(dp),dimension(n_fit),intent(out)::fit_params
    !     real(dp),dimension(n_fit, 2),intent(out)::fit_par_minmax
        
    !     ! Local
    !     real(dp),dimension(:),allocatable::fitting_parameters

    !     call init_param(system_parameters, fitting_parameters)
    !     fit_params = fitting_parameters
    !     deallocate(fitting_parameters)
    !     fit_par_minmax(:, 1) = minpar
    !     fit_par_minmax(:, 2) = maxpar
    
    !     return
    ! end subroutine get_default_fitting_parameters_and_minmax

end module f90trades
