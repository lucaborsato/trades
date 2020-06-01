module custom_type
  use constants,only:dp,zero,one
  implicit none

! ==============================================================================

  ! OBS DATA AS STRUCTURE
  type dataRV
    sequence
    integer::nRVset
    ! nRV will be used as the number of RV datapoints for the observed one,
    ! and as a counter for the simulated one
    integer::nRV=0
    real(dp),dimension(:),allocatable::jd,RV,eRV,trend
    integer,dimension(:),allocatable::RV_stat
    real(dp),dimension(:,:),allocatable::gamma
    integer,dimension(:),allocatable::RVsetID,nRVsingle
  end type dataRV
  
! ==============================================================================

  ! T0 data for one body!
  type dataT0
    sequence
    ! nT0 will be used as the number of TTs for the observed one,
    ! and as a counter for the simulated one
    integer::nT0=0,nDur=0
    integer,dimension(:),allocatable::epo
    real(dp),dimension(:),allocatable::T0,eT0,dur,edur,oc
    integer,dimension(:),allocatable::T0_stat,dur_stat
    real(dp)::Tephem=zero,eTephem=zero,Pephem=zero,ePephem=zero
  end type dataT0

! ==============================================================================

  ! full data type
  type dataObs
    sequence ! needed to store variables contiguous
    integer::ndata,nfree,dof=1
    real(dp)::inv_dof
    ! RV data
    type(dataRV)::obsRV
    ! T0 and duration data
    ! and as a counter for the simulated ones
    integer::nTTs=0,nDurs=0 ! default set to zero
    type(dataT0),dimension(:),allocatable::obsT0
  end type dataObs

! ==============================================================================

  ! define new data type for grid search
  type parameter_grid
    sequence
    character(15)::name ! name/id of the parameter
    real(dp),dimension(3)::input_values=zero ! min, max, step as read from input file
    character(2)::step_type='rn' ! input step type as read from the input file: 'ss/rn/sn'
    real(dp)::step_grid=one ! calculated step size of the parameter
    integer::n_steps=1 ! number of steps calculated, it take into account the min and max values
    real(dp),dimension(:),allocatable::grid_values ! values of the grid parameter with dimension n_steps
  end type parameter_grid
 
! ==============================================================================
 
  contains
  
! ==============================================================================

  subroutine init_dataRV(nRV,RV)
    integer,intent(in)::nRV
    type(dataRV),intent(inout)::RV
    
    RV%nRV=nRV
    allocate(RV%jd(nRV),RV%RV(nRV),RV%eRV(nRV),RV%trend(nRV))
    allocate(RV%RV_stat(nRV),RV%RVsetID(nRV))
    RV%jd=zero
    RV%RV=zero
    RV%eRV=zero
    RV%trend=zero
    RV%RV_stat=0
    RV%RVsetID=0
    
    !!
    ! nRVset & RVsetID & nRVsingle HAVE TO BE SET BY 'HAND'
    !!
  
    return
  end subroutine init_dataRV
  
! ==============================================================================

  subroutine deallocate_dataRV(RV)
    type(dataRV),intent(inout)::RV
    
    if(allocated(RV%jd)) deallocate(RV%jd)
    if(allocated(RV%RV)) deallocate(RV%RV)
    if(allocated(RV%eRV)) deallocate(RV%eRV)
    if(allocated(RV%trend)) deallocate(RV%trend)
    if(allocated(RV%RV_stat)) deallocate(RV%RV_stat)
    if(allocated(RV%gamma)) deallocate(RV%gamma)
    if(allocated(RV%RVsetID)) deallocate(RV%RVsetID)
    if(allocated(RV%nRVsingle)) deallocate(RV%nRVsingle)
    RV%nRV=0
    
    return
  end subroutine deallocate_dataRV

! ==============================================================================
  
  subroutine init_dataT0(nT0,T0,dur_check)
    integer,intent(in)::nT0
    type(dataT0),intent(inout)::T0
    integer,intent(in)::dur_check
    
    T0%nT0=nT0
    allocate(T0%epo(nT0),T0%T0(nT0),T0%eT0(nT0),T0%oc(nT0),T0%T0_stat(nT0))
    T0%epo=0
    T0%T0=zero
    T0%eT0=zero
    T0%oc=zero
    T0%T0_stat=0
    ! duration
    allocate(T0%dur(nT0),T0%edur(nT0),T0%dur_stat(nT0))
    T0%dur=zero
    T0%edur=zero
    T0%dur_stat=0
    if(dur_check.eq.1) T0%nDur=nT0

    return
  end subroutine init_dataT0
  
! ==============================================================================

  subroutine deallocate_dataT0(T0)
    type(dataT0),intent(inout)::T0
    
    if(allocated(T0%epo)) deallocate(T0%epo)
    if(allocated(T0%T0)) deallocate(T0%T0)
    if(allocated(T0%eT0)) deallocate(T0%eT0)
    if(allocated(T0%oc)) deallocate(T0%oc)
    if(allocated(T0%T0_stat)) deallocate(T0%T0_stat)
    if(allocated(T0%dur)) deallocate(T0%dur)
    if(allocated(T0%edur)) deallocate(T0%edur)
    if(allocated(T0%dur_stat)) deallocate(T0%dur_stat)
    if(allocated(T0%oc)) deallocate(T0%oc)
    T0%nT0=0
    T0%nDur=0
    T0%Tephem=zero
    T0%eTephem=zero
    T0%Pephem=zero
    T0%ePephem=zero
  
    return
  end subroutine deallocate_dataT0
  
! ==============================================================================

  subroutine deallocate_dataObs(oData)
    type(dataObs),intent(inout)::oData
    
    integer::n,i
    
    call deallocate_dataRV(oData%obsRV)
    n=size(oData%obsT0)
    if(oData%ntts.gt.0)then
      do i=1,n
        call deallocate_dataT0(oData%obsT0(i))
      end do
    end if
  
  end subroutine deallocate_dataObs

! ==============================================================================


end module custom_type
