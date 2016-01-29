module derived_parameters_mod
  use constants
  use parameters
  use init_trades,only:get_unit,get_rows
  implicit none
  
  ! modify the init_derived_parameters, the compute_derived_parameters and the check_derived_parameters
  ! subroutines...
  ! or you can rename them as you want, adapt yours and then use interface so the code will call the same
  ! subroutine/function names without problems
  ! please use the same variables declared in parameters.f90 module in the the field ! derived parameters
  
  interface init_derived_parameters
    module procedure init_derived_phase
!     module procedure init_derived_ttvmax ! use the phase one
  end interface init_derived_parameters
  
  interface compute_derived_parameters
!     module procedure compute_derived_phase
    module procedure compute_derived_ttvmax
  end interface compute_derived_parameters
  
  interface check_derived_parameters
!     module procedure check_derived_phase
    module procedure check_derived_ttvmax
  end interface check_derived_parameters
  

contains

!! ---
!! --- RV phase by Luca 'Ranger' Malavolta ---
!! ---

  ! check if exist derived_boundaries.dat file, if True it initialises the
  ! min,max values of the derived parameters, otherwise it set to False flag
  ! derived_boundaries.dat has to be in this kind of form:
  !# name derived_min derived_max
  !phi2 30. 270.
  !phi3 40. 95.
  subroutine init_derived_phase(cpuid,path_in)
    integer,intent(in)::cpuid
    character(512),intent(in)::path_in
    character(512)::read_file,row
    integer::uread,irow,io_stat
    logical::fstat
    
!     write(*,*)trim(path_in)
    read_file=trim(adjustl(trim(path_in))//"derived_boundaries.dat")
!     write(*,*)trim(read_file)
    inquire(file=trim(read_file),exist=fstat)
    if(fstat)then
      uread=get_unit(cpuid)
      open(uread,file=read_file,status='OLD')
      n_derived=get_rows(uread)
!       write(*,*)n_derived
      allocate(derived_names(n_derived), derived_boundaries(2,n_derived))
!       write(*,*)'allocated derived_names and derived_boundaries'
      irow=0
      readdo: do
        read(uread,'(a512)',IOSTAT=io_stat)row
        if(IS_IOSTAT_END(io_stat)) exit readdo
        row=trim(adjustl(row))
        if(row(1:1).ne."#")then
          irow=irow+1
          read(row,*)derived_names(irow),derived_boundaries(1,irow),derived_boundaries(2,irow)
          write(*,*)derived_names(irow),derived_boundaries(1,irow),derived_boundaries(2,irow)
        end if
      end do readdo
      close(uread)
      check_derived = .true.
    else
      n_derived=0
      check_derived = .false.
    end if
!     write(*,*)' check_derived =',check_derived
  
    return
  end subroutine init_derived_phase

  ! SUBROUTINE TO MODIFY AT YOUR NEEDS!
  ! in my case:
  ! phi_pl = mean_anomaly_pl - argument_pericenter_pl
  subroutine compute_derived_phase(fitting_parameters,derived_parameters)
    real(dp),dimension(:),intent(in)::fitting_parameters
    real(dp),dimension(:),intent(out),allocatable::derived_parameters
    real(dp)::temp_par
    integer::ii,i_der
    ! global variables: parid,nfit
    
    allocate(derived_parameters(n_derived))
    derived_parameters=zero
    i_der=0
    do ii=1,nfit
      if(parid(ii)(1:5).eq.'ecosw')then
!         write(*,*)trim(parid(ii))
        if(parid(ii+1)(1:5).eq.'esinw')then
!           write(*,*)trim(parid(ii+1))
          temp_par=zero
          temp_par=mod(atan2(fitting_parameters(ii+1),fitting_parameters(ii))*rad2deg+360._dp,360._dp)
          if(parid(ii+2)(1:2).eq.'mA')then
            i_der=i_der+1
!             write(*,*)trim(parid(ii+2)),i_der
!             derived_parameters(i_der)=mod(fitting_parameters(ii+2)-temp_par+360._dp,360._dp) ! WRONG
            derived_parameters(i_der)=mod(fitting_parameters(ii+2)+temp_par+360._dp,360._dp)
           end if
        end if
      end if
    end do
  
    return
  end subroutine compute_derived_phase
  
  ! check if each derived parameter is within the provided boundaries
  ! if at least one is not within the range it return a False logical
  function check_derived_phase(fitting_parameters) result(check_status)
    logical::check_status
    real(dp),dimension(:),intent(in)::fitting_parameters
    integer::ii
    real(dp),dimension(:),allocatable::derived_parameters
    
    check_status=.true.
    call compute_derived_parameters(fitting_parameters,derived_parameters)
    check_der:do ii=1,n_derived
      if(derived_parameters(ii).ge.derived_boundaries(1,ii))then ! derived_par >= derived_min
        if(derived_parameters(ii).le.derived_boundaries(2,ii))then ! derived_par <= derived_max
          check_status=.true.
        else
          check_status=.false.
          exit check_der
        end if
      else
        check_status=.false.
        exit check_der
      end if
    end do check_der
    if(allocated(derived_parameters)) deallocate(derived_parameters)
  
    return
  end function check_derived_phase

  
! -------------------------------------------

!! ---
!! --- TTVmax by Luigi Mancini ---
!! ---
  subroutine compute_derived_ttvmax(fitting_parameters,derived_parameters)
    use celestial_mechanics,only:par2kel_fit
    real(dp),dimension(:),intent(in)::fitting_parameters
    real(dp),dimension(:),intent(out),allocatable::derived_parameters
    real(dp),dimension(:),allocatable::mass,Rad,Per,sma,ecc,argper,meanAnom,inc,lNode
    logical::checkpar
    
    
!     TTVmax =  \mu*e3*(a2/a3)^3*P3;
!     dove \mu è la massa ridotta (m3/Mstar), e3 è l'eccentricità del
!     pianeta c, a2 è il semi-asse maggiore del pianeta b, a3 è il semi-asse
!     maggiore del pianeta c, P3 è il periodo orbitale del pianeta c.
    allocate(derived_parameters(n_derived))
    derived_parameters=zero
    allocate(mass(NB),Rad(NB),Per(NB),sma(NB),ecc(NB),argper(NB),meanAnom(NB),inc(NB),lNode(NB))
    checkpar=.true.
    call par2kel_fit(system_parameters,fitting_parameters,mass,Rad,Per,sma,ecc,argper,meanAnom,inc,lNode,checkpar)
    if(checkpar)then
      derived_parameters = (mass(3)/mass(1))*ecc(3)*Per(3)*((sma(2)/sma(3))**3)*1440.
    else
      derived_parameters=-one
    end if
    deallocate(mass,Rad,Per,sma,ecc,argper,meanAnom,inc,lNode)
  
    return
  end subroutine compute_derived_ttvmax
  
    
  function check_derived_ttvmax(fitting_parameters) result(check_status)
    logical::check_status
    real(dp),dimension(:),intent(in)::fitting_parameters
    integer::ii
    real(dp),dimension(:),allocatable::derived_parameters
  
    check_status=.true.
    call compute_derived_parameters(fitting_parameters,derived_parameters)
    check_der:do ii=1,n_derived
      if(derived_parameters(ii).ge.derived_boundaries(1,ii))then ! derived_par >= derived_min
        if(derived_parameters(ii).le.derived_boundaries(2,ii))then ! derived_par <= derived_max
          check_status=.true.
        else
          check_status=.false.
          exit check_der
        end if
      else
        check_status=.false.
        exit check_der
      end if
    end do check_der
    if(allocated(derived_parameters)) deallocate(derived_parameters)
  
    return
  end function check_derived_ttvmax
end module derived_parameters_mod

