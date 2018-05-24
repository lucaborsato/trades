module linear_ephem
  use constants
  use parameters
  use lin_fit,only:linfit
  implicit none

  interface set_ephem
    module procedure set_ephem_noinput,set_ephem_winput
  end interface set_ephem

    contains

  ! ------------------------------------------------------------------ !
  ! given the T0 data it does a linear fit to the data and it finds the
  ! ephemeris T and P: tn = Tref + Pref*n
  subroutine set_ephem_noinput()
    integer,dimension(:),allocatable::x
    real(dp),dimension(:),allocatable::y,ey
    integer::j
    character(80)::fmt

    write(*,'(a)')" COMPUTING LINEAR EPHEMERIS OF: "
    if(.not.allocated(Tephem))&
        &allocate(Tephem(NB),Pephem(NB),eTephem(NB),ePephem(NB))
    Tephem=zero
    Pephem=zero
    eTephem=zero
    ePephem=zero
    do j=2,NB
      if(nT0(j).gt.0)then
        allocate(x(nT0(j)),y(nT0(j)),ey(nT0(j)))
        x=epoT0obs(1:nT0(j),j)
        y=T0obs(1:nT0(j),j)
        ey=eT0obs(1:nT0(j),j)
        call linfit(x,y,ey,Pephem(j),ePephem(j),Tephem(j),eTephem(j))
        fmt=adjustl("(a,i3,a,4("//trim(sprec)//",a))")
        write(*,trim(fmt))" body ",j,": t_N = (",&
            &Tephem(j),"+/-",eTephem(j),") + (",&
            &Pephem(j),"+/-",ePephem(j),") x N"
        deallocate(x,y,ey)
      end if
    end do

    return
  end subroutine set_ephem_noinput
  
  subroutine set_ephem_winput(n_body,n_t0,t0_num,t0_obs,et0_obs)
    integer,intent(in)::n_body
    integer,dimension(:),intent(in)::n_t0
    integer,dimension(:,:),intent(in)::t0_num
    real(dp),dimension(:,:),intent(in)::t0_obs,et0_obs
    
    integer,dimension(:),allocatable::x
    real(dp),dimension(:),allocatable::y,ey
    integer::j

    if(.not.allocated(Tephem))&
        &allocate(Tephem(n_body),Pephem(n_body),eTephem(n_body),ePephem(n_body))
    Tephem=zero
    Pephem=zero
    eTephem=zero
    ePephem=zero
    do j=2,n_body
      if(n_T0(j).gt.0)then
        allocate(x(n_T0(j)),y(n_T0(j)),ey(n_T0(j)))
        x=t0_num(1:n_t0(j),j)
        y=t0_obs(1:n_t0(j),j)
        ey=et0_obs(1:n_t0(j),j)
        call linfit(x,y,ey,Pephem(j),ePephem(j),Tephem(j),eTephem(j))
        deallocate(x,y,ey)
      end if
    end do

    return
  end subroutine set_ephem_winput


end module linear_ephem
