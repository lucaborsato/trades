! ********************* !
! MODULE WITH CONSTANTS !
! ********************* !

module constants
  implicit none
  ! precision parameters
  public
  integer,parameter::sp=selected_real_kind(6) ! single precision
  integer,parameter::dp=selected_real_kind(8) !define the KIND, the precision of all the constants
  integer,parameter::prec=precision(0._dp)
  character(8)::sprec

  !radiants, degrees conversions etc.
  real(dp),parameter::pi=4._dp*atan(1._dp)
  real(dp),parameter::dpi=2._dp*pi
  real(dp),parameter::deg2rad=pi/180._dp
  real(dp),parameter::rad2deg=180._dp/pi

  !various
  real(dp),parameter::s24h=86400._dp !seconds in a day = 24h = 86400s
  real(dp),parameter::zero=0._dp,one=1._dp,two=2._dp,three=3._dp
  real(dp),parameter::TOLERANCE=epsilon(zero)
  real(dp),parameter::sqrt2 = sqrt(2._dp)
  real(dp),parameter::sqrt_12 = sqrt(12._dp) ! = 2*sqrt(3) needed for stability criterion with the mutual R_Hill
  real(dp),parameter::onethird=one/three
  
  ! Constants from USNO 2013
  ! http://asa.usno.navy.mil/
  ! http://maia.usno.navy.mil/NSFA/

  !masses conversions
  real(dp),parameter::Msmer=6.0236e6_dp ! Msun to Mmer
  real(dp),parameter::Mmers=1._dp/Msmer ! Mmer to Msun
  real(dp),parameter::Msven=4.08523719e5_dp ! Msun to Mven
  real(dp),parameter::Mvens=1._dp/Msven     ! Mven to Msun
  real(dp),parameter::Msear=332946.0487_dp ! Msun to Mear
  real(dp),parameter::Mears=1._dp/Msear    ! Mear to Msun
  real(dp),parameter::Msmar=3.09870359e6_dp ! Msun to Mmar
  real(dp),parameter::Mmars=1._dp/Msmar     ! Mmar to Msun
  real(dp),parameter::Msunj=1.047348644e3_dp ! Msun to Mjup
  real(dp),parameter::Mjups=1._dp/Msunj      ! Mjup to Msun
  real(dp),parameter::Mssat=3.4979018e3_dp ! Msun to Msat
  real(dp),parameter::Msats=1._dp/Mssat    ! Msat to Msun
  real(dp),parameter::Msura=2.290298e4_dp ! Msun to Mura
  real(dp),parameter::Muras=1._dp/Msura   ! Mura to Msun
  real(dp),parameter::Msnep=1.941226e4_dp ! Msun to Mnep
  real(dp),parameter::Mneps=1._dp/Msnep   ! Mnep to Msun

  !masses of Solar System objects
  real(dp),parameter::Msun=1.9884e30_dp ! Sun mass in kg
  real(dp),parameter::Mmer=Msun*Mmers   ! Mercury mass in kg
  real(dp),parameter::Mven=Msun*Mvens   ! Venus mass in kg
  real(dp),parameter::Mear=5.9722e24_dp ! Earth mass in kg
  real(dp),parameter::Mmar=Msun*Mmars   ! Mars mass in kg
  real(dp),parameter::Mjup=Msun*Mjups   ! Jupiter mass in kg
  real(dp),parameter::Msat=Msun*Msats   ! Saturn mass in kg
  real(dp),parameter::Mura=Msun*Muras   ! Uranus mass in kg
  real(dp),parameter::Mnep=Msun*Mneps   ! Neptune mass in kg

  !radii of Solar System objects
  real(dp),parameter::Rsun=696000._dp   ! Sun radius in km
  real(dp),parameter::Rmer=2439.7_dp    ! Mercury radius in km
  real(dp),parameter::Rven=6051.8_dp    ! Venus radius in km
  real(dp),parameter::Rear=6378.1366_dp ! Earth radius in km
  real(dp),parameter::Rmar=3396.19_dp   ! Mars radius in km
  real(dp),parameter::Rjup=71492._dp    ! Jupiter radius in km
  real(dp),parameter::Rsat=60268._dp    ! Saturn radius in km
  real(dp),parameter::Rura=25559._dp    ! Uranus radius in km
  real(dp),parameter::Rnep=24764._dp    ! Neptune radius in km
  real(dp),parameter::Rplu=1195._dp     ! Pluto radius in km
  !
  ! Radius conversion
  real(dp),parameter::Rsunj=Rsun/Rjup   ! Rsun to Rjup
  real(dp),parameter::Rjups=Rjup/Rsun   ! Rjup to Rsun

  !astronomical constants
  real(dp),parameter::AU=149597870700._dp !Astronomical Unit in meters
  real(dp),parameter::kappa=0.01720209895_dp ! Gaussian gravitational constant
  real(dp),parameter::Giau=kappa*kappa ! G [AU^3/Msun/d^2]
  real(dp),parameter::Gsi=6.67428e-11_dp !Gravitational Constants in SI system [m^3/kg/s^2]
  real(dp),parameter::Gaumjd=Gsi*s24h*s24h*Mjup/(AU**3) ! G in [AU,Mjup,day]
  real(dp),parameter::speed=299792458._dp ! speed of light (c) in [m/s]
  real(dp),parameter::speedaud=speed*s24h/AU ! speed of light in [AU/d]

  !others
  real(dp),parameter::RsunAU=(Rsun*1.e3_dp)/AU !Sun radius in AU
  real(dp),parameter::RjupAU=(Rjup*1.e3_dp)/AU !Jupiter radius in AU

  ! for parallel openMP writing files
  integer,parameter::nfiles=90
  integer::ncpu
  integer,dimension(:,:),allocatable::unit2

  ! define new data type for grid search
  type parameter_grid
    character(15)::name ! name/id of the parameter
    real(dp),dimension(3)::input_values=zero ! min, max, step as read from input file
    character(2)::step_type='rn' ! input step type as read from the input file: 'ss/rn/sn'
    real(dp)::step_grid=one ! calculated step size of the parameter
    integer::n_steps=1 ! number of steps calculated, it take into account the min and max values
    real(dp),dimension(:),allocatable::grid_values ! values of the grid parameter with dimension n_steps
  end type parameter_grid

  
  
end module constants
