! ********************* !
! MODULE WITH CONSTANTS !
! ********************* !

module constants
    implicit none
    ! precision parameters
    public
    integer, parameter::sp = selected_real_kind(6) ! single precision
    integer, parameter::dp = selected_real_kind(8) !define the KIND, the precision of all the constants

    ! FOR TESTING: QUAD PRECISION
    integer, parameter::qp = selected_real_kind(32) !define the KIND, the precision of all the constants
!   integer,parameter::dp=selected_real_kind(32) !define the KIND, the precision of all the constants

    integer, parameter::prec = precision(0.0_dp)
    character(8)::sprec = 'es23.16'

    !various
    real(dp), parameter::s24h = 86400.0_dp !seconds in a day = 24h = 86400s
    real(dp), parameter::m24h = 1440.0_dp !minutes in a day = 24h = 1440m
    real(dp), parameter::zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, three = 3.0_dp
    real(dp), parameter::ten = 10.0_dp
    real(dp), parameter::TOL_sp = epsilon(0.0_sp)
    real(dp), parameter::TOL_dp = epsilon(zero)
    real(dp), parameter::TOL_qp = epsilon(0.0_qp)
    real(dp), parameter::TOLERANCE = 1.0e-9_dp
    real(dp), parameter::half = 0.5_dp
    real(dp), parameter::sqrt_half = sqrt(half)
    real(dp), parameter::sqrt2 = sqrt(two)
    real(dp), parameter::sqrt_12 = sqrt(12.0_dp) ! = 2*sqrt(3) needed for stability criterion with the mutual R_Hill
    real(dp), parameter::onethird = one/three

    !radiants, degrees conversions etc.
    real(dp), parameter::pi = 4.0_dp*atan(one)
    real(dp), parameter::dpi = two*pi
    real(dp), parameter::deg2rad = pi/180.0_dp
    real(dp), parameter::rad2deg = 180.0_dp/pi
    real(dp), parameter::circ = 360.0_dp

    ! Constants from USNO 2013
    ! http://asa.usno.navy.mil/
    ! http://maia.usno.navy.mil/NSFA/
    ! TO UPDATE WITH
    ! http://asa.hmnao.com/SecK/Constants.html

    !masses conversions
    real(dp), parameter::Msmer = 6.0236e6_dp ! Msun to Mmer
    real(dp), parameter::Mmers = one/Msmer ! Mmer to Msun
    real(dp), parameter::Msven = 4.08523719e5_dp ! Msun to Mven
    real(dp), parameter::Mvens = one/Msven     ! Mven to Msun
    real(dp), parameter::Msear = 332946.0487_dp ! Msun to Mear
    real(dp), parameter::Mears = one/Msear    ! Mear to Msun
    real(dp), parameter::Msmar = 3.09870359e6_dp ! Msun to Mmar
    real(dp), parameter::Mmars = one/Msmar     ! Mmar to Msun
    real(dp), parameter::Msunj = 1.047348644e3_dp ! Msun to Mjup
    real(dp), parameter::Mjups = one/Msunj      ! Mjup to Msun
    real(dp), parameter::Mssat = 3.4979018e3_dp ! Msun to Msat
    real(dp), parameter::Msats = one/Mssat    ! Msat to Msun
    real(dp), parameter::Msura = 2.290298e4_dp ! Msun to Mura
    real(dp), parameter::Muras = one/Msura   ! Mura to Msun
    real(dp), parameter::Msnep = 1.941226e4_dp ! Msun to Mnep
    real(dp), parameter::Mneps = one/Msnep   ! Mnep to Msun

    !masses of Solar System objects
    real(dp), parameter::Msun = 1.9884e30_dp ! Sun mass in kg
    real(dp), parameter::Mmer = Msun*Mmers   ! Mercury mass in kg
    real(dp), parameter::Mven = Msun*Mvens   ! Venus mass in kg
    real(dp), parameter::Mear = 5.9722e24_dp ! Earth mass in kg
    real(dp), parameter::Mmar = Msun*Mmars   ! Mars mass in kg
    real(dp), parameter::Mjup = Msun*Mjups   ! Jupiter mass in kg
    real(dp), parameter::Msat = Msun*Msats   ! Saturn mass in kg
    real(dp), parameter::Mura = Msun*Muras   ! Uranus mass in kg
    real(dp), parameter::Mnep = Msun*Mneps   ! Neptune mass in kg

    !radii of Solar System objects
    real(dp), parameter::Rsun = 696000.0_dp   ! Sun radius in km
    real(dp), parameter::Rmer = 2439.7_dp    ! Mercury radius in km
    real(dp), parameter::Rven = 6051.8_dp    ! Venus radius in km
    real(dp), parameter::Rear = 6378.1366_dp ! Earth radius in km
    real(dp), parameter::Rmar = 3396.19_dp   ! Mars radius in km
    real(dp), parameter::Rjup = 71492.0_dp    ! Jupiter radius in km
    real(dp), parameter::Rsat = 60268.0_dp    ! Saturn radius in km
    real(dp), parameter::Rura = 25559.0_dp    ! Uranus radius in km
    real(dp), parameter::Rnep = 24764.0_dp    ! Neptune radius in km
    real(dp), parameter::Rplu = 1195.0_dp     ! Pluto radius in km
    !
    ! Radius conversion
    real(dp), parameter::Rsunj = Rsun/Rjup   ! Rsun to Rjup
    real(dp), parameter::Rjups = Rjup/Rsun   ! Rjup to Rsun

    !astronomical constants
    real(dp), parameter::AU = 149597870700.0_dp !Astronomical Unit in meters
    real(dp), parameter::m_au = one/AU ! = 6.684587122268445e-12 au: 1 meter in au, needed as smoothing parameters in eq. of motion
    real(dp), parameter::cm_au = one/(AU*100.0_dp) ! = 6.684587122268446e-14 au: 1 cm in au, needed as smoothing parameters in eq. of motion
    real(dp), parameter::cm_au_2 = cm_au*cm_au
    real(dp), parameter::kappa = 0.01720209895_dp ! Gaussian gravitational constant
    real(dp), parameter::Giau = kappa*kappa ! G [AU^3/Msun/d^2]
    real(dp), parameter::Gsi = 6.67428e-11_dp !Gravitational Constants in SI system [m^3/kg/s^2]
    real(dp), parameter::Gaumjd = Gsi*s24h*s24h*Mjup/(AU**3) ! G in [AU,Mjup,day]
    real(dp), parameter::speed = 299792458._dp ! speed of light (c) in [m/s]
    real(dp), parameter::speedaud = speed*s24h/AU ! speed of light in [AU/d]

    !others
    real(dp), parameter::RsunAU = (Rsun*1.e3_dp)/AU !Sun radius in AU
    real(dp), parameter::RjupAU = (Rjup*1.e3_dp)/AU !Jupiter radius in AU

    ! for parallel openMP writing files
    integer, parameter::nfiles = 90
    integer::ncpu = 1
    integer, dimension(:, :), allocatable::unit2

end module constants
