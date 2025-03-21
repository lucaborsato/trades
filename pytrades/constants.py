#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import numpy as np  # array


# radiants, degrees conversions etc.
pi = 4.0 * np.arctan(1.0)
dpi = 2.0 * pi
deg2rad = pi / 180.0
rad2deg = 180.0 / pi

# various
TOLERANCE = np.finfo(float(1.0)).eps
day2sec = 86400.0  # seconds in a day  =  24h  =  86400 s
day2min = 1440.0  # min in a day  =  1440. min
day2hour = 24.0
sec2day = 1.0 / 86400.0
min2day = 1.0 / 1440.0
hour2day = 1.0 / 24.0


# masses conversions
Msmer = 6.0236e6  # Msun to Mmer
Mmers = 1.0 / Msmer  # Mmer to Msun
Msven = 4.08523719e5  # Msun to Mven
Mvens = 1.0 / Msven  #  Mven to Msun
Msear = 332946.0487  # Msun to Mear
Mears = 1.0 / Msear  #  Mear to Msun
Msmar = 3.09870359e6  # Msun to Mmar
Mmars = 1.0 / Msmar  #  Mmar to Msun
Msjup = 1.047348644e3  # Msun to Mjup
Mjups = 1.0 / Msjup  #  Mjup to Msun
Mssat = 3.4979018e3  # Msun to Msat
Msats = 1.0 / Mssat  #  Msat to Msun
Msura = 2.290298e4  # Msun to Mura
Muras = 1.0 / Msura  #  Mura to Msun
Msnep = 1.941226e4  # Msun to Mnep
Mneps = 1.0 / Msnep  #  Mnep to Msun

Mejup = Mears * Msjup  # Mear to Mjup
Mjear = Mjups * Msear  # Mjup to Mear

# masses of Solar System objects
Msun = 1.9884e30  # Sun mass in kg
Mmer = Msun * Mmers  #  Mercury mass in kg
Mven = Msun * Mvens  #  Venus mass in kg
Mear = 5.9722e24  # Earth mass in kg
Mmar = Msun * Mmars  #  Mars mass in kg
Mjup = Msun * Mjups  #  Jupiter mass in kg
Msat = Msun * Msats  #  Saturn mass in kg
Mura = Msun * Muras  #  Uranus mass in kg
Mnep = Msun * Mneps  #  Neptune mass in kg

# radii of Solar System objects
Rsun = 696000.0  #  Sun radius in km
Rmer = 2439.7  #  Mercury radius in km
Rven = 6051.8  #  Venus radius in km
Rear = 6378.1366  # Earth radius in km
Rmar = 3396.19  #  Mars radius in km
Rjup = 71492.0  #  Jupiter radius in km
Rsat = 60268.0  #  Saturn radius in km
Rura = 25559.0  #  Uranus radius in km
Rnep = 24764.0  #  Neptune radius in km
Rplu = 1195.0  #  Pluto radius in km
#
Rsjup = Rsun / Rjup  # Rsun to Rjup
Rjups = Rjup / Rsun  # Rjup to Rsun
Rsear = Rsun / Rear  # Rsun to Rear
Rears = Rear / Rsun  # Rear to Rsun
Rsmar = Rsun / Rmar  # Rsun to Rnep
Rmars = Rmar / Rsun  # Rnep to Rsun
Rsnep = Rsun / Rnep  # Rsun to Rnep
Rneps = Rnep / Rsun  # Rnep to Rsun
#
Rejup = Rear / Rjup  # Rearth to Rjupiter
Rjear = Rjup / Rear  # Rjupiter to Rearth


# astronomical constants
AU = 149597870700.0  # Astronomical Unit in meters
kappa = 0.01720209895  # Gaussian gravitational constant
Giau = kappa * kappa  # G [AU^3/Msun/d^2]
Gsi = 6.67428e-11  # Gravitational Constants in SI system [m^3/kg/s^2]
Gaumjd = Gsi * day2sec * day2sec * Mjup / (AU**3)  # G in [AU,Mjup,day]
speed = 299792458.0  # speed of light (c) in [m/s]
speedaud = speed * day2sec / AU  # speed of light in [AU/d]
pc2AU = 206264.806


# others
RsunAU = (Rsun * 1.0e3) / AU  # Sun radius in AU
RjupAU = (Rjup * 1.0e3) / AU  # Jupiter radius in AU

sphere_earth = 4.0*np.pi*((Rear*100000.0)**3)/3.0 # Rear in km to cm
rho_earth = (Mear*1000.0)/sphere_earth # gr/cm^3

MJD = 2400000.5  # MJD ref time to convert to JD
btjd = 2457000.0  # TESS BJD ref time

# others
onethird = 1.0 / 3.0
huge = np.finfo(0.0).max
