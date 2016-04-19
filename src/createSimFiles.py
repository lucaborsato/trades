#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os #  os: operating system
import glob # glob: globbing file...loading multiple files as *.pippa
import sys # sys: system...boh
import numpy as np # array
#sys.path.append("/media/Data/Dropbox/Dottorato/pyscript")
#sys.path.append("../")
#from path_check import checkpath
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# THIS SCRIPT CREATES THE FILES NEEDED BY TRADES TO RUN
# IN A SPECIFIC FOLDER
#
# THE STANDARD WILL BE KEPLER-9 HOLMAN 2010 (DISCOVERY PAPER)
# 3 BODIES, STAR+B+C + 6RV
#
# --------------------------------------------------------------------------

def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-p', action='store', dest='fpath', required=True, help='Folder path')
  cli = parser.parse_args()
  return os.path.abspath(cli.fpath)

# check if directory does not exist, then creates iteration
def createFolder(fpath):
  if not os.path.isdir(fpath):
    print " CREATING FOLDER " + fpath
    os.makedirs(fpath)

#
# STANDARD arg.in FILE

def createARG(fpath):
  farg = os.path.join(fpath, "arg.in")
  print " CREATING: " + farg
  print ""

  ofarg = open(farg, 'w')

  linel = "progtype = 0" # progtype
  liner = "# 0=only integration, 1=grid search, 2=Levenberg-Marquardt, 3=PIKAIA (GA), 4=PSO, 5=PolyChord."
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "nboot = 0" # nboot
  liner = "# bootstrap: <=0 no bootstrap, >0 yes bootstrap (Nboot set to 100 if <100)"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line
  
  linel = "bootstrap_scaling = T" # nboot
  liner = "# bootstrap_scaling = .true. or T, .false. or F, default .false."
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "tepoch = 2455088.212" # epoch
  liner = "# epoch of the elements [JD]."
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "tstart = 2454965." # time start
  liner = "# time start of the integration [JD]. Set this to a value > 9.e7 if you want to use the tepoch (next parameter) "
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "tint = 500." # integration time
  liner = "# time duration of the integration in days"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "step = 1.e-3" # initial step size
  liner = "# initial time step size in days"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "wrttime = 0.04167" # write interval time
  liner = "# time interval in days of write data in files (if < stepsize the program will write every step)"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "NB = 3" # Number of bodies
  liner = "# number of bodies to use (from 2 to N, where the max N is the number of files in bodies.lst)"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "idtra = 1" # T_0 fit, duration fit
  liner = "# number of the body to check if it transits (from 2 to N, 1 for everyone, 0 for no check)."
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "durcheck = 0" # T_0 fit, duration fit
  liner = "# 0/1=no/yes duration fit"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print " " + line

  linel = "tol_int = 1.e-13" # TOLERANCE
  liner = "# tolerance in the integration (stepsize selection etc)"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "wrtorb = 1" # write orbit
  liner = "# write orbit condition: 1[write], 0[do not write]"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "wrtconst = 1" # write const of motion
  liner = "# write constants condition: 1[write], 0[do not write]"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "wrtel = 1" # write orb elements
  liner = "# write orbital elements condition: 1[write], 0[do not write]"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "rvcheck = 1" # check RV
  liner = "# check of Radial Velocities condition: 1[check, read from obsRV.dat], 0[do not check]"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "idpert = 3" # grid
  liner = "# grid option: id of the perturber body [integer >1, <= tot bodies; else no perturber]"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  linel = "lmon = 0" # lmon
  liner = "# lmon: Levenberg-Marquardt off = 0[no LM], on=1[yes LM]. Default lmon = 0"
  line = liner + "\n" + linel
  ofarg.write(line + "\n")
  print line

  line = "# weight_chi_square: parameter needed to weight the fitness function, such that\n# fitness = Chi^2_r * weight_chi_square + Chi^2_wr * (1. - weight_chi_square),\n# where Chi^2_wr is the reduced Chi Square but weighted by the number of data for each data set (Chi2r * Ndata/Nset for all dataset).\n# Default weight_chi_square = 1.\nweight_chi_square = 1."
  ofarg.write(line + "\n")
  print line

  line = "# secondary_parameters: define if the program has to check only boundarie for derived parameters (1) or if it has also to fix values due to derived parameters (2) or do nothing (0) with derived parameters.\n# Default secondary_parameters = 0.\nsecondary_parameters = 0"
  ofarg.write(line + "\n")
  print line
  ofarg.close()

# CREATES bodies.lst file
# star.dat 0 0               #filename Mass Radius [1=to fit, 0=fixed]
# b.dat 0 0 1 1 1 1 0 0  #filename Mass Radius Period eccentricity Arg.ofPericenter MeanAnomaly inclination Long.ofNodes
# c.dat 0 0 1 1 1 1 0 1  #filename Mass Radius Period eccentricity Arg.ofPericenter MeanAnomaly inclination Long.ofNodes
def createBODIESLST(fpath):
  fbodies = os.path.join(fpath, "bodies.lst")
  print " CREATING: " + fbodies
  print ""

  ofbodies = open(fbodies, 'w')

  linel = "star.dat 0 0"
  liner = "#filename Mass Radius [1=to fit, 0=fixed]"
  line = linel + "          " + liner
  ofbodies.write(line + "\n")
  print " " + line

  linel = "b.dat 0 0 1 1 1 1 0 0"
  liner = "#filename Mass Radius Period eccentricity Arg.ofPericenter MeanAnomaly inclination Long.ofNodes"
  line = linel + "          " + liner
  ofbodies.write(line + "\n")
  print " " + line

  linel = "c.dat 0 0 1 1 1 1 0 1"
  liner = "#filename Mass Radius Period eccentricity Arg.ofPericenter MeanAnomaly inclination Long.ofNodes"
  line = linel + "          " + liner
  ofbodies.write(line + "\n")
  print " " + line

  ofbodies.close()

# CREATE star.dat
# 1.00    # Mstar [Msun]
# 1.10    # Rstar [Rsun]
def createSTAR(fpath):
  fstar = os.path.join(fpath, "star.dat")
  print " CREATING: " + fstar
  print ""

  ofstar = open(fstar, 'w')

  linel = "1.00"
  liner = "# Mstar [Msun]"
  line = linel + "          " + liner
  ofstar.write(line + "\n")
  print " " + line

  linel = "1.10"
  liner = "# Rstar [Rsun]"
  line = linel + "          " + liner
  ofstar.write(line + "\n")
  print " " + line
  
  ofstar.close()
  
# CREATES PLANETS FILES: b.dat
# 0.136985 0. 0.01 ss          # Mmin Mmax X ss/rn/sn=StepSize/RandomNumber/StepNumber [Mjup]
# 0.842 0. 0. rn               #Radius of Planet [Rjup]
# 19.238756 0. 0. rn           #Period [day] - set it > 9000000. to calculate from semi-major axis
# 999. 0. 0. rn                #semi major axis [AU] - set it to 999.0 if there is no value available
# 0.058152 0. 0. rn            #eccentricity
# 356.056648 0. 0. rn          #argument of the pericenter [deg]
# 3.771733 0. 0. rn            #mean anomaly [deg]
# 9.e8 0. 0. rn                # time pericenter[JD]
# 88.55 0. 0. rn               #orbit inclination [deg] #88.55 !!TASK8 88.2781373
# 0. 0. 0. rn                  #longitude of the ascending node [deg]
  
def createPLANETB(fpath):
  fb = os.path.join(fpath, "b.dat")
  print " CREATING: " + fb
  print ""

  ofb = open(fb, 'w')

  linel = "0.136985 0. 0. rn"
  liner = "# Mmin Mmax X ss/rn/sn=StepSize/RandomNumber/StepNumber [Mjup]"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line

  linel = "0.842 0. 0. rn"
  liner = "# Radius of Planet [Rjup]"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line
  
  linel = "19.238756 0. 0. rn"
  liner = "# Period [day] - set it > 9000000. to calculate from semi-major axis"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line

  linel = "999. 0. 0. rn"
  liner = "# semi major axis [AU] - set it to 999.0 if there is no value available"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line

  linel = "0.058152 0. 0. rn"
  liner = "# eccentricity"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line
  
  linel = "356.056648 0.6 0. rn"
  liner = "# argument of the pericenter [deg]"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line

  linel = "3.771733 0. 0. rn"
  liner = "# mean anomaly [deg]: if set to >= 999. the time of pericenter passage (next row) will be used"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line

  linel = "9.e8 0. 0. rn"
  liner = "# time of pericenter passage [JD]: if set to >= 9.e8 the mean anomaly (previous row) will be used"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line

  linel = "88.55 0. 0. rn"
  liner = "# orbit inclination [deg]"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line
  
  linel = "0. 0. 0. rn"
  liner = "# longitude of the ascending node [deg]"
  line = linel + "          " + liner
  ofb.write(line + "\n")
  print " " + line

  ofb.close()

# CREATES PLANETS FILES: c.dat
# 0.094348 0.094349 0.01 ss    # Mmin Mmax X ss/rn/sn=StepSize/RandomNumber/StepNumber [Mjup]
# 0.823 0. 0. rn               # Radius of Planet [Rjup]
# 38.986097 0. 0. rn           # Period [day] - set it > 9000000. to calculate from semi-major axis
# 999. 0. 0. rn                # semi major axis [AU] - set it to 999.0 if there is no value available
# 0.067645 0. 0. rn            # eccentricity
# 167.573369 0. 0. rn          # argument of the pericenter [deg]
# 307.432620 0. 0. rn          # mean anomaly[deg]
# 9.e8 0. 0. rn                # time pericenter[JD]
# 89.12 0. rn                  # orbit inclination [deg]
# 359.889907 0. rn             # longitude of the ascending node [deg]

def createPLANETC(fpath):
  fc = os.path.join(fpath, "c.dat")
  print " CREATING: " + fc
  print ""

  ofc = open(fc, 'w')

  linel = "0.094348 0.094349 0.01 ss"
  liner = "# Mmin Mmax X ss/rn/sn=StepSize/RandomNumber/StepNumber [Mjup]"
  line = linel + "    " + liner
  ofc.write(line + "\n")
  print " " + line

  linel = "0.823 0. 0. rn"
  liner = "# Radius of Planet [Rjup]"
  line = linel + "    " + liner
  ofc.write(line + "\n")
  print " " + line
  
  linel = "38.986097 0. 0. rn"
  liner = "# Period [day] - set it > 9.e6 to calculate from semi-major axis"
  line = linel + "    " + liner
  ofc.write(line + "\n")
  print " " + line

  linel = "999. 0. 0. rn"
  liner = "# semi major axis [AU] - set it to 999.0 if there is no value available"
  line = linel + "    " + liner
  ofc.write(line + "\n")
  print " " + line

  linel = "0.067645 0. 0. rn"
  liner = "# eccentricity"
  line = linel + "    " + liner
  ofc.write(line + "\n")
  print " " + line
  
  linel = "167.573369 0. 0. rn"
  liner = "# argument of the pericenter [deg]"
  line = linel + "    " + liner
  ofc.write(line + "\n")
  print " " + line

  linel = "307.432620 0. 0. rn"
  liner = "# mean anomaly [deg]: if set to >= 999. the time of pericenter passage (next row) will be used"
  line = linel + "    " + liner
  ofc.write(line + "\n")
  print " " + line

  linel = "9.e8 0. 0. rn"
  liner = "# time of pericenter passage [JD]: if set to >= 9.e8 the mean anomaly (previous row) will be used"
  line = linel + "          " + liner
  ofc.write(line + "\n")
  print " " + line

  linel = "89.12 0. 0. rn"
  liner = "# orbit inclination [deg]"
  line = linel + "    " + liner
  ofc.write(line + "\n")
  print " " + line
  
  linel = "359.889907 0. 0. rn"
  liner = "# longitude of the ascending node [deg]"
  line = linel + "    " + liner
  ofc.write(line + "\n")
  print " " + line

  ofc.close()


# CREATES lm.opt
# -1        # max function evaluation (if <=0 it will be set to 200*(n+1), where n=number of parameters to be fitted)
# 1.e-8  # ftol (must be > 0., if <=0 it will be setted with an internal TOLERANCE) DETERMINED IN THE PROGRAM
# 1.e-8  # xtol (must be > 0., if <=0 it will be setted with an internal TOLERANCE) DETERMINED IN THE PROGRAM
# 1.e-8  # gtol (must be >= 0., if <0 it will be setted to zero)
# -1        # value that multiply the parameters to detect a parameter step in lmdif (usually to be setted to sqrt( machine precision ), if <0 it will be set to this value) DETERMINED IN THE PROGRAM
# 0         # nprint: it prints some results every nprint function evaluations (if <0 it will be set to 0)

def createLM(fpath):
  flm = os.path.join(fpath, "lm.opt")
  print " CREATING: " + flm
  print ""

  oflm = open(flm, 'w')

  linel = "-1"
  liner = "# max function evaluation (if <=0 it will be set to 200*(n+1), where n=number of parameters to be fitted)"
  line = linel + "    " + liner
  oflm.write(line + "\n")
  print " " + line

  linel = "1.e-8"
  liner = "# ftol (must be > 0., if <=0 it will be setted with an internal TOLERANCE) DETERMINED IN THE PROGRAM"
  line = linel + "    " + liner
  oflm.write(line + "\n")
  print " " + line

  linel = "1.e-8"
  liner = "# xtol (must be > 0., if <=0 it will be setted with an internal TOLERANCE) DETERMINED IN THE PROGRAM"
  line = linel + "    " + liner
  oflm.write(line + "\n")
  print " " + line

  linel = "1.e-8"
  liner = "# gtol (must be >= 0., if <0 it will be setted to zero)"
  line = linel + "    " + liner
  oflm.write(line + "\n")
  print " " + line

  linel = "-1."
  liner = "# value that multiply the parameters to detect a parameter step in lmdif (usually to be setted to sqrt( machine precision ), if <0 it will be set to this value) DETERMINED IN THE PROGRAM"
  line = linel + "    " + liner
  oflm.write(line + "\n")
  print " " + line

  linel = "0"
  liner = "# nprint: it prints some results every nprint function evaluations (if <0 it will be set to 0)"
  line = linel + "    " + liner
  oflm.write(line + "\n")
  print " " + line

  oflm.close()

# CREATES pikaia.opt
# 10.     # ctrl(1) = number of individuals
# 10.     # ctrl(2) = number of generations
# 7.      # ctrl(3) = number of significant digits
# 0.85    # ctrl(4) = crossover probability; must be  <= 1.0 (default is 0.85)
# 2.      # ctrl(5) = mutation mode; 1/2=steady/variable (default is 2)
# 0.02    # ctrl(6) = initial mutation rate; should be small (default is 0.005) (Note: the mutation rate is the probability that any one gene locus will mutate in any one generation.)
# 0.0005  # ctrl(7) = minimum mutation rate; must be >= 0.0 (default is 0.0005)
# 0.3     # ctrl(8) = maximum mutation rate; must be <= 1.0 (default is 0.25)
# 1.      # ctrl(9) = relative fitness differential; range from 0 (none) to 1 (maximum).  (default is 1.)
# 1       # ctrl(10) = reproduction plan; 1/2/3=Full generational replacement/Steady-state-replace-random/Steady-state-replace-worst (default is 3)
# 0       # ctrl(11) = elitism flag; 0/1=off/on (default is 0) (Applies only to reproduction plans 1 and 2)
# 0       # ctrl(12) = printed output 0/1/2=None/Minimal/Verbose (default is 0)
# 733625  # seed
# 0       # wrtAll = 0 [not writing all individuals for each iteration] 1 [writing all individuals for each iteration]
# 1       # nGlobal = number of global search, number of times that PIKAIA will be run, i.e., Nindividual x Ngeneration x nGlobal

def createPIKAIA(fpath):
  fpik = os.path.join(fpath, "pikaia.opt")
  print " CREATING: " + fpik
  print ""

  ofpik = open(fpik, 'w')
  
  #line = []
  #line.append("10.     # ctrl(1) = number of individuals\n")
  #line.append("10.     # ctrl(2) = number of generations\n")
  #line.append("7.      # ctrl(3) = number of significant digits\n")
  #line.append("0.85    # ctrl(4) = crossover probability; must be  <= 1.0 (default is 0.85)\n")
  #line.append("2.      # ctrl(5) = mutation mode; 1/2=steady/variable (default is 2)\n")
  #line.append("0.02    # ctrl(6) = initial mutation rate; should be small (default is 0.005) (Note: the mutation rate is the probability that any one gene locus will mutate in any one generation.)\n")
  #line.append("0.0005  # ctrl(7) = minimum mutation rate; must be >= 0.0 (default is 0.0005)\n")
  #line.append("0.3     # ctrl(8) = maximum mutation rate; must be <= 1.0 (default is 0.25)\n")
  #line.append("1.      # ctrl(9) = relative fitness differential; range from 0 (none) to 1 (maximum).  (default is 1.)\n")
  #line.append("1       # ctrl(10) = reproduction plan; 1/2/3=Full generational replacement/Steady-state-replace-random/Steady-state-replace-worst (default is 3)\n")
  #line.append("0       # ctrl(11) = elitism flag; 0/1=off/on (default is 0) (Applies only to reproduction plans 1 and 2)\n")
  #line.append("0       # ctrl(12) = printed output 0/1/2=None/Minimal/Verbose (default is 0)\n")
  #line.append("123456  # seed\n")
  #line.append("0       # wrtAll = 0 [not writing all individuals for each iteration] 1 [writing all individuals for each iteration]\n")
  #line.append("1       # nGlobal = number of global search, number of times that GA will be run, i.e., Nindividual x Ngeneration x nGlobal\n")
  #for ll in line:
    #ofpik.write(ll)
    #print ll.strip('\n')
  #ofpik.close()

  ofpik.write("10.     # ctrl(1) = number of individuals\n")
  ofpik.write("10.     # ctrl(2) = number of generations\n")
  ofpik.write("7.      # ctrl(3) = number of significant digits\n")
  ofpik.write("0.85    # ctrl(4) = crossover probability; must be  <= 1.0 (default is 0.85)\n")
  ofpik.write("2.      # ctrl(5) = mutation mode; 1/2=steady/variable (default is 2)\n")
  ofpik.write("0.02    # ctrl(6) = initial mutation rate; should be small (default is 0.005) (Note: the mutation rate is the probability that any one gene locus will mutate in any one generation.)\n")
  ofpik.write("0.0005  # ctrl(7) = minimum mutation rate; must be >= 0.0 (default is 0.0005)\n")
  ofpik.write("0.3     # ctrl(8) = maximum mutation rate; must be <= 1.0 (default is 0.25)\n")
  ofpik.write("1.      # ctrl(9) = relative fitness differential; range from 0 (none) to 1 (maximum).  (default is 1.)\n")
  ofpik.write("1       # ctrl(10) = reproduction plan; 1/2/3=Full generational replacement/Steady-state-replace-random/Steady-state-replace-worst (default is 3)\n")
  ofpik.write("0       # ctrl(11) = elitism flag; 0/1=off/on (default is 0) (Applies only to reproduction plans 1 and 2)\n")
  ofpik.write("0       # ctrl(12) = printed output 0/1/2=None/Minimal/Verbose (default is 0)\n")
  ofpik.write("123456  # seed\n")
  ofpik.write("0       # wrtAll = 0 [not writing all individuals for each iteration] 1 [writing all individuals for each iteration]\n")
  ofpik.write("1       # nGlobal = number of global search, number of times that GA will be run, i.e., Nindividual x Ngeneration x nGlobal\n")
  ofpik.close()

# CREATES pso.opt
# 10     # number of particle
# 10     # number of iterations to do
# 10     # counter to write to screen PSO run...if 0 = no write
# 0      # wrtAll = 0 [not writing all individuals for each iteration] 1 [writing all individuals for each iteration]
# 1      # nGlobal = number of global search, number of times that PIKAIA will be run, i.e., Npart x Niter x nGlobal

def createPSO(fpath):
  fpso = os.path.join(fpath, "pso.opt")
  print " CREATING: " + fpso
  print ""

  ofpso = open(fpso, 'w')
  
  #line = []
  #line.append("10     # number of particle\n")
  #line.append("10     # number of iterations to do\n")
  #line.append("0      # counter to write to screen PSO run...if 0 = no write\n")
  #line.append("0      # wrtAll = 0 [not writing all individuals for each iteration] 1 [writing all individuals for each iteration]\n")
  #line.append("1      # nGlobal = number of global search, number of times that PSO will be run, i.e., Npart x Niter x nGlobal\n")
  #for ll in line:
    #ofpso.write(ll)
    #print ll.strip('\n')
  #ofpso.close()

  ofpso.write("10     # number of particle\n")
  ofpso.write("10     # number of iterations to do\n")
  ofpso.write("0      # counter to write to screen PSO run...if 0 = no write\n")
  ofpso.write("0      # wrtAll = 0 [not writing all individuals for each iteration] 1 [writing all individuals for each iteration]\n")
  ofpso.write("1      # nGlobal = number of global search, number of times that PSO will be run, i.e., Npart x Niter x nGlobal\n")
  ofpso.write("123456  # seed\n")
  ofpso.close()

# create PolyChord configuration file
def createPC(fpath):
  fpc = os.path.join(fpath, "PolyChord.opt")
  print " CREATING: " + fpc
  print ""

  ofpc = open(fpc, 'w')
  
  ofpc.write("# The number of live points: nlive = [integer]\n")
  ofpc.write("nlive = 500\n")
  
  ofpc.write("# The number of chords to draw: num_repeats = [integer]\n")
  ofpc.write("num_repeats = 1\n")
  
  ofpc.write("# whether or not to do clustering: do_clustering = [logical] .false. or F, .true. or T\n")
  ofpc.write("do_clustering = T\n")
  
  #ofpc.write("# The maximum number of active clusters (Memory storage considerations): ncluster = [integer]\n")
  #ofpc.write("ncluster = 30\n")
  
  ofpc.write("# The degree of precision in the final answer: precision_criterion = [real]\n")
  ofpc.write("precision_criterion = 1.e-3\n")
  
  ofpc.write("# The maximum number of dead points/samples: max_ndead = [integer], set to -1 for no maximum number\n")
  ofpc.write("max_ndead = -1\n")
  
  ofpc.write("# The degree of feedback to provide: feedback = 0 [no feedback], 1 [normal feedback], 2 [high feedback]\n")
  ofpc.write("feedback = 1.\n")
  
  ofpc.write("# Whether or not to calculate weighted posteriors: posteriors = [logical] .false. or F, .true. or T\n")
  ofpc.write("posteriors = T\n")
  
  ofpc.write("# Whether or not to calculate equally weighted posteriors: equals = [logical] .false. or F, .true. or T\n")
  ofpc.write("equals = T\n")
  
  ofpc.write("# Whether or not to calculate clustered posteriors: cluster_posteriors = [logical] .false. or F, .true. or T\n")
  ofpc.write("cluster_posteriors = T\n")
  
  ofpc.write("# Whether or not to calculate weighted posterior: updated_posterior = [integer]\n")
  ofpc.write("update_posterior = 500\n")
  
  ofpc.write("# What factor should we bulk up the posterior points by (usinginter-chain points): boost_posterior = [real], <=0 to use all of them \n")
  ofpc.write("boost_posterior = 5.\n")
  
  ofpc.write("# Whether or not to resume from file: read_resume = [logical] .false. or F, .true. or T\n")
  ofpc.write("read_resume = F\n")
  
  ofpc.write("# Whether or not to write resume files: write_resume = [logical] .false. or F, .true. or T\n")
  ofpc.write("write_resume = T\n")
  
  #ofpc.write("# How often to update the resume file: update_resume = [integer], if omitted it will be set equal to nlive\n")
  #ofpc.write("update_resume = 500\n")
  
  ofpc.write("# Whether or not to write phys_live points: write_live = [logical] .false. or F, .true. or T\n")
  ofpc.write("write_live = T\n")
  
  ofpc.write("# Whether or not to write stats file: write_stats = [logical] .false. or F, .true. or T\n")
  ofpc.write("write_stats = T\n")
  
  ofpc.write("# The directory to put outputs in: base_dir = [character(512)], if omitted it will be set to 'chains'\n")
  ofpc.write("base_dir = 'chains'\n")
  
  ofpc.write("# The file root for outputs: file_root = [character(512)], if omitted it will be set to 'trades_PolyChord'\n")
  ofpc.write("file_root = 'trades_PolyChord'\n")
  
  ofpc.write("# The fraction of time spent in each grade: grade_frac = [real] or list of [real]s, e.g. 1. or 0.3 0.3 0.4\n")
  ofpc.write("grade_frac = 1.\n")
  
  ofpc.close()

# CREATES obsRV.dat
# JD RVobs eRVobs RVsetID
#2455342.947836 17.15 2.51 1
#2455344.061419 2.23 2.74 1
#2455351.037415 -7.89 2.36 1
#2455352.011289 -4.61 2.21 1
#2455367.083543 -25.45 2.49 1
#2455376.930532 17.39 2.61 1
def createOBSRV(fpath):
  fobsRV = os.path.join(fpath, "obsRV.dat")
  print " CREATING: " + fobsRV
  print ""

  ofobsRV = open(fobsRV, 'w')
  
  #line = []
  #line.append("# JD RVobs eRVobs RVsetID\n")
  #line.append("2455342.947836  17.15  2.51 1\n")
  #line.append("2455344.061419   2.23  2.74 1\n")
  #line.append("2455351.037415  -7.89  2.36 1\n")
  #line.append("2455352.011289  -4.61  2.21 1\n")
  #line.append("2455367.083543 -25.45  2.49 1\n")
  #line.append("2455376.930532  17.39  2.61 1\n")
  #for ll in line:
    #ofobsRV.write(ll)
    #print ll.strip('\n')
  #ofobsRV.close()

  ofobsRV.write("# JD RVobs eRVobs RVsetID\n")
  ofobsRV.write("2455342.947836  17.15  2.51 1\n")
  ofobsRV.write("2455344.061419   2.23  2.74 1\n")
  ofobsRV.write("2455351.037415  -7.89  2.36 1\n")
  ofobsRV.write("2455352.011289  -4.61  2.21 1\n")
  ofobsRV.write("2455367.083543 -25.45  2.49 1\n")
  ofobsRV.write("2455376.930532  17.39  2.61 1\n")
  ofobsRV.close()


# CREATES NB2_observations.dat
# epoch T_0 eT_0
#  -5  2454977.24875  0.00059
#  -4  2454996.48240  0.00062
#  -2  2455034.95437  0.00052
#  -1  2455054.19058  0.00052
#   0  2455073.43412  0.00072
#   2  2455111.92589  0.00050
#   3  2455131.17167  0.00048
#   4  2455150.42951  0.00048
#   5  2455169.68103  0.00046
def createNB2OBS(fpath):
  fNB2obs = os.path.join(fpath, "NB2_observations.dat")
  print " CREATING: " + fNB2obs
  print ""

  ofNB2obs = open(fNB2obs, 'w')

  #line = []
  #line.append("# epoch T_0 eT_0\n"
  #line.append(" -5  2454977.24875  0.00059\n")
  #line.append(" -4  2454996.48240  0.00062\n")
  #line.append(" -2  2455034.95437  0.00052\n")
  #line.append(" -1  2455054.19058  0.00052\n")
  #line.append("  0  2455073.43412  0.00072\n")
  #line.append("  2  2455111.92589  0.00050\n")
  #line.append("  3  2455131.17167  0.00048\n")
  #line.append("  4  2455150.42951  0.00048\n")
  #line.append("  5  2455169.68103  0.00046\n")
  #for ll in line:
    #ofNB2obs.write(ll)
    #print ll.strip('\n')
  #ofNB2obs.close()

  ofNB2obs.write("# epoch T_0 eT_0\n")
  ofNB2obs.write(" -5  2454977.24875  0.00059\n")
  ofNB2obs.write(" -4  2454996.48240  0.00062\n")
  ofNB2obs.write(" -2  2455034.95437  0.00052\n")
  ofNB2obs.write(" -1  2455054.19058  0.00052\n")
  ofNB2obs.write("  0  2455073.43412  0.00072\n")
  ofNB2obs.write("  2  2455111.92589  0.00050\n")
  ofNB2obs.write("  3  2455131.17167  0.00048\n")
  ofNB2obs.write("  4  2455150.42951  0.00048\n")
  ofNB2obs.write("  5  2455169.68103  0.00046\n")
  ofNB2obs.close()

# CREATES NB3_observations.dat
# epoch T_0 eT_0
# -3  2454969.30577  0.00070
# -2  2455008.33086  0.00061
# -1  2455047.33560  0.00058
#  0  2455086.31251  0.00059
#  1  2455125.26284  0.00053
#  2  2455164.18168  0.00055
def createNB3OBS(fpath):
  fNB3obs = os.path.join(fpath, "NB3_observations.dat")
  print " CREATING: " + fNB3obs
  print ""

  ofNB3obs = open(fNB3obs, 'w')
  
  #line = []
  #line.append("# epoch T_0 eT_0\n")
  #line.append(" -3  2454969.30577  0.00070\n")
  #line.append(" -2  2455008.33086  0.00061\n")
  #line.append(" -1  2455047.33560  0.00058\n")
  #line.append("  0  2455086.31251  0.00059\n")
  #line.append("  1  2455125.26284  0.00053\n")
  #line.append("  2  2455164.18168  0.00055\n")
  #for ll in line:
    #ofNB3obs.write(ll)
    #print ll.strip('\n')
  #ofNB3obs.close()

  ofNB3obs.write("# epoch T_0 eT_0\n")
  ofNB3obs.write(" -3  2454969.30577  0.00070\n")
  ofNB3obs.write(" -2  2455008.33086  0.00061\n")
  ofNB3obs.write(" -1  2455047.33560  0.00058\n")
  ofNB3obs.write("  0  2455086.31251  0.00059\n")
  ofNB3obs.write("  1  2455125.26284  0.00053\n")
  ofNB3obs.write("  2  2455164.18168  0.00055\n")
  ofNB3obs.close()


# --------------------------------------------------------------------------
# MAIN


fpath = get_args()
createFolder(fpath)
print ""
createARG(fpath)
print ""
createBODIESLST(fpath)
print ""
createSTAR(fpath)
print ""
createPLANETB(fpath)
print ""
createPLANETC(fpath)
print ""
createLM(fpath)
print ""
createPIKAIA(fpath)
print ""
createPSO(fpath)
print ""
createPC(fpath)
print ""
createOBSRV(fpath)
print ""
createNB2OBS(fpath)
print ""
createNB3OBS(fpath)
print ""
