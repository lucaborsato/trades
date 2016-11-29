#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os #  os: operating system
#import time # time: execution time
import glob # glob: globbing file...loading multiple files as *.pippa
import sys # sys: system...boh
#import datetime # datetime: day and so...
import numpy as np # array
#import veusz.embed as vsz
#sys.path.append("/media/Data/Dropbox/Dottorato/pyscript")
sys.path.append("../")
import constants as cst
#from path_check import checkpath
#from shapes import circle
# --------------------------------------------------------------------------

# WARNING
# TO FIX

# --------------------------------------------------------------------------
# THIS SCRIPT READ #_rotorbit.dat WITH ORBIT DATA and #_initialElements.dat.
# IT PREPARES DATA TO BE PLOTTED BY ANOTHER SCRIPT IN ../veusz/
# NEEDED TO PLOT SKY-PLANE ORBITS ZOOMED ON THE STAR, ORBIT PROJECTIONS OF THE ORBITS,
# AND THE RADIAL VELOCITY MODEL.
# THE MARKER SIZE IS SCALED ON THE MASSES OF THE PLANET
# --------------------------------------------------------------------------

# read command arguments: full_path_folder, id_sim, Rstar[Rsun]
def readArgs():
  args = sys.argv
  del args[0]
  nargs = len(sys.argv)
  if (nargs < 2):
    sys.exit('Enter 4 arguments: full_path_folder id_sim LMfitted[0=no,1=yes] zoom/nozoom')
  elif (nargs == 2):
    args.append('0')
    args.append('zoom')
  elif (nargs == 3):
    args.append('zoom')
  return args

# convert command line arguments into variables
def getArgs():
  args = readArgs()
  fpath = os.path.abspath(args[0])
  idsim = args[1]
  lmf = args[2]
  zoom = args[3]
  print 
  print " FOLDER PATH: " + fpath
  print "       IDsim: " + idsim
  print "      LMflag: " + lmf
  print "    plotting: " + zoom
  print 
  return fpath, idsim, lmf, zoom


# read first row of bodies.lst file and then read star file with Mass and Radius in Solar values
def readStar(fpath):
  fbd = os.path.join(fpath, "bodies.lst")
  ofbd = open(fbd, 'r')
  line = ofbd.readline()
  ofbd.close()
  fstar = os.path.join(fpath, line.split(' ')[0])
  ofstar = open(fstar, 'r')
  Mstar = float( ofstar.readline().split('\t')[0].split(' ')[0] )
  Rstar = float( ofstar.readline().split('\t')[0].split(' ')[0] )
  ofstar.close()
  print " M_star = %.4f M_sun --- R_star = %.4f R_sun = %.12f R_au" %(Mstar, Rstar, Rstar*cst.RsunAU)
  return Mstar, Rstar

# read #_rotorbit.dat
def readOrbit(fpath, idsim, lmf):
  orbnm = str(idsim) + "_" + str(lmf) + "_rotorbit.dat"
  forb = os.path.join(fpath, orbnm)
  data = np.genfromtxt(forb)
  print " read file " + forb
  data0 = data[0].copy()
  print " selected initial state vector "
  idx = np.argsort(data[:,0])
  #print data.shape, data[0,0]
  datasorted = data[idx,:]
  print " sorted orbit data "
  #print datasorted.shape, datasorted[0,0]
  ncols = data.shape[1]
  nb = int((ncols - 3) / 6)
  print " determined NB = ",nb
  #data0 = initial state vector
  #datasorted = orbit data sorted in time
  #nb = number of bodies of the system
  return data0, datasorted, nb

# select the BJD epoch = BJDref,
# the timing of the integration BJD,
# the BJD respect to the BJDref (BJDepo),
# the BJD respect to the lower BJD (BJDint)
def getBJD(iState, orbit):
  BJDref = iState[0]
  BJD = orbit[:,0]
  BJDepo = BJD - BJDref
  BJDint = BJD - BJD[0]
  intTime = BJD[-1] - BJD[0]
  #print " epoch BJD = %25.12f " %(BJDref)
  #print "   min BJD = %25.12f ==> BJD - BJDref = %25.12f" %(BJD[0], BJDepo[0])
  #print "   max BJD = %25.12f ==> BJD - BJDref = %25.12f" %(BJD[-1], BJDepo[-1])
  #print " integration time = %.4f days" %(intTime)
  return BJDref, BJDint, BJDepo, BJD, intTime

# it returns the column indices of the orbit array
def idxXYZ(orbit, nb):
  ibd = np.arange(2,nb+1)
  ix = 2 + (ibd - 1) * 6
  iy = ix + 1
  iz = ix + 2
  return ix, iy, iz

# read state vector for body number ibd and it gives x,y,z as output
def orb2xyz(state, ibd):
  vtype = np.array(state.shape).shape[0]
  ix = 2 + (ibd - 1) * 6
  iy = ix + 1
  iz = ix + 2
  if vtype == 1:
    x = state[ix].copy()
    y = state[iy].copy()
    z = state[iz].copy()
  else:
    x = state[:,ix].copy()
    y = state[:,iy].copy()
    z = state[:,iz].copy()
  return x,y,z

# read initial orbital elements from the #_initialElements.dat
def readElements(path, idsim, lmf):
  elef = str(idsim) + "_" + str(lmf) + "_initialElements.dat"
  fele = os.path.join(path, elef)
  elements = np.genfromtxt(fele)
  print " read elements file: " + elef
  return elements

# select X, Y, Z for Z values higher than selZ
def selZhigher1(idb, Xall, Yall, Zall, selZ):
  iib = idb - 2
  Xx = Xall[:,iib]
  Yy = Yall[:,iib]
  Zz = Zall[:,iib]
  Xok = Xx[Zz>=selZ]
  Yok = Yy[Zz>=selZ]
  Zok = Zz[Zz>=selZ]
  return Xok, Yok, Zok

# select X, Y, Z for Z values higher than selZ
def selZhigher2(idb, BJDall, Xall, Yall, Zall, selZ):
  iib = idb - 2
  Xx = Xall[:,iib]
  Yy = Yall[:,iib]
  Zz = Zall[:,iib]
  Xok = Xx[Zz>=selZ]
  Yok = Yy[Zz>=selZ]
  Zok = Zz[Zz>=selZ]
  BJDok = BJDall[Zz>=selZ]
  return BJDok, Xok, Yok, Zok


# select X, Y, Z for Z values lower than selZ
def selZlower1(idb, Xall, Yall, Zall, selZ):
  iib = idb - 2
  Xx = Xall[:,iib]
  Yy = Yall[:,iib]
  Zz = Zall[:,iib]
  Xok = Xx[Zz<selZ]
  Yok = Yy[Zz,selZ]
  Zok = Zz[Zz<selZ]
  return Xok, Yok, Zok

# select X, Y, Z for Z values lower than selZ
def selZlower2(idb, BJDall, Xall, Yall, Zall, selZ):
  iib = idb - 2
  Xx = Xall[:,iib]
  Yy = Yall[:,iib]
  Zz = Zall[:,iib]
  Xok = Xx[Zz<selZ]
  Yok = Yy[Zz<selZ]
  Zok = Zz[Zz<selZ]
  BJDok = BJDall[Zz<selZ]
  return BJDok, Xok, Yok, Zok

# select the axis (X=0, Y=1) to sort the orbit vector of one body
def selXYSort(idb, elem):
  if idb == 2:
    ii = 0
  else:
    ii = idb - 2
  nRowEl = len(elem.shape)
  if nRowEl == 1:
    lN = elem[-1]
  else:
    lN = elem[ii,-1]
  lN1 = lN%360.
  if lN1 < 0.:
    lN1 = lN1 + 360.
  selAx = 0 # X axis to check
  if not( (lN1 >= 0. and lN1 <= 45.) or (lN1 > 135. and lN1 <= 225.) or (lN1 > 315. and lN1 < 360.) ):
    selAx = 1
  return selAx

# call selXYSort for all the bodies and create vector con axis selection
def getSelXYcheck(nb, elem):
  selAx = []
  for i in range(2, nb+1):
    selAx.append(selXYSort(i, elem))
  selAx = np.array(selAx)
  return selAx

# rename to variable the RVmodel data
def getRVmodel(orbit):
  RVmod = orbit[:,-1]
  return RVmod

# creates if the dir for plots and files if it does not exist
def createPlotsDir(fpath):
  plotFolder = os.path.join(fpath, "plots")
  if not os.path.isdir(plotFolder):
    createDir = "mkdir -p " + plotFolder
    os.system(createDir)
  return plotFolder

# given an orbit it take a mid time of the integration
# finds a X or Y variation and takes just +/- P/2 from that point
def selOnePeriodOrbit(idb, elem, intTime, BJD, Xa, Ya, selAx):
  nRowEl = len(elem.shape)
  idj = np.argsort(BJD)
  BJD1 = BJD[idj]
  Xa1 = Xa[idj]
  Ya1 = Ya[idj]
  nsteps = BJD1.shape[0]
  if nRowEl == 1:
    Ph = elem[2] * 0.5
    Pq = elem[2] * 0.33
  else:
    Ph = elem[idb-2,2] * 0.5
    Pq = elem[idb-2,2] * 0.33
  if Ph < intTime*0.5:
    nmid = int(nsteps * 0.5)
  else:
    nmid = 0
  print " nmid = ", nmid, " nsteps = ", nsteps
  ntra = 0
  for i in range(nmid, nsteps):
    if selAx[idb-2] == 0:
      if Xa1[i-1]*Xa1[i] <= 0.:
        ntra = i
        break
    else:
      if Ya1[i-1]*Ya1[i] <= 0.:
        ntra = i
        break
  BJDtra = BJD1[ntra]
  #BJDmin = BJDtra - Ph
  #BJDmax = BJDtra + Ph
  BJDmin = BJDtra - Pq
  BJDmax = BJDtra + Pq
  BJDt = BJD1[BJD1>=BJDmin]
  Xt = Xa1[BJD1>=BJDmin]
  Yt = Ya1[BJD1>=BJDmin]
  print " BJDtra = ", BJDtra
  print " BJDmin = ", BJDmin, " BJDt[0] = ", BJDt[0]
  BJDok = BJDt[BJDt<=BJDmax]
  Xok = Xt[BJDt<=BJDmax]
  Yok = Yt[BJDt<=BJDmax]
  print " BJDmax = ", BJDmax, " BJDok[-1] = ", BJDok[-1]
  return BJDok, Xok, Yok

# saves single body data for sky-plane plot
def saveSkyFile(plotFolder, idsim, lmf, idb, BJD, Xa, Ya, selAx):
  skyF = str(idsim) + "_" + str(lmf) + "_NB" + str(idb) + "_skyplane.dat"
  skyP = os.path.join(plotFolder, skyF)
  print " ready to write file " + skyP
  oSky = open(skyP, 'w')
  ids = np.argsort(BJD)
  #if selAx[idb-2] == 0:
    #ids = np.argsort(Xa)
  #else:
    #ids = np.argsort(Ya)
  Xs = Xa[ids]
  Ys = Ya[ids]
  line = "# BJD_%d X_%d Y_%d\n" %(idb, idb, idb)
  oSky.write(line)
  nline = Xs.shape[0]
  for i in range(0, nline):
    line = "%25.15f %25.15f %25.15f\n" %(BJD[i], Xs[i], Ys[i])
    oSky.write(line)
  oSky.close()

# saves all data for sky-plane plot
def saveAllSky(plotFolder, idsim, lmf, nb, elem, iX, iY, iZ, intTime, BJD, Xall, Yall, Zall, selAx):
  for idb in range(2, nb+1):
    print " body number " + str(idb)
    BJDa, Xa, Ya, Za = selZhigher2(idb, BJD, Xall, Yall, Zall, 0.)
    print " number of steps with Z > 0: %d" %(BJDa.shape[0])
    if BJDa.shape[0] > 0:
      BJDok, Xok, Yok = selOnePeriodOrbit(idb, elem, intTime, BJDa, Xa, Ya, selAx)
      saveSkyFile(plotFolder, idsim, lmf, idb, BJDok, Xok, Yok, selAx)

# ===== Function Main
def orbitPrepare(fpath, idsim, lmf):
  # read Star data
  Mstar, Rstar = readStar(fpath)
  Rsau = Rstar * cst.RsunAU
  elem = readElements(fpath, idsim, lmf)
  iState, orbit, nb = readOrbit(fpath, idsim, lmf)
  idX, idY, idZ = idxXYZ(orbit, nb)
  BJDref, BJDint, BJDepo, BJD, intTime = getBJD(iState, orbit)
  Xall = orbit[:,idX].copy()
  Yall = orbit[:,idY].copy()
  Zall = orbit[:,idZ].copy()
  #print "Xall: ", Xall.shape, Xall[0]
  #print "Yall: ", Yall.shape, Yall[0]
  RVmod = getRVmodel(orbit)
  selAx = getSelXYcheck(nb, elem)
  plotFolder = createPlotsDir(fpath)
  saveAllSky(plotFolder, idsim, lmf, nb, elem, idX, idY, idZ, intTime, BJD, Xall, Yall, Zall, selAx)
  return nb, Mstar, Rsau, elem, BJDref, BJD, BJDepo, iState, Xall, Yall, Zall, idX, idY, idZ, RVmod, plotFolder



