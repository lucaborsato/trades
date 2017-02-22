#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os #  os: operating system
import time # time: execution time
import glob # glob: globbing file...loading multiple files as *.pippa
import sys # sys: system...boh
import numpy as np # array
import veusz.embed as vsz
import constants as cst
#import numpy.polynomial.polynomial as poly

# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# THIS SCRIPT PLOTS THE O-C DIAGRAMS OF SELECTED BODIES (IN RANGE) AND
# IT PLOTS THE COMPARISON OBS-SIM RV AT THE BOTTOM.
# --------------------------------------------------------------------------

def get_bool(bool_in):
  true_list = 'true tru tr t yes ye y 1'.split()
  if(str(bool_in).lower() in true_list):
    bool_out = True
  else:
    bool_out = False
  return bool_out

def get_args():
  print " EMBED VEUSZ VERSION ", vsz.API_VERSION
  parser = argparse.ArgumentParser()
  parser.add_argument('-p', action='store', dest='full_path', required=True, help='Folder path')
  parser.add_argument('-s', action='store', dest='idsim', required=True, help='Simulation ID number')
  parser.add_argument('-lm', action='store', dest='lmflag', default='0', help='LM flag: 0 = not used LM (default), 1 = used LM')
  parser.add_argument('-fb', action='store', dest='nbI', default='2', help='First body to plot: e.g., body = 2 is the first planet (default 2)')
  parser.add_argument('-lb', action='store', dest='nbF', default='2', help='Last body to plot: e.g., body = 3 is the first planet (default 2)')
  parser.add_argument('-rv', action='store', dest='RVset', default='0', help='RVsetlag: 0/No/False = do not use RV data, 1/Yes/True = use RV data')
  parser.add_argument('-xs', '--x-scale', action='store', dest='xscale', default=None, help='Value to be subtract to the x-axis (time) in the plots. Default is None that means 2440000.5 will be subtract.')
  parser.add_argument('-oc', '--oc-fit', action='store', dest='oc_fit', default='False', help='Fit obs and sim linear ephem (set to False) or only the obs ephemem (set to True)')
  cli = parser.parse_args()
  full_path = os.path.abspath(cli.full_path)
  cli.idsim = int(cli.idsim)
  cli.lmflag = int(cli.lmflag)
  cli.nbI = int(cli.nbI)
  cli.nbF = int(cli.nbF)
  cli.oc_fit = get_bool(cli.oc_fit)
  cli.RVset = get_bool(cli.RVset)
  
  #return full_path, cli.idsim, cli.lmflag, int(cli.nbI), int(cli.nbF), cli.RVset, cli.xscale, cli.oc_fit
  return cli

def intro(fullpath, nbI, nbF):
  print ""
  print "--- O-C plot ---"
  print ""
  print "Folder:"
  print fullpath
  print "Bodies to plot, from %d to %d"%(nbI,nbF)
  print ""

# linear fit function with errors on y-axis.
# It returns also the errors.
# y = b_ls + m_ls * x
def lstsq_fit(x, y, yerr):
  A = np.vstack((np.ones_like(x), x)).T
  C = np.diag(yerr * yerr)
  cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
  b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))
  coeff_err = np.sqrt(np.diag(cov))
  return b_ls, m_ls, coeff_err


def read_simT0(path, ssim, lmf, ibd, t_subtract=0., oc_fit=False):
  if(ibd >= 2):
    #simT0 = ssim + "_" + str(lmf) + "_NB" + str(ibd) + "_simT0.dat"
    simT0 = '%d_%d_NB%d_simT0.dat' %(ssim, lmf, ibd)
    fT0 = os.path.join(path, simT0)
    print " OPENING FILE " + fT0
    try:
      with open(fT0):

        print "Reading file %s" %(simT0)
        data = np.genfromtxt(fT0) # epoT0obs T0obs eT0obs T0_sim T0obs-T0_sim T0_stat

        epo = data[:,0].astype(int)
        nt = epo.shape[0]
    
        t_obs_subtract = data[:,1]-t_subtract
    
        err = np.zeros((nt,3))
        err[:,0] = data[:,2].copy()
        err[:,1] = err[:,0] * 1440.
        err[:,2] = err[:,0] * 86400.

        res = np.zeros((nt,3))

        if(not oc_fit):
          res[:,0] = data[:,4].copy()
          res[:,1] = res[:,0] * 1440.
          res[:,2] = res[:,0] * 86400.

          #Teph, Peph = poly.polyfit(epo, data[:,1], deg=1, w=1./err[:,0])
          Teph, Peph, TP_err = lstsq_fit(epo, data[:,1], err[:,0])
          tc = Teph + Peph * epo
          print "Linear ephemeris:"
          print "T_ephem [JD] = %20.7f +/- %9.7f" %(Teph, TP_err[0])
          print "P_ephem [d]  = %20.7f +/- %9.7f" %(Peph, TP_err[1])

          oc = np.zeros((nt,6))
          oc[:,0] = data[:,1] - tc
          oc[:,1] = oc[:,0] * 1440.
          oc[:,2] = oc[:,0] * 86400.
          oc[:,3] = data[:,3] - tc
          oc[:,4] = oc[:,3] * 1440.
          oc[:,5] = oc[:,3] * 86400.

        else:
          #Teph_obs, Peph_obs = poly.polyfit(epo, data[:,1], deg=1, w=1./err[:,0])
          Teph_obs, Peph_obs, TP_err_obs = lstsq_fit(epo, data[:,1], err[:,0])
          tc_o = Teph_obs + Peph_obs * epo
          print "Linear ephemeris observed:"
          print "T_ephem [JD] = %20.7f +/- %9.7f" %(Teph_obs, TP_err_obs[0])
          print "P_ephem [d]  = %20.7f +/- %9.7f" %(Peph_obs, TP_err_obs[1])
          
          #Teph_sim, Peph_sim = poly.polyfit(epo, data[:,3], deg=1, w=np.ones((nt))/err[:,0].mean())
          Teph_sim, Peph_sim, TP_err_sim = lstsq_fit(epo, data[:,3], np.ones((nt))*err[:,0].mean())
          tc_s = Teph_sim + Peph_sim * epo
          print "Linear ephemeris simulated:"
          print "T_ephem [JD] = %20.7f +/- %9.7f" %(Teph_sim, TP_err_sim[0])
          print "P_ephem [d]  = %20.7f +/- %9.7f" %(Peph_sim, TP_err_sim[1])

          oc = np.zeros((nt,6))
          oc[:,0] = data[:,1] - tc_o
          oc[:,1] = oc[:,0] * 1440.
          oc[:,2] = oc[:,0] * 86400.
          oc[:,3] = data[:,3] - tc_s
          oc[:,4] = oc[:,3] * 1440.
          oc[:,5] = oc[:,3] * 86400.
          
          res[:,0] = oc[:,0] - oc[:,3]
          res[:,1] = oc[:,1] - oc[:,4]
          res[:,2] = oc[:,2] - oc[:,5]

    except IOError:
      print '%s does not exist ... closing script' %(simT0)
      sys.exit()

  return epo, t_obs_subtract, oc, res, err

# define the scaling of the rows
def scalingRows(gRows):
  scRows = []
  for i in range(0, gRows, 2):
    scRows.append(1.)
    scRows.append(0.90)
  return scRows


def lims(v1, v2, inc=0.05):
  min1 = min(v1)
  min2 = min(v2)
  max1 = max(v1)
  max2 = max(v2)
  min3 = min(min1,min2)
  max3 = max(max1,max2)
  d3 = max3 - min3
  #inc = 0.05
  min4 = float(min3 - abs(d3)*inc) # use float because veusz want float and not np.float64 ...
  max4 = float(max3 + abs(d3)*inc) # use float because veusz want float and not np.float64 ...
  return min4, max4

def read_obsRV(fullpath):
  obsRV = os.path.join(fullpath, "obsRV.dat")
  RV_o = np.zeros(1)
  stat = False
  if os.path.exists(obsRV):
    RV_o = np.genfromtxt(obsRV)
    stat = True
  return RV_o, stat

# read #_simRV_1.dat with columns:
# jd RVobs eRVobs rv_sim RV_sim gamma e_gamma RVsetID RV_stat
def read_simRV(fullpath, sim, lmf):
  fRV = os.path.join(fullpath, '%d_%d_simRV.dat' %(sim, lmf))
  RV = np.zeros(1)
  stat = False
  gamma = np.zeros(2)
  nRV_set = 1
  if (os.path.exists(fRV)):
    RV = np.genfromtxt(fRV)
    gamma = RV[:,5:6]
    n_rv = RV.shape[0]
    for irv in range(1, n_rv):
      if (RV[irv,7] != RV[irv-1,7]):
        nRV_set += 1
    stat = True
  return RV, gamma, nRV_set, stat

def read_orbit(fullpath, sim, lmf):
  #orbit = sim + "_" + str(lmf) + "_rotorbit.dat"
  orbit = '%d_%d_rotorbit.dat' %(sim, lmf)
  forbit = os.path.join(fullpath, orbit)
  jd_lte = np.zeros(1)
  rv = np.zeros(1)
  stat = False
  if (os.path.exists(forbit)):
    data = np.genfromtxt(forbit, usecols = (0,1,-1))
    idx = np.argsort(data[:,0])
    jd_lte = data[idx,0] + data[idx,1]
    rv = data[idx,2]
    stat = True
  return jd_lte, rv, stat

# x-axis of O-C
def xOC(ocPlot, xname, xlab, min_time, max_time, font, sizeLab, LabHide, sizeTL):
  x_oc = ocPlot.Add('axis', name = xname, autoadd = False,
                    label = xlab,
                    min = float(min_time),
                    max = float(max_time),
                    log = False,
                    mode = "numeric",
                    autoExtend = False,
                    autoExtendZero = False,
                    autoMirror = True,
                    Label__font = font,
                    Label__size = sizeLab,
                    Label__hide = LabHide,
                    Label__offset = '0pt',
                    TickLabels__font = font,
                    TickLabels__size = sizeTL,
                    TickLabels__format = '%.1f'
                    )
  return x_oc

# x-axis of RV
def xRV(ocPlot, xname, xlab, min_time, max_time, font, sizeLab, LabHide, sizeTL):
  
  x_oc = ocPlot.Add('axis', name = xname, autoadd = False,
                    label = xlab,
                    min = float(min_time),
                    max = float(max_time),
                    log = False,
                    mode = "numeric",
                    autoExtend = False,
                    autoExtendZero = False,
                    autoMirror = True,
                    Label__font = font,
                    Label__size = sizeLab,
                    Label__hide = LabHide,
                    Label__offset = '0pt',
                    TickLabels__font = font,
                    TickLabels__size = sizeTL,
                    TickLabels__format = '%.1f'
                    )
  return x_oc


# y-axis of O-C ... position: 0=left, 1=right
def yOC(ocPlot, yname, ylab, oPosition, rotation, sizeLab, minOC, maxOC, lPos, uPos, font, sizeTL):
  if oPosition == 0:
    formatTickLabels = "%.3f"
  elif oPosition == 1:
    formatTickLabels = "%.1f"
  y_oc = ocPlot.Add('axis', name = yname, autoadd = False,
                    label = ylab,
                    min = minOC,
                    max = maxOC,
                    log = False,
                    direction = 'vertical',
                    mode = "numeric",
                    autoExtend = False,
                    autoExtendZero = False,
                    autoMirror = False,
                    otherPosition = oPosition,
                    Label__rotate = rotation,
                    Label__offset = '0pt',
                    Label__size = sizeLab,
                    Label__atEdge = True,
                    lowerPosition = lPos,
                    upperPosition = uPos,
                    Label__hide = False,
                    Label__font = font,
                    TickLabels__size = sizeTL,
                    TickLabels__font = font,
                    TickLabels__format = formatTickLabels,
                    MajorTicks__number = 5,
                    MinorTicks__number = 10,
                    )
  return y_oc

def yRES(plot, yname, ylab, oPosition, rotation, sizeLab, minOC, maxOC, lPos, uPos, font, sizeTL):
  if oPosition == 0:
    formatTickLabels = "%.3f"
  elif oPosition == 1:
    formatTickLabels = "%.1f"
  y_res = plot.Add('axis', name = yname, autoadd = False,
                    label = ylab,
                    min = minOC,
                    max = maxOC,
                    log = False,
                    direction = 'vertical',
                    mode = "numeric",
                    autoExtend = False,
                    autoExtendZero = False,
                    autoMirror = False,
                    otherPosition = oPosition,
                    Label__rotate = rotation,
                    Label__offset = '0pt',
                    Label__size = sizeLab,
                    Label__atEdge = True,
                    lowerPosition = lPos,
                    upperPosition = uPos,
                    Label__hide = False,
                    Label__font = font,
                    TickLabels__size = sizeTL,
                    TickLabels__font = font,
                    TickLabels__format = formatTickLabels,
                    MajorTicks__number = 5,
                    MinorTicks__hide = True,
                    )
  return y_res


# --------------------------------------------------------------------------


def main():
  #cli.full_path, s_sim, lmf, nbI, nbF, RVset, xscale, oc_fit = get_args()
  cli = get_args()

  # write to screen the intro
  intro(cli.full_path, cli.nbI, cli.nbF)

  # prepare stuff for plots
  #font = "Times New Roman"
  #font = "Utopia"
  font = 'URW Palladio L'
  markers = ["circle", "square", "triangle", "diamond", "star", "pus"]
  sizeO = '0.8pt' # size of Obsevation markers
  sizeS = '1.5pt' # size of Simulation markers
  sizeTL = '8pt' # size of ticklabels
  #sizeLab = '10pt'
  sizeLab = '8pt' # size of labels

  #gRows = (cli.nbF + 1) - cli.nbI + 2 # added another +1 to allow legend
  nbPlt = (cli.nbF - cli.nbI) + 1
  print " First body %d , last body %d ==> bodies to plot: %d" %(cli.nbI, cli.nbF, nbPlt)
  gRows = nbPlt * 2 + 1 # O-C + res x nbPlt + 1 legend
  #if cli.nbI == 0:
  if (cli.nbI <= 1):
    gRows = 1 # if you do not want to plot any O-C, only legend

  if cli.RVset != '0':
    RV, gamma, nRV_set, statRV = read_simRV(cli.full_path, cli.idsim, cli.lmflag)
    #read_obsRV(fullpath)
    if statRV:
      gRows = gRows + 2 # added +2 rows to allow RV plot

  print "Number of rows of the grid: %d" %(gRows)

  gCols = 1
  # define height of the page in cm
  hPage = "15cm"
  #if (cli.nbF-cli.nbI+1) > 3:
  if gRows > 7:
    hVal = (15. * float(gRows)/3.)
    if hVal > 23.:
      hVal = 23.
    hPage = "%4.2fcm"%(hVal)
  scRows = scalingRows(gRows)
  print 'scaling rows: ', scRows

  letters = ["b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"]

  ocWin = vsz.Embedded("O-C Plot")
  ocWin.ResizeWindow(800, 800)
  ocPage = ocWin.Root.Add('page',
                          height = hPage, width = hPage # squared page
                          )
  gridOC = ocPage.Add('grid', name = 'gridOC', autoadd = False,
                      rows = gRows,
                      columns = gCols,
                      scaleRows = scRows,
                      leftMargin = '0.2cm',
                      rightMargin = '0.2cm',
                      topMargin = '0.5cm',
                      bottomMargin = '0.1cm',
                      )

  # prepare margins for each plot
  # O-C plot
  loc = "1.3cm" # left margin
  roc = "1.3cm" # right margin
  toc = "0.0cm" # top margin
  boc = "0.0cm" # bottom margin
  # res plot
  lres = "1.3cm" # left margin
  rres = "1.3cm" # right margin
  tres = "0.1cm" # top margin
  bres = "1.cm" # bottom margin
  # define Min and Max Position of the plots from axis. Range [0.-1.]
  lPos = 0.03 # lowerPosition
  uPos = 1. - lPos # upperPosition

  if(cli.xscale is None):
    t_subtract=2440000.5
  else:
    t_subtract=np.float64(cli.xscale)
    
  xoc_list = [] # xaxis of O-C and res plots <-> T0

  if (cli.nbI > 1):
    min_max_time = np.zeros((nbPlt,2))
    for inb in range(cli.nbI,cli.nbF+1):
      #epo, oc, res, err = read_simT0(cli.full_path, cli.idsim, cli.lmflag, inb)
      epo, t_obs_base, oc, res, err = read_simT0(cli.full_path, cli.idsim, cli.lmflag, inb, t_subtract, cli.oc_fit)
      print ""

      keyObs = ""
      keySim = ""
      if (inb == cli.nbI):
        keyObs = "observations"
        keySim = "simulations"

      #minEpo = int(np.min(epo)) - 1
      #maxEpo = int(np.max(epo)) + 1
      min_time, max_time = lims(t_obs_base, t_obs_base)
      min_max_time[inb-cli.nbI,0] = min_time
      min_max_time[inb-cli.nbI,1] = max_time

      if(t_subtract != 0.):
        xlab = "BJD_{TDB} - %.1f" %(t_subtract)
      else:
        xlab = "BJD_{TDB}"

      # O-C

      #minOCd, maxOCd = lims(oc[:,0], oc[:,3])
      #minOCm, maxOCm = lims(oc[:,1], oc[:,4])

      min0, max0 = lims((oc[:,0]-err[:,0]), (oc[:,0]+err[:,0]))
      minOCd, maxOCd = lims((min0, max0), oc[:,3])
      minOCm = minOCd * 1440.
      maxOCm = maxOCd * 1440.

      ylab_l = "O-C [d]"
      ylab_r = "O-C [m]"

      ocname = "oc_%d_%s" %(inb, letters[inb-2])
      ocPlot = gridOC.Add('graph', name = ocname, autoadd = False,
                          leftMargin = loc,
                          rightMargin = roc,
                          topMargin = toc,
                          bottomMargin = boc,
                          )

      xname = 'x%d' %(inb)
      #x_oc = xOC(ocPlot, xname, xlab, minEpo, maxEpo, font, sizeLab, True, sizeTL)
      x_oc = xOC(ocPlot, xname, xlab, min_time, max_time, font, sizeLab, True, sizeTL)
      xoc_list.append(x_oc)

      yname = 'y%d_l' %(inb)
      y_l_oc = yOC(ocPlot, yname, ylab_l, 0, '0', sizeLab, minOCd, maxOCd, lPos, uPos, font, sizeTL)

      xyname = "xy_d_obs_%s" %(letters[inb-2])
      xy_d_obs = ocPlot.Add('xy', name = xyname, autoadd = False,
                            #xData = epo,
                            #yData = oc[:,0],
                            xAxis = xname,
                            yAxis = yname,
                            PlotLine__hide = True,
                            marker = markers[0],
                            markerSize = sizeO,
                            MarkerFill__color = 'black',
                            MarkerFill__hide = False,
                            MarkerLine__color = 'black',
                            errorStyle = 'bar',
                            ErrorBarLine__color = 'black',
                            ErrorBarLine__hideHorz = True,
                            FillBelow__color = 'black',
                            FillBelow__hide = True,
                            FillBelow__hideerror = False,
                            FillAbove__color = 'black',
                            key = keyObs,
                            )
      #epoName = 'epo_' + letters[inb-2]
      #ocWin.SetData(epoName, epo)
      time_name = 'xtime_%s' %(letters[inb-2])
      ocWin.SetData(time_name, t_obs_base)
      dataYname = "oc_d_obs_%s" %(letters[inb-2])
      ocWin.SetData(dataYname, oc[:,0], symerr = err[:,0])
      #xy_d_obs.xData.val = epoName
      xy_d_obs.xData.val = time_name
      xy_d_obs.yData.val = dataYname
      #xy_d_obs.yData.val = oc[:,0]

      if (inb == cli.nbI):
        legend = ocPlot.Add('key', name = 'legend', autoadd = False,
                            #Text__size = '12pt',
                            Text__font = font,
                            Text__size = sizeLab,
                            Border__hide = True,
                            horzPosn = 'centre',
                            vertPosn = 'manual',
                            keyLength = '0.35cm',
                            keyAlign = 'top',
                            #horzManual = 0,
                            vertManual = 1.001,
                            marginSize = 0.3,
                            columns = 2,
                            )

      xyname = "xy_d_sim_%s" %(letters[inb-2])
      xy_d_sim = ocPlot.Add('xy', name = xyname, autoadd = False,
                            #xData = epo,
                            #yData = oc[:,3],
                            xAxis = xname,
                            yAxis = yname,
                            PlotLine__hide = True,
                            marker = markers[0],
                            markerSize = sizeS,
                            MarkerFill__color = 'white',
                            MarkerFill__hide = False,
                            MarkerLine__color = 'blue',
                            key = keySim,
                            )
      #epoName = 'epo_' + letters[inb-2]
      #ocWin.SetData(epoName, epo)
      time_name = 'xtime_%s' %(letters[inb-2])
      ocWin.SetData(time_name, t_obs_base)
      dataYname = "oc_d_sim_%s" %(letters[inb-2])
      ocWin.SetData(dataYname, oc[:,3])
      #xy_d_sim.xData.val = epoName
      xy_d_sim.xData.val = time_name
      xy_d_sim.yData.val = dataYname

      zeroLine = ocPlot.Add('function', name = 'zeroLine', autoadd = False,
                            function = '0.',
                            Line__color = 'black',
                            Line__style = 'solid',
                            Line__width = '0.5pt',
                            xAxis = xname,
                            yAxis = yname,
                            )

      yname = 'y%d_r' %(inb)
      y_r_oc = yOC(ocPlot, yname, ylab_r, 1, '180', sizeLab, minOCm, maxOCm, lPos, uPos, font, sizeTL)

      # Residuals

      #minResd, maxResd = lims(res[:,0], res[:,0])
      #minResm, maxResm = lims(res[:,1], res[:,1])

      minResd, maxResd = lims((res[:,0]-err[:,0]), (res[:,0]+err[:,0]))
      minResm = minResd * 1440.
      maxResm = maxResd * 1440.

      ylab_l = "res [d]"
      ylab_r = "res [m]"
      #print xlab, ylab_l, ylab_r

      resname = "res_%d_%s" %(inb, letters[inb-2])
      resPlot = gridOC.Add('graph', name = resname, autoadd = False,
                            leftMargin = lres,
                            rightMargin = rres,
                            topMargin = tres,
                            bottomMargin = bres,
                            )

      xname = 'x%d' %(inb)
      #x_res = xOC(resPlot, xname, xlab, minEpo, maxEpo, font, sizeLab, False, sizeTL)
      x_res = xOC(resPlot, xname, xlab, min_time, max_time, font, sizeLab, False, sizeTL)
      xoc_list.append(x_res)

      yname = 'y%d_l' %(inb)
      y_l_res = yRES(resPlot, yname, ylab_l, 0., '0', sizeLab, minResd, maxResd, lPos, uPos, font, sizeTL)

      xyname = "res_d_%s" %(letters[inb-2])
      res_d = resPlot.Add('xy', name = xyname, autoadd = False,
                          xAxis = xname,
                          yAxis = yname,
                          PlotLine__hide = True,
                          marker = markers[0],
                          markerSize = sizeS,
                          MarkerFill__color = 'blue',
                          MarkerFill__hide = True,
                          MarkerLine__color = 'blue',
                          errorStyle = 'bar',
                          ErrorBarLine__color = 'blue',
                          ErrorBarLine__hideHorz = True,
                          FillBelow__color = 'blue',
                          FillBelow__hide = True,
                          FillBelow__hideerror = False,
                          FillAbove__color = 'blue',
                          )
      #epoName = 'epo_' + letters[inb-2]
      #ocWin.SetData(epoName, epo)
      time_name = 'xtime_%s' %(letters[inb-2])
      ocWin.SetData(time_name, t_obs_base)
      dataYname = "res_d_%s" %(letters[inb-2])
      ocWin.SetData(dataYname, res[:,0], symerr = err[:,0])
      #res_d.xData.val = epoName
      res_d.xData.val = time_name
      res_d.yData.val = dataYname


      zeroLine = resPlot.Add('function', name = 'zeroLine', autoadd = False,
                              function = '0.',
                              Line__color = 'black',
                              Line__style = 'solid',
                              Line__width = '0.5pt',
                              xAxis = xname,
                              yAxis = yname,
                              )

      yname = 'y%d_r' %(inb)
      y_r_res = yRES(resPlot, yname, ylab_r, 1, '180', sizeLab, minResm, maxResm, lPos, uPos, font, sizeTL)

      # print T0 residuals rms
      rms_oc = np.sqrt(np.mean(oc*oc, axis=0))
      rms_resT0 = np.sqrt(np.mean(res*res, axis=0))
      print 'rms(OC)obs = %.8f days = %.6f min = %.4f sec' %(rms_oc[0], rms_oc[1], rms_oc[2])
      print 'rms(OC)sim = %.8f days = %.6f min = %.4f sec' %(rms_oc[3], rms_oc[4], rms_oc[5])
      print 'rms(T0)    = %.8f days = %.6f min = %.4f sec' %(rms_resT0[0], rms_resT0[1], rms_resT0[2])
      print ''

    # set properly the xaxis of all the O-C and res plots
    min_time = float(np.min(min_max_time))
    max_time = float(np.max(min_max_time))
    print 'min time = %.5f max time = %.5f' %(min_time, max_time)
    for ii in range(len(xoc_list)):
      xoc_list[ii].min.val=min_time
      xoc_list[ii].max.val=max_time


  #if cli.RVset != '0':
  if (cli.RVset):
    if statRV:
      #print " Velocity of the system "
      #print " RV gamma = %.7f +/- %.7f" %(gamma[0], gamma[1])
      RVbkc = RV.copy()
      RV[:,1] = RV[:,1] - gamma[:,0]
      RV[:,4] = RV[:,4] - gamma[:,0]
      min0, max0 = lims((RV[:,1]-RV[:,2]), (RV[:,1]+RV[:,2]))
      minRV, maxRV = lims( (min0, max0), RV[:,4] )
      # insert model RV from #_rotorbit.dat file
      jd_lte, RVmod, statMod = read_orbit(cli.full_path, cli.idsim, cli.lmflag)
      jd_lte_bck = jd_lte.copy()
      jd_lte = jd_lte - t_subtract
      jd = RV[:,0] - t_subtract
      #minJD, maxJD = lims(jd, jd_lte)
      minJD, maxJD = lims(jd, jd)
      print " limits of BJD = ", minJD, maxJD
      minRV, maxRV = lims( (minRV, maxRV), RVmod)
      #print minRV, maxRV
      dRV = RV[:,1] - RV[:,4]
      RV_rms = np.sqrt(np.mean(dRV*dRV))
      print ' RV summary '
      print ' RMS( RV_obs - RV_sim ) = %.7f m/s' %(RV_rms)
      print_gamma = np.zeros((nRV_set,2)) # gamma, err_gamma
      print_gamma[0,:] = RV[0,5:7]
      for i in range(1,nRV_set):
        for j in range(1, RV.shape[0]):
          if (RV[j,5] != RV[j-1,5]):
            #print i,j,RV[j,5],RV[j-1,5]
            print_gamma[i,:] = RV[j,5:7]
      print 'gamma'.rjust(15),'+/-', 'err_gamma'.rjust(15)
      for i in range(0, nRV_set):
        print '%15.6f +/- %15.6f' %(print_gamma[i,0], print_gamma[i,1])

      # RVplot
      RVname = 'RVPlot'
      RVPlot = gridOC.Add('graph', name = RVname, autoadd = False,
                          leftMargin = loc,
                          rightMargin = roc,
                          topMargin = toc,
                          bottomMargin = boc,
                          )

      if(t_subtract != 0.):
        xlab = "BJD_{TDB} - %.1f" %(t_subtract)
      else:
        xlab = "BJD_{TDB}"
      xname = 'x_rv'
      x_RV = xRV(RVPlot, xname, xlab, minJD, maxJD, font, sizeLab, True, sizeTL)

      ylab = "RV [m/s]"
      yname = 'y_rv'
      y_RV = RVPlot.Add('axis', name = yname, autoadd = False,
                        label = ylab,
                        direction = 'vertical',
                        mode = "numeric",
                        autoExtend = False,
                        autoExtendZero = False,
                        otherPosition = 0,
                        Label__rotate = 0.,
                        Label__offset = '0pt',
                        log = False,
                        #min = 'Auto',
                        #max = 'Auto',
                        min = minRV,
                        max = maxRV,
                        lowerPosition = lPos,
                        upperPosition = uPos,
                        autoMirror = True,
                        Label__hide = False,
                        Label__font = font,
                        Label__size = sizeLab,
                        Label__atEdge = True,
                        TickLabels__size = sizeTL,
                        TickLabels__font = font,
                        MajorTicks__number = 5,
                        TickLabels__format = "%.2f",
                        )

      # plot with different colors and markers
      #xyname = 'obsRV'
      #xy_RVobs = RVPlot.Add('xy', name = xyname, autoadd = False,
                #xAxis = xname,
                #yAxis = yname,
                #PlotLine__hide = True,
                #marker = markers[0],
                #markerSize = sizeO,
                #MarkerFill__color = 'black',
                #MarkerFill__hide = False,
                #MarkerLine__color = 'black',
                #errorStyle = 'bar',
                #ErrorBarLine__color = 'black',
                #ErrorBarLine__hideHorz = True,
                #FillBelow__color = 'black',
                #FillBelow__hide = True,
                #FillBelow__hideerror = False,
                #FillAbove__color = 'black',
                #)
      #jdName = "jdRV"
      #ocWin.SetData(jdName, jd)
      #dataYname = "RVo"
      #ocWin.SetData(dataYname, RV[:,1], symerr = RV[:,2])
      #xy_RVobs.xData.val = jdName
      #xy_RVobs.yData.val = dataYname
      
      rv_colors = ['red', 'green', 'black', 'blue', 'yellow', 'cyan']
      rv_obs_size = '1.25pt'
      for i_set in range(0, nRV_set):
        jd_sel   = jd[ RV[:,7] == i_set+1   ]
        RVo_sel  = RV[ RV[:,7] == i_set+1,1 ]
        eRVo_sel = RV[ RV[:,7] == i_set+1,2 ]
        xyname = 'obsRV_%03d' %(i_set+1)
        label_key = "RVset#%02d" %(i_set+1)
        xy_RVobs = RVPlot.Add('xy', name = xyname, autoadd = False,
                              xAxis = xname,
                              yAxis = yname,
                              PlotLine__hide = True,
                              marker = markers[i_set],
                              markerSize = rv_obs_size,
                              MarkerFill__color = rv_colors[i_set],
                              MarkerFill__hide = False,
                              MarkerLine__color = 'black',
                              errorStyle = 'bar',
                              ErrorBarLine__color = 'black',
                              ErrorBarLine__hideHorz = True,
                              FillBelow__color = 'black',
                              FillBelow__hide = True,
                              FillBelow__hideerror = False,
                              FillAbove__color = 'black',
                              key = label_key
                              )
        jdName = "jdRV_%03d" %(i_set+1)
        ocWin.SetData(jdName, jd_sel)
        dataYname = "RVo_%03d" %(i_set+1)
        ocWin.SetData(dataYname, RVo_sel, symerr = eRVo_sel)
        xy_RVobs.xData.val = jdName
        xy_RVobs.yData.val = dataYname

      xyname = 'simRV'
      xy_RVsim = RVPlot.Add('xy', name = xyname, autoadd = False,
                            xAxis = xname,
                            yAxis = yname,
                            PlotLine__hide = True,
                            marker = markers[0],
                            markerSize = sizeS,
                            MarkerFill__color = 'white',
                            MarkerFill__hide = False,
                            MarkerLine__color = 'blue',
                            )
      jdName = "jdRV"
      ocWin.SetData(jdName, jd)
      dataYname = "RVs"
      ocWin.SetData(dataYname, RV[:,4])
      xy_RVsim.xData.val = jdName
      xy_RVsim.yData.val = dataYname


      if (statMod):
        jdmod = jd_lte# - 2454000.
        xyname = 'modRV'
        xy_modRV = RVPlot.Add('xy', name = xyname, autoadd = False,
                              xAxis = xname,
                              yAxis = yname,
                              PlotLine__hide = False,
                              #PlotLine__color = 'blue',
                              #PlotLine__color = '#33ADFF',
                              PlotLine__color = '#B2B2B2',
                              PlotLine__width = '1pt',
                              PlotLine__style = 'dot1',
                              marker = 'none',
                              key = "RV model"
                              )
        jdName = "jdMod"
        ocWin.SetData(jdName, jdmod)
        dataYname = "RVmod"
        ocWin.SetData(dataYname, RVmod)
        xy_modRV.xData.val = jdName
        xy_modRV.yData.val = dataYname

        legend = RVPlot.Add('key', name = 'legendRV', autoadd = False,
                            #Text__font ='URW Palladio L',
                            #Text__size = '8pt',
                            Text__font = font,
                            Text__size = '6pt',
                            Background__color = 'white',
                            Background__hide = True,
                            Border__hide = True,
                            horzPosn = 'manual',
                            vertPosn = 'manual',
                            keyLength = '0.2cm',
                            keyAlign = 'centre',
                            horzManual = 1.0025,
                            #vertManual = -0.7,
                            vertManual = 0.8,
                            marginSize = 0.,
                            columns = 1,
                            symbolswap = True,
                            )

      else:
        print "I did not find the orbit file to plot RV model. Skipping."

      zeroLine = RVPlot.Add('function', name = 'zeroLine', autoadd = False,
                            function = '0.',
                            Line__color = 'black',
                            Line__style = 'solid',
                            Line__width = '0.5pt',
                            xAxis = xname,
                            yAxis = yname,
                            )

      # RVresiduals

      resRV = RV[:,1] - RV[:,4]
      minRVres, maxRVres = lims( (resRV-RV[:,2]), (resRV+RV[:,2]))

      resname = "res_RV"
      resRVplot = gridOC.Add('graph', name = resname, autoadd = False,
                              leftMargin = lres,
                              rightMargin = rres,
                              topMargin = tres,
                              bottomMargin = bres,
                              )

      xname = 'x_RVres'
      x_RVres = xOC(resRVplot, xname, xlab, minJD, maxJD, font, sizeLab, False, sizeTL)

      ylab = "res [m/s]"
      yname = 'y_RVres'
      y_RVres = resRVplot.Add('axis', name = yname, autoadd = False,
                              label = ylab,
                              direction = 'vertical',
                              mode = "numeric",
                              autoExtend = False,
                              autoExtendZero = False,
                              otherPosition = 0,
                              Label__rotate = 0.,
                              Label__offset = '0pt',
                              log = False,
                              min = minRVres,
                              max = maxRVres,
                              lowerPosition = lPos,
                              upperPosition = uPos,
                              autoMirror = True,
                              Label__hide = False,
                              Label__font = font,
                              Label__size = sizeLab,
                              Label__atEdge = True,
                              TickLabels__size = sizeTL,
                              TickLabels__format = "%.2f",
                              MajorTicks__number = 3,
                              MinorTicks__hide = True,
                              )

      xyname = "res_RV"
      res_RV = resRVplot.Add('xy', name = xyname, autoadd = False,
                              xAxis = xname,
                              yAxis = yname,
                              PlotLine__hide = True,
                              marker = markers[0],
                              markerSize = sizeS,
                              MarkerFill__color = 'blue',
                              MarkerFill__hide = True,
                              MarkerLine__color = 'blue',
                              errorStyle = 'bar',
                              ErrorBarLine__color = 'blue',
                              ErrorBarLine__hideHorz = True,
                              FillBelow__color = 'blue',
                              FillBelow__hide = True,
                              FillBelow__hideerror = False,
                              FillAbove__color = 'blue',
                              )
      jdName = "jdRV"
      ocWin.SetData(jdName, jd)
      dataYname = "res_RV"
      ocWin.SetData(dataYname, resRV, symerr = RV[:,2])
      res_RV.xData.val = jdName
      res_RV.yData.val = dataYname

      zeroLine = resRVplot.Add('function', name = 'zeroLine', autoadd = False,
                                function = '0.',
                                Line__color = 'black',
                                Line__style = 'solid',
                                Line__width = '0.5pt',
                                xAxis = xname,
                                yAxis = yname,
                                )



  # ---------
  # save plot
  s_nbIF = "_%d-%d" %(cli.nbI,cli.nbF)
  plotFolder = os.path.join(cli.full_path, "plots")
  if not os.path.isdir(plotFolder):
    createDir = "mkdir -p " + plotFolder
    os.system(createDir)
  fveusz = os.path.join(plotFolder, '%d_%d_oc%s' %(cli.idsim, cli.lmflag, s_nbIF))
  ocWin.Save(fveusz + ".vsz")
  ocWin.Export(fveusz + ".eps")
  ocWin.Export(fveusz + ".pdf")
  ocWin.Export(fveusz + "_150dpi.png", dpi=150, antialias=True, quality=100, backcolor='white')
  ocWin.Export(fveusz + "_300dpi.png", dpi=300, antialias=True, quality=100, backcolor='white')
  ocWin.Export(fveusz + ".svg")
  ocWin.Close()


if __name__ == '__main__':
  main()
