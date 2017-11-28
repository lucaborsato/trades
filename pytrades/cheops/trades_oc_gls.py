#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# IMPORT

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import numpy as np # array
import os
import sys
import glob
import astropy.time as atime

import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt

# custom modules
script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), '../'))
sys.path.append(module_path)

import gls


# ==============================================================================

def get_args():
  
  parser = argparse.ArgumentParser(description='CHEOPS TTV LEVEL ANALYSIS')
  
  parser.add_argument('-f', '--f', '-folder-path','--folder-path',  
                      action='store', dest='folder_path', required=True, 
                      help='The path (absolute or relative) with simulation files for TRADES.')
  
  parser.add_argument('-s', '--s', '-sim-id','--sim-id',  
                      action='store', dest='sim_id', default='1', 
                      help='Number ID of the simulation. Default is 1.')
  
  parser.add_argument('-b', '--b', '-body-id','--body-id',  
                      action='store', dest='body_id', default='2', 
                      help='Number ID of the body: star is 1, planets from 2 to NB. Default is 2.')
  
  parser.add_argument('-t', '--t', '-t-ref','--t-ref',  
                      action='store', dest='tref', default=None, 
                      help='Transit time of reference for/from linear ephemeris. Default is None, and it will be computed from a linear ephemeris with guess the first transit time.')
  
  parser.add_argument('-p', '--p', '-p-ref','--p-ref',  
                      action='store', dest='pref', default=None, 
                      help='Period (in days) of reference for/from linear ephemeris. Default is None, and it will computed from the median of the differences between each transit time.')
  
  cli = parser.parse_args()
  
  cli.folder_path = os.path.abspath(cli.folder_path)
  
  try:
    cli.sim_id = int(cli.sim_id)
    if(cli.sim_id <= 0): cli.sim_id = 1
  except:
    cli.sim_id = 1
  
  try:
    cli.body_id = int(cli.body_id)
    if(cli.body_id <= 1): cli.body_id = 2
  except:
    cli.body_id = 2
  
  try:
    cli.tref = np.float64(cli.tref)
  except:
    cli.tref = None
  
  try:
    cli.pref = np.float64(cli.pref)
  except:
    cli.pref = None
  
  return cli

# ==============================================================================

def printlog(line, log=None):
  
  print '%s' %(line)
  if(log is not None):
    log.write('%s\n' %(line))
  
  return

# ==============================================================================

def read_file(full_file):
  
  ttra_lte = np.genfromtxt(full_file, usecols=(0,1))
  ttra = ttra_lte[:,0] + ttra_lte[:,1]
  nttra = np.shape(ttra)[0]
  if(nttra == 0): ttra = None
  
  return ttra, nttra

# ==============================================================================

# Linear fit function with errors on y-axis.
# It returns also the errors.
# y = b_ls + m_ls * x

# ==============================================================================

def lstsq_fit(x, y, yerr):
  
  A = np.vstack((np.ones_like(x), x)).T
  C = np.diag(yerr * yerr)
  cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
  b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))
  coeff_err = np.sqrt(np.diag(cov))
  
  return b_ls, m_ls, coeff_err

# ==============================================================================

def refine_linephem(ttra, tref0, pref0):
  
  nttra = np.shape(ttra)[0]
  tref, pref = tref0, pref0
  if(tref is None):
    tref = ttra[0]
  if(pref is None):
    pref = np.median([ttra[ii+1] - ttra[ii] for ii in range(nttra-1)])
  epo = np.rint((ttra-tref)/pref).astype(int)
  tref, pref, _ = lstsq_fit(epo, ttra, np.zeros((nttra))+1./86400.)
  
  return tref, pref, epo

# ==============================================================================

class Transit:
  
  def __init__(self, full_file, sim_id, body_id, tref, pref, tstart, tend):
    self.full_file = full_file
    self.dirname = os.path.dirname(full_file)
    self.filename = os.path.basename(full_file)
    self.sim_id = sim_id
    self.body_id = body_id
    self.ttra, self.nttra = read_file(self.full_file)
    self.tref = tref
    self.pref = pref
    self.ttra_sel = self.ttra[np.logical_and(self.ttra >= tstart, self.ttra <= tend)]
    self.nttra_sel = np.shape(self.ttra_sel)[0]

# ==============================================================================

def main():
  
  # get folder name (and sim id?)
  cli = get_args()
  
  # determines number of files with transit times
  full_file = os.path.join(cli.folder_path, '%d_0_NB%d_tra.dat' %(cli.sim_id, cli.body_id))
  target = os.path.basename(cli.folder_path).upper()
  
  olog = open(os.path.join(os.path.dirname(full_file), 'log_oc_%d_NB%d.txt' %(cli.sim_id, cli.body_id)), 
              'w')
  
  printlog('\n%s\nREAD FILE: %s' %(target, full_file), olog)
  
  # cheops launch date
  #t_launch = atime.Time('2018-06-30T06:00:00.0', format='isot', scale='tdb')
  t_launch = atime.Time('2019-01-01T00:00:00.0', format='isot', scale='tdb')
  # cheops end nominal mission
  t_end = t_launch.jd + 365.25*3.5
  cheops_dur = t_end - t_launch.jd
  
  printlog('TIME SELECTION BJD: %.6f - %.6f' %(t_launch.jd, t_end), olog)
  
  # find transit times in file
  transits = Transit(full_file, cli.sim_id, cli.body_id, 
                     cli.tref, cli.pref, 
                     t_launch.jd, t_end)
  
  #transits.tref, transits.pref, epo_sel = refine_linephem(transits.ttra_sel,
                                                 #transits.tref, transits.pref)
  epo_sel, transits.tref, transits.pref = anc.compute_lin_ephem(transits.ttra_sel)
  printlog('LINEAR EPHEMERIS: %.6f + N x %.6f' %(transits.tref, transits.pref), olog)
  printlog('FOUND %d TRANSIT DURING CHEOPS ON A TOTAL OF %d' %(transits.nttra_sel, transits.nttra), olog)
  
  tlin = transits.tref + epo_sel * transits.pref
  oc_d = transits.ttra_sel - tlin
  eoc_d = np.zeros((transits.nttra_sel))+1./86400.
  
  tt = transits.ttra_sel - transits.tref
  gls_o = gls.Gls((tt, oc_d, eoc_d), Pbeg=3.*transits.pref, Pend=2.*cheops_dur, ofac=10,
                  verbose = False)
  gls_o.info()
  amp_ttv = gls_o.hpstat['amp']
  p_ttv = 1./gls_o.hpstat['fbest']
  printlog('amp_ttv = %.6f d P_ttv = %.6f' %(amp_ttv, p_ttv), olog)
  
  oc_file = os.path.join(cli.folder_path, 'oc_%d_NB%d.dat' %(cli.sim_id, cli.body_id))
  
  np.savetxt(oc_file, np.column_stack((epo_sel, transits.ttra_sel, oc_d)), 
                      fmt='%5d %23.16e %23.16e', 
                      header='Tlin = %23.16e + N x %23.16e\namp_TTV = %.12f min P_TTV = %.12f d\nepoch TT_d OC_d' 
                              %(transits.tref,
                                transits.pref,
                                amp_ttv*1440.,
                                p_ttv)
            )
  
  plot_folder = os.path.join(os.path.dirname(full_file), 'plots')
  if(not os.path.isdir(plot_folder)):
    os.makedirs(plot_folder)
  plot_file = os.path.join(plot_folder, 'oc_%d_NB%d' %(cli.sim_id, cli.body_id))
  
  printlog('SAVING O-C PLOT %s.png\n' %(plot_file), olog)
  
  fig = plt.figure(figsize=(12,12))

  plt.axhline(0., color='black', ls='-', lw=0.33, alpha=0.77)
  plt.plot(tt, oc_d*1440., color='lightgray', marker='o', ms=3, mfc='#1f77b4', ls='-', lw=0.45)

  plt.title('%s (O-C) body %d' %(target, cli.body_id,), loc='left')
  plt.title('Tlin = %.6f + N x %.6f' %(transits.tref, transits.pref), loc='center', fontsize=8)
  plt.title(r'$\sigma(TTV)=%.3f$ min P(TTV)$=%.4f$ d' %(amp_ttv*1440., p_ttv), loc='right', fontsize=8)

  plt.xlabel('BJD - %.4f' %(transits.tref))
  plt.ylabel('O-C (min)')
  
  fig.savefig('%s.png' %(plot_file), bbox_inches='tight', dpi=300)
  plt.close(fig)
  
  olog.close()
  
  return

# ==============================================================================
# ==============================================================================
if __name__ == "__main__":
  main()


