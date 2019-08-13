#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import sys
import argparse
import os
import numpy as np # array
import time
#import h5py
#import random
#import constants as cst # local constants module
#from scipy.stats import norm as scipy_norm
import ancillary as anc
from matplotlib import use as mpluse
mpluse("Agg")
#mpluse("Qt4Agg")
import matplotlib.pyplot as plt
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
#from matplotlib import rcParams
#rcParams['text.latex.unicode']=True

import logging
import warnings

def main():
  # ---
  # initialize logger

  logger = logging.getLogger("Main_log")
  logger.setLevel(logging.DEBUG)
  formatter = logging.Formatter("%(asctime)s - %(message)s")

  # read cli arguments
  cli = anc.get_args()

  #plot_folder = prepare_plot_folder(working_path)
  emcee_plots = anc.prepare_emcee_plot_folder(cli.full_path)
  log_file = os.path.join(emcee_plots, 'GelmanRubin_log.txt')
  flog = logging.FileHandler(log_file, 'w')
  flog.setLevel(logging.DEBUG)
  flog.setFormatter(formatter)
  logger.addHandler(flog)
  # log screen
  slog = logging.StreamHandler()
  slog.setLevel(logging.DEBUG)
  slog.setFormatter(formatter)
  logger.addHandler(slog)

  # computes mass conversion factor
  #m_factor = anc.mass_conversion_factor(cli.m_type)
  m_factor, m_unit = anc.mass_type_factor(1., cli.m_type, False)

  # set emcee and trades folder
  emcee_folder = cli.full_path
  trades_folder = os.path.join(os.path.dirname(cli.full_path), '')
  # and best folder
  emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(emcee_folder, cli.temp_status)

  # get data from the hdf5 file
  parameter_names_emcee, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)
  # print Memory occupation of ...
  anc.print_memory_usage(chains)

  nfit, nwalkers, nruns, nburnin, nruns_sel = anc.get_emcee_parameters(chains, cli.temp_status, cli.nburnin, completed_steps)
  logger.info('nfit(%d), nwalkers(%d), nruns(%d), nburnin(%d), nruns_sel(%d)' %(nfit, nwalkers, nruns, nburnin, nruns_sel))

  # set label and legend names
  kel_labels = anc.keplerian_legend(parameter_names_emcee, cli.m_type)

  chains_T, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, nburnin, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries, chains)

  if(cli.temp_status):
    n_steps = completed_steps
  else:
    n_steps = nruns
  sel_steps = int(cli.sel_steps)
  if (sel_steps == 0):
    sel_steps = n_steps

  steps = np.linspace(start=0, stop=n_steps, num=sel_steps, endpoint=True, dtype=np.int)
  steps[0] = 10
  #print steps
  sel_steps = steps.shape[0]

  gr_Rc_2 = np.ones((sel_steps, nfit)) + 99.

  for ifit in range(0, nfit):
    logger.info('Parameter: %13s' %(parameter_names_emcee[ifit]))
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot2grid((1, 1), (0,0))
    
    for istep in range(0,sel_steps):
      
      time0 = time.time()
      gr_Rc_2[istep,ifit] = anc.GelmanRubin(chains_T[:steps[istep], :, ifit])
      if(istep == sel_steps-1):
        LBo_d, LBo_h, LBo_m, LBo_s = anc.computation_time(time.time()-time0)
        logger.info('steps = %6d for %13s ==> Gelman-Rubin test:  LBo time = %2d m %6.3f s' %(steps[istep], parameter_names_emcee[ifit], LBo_m, LBo_s))
      
    ax.axhline(1.01, color='gray')
    ax.plot(steps, gr_Rc_2[:,ifit], '-', color='k', lw=1.3, label='LBo 2')
    ax.set_ylim(0.95, 2.3)
    ax.set_xlabel('steps (%s)' %(parameter_names_emcee[ifit].strip()))
    ax.legend(loc='center left', fontsize=9, bbox_to_anchor=(1, 0.5))
    fig.savefig(os.path.join(emcee_plots, 'GR_%03d_%s.png' %(ifit+1, parameter_names_emcee[ifit])), bbox_inches='tight', dpi=200)
    plt.close(fig)
    logger.info('saved plot %s' %(os.path.join(emcee_plots, 'GRtrace_pam_%s.png' %(parameter_names_emcee[ifit]))))
    
  logger.info('')
  
  return

if __name__ == "__main__":
  main()

