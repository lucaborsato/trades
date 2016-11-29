#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import sys
import argparse
import os
import numpy as np # array
import h5py
import random
import constants as cst # local constants module
from scipy.stats import norm as scipy_norm
import ancillary as anc
from matplotlib import use as mpluse
mpluse("Agg")
#mpluse("Qt4Agg")
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)
#from matplotlib import rcParams
#rcParams['text.latex.unicode']=True
#import corner


def main():

  print 
  print ' ======================== '
  print ' TRADES+EMCEE CHAIN PLOTS'
  print ' ======================== '
  print

  # read cli arguments
  cli = anc.get_args()
  # computes mass conversion factor
  #m_factor, m_unit = anc.mass_conversion_factor_and_unit(cli.m_type)
  m_factor, m_unit = anc.mass_type_factor(1., cli.m_type, False)

  # set emcee and trades folder
  emcee_folder = cli.full_path
  trades_folder = os.path.join(os.path.dirname(cli.full_path), '')
  # and best folder
  emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(emcee_folder, cli.temp_status)

  parameter_names_emcee, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)

  # set label and legend names
  kel_labels = anc.keplerian_legend(parameter_names_emcee, cli.m_type)

  nfit, nwalkers, nruns, npost, nruns_sel = anc.get_emcee_parameters(chains, cli.temp_status, cli.npost, completed_steps)

  anc.print_memory_usage(chains)

  chains_T, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries, chains)

  try:
    n_acor = autocor_time.shape[0]
  except:
    n_acor = len(autocor_time)
  if(n_acor == 0):
    #autocor_time = anc.compute_autocor_time(flatchain_posterior_0)
    autocor_time = anc.compute_autocor_time(chains_T)
  print autocor_time
  thin_steps = np.rint(np.mean(np.array(autocor_time, dtype=np.float64))).astype(int)
  print 'thin_steps = ',thin_steps
  #sys.exit()


  # create a flat array of the posterior: from (nruns_sel, nwalkers, nfit) -> (nruns_sel * nwalkers, nfit)
  flatchain_posterior_0 = chains_T[:,:,:].reshape((-1, nfit))
  flatchain_posterior_1 = flatchain_posterior_0

  if(cli.boot_id > 0):
    flatchain_posterior_msun = anc.posterior_back_to_msun(m_factor,parameter_names_emcee,flatchain_posterior_0)
    boot_file = anc.save_bootstrap_like(emcee_folder, cli.boot_id, parameter_names_emcee, flatchain_posterior_msun)
    logger.info('saved bootstrap like file: %s' %(boot_file))
    del flatchain_posterior_msun

  nburnin=npost

  median_parameters, median_perc68, median_confint = anc.get_median_parameters(flatchain_posterior_0)

  k = np.ceil(2. * flatchain_posterior_0.shape[0]**(1./3.)).astype(int)
  #if(k>11): k=11
  if(k>50): k=50
  mode_bin, mode_parameters, mode_perc68, mode_confint = anc.get_mode_parameters_full(flatchain_posterior_0, k)
  print 

  emcee_plots = os.path.join(cli.full_path,'plots')
  if (not os.path.isdir(emcee_plots)):
    os.makedirs(emcee_plots)

  for i in range(0, nfit):
    emcee_fig_file = os.path.join(emcee_plots, 'chain_%s.png' %(parameter_names_emcee[i].strip()))
    print ' %s' %(emcee_fig_file),
    fig, (axChain, axHist) = plt.subplots(nrows=1, ncols=2, figsize=(12,12))

    
    (counts, bins_val, patches) = axHist.hist(flatchain_posterior_1[:,i], bins=k, range=(flatchain_posterior_1[:,i].min(), flatchain_posterior_1[:,i].max()), orientation='horizontal', normed=True, stacked=True, histtype='stepfilled', color='darkgrey', edgecolor='lightgray', align='mid')

    xpdf = scipy_norm.pdf(flatchain_posterior_1[:,i], loc = flatchain_posterior_1[:,i].mean(), scale = flatchain_posterior_1[:,i].std())
    idx = np.argsort(flatchain_posterior_1[:,i])
    axHist.plot(xpdf[idx], flatchain_posterior_1[idx,i], color='black', marker='None', ls='-.', lw=1.5, label='pdf')

    axChain.plot(chains_T[:,:,i], '-', alpha=0.3)

    post_15th = np.percentile(flatchain_posterior_1[:,i], anc.percentile_val[2], interpolation='midpoint')
    post_84th = np.percentile(flatchain_posterior_1[:,i], anc.percentile_val[3], interpolation='midpoint')

    median  = median_parameters[i]
    lower = median_parameters[i]+median_confint[0,i]
    upper = median_parameters[i]+median_confint[3,i]

    # plot of mode (mean of higher peak/bin)
    axChain.axhline(mode_parameters[i], color='red', ls='-', lw=2.1, alpha=1, label='mode')
    axChain.axhline(mode_parameters[i]+mode_confint[0,i], color='red', ls='--', lw=1.3, alpha=0.5, label='lower mode')
    axChain.axhline(mode_parameters[i]+mode_confint[3,i], color='red', ls='-.', lw=1.3, alpha=0.5, label='upper mode')
    
    # plot of median
    #axChain.axhline(max_lnprob_parameters[i], marker='None', c='cyan',ls='-', lw=1.5, label='max lnprob fit', alpha=0.65)
    axChain.axhline(median, marker='None', c='blue',ls='-', lw=2.1, alpha=1.0, label='median fit')
    axChain.axhline(lower, marker='None', c='blue',ls='-.', lw=1.1, alpha=0.5, label='lower median')
    axChain.axhline(upper, marker='None', c='blue',ls='--', lw=1.1, alpha=0.5, label='upper median')
    
    
    axChain.ticklabel_format(useOffset=False)
    axChain.set_xlabel('$N_\mathrm{steps}$')
    axChain.set_ylabel(kel_labels[i])
    
    y_min, y_max = anc.compute_limits(flatchain_posterior_1[:,i], 0.05)
    if(y_min == y_max):
      y_min = parameter_boundaries[i,0]
      y_max = parameter_boundaries[i,1]
    
    axChain.set_ylim([y_min, y_max])
    axChain.set_title('Full chain %s:=[%.3f , %.3f]' %(kel_labels[i], parameter_boundaries[i,0], parameter_boundaries[i,1]))
    plt.draw()

    axHist.ticklabel_format(useOffset=False)
    axHist.set_ylim([y_min, y_max])

    # plot mode
    axHist.axhline(mode_parameters[i], color='red', ls='-', lw=2.1, alpha=1, label='mode')
    axHist.axhline(mode_parameters[i]+mode_confint[0,i], color='red', ls='--', lw=1.3, alpha=0.5, label='lower mode')
    axHist.axhline(mode_parameters[i]+mode_confint[3,i], color='red', ls='-.', lw=1.3, alpha=0.5, label='upper mode')
    # plot median
    axHist.axhline(median, marker='None', c='blue',ls='-', lw=2.1, alpha=1.0, label='median fit')
    axHist.axhline(lower, marker='None', c='blue',ls='-.', lw=1.1, alpha=0.5, label='lower median')
    axHist.axhline(upper, marker='None', c='blue',ls='--', lw=1.1, alpha=0.5, label='upper median')
    # plot CI
    axHist.axhline(post_15th, marker='None', color='lightgreen', ls='-', lw=1.5, alpha=0.4, label='posterior[15.86th]')
    axHist.axhline(post_84th, marker='None', color='darkgreen', ls='-', lw=1.5, alpha=0.4, label='posterior[84.14th]')
    
    axHist.set_title('Distribution of posterior chain')
    axHist.legend(loc='center left', fontsize=9, bbox_to_anchor=(1, 0.5))
    plt.draw()

    fig.savefig(emcee_fig_file, bbox_inches='tight', dpi=150)
    print 'saved'
    print '   median = ', median,             ' lower = ', lower,                                ' upper = ', upper
    print '     mode = ', mode_parameters[i], ' lower = ', mode_parameters[i]+mode_confint[0,i], ' upper = ', mode_parameters[i]+mode_confint[3,i]
    print 'post 15th = ', post_15th,          '  84th = ', post_84th
    print

  if(cli.temp_status):
    lnprob_burnin = lnprobability[:,nburnin:completed_steps]
  else:
    lnprob_burnin = lnprobability[:,nburnin:]
  fig = plt.figure(figsize=(12,12))
  plt.plot(lnprob_burnin.T, '-', alpha=0.8)
  plt.xlabel('$N_\mathrm{steps}$')
  plt.ylabel('logprob')
  min_lnp = np.min(lnprob_burnin.T, axis=0).min()
  max_lnp = np.max(lnprob_burnin.T, axis=0).max()
  y_min, y_max = anc.compute_limits(np.asarray([min_lnp, max_lnp]), 0.05)
  plt.ylim((y_min, y_max))
  plt.draw()
  fig.savefig(os.path.join(emcee_plots, 'emcee_lnprobability.png'), bbox_inches='tight', dpi=150)
  print ' %s saved' %(os.path.join(emcee_plots, 'emcee_lnprobability.png'))

  return

if __name__ == "__main__":
  main()




