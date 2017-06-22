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

  chains_T_full, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries, chains)
  
  if(cli.use_thin):
    chains_T, flatchain_posterior_0, lnprob_burnin, thin_steps, chains_T_full_thinned = anc.thin_the_chains(cli.use_thin, npost, nruns, nruns_sel, autocor_time, chains_T_full, lnprobability, burnin_done=False, full_chains_thinned=True)
    npost_plt = np.rint(npost / thin_steps).astype(int)
    nend = np.rint(nruns / thin_steps).astype(int)
  
  else:
    chains_T, flatchain_posterior_0, lnprob_burnin, thin_steps = anc.thin_the_chains(cli.use_thin, npost, nruns, nruns_sel, autocor_time, chains_T_full, lnprobability, burnin_done=False, full_chains_thinned=False)
    npost_plt = npost
    nend = nruns
  
  #name_par, name_excluded = anc.get_sample_list(cli.sample_str, parameter_names_emcee)
  #sample_parameters, idx_sample = anc.pick_sample_parameters(flatchain_posterior_0, parameter_names_emcee, name_par = name_par, name_excluded = name_excluded)

  #flatchain_posterior_1 = flatchain_posterior_0

  if(cli.boot_id > 0):
    flatchain_posterior_msun = anc.posterior_back_to_msun(m_factor,parameter_names_emcee,flatchain_posterior_0)
    boot_file = anc.save_bootstrap_like(emcee_folder, cli.boot_id, parameter_names_emcee, flatchain_posterior_msun)
    logger.info('saved bootstrap like file: %s' %(boot_file))
    del flatchain_posterior_msun

  #median_parameters, median_perc68, median_confint = anc.get_median_parameters(flatchain_posterior_0)

  #k = np.ceil(2. * flatchain_posterior_0.shape[0]**(1./3.)).astype(int)
  #if(k>50): k=50
  k = anc.get_bins(flatchain_posterior_0, rule='doane')
  #mode_bin, mode_parameters, mode_perc68, mode_confint = anc.get_mode_parameters_full(flatchain_posterior_0, k)
  
  # max lnprob
  #max_lnprob, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint = anc.get_maxlnprob_parameters(lnprob_burnin, chains_T, flatchain_posterior_0)
  #print 
  
  try:
    overplot = int(cli.overplot)
  except:
    overplot = None
  #print overplot
  #sys.exit()
  
  ## OPEN summary_parameters.hdf5 FILE
  s_h5f = h5py.File(os.path.join(cli.full_path, 'summary_parameters.hdf5'), 'r')
  if(overplot is not None):
    # sample_parameters
    ci_fitted = s_h5f['confidence_intervals/fitted/ci'][...]
    sample_parameters = s_h5f['parameters/0666/fitted/parameters'][...]
    sample_lgllhd = s_h5f['parameters/0666'].attrs['lgllhd']
    
    try:
      sample2_parameters = s_h5f['parameters/0667/fitted/parameters'][...]
      sample2_lgllhd = s_h5f['parameters/0667'].attrs['lgllhd']
    except:
      sample2_parameters = None
      sample2_lgllhd = None
      
    try:
      sample3_parameters = s_h5f['parameters/0668/fitted/parameters'][...]
      sample3_lgllhd = s_h5f['parameters/0668'].attrs['lgllhd']
    except:
      sample3_parameters = None
      sample3_lgllhd = None
    
    median_parameters = s_h5f['parameters/1051/fitted/parameters'][...]
    median_lgllhd = s_h5f['parameters/1051'].attrs['lgllhd']
    
    max_lnprob_parameters = s_h5f['parameters/2050/fitted/parameters'][...]
    max_lgllhd = s_h5f['parameters/2050'].attrs['lgllhd']
    
    try:
      mode_parameters = s_h5f['parameters/3051/fitted/parameters'][...]
      mode_lgllhd = s_h5f['parameters/3051'].attrs['lgllhd']
    except:
      mode_parameters = None
      mode_lgllhd = None
      
    overp_par = s_h5f['parameters/%04d/fitted/parameters' %(overplot)][...]
    overp_lgllhd = s_h5f['parameters/%04d' %(overplot)].attrs['lgllhd']
  
  #nfit = s_h5f['confidence_intervals/fitted'].attrs['nfit']
  ndata = s_h5f['confidence_intervals/fitted'].attrs['ndata']
  dof = s_h5f['confidence_intervals/fitted'].attrs['dof']
  
  s_h5f.close()

  emcee_plots = os.path.join(cli.full_path,'plots')
  if (not os.path.isdir(emcee_plots)):
    os.makedirs(emcee_plots)

  for i in range(0, nfit):
    if('Ms' in parameter_names_emcee[i]):
      conv_plot = m_factor
    else:
      conv_plot = 1.
    
    emcee_fig_file = os.path.join(emcee_plots, 'chain_%03d_%s.png' %(i+1, parameter_names_emcee[i].strip()))
    print ' %s' %(emcee_fig_file),
    fig, (axChain, axHist) = plt.subplots(nrows=1, ncols=2, figsize=(12,12))

    
    (counts, bins_val, patches) = axHist.hist(flatchain_posterior_0[:,i], bins=k, range=(flatchain_posterior_0[:,i].min(), flatchain_posterior_0[:,i].max()), orientation='horizontal', normed=True, stacked=True, histtype='stepfilled', color='darkgrey', edgecolor='lightgray', align='mid')

    xpdf = scipy_norm.pdf(flatchain_posterior_0[:,i], loc = flatchain_posterior_0[:,i].mean(), scale = flatchain_posterior_0[:,i].std())
    idx = np.argsort(flatchain_posterior_0[:,i])
    axHist.plot(xpdf[idx], flatchain_posterior_0[idx,i], color='black', marker='None', ls='-.', lw=1.5, label='pdf')

    # chains after burn-in
    #axChain.plot(chains_T[:,:,i], '-', alpha=0.3)
    # chains with the burn-in
    if(cli.use_thin):
      axChain.plot(chains_T_full_thinned[:,:,i], '-', alpha=0.3)
    else:
      axChain.plot(chains_T_full[:,:,i], '-', alpha=0.3)

    axChain.axvspan(0, npost_plt, color='gray', alpha=0.45)
    axChain.axvline(npost_plt, color='gray', ls='-', lw=1.5)


    if(overplot is not None):
      if(mode_parameters is not None):
        # plot of mode (mean of higher peak/bin)
        axChain.axhline(mode_parameters[i]*conv_plot, color='red', ls='-', lw=2.1, alpha=1, label='mode')
      
      # plot of median
      axChain.axhline(median_parameters[i]*conv_plot, marker='None', c='blue',ls='-', lw=2.1, alpha=1.0, label='median fit')
      
      # plot of max_lnprob
      axChain.axhline(max_lnprob_parameters[i]*conv_plot, marker='None', c='black',ls='-', lw=1.1, alpha=1.0, label='max lnprob')
      
      if(sample_parameters is not None):
        # plot of sample_parameters
        axChain.axhline(sample_parameters[i]*conv_plot, marker='None', c='orange',ls='--', lw=2.3, alpha=0.77, label='picked: %12.7f' %(sample_parameters[i]))
      
      if(sample2_parameters is not None):
        # plot of sample2_parameters
        axChain.axhline(sample2_parameters[i]*conv_plot, marker='None', c='cyan',ls=':', lw=2.7, alpha=0.77, label='close lgllhd: %12.7f' %(sample2_parameters[i]))
      
      if(sample3_parameters is not None):
        # plot of sample3_parameters
        axChain.axhline(sample3_parameters[i]*conv_plot, marker='None', c='yellow',ls='-', lw=3.1, alpha=0.66, label='close lgllhd: %12.7f' %(sample3_parameters[i]))
      
      if(overplot not in [1050, 1051, 2050, 3050, 3051]):
        axChain.axhline(overp_par[i]*conv_plot, marker='None', c='black',ls='--', lw=2.5, alpha=0.6, label='overplot %d' %(overplot))
      
      
      # plot ci
      axChain.axhline(ci_fitted[i,0]*conv_plot, marker='None', c='forestgreen',ls='-', lw=2.1, alpha=1.0, label='CI 15.865th (%.5f)' %(ci_fitted[i,0]*conv_plot))
      axChain.axhline(ci_fitted[i,1]*conv_plot, marker='None', c='forestgreen',ls='-', lw=2.1, alpha=1.0, label='CI 84.135th (%.5f)' %(ci_fitted[i,1]*conv_plot))
    
    axChain.ticklabel_format(useOffset=False)
    xlabel = '$N_\mathrm{steps}$'
    if(cli.use_thin):
      xlabel = '$N_\mathrm{steps} \\times %d$' %(thin_steps)
    axChain.set_xlabel(xlabel)
    axChain.set_xlim([0, nend])
    axChain.set_ylabel(kel_labels[i])
    
    #y_min, y_max = anc.compute_limits(flatchain_posterior_0[:,i], 0.05)
    #if(y_min == y_max):
      #y_min = parameter_boundaries[i,0]
      #y_max = parameter_boundaries[i,1]
    y_min = flatchain_posterior_0[:,i].min()
    y_max = flatchain_posterior_0[:,i].max()
    
    axChain.set_ylim([y_min, y_max])
    axChain.set_title('Full chain %s:=[%.3f , %.3f]' %(kel_labels[i], parameter_boundaries[i,0], parameter_boundaries[i,1]))
    plt.draw()

    axHist.ticklabel_format(useOffset=False)
    axHist.set_ylim([y_min, y_max])

    if(overplot is not None):
      if(mode_parameters is not None):
        # plot mode
        axHist.axhline(mode_parameters[i]*conv_plot, color='red', ls='-', lw=2.1, alpha=1, label='mode')
        
      # plot median
      axHist.axhline(median_parameters[i]*conv_plot, marker='None', c='blue',ls='-', lw=2.1, alpha=1.0, label='median fit')
      
      # plot of max_lnprob
      axHist.axhline(max_lnprob_parameters[i]*conv_plot, marker='None', c='black',ls='-', lw=1.1, alpha=1.0, label='max lnprob')
      
      if(sample_parameters is not None):
        # plot of sample_parameters
        axHist.axhline(sample_parameters[i]*conv_plot, marker='None', c='orange',ls='--', lw=2.3, alpha=0.77, label='picked: %12.7f' %(sample_parameters[i]*conv_plot))
      
      if(sample2_parameters is not None):
        # plot of sample2_parameters
        axHist.axhline(sample2_parameters[i]*conv_plot, marker='None', c='cyan',ls=':', lw=2.7, alpha=0.77, label='close lgllhd: %12.7f' %(sample2_parameters[i]))
      
      if(sample3_parameters is not None):
        # plot of sample3_parameters
        axHist.axhline(sample3_parameters[i]*conv_plot, marker='None', c='yellow',ls='-', lw=3.1, alpha=0.66, label='close lgllhd: %12.7f' %(sample3_parameters[i]))
      
      if(overplot not in [1050, 1051, 2050, 3050, 3051]):
        axHist.axhline(overp_par[i]*conv_plot, marker='None', c='black',ls='--', lw=2.5, alpha=0.8, label='overplot %d' %(overplot))
      
      # plot ci
      axHist.axhline(ci_fitted[i,0]*conv_plot, marker='None', c='forestgreen',ls='-', lw=2.1, alpha=1.0, label='CI 15.865th (%.5f)' %(ci_fitted[i,0]*conv_plot))
      axHist.axhline(ci_fitted[i,1]*conv_plot, marker='None', c='forestgreen',ls='-', lw=2.1, alpha=1.0, label='CI 84.135th (%.5f)' %(ci_fitted[i,1]*conv_plot))
    
    axHist.set_title('Distribution of posterior chain')
    axHist.legend(loc='center left', fontsize=9, bbox_to_anchor=(1, 0.5))
    plt.draw()

    fig.savefig(emcee_fig_file, bbox_inches='tight', dpi=150)
    print ' saved'
    print

  fig = plt.figure(figsize=(12,12))
  
  # lnprob
  xlabel = '$N_\mathrm{steps}$'
  if(cli.use_thin):
    xlabel = '$N_\mathrm{steps} \\times %d$' %(thin_steps)
  
  ax = plt.subplot2grid((2,1), (0,0))
  ax.plot(lnprob_burnin.T, '-', alpha=0.3)

  if(overplot is not None):
    posterior_msun = anc.posterior_back_to_msun(m_factor,parameter_names_emcee,flatchain_posterior_0)
    post_sel, lnprob_sel = anc.select_within_all_ci(posterior_msun, ci_fitted[:,0:2], lnprob_burnin.T.reshape(-1))
    #lnprob_sel = lnprob_burnin.T.reshape((-1))
    lgllhd_med = np.percentile(lnprob_burnin.T.reshape(-1), 50., interpolation='midpoint')
    abs_dlg = np.abs(lnprob_sel - lgllhd_med)
    lgllhd_mad = np.percentile(abs_dlg, 50., interpolation='midpoint')
    
    #lnp_min = np.min(lnprob_sel)
    #lnp_max = np.max(lnprob_sel)
    lnp_min = lgllhd_med - lgllhd_mad
    lnp_max = lgllhd_med + lgllhd_mad
    print ' lgllhd_med & mad = ',lgllhd_med, lgllhd_mad
    print ' lnp_min = ',lnp_min, ' lnp_max = ',lnp_max
    print ' lnl_668 = ',sample3_lgllhd
  
    ax.axhline(lgllhd_med, color='black', ls='-', lw=1.6, alpha=0.77)
  
  #if(sample2_lgllhd is not None):
    #ax.axhline(sample2_lgllhd, marker='None', c='cyan',ls=':', lw=2.7, alpha=0.9)
  
    if(sample3_lgllhd is not None):
      ax.axhline(sample3_lgllhd, marker='None', c='yellow',ls='-', lw=3.1, alpha=0.9)
    
    ax.axhspan(lnp_min, lnp_max, color='gray', alpha=0.77)
    ax.axhline(lnp_min, color='black', ls='--', lw=1.6, alpha=0.77)
    ax.axhline(lnp_max, color='black', ls='--', lw=1.6, alpha=0.77)
  
  min_lnp = np.min(lnprob_burnin.T, axis=0).min()
  max_lnp = np.max(lnprob_burnin.T, axis=0).max()
  y_min, y_max = anc.compute_limits(np.asarray([min_lnp, max_lnp]), 0.05)
  ax.set_ylim((y_min, y_max))
  ax.set_ylabel('lnprob')
  #ax.get_xaxis().set_visible(False)
  ax.set_xlabel(xlabel)
  
  # chi2r
  chi2r = -2.*(lnprob_burnin.T-ln_err_const)/np.float64(dof)
  
  ax = plt.subplot2grid((2,1), (1,0))
  ax.axhline(1.0, color='gray', ls='-')
  ax.plot(chi2r, '-', alpha=0.3)
  
  if(overplot is not None):
    c2r_med = -(2.*(lgllhd_med - ln_err_const))/np.float64(dof)
    c2r_smax = -(2.*(lnp_min - ln_err_const))/np.float64(dof)
    c2r_smin = -(2.*(lnp_max - ln_err_const))/np.float64(dof)
  
    print ' c2r_med = ',c2r_med
    print ' c2r_smin = ',c2r_smin, ' c2r_smax = ', c2r_smax
    
    ax.axhline(c2r_med, color='black', ls='-', lw=1.6, alpha=0.77)
    ax.axhspan(c2r_smin, c2r_smax, color='gray', alpha=0.77)
    ax.axhline(c2r_smin, color='black', ls='--', lw=1.6, alpha=0.77)
    ax.axhline(c2r_smax, color='black', ls='--', lw=1.6, alpha=0.77)
    #if(sample2_lgllhd is not None):
      #c2r_sample2 = -2.*(sample2_lgllhd - ln_err_const)/np.float64(dof)
      #ax.axhline(c2r_sample2, marker='None', c='cyan',ls=':', lw=2.7, alpha=0.9)
    if(sample3_lgllhd is not None):
      c2r_sample3 = -2.*(sample3_lgllhd - ln_err_const)/np.float64(dof)
      ax.axhline(c2r_sample3, marker='None', c='yellow',ls='-', lw=3.1, alpha=0.9)
  
  c2r_min = -2.*(y_max - ln_err_const)/np.float64(dof)
  c2r_max = -2.*(y_min - ln_err_const)/np.float64(dof)
  ax.set_ylim((c2r_min, c2r_max))
  ax.set_ylabel('$\chi^{2}/\mathrm{dof}$')
  #ax.get_xaxis().set_visible(True)
  ax.set_xlabel(xlabel)
  
  fig.savefig(os.path.join(emcee_plots, 'emcee_lnprobability.png'), bbox_inches='tight', dpi=150)
  print ' %s saved' %(os.path.join(emcee_plots, 'emcee_lnprobability.png'))


  return

if __name__ == "__main__":
  main()




