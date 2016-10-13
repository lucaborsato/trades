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
from ancillary import *
from matplotlib import use as mpluse
mpluse("Agg")
#mpluse("Qt4Agg")
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)
#from matplotlib import rcParams
#rcParams['text.latex.unicode']=True
#import corner


## common variables needed to create labels and parameter names
#kel_id = ['$M$', '$R$', '$P$', '$e\cos\omega$', '$e\sin\omega$', '$\\nu$', '$i$', '$\Omega$']
#kel_id_2 = ['M', 'R', 'P', 'e\cos\omega', 'e\sin\omega', '\\nu', 'i', '\Omega']
#kel_units = ['$[M_{\oplus}]$', '$[R_\mathrm{Jup}]$', '[days]', '', '', '[deg]', '[deg]', '[deg]']
#kel_fmt = ['%.3f', '%.4f', '%.3f', '%.1f', '%.1f', '%.1f', '%.1f', '%.1f']
#nelem = len(kel_id)
#letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'l', 'm', 'n', 'o', 'p']

#def set_bool_argument(arg):
  #if (arg != False):
    #if (arg.lower() in ['t', 'tr', 'tru', 'true']):
      #arg = True
    #else:
      #arg = False
  #return arg

# read command line (cli) arguments
#def get_args():
  #parser = argparse.ArgumentParser(description='TRADES+EMCEE PLOT')
  #parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')
  #parser.add_argument('-nb', '--nburn', '-np', '--npost', action='store', dest='npost', required=True, help='The number of posterior/burn in steps to discard at the beginning of each chain. It has to be > 0')
  #parser.add_argument('-m', '--mtype', '--mass-type', action='store', dest='m_type', default='e', help='Mass type: j = Jupiter, e = Earth, s = Sun. Default is Earth = e.')
  #parser.add_argument('-g', '--good-parameters', action='store', dest='good', default=False, help='If you want to use a previous solution, set it to True and it will search for a good_parameters.dat file with first column the name of the parameter and the value in the second. Mass paramenters in Jupiter mass. Default is False')
  #parser.add_argument('-t', '--temp-file', action='store', dest='temp_status', default=False, help='If you want to read temporary emcee_temp.hdf5 file, because simulation is not finished yet. Default is False')
  
  #cli = parser.parse_args()
  #cli.full_path = os.path.abspath(cli.full_path)
  #cli.m_type = cli.m_type.lower()
  
  ##if (cli.good.lower() in ['t', 'tr', 'tru', 'true']):
    ##cli.good = True
  ##else:
    ##cli.good = False
  ##cli.good = set_bool_argument(cli.good)
  
  ##if (cli.temp_status != False):
    ##if (cli.temp_status.lower() in ['t', 'tr', 'tru', 'true']):
      ##cli.temp_status = True
    ##else:
      ##cli.temp_status = False
  #cli.temp_status = set_bool_argument(cli.temp_status)
  
  #return cli

print 
print ' ================== '
print ' TRADES+EMCEE PLOTS'
print ' ================== '
print

# read cli arguments
cli = get_args()
# computes mass conversion factor
m_factor = mass_conversion_factor(cli.m_type)

# set emcee and trades folder
emcee_folder = cli.full_path
trades_folder = os.path.join(os.path.dirname(cli.full_path), '')
# and best folder
emcee_file, emcee_best, folder_best = get_emcee_file_and_best(emcee_folder, cli.temp_status)

parameter_names_emcee, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = get_data(emcee_file, cli.temp_status)

# set label and legend names
kel_legends, labels_list = keplerian_legend(parameter_names_emcee, cli.m_type)

nfit, nwalkers, nruns, npost, nruns_sel = get_emcee_parameters(chains, cli.temp_status, cli.npost, completed_steps)

print_memory_usage(chains)

chains_T, parameter_boundaries = select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries, chains)

# create a flat array of the posterior: from (nruns_sel, nwalkers, nfit) -> (nruns_sel * nwalkers, nfit)
flatchain_posterior_0 = chains_T[:,:,:].reshape((-1, nfit))
#flatchain_posterior_1 = flatchain_posterior_0.copy()
flatchain_posterior_1 = flatchain_posterior_0

derived_names, derived_chains_T, derived_posterior = get_derived_posterior_parameters(parameter_names_emcee, chains_T, flatchain_posterior_0)
nder = len(derived_names)

if(cli.boot_id > 0):
  flatchain_posterior_msun = posterior_back_to_msun(m_factor,parameter_names_emcee,flatchain_posterior_0)
  boot_file = save_bootstrap_like(emcee_folder, cli.boot_id, parameter_names_emcee, flatchain_posterior_msun)
  logger.info('saved bootstrap like file: %s' %(boot_file))
  del flatchain_posterior_msun

nburnin=npost

median_parameters, median_perc68, median_confint = get_median_parameters(flatchain_posterior_0)
median_parameters_der, median_perc68_der, median_confint_der = get_median_parameters(derived_posterior)
print '# MEDIAN PARAMETER VALUES -> 1050'
print_parameters(parameter_names_emcee, median_parameters, median_perc68, median_confint, 'median')
print_parameters(derived_names, median_parameters_der, median_perc68_der, median_confint_der, 'median_der')
print '# ' + '-'*220

k = np.ceil(2. * flatchain_posterior_0.shape[0]**(1./3.)).astype(int)
#if(k>11): k=11
if(k>50): k=50
mode_bin, mode_parameters, mode_perc68, mode_confint = get_mode_parameters(flatchain_posterior_0, k)
mode_bin_der, mode_parameters_der, mode_perc68_der, mode_confint_der = get_mode_parameters(derived_posterior, k)
print '# MODE PARAMETER VALUES -> 3050'
print_parameters(parameter_names_emcee, mode_parameters, mode_perc68, mode_confint, 'mode')
print_parameters(derived_names, mode_parameters_der, mode_perc68_der, mode_confint_der, 'mode_der')
print '# ' + '-'*220

# read good parameters file
if (cli.good):
  good_file = os.path.join(cli.full_path, 'good_parameters.dat')
  if (os.path.exists(good_file)):
    good_parameters = np.genfromtxt(good_file, usecols=(1), dtype=np.float64)
    good_id = np.genfromtxt(good_file, usecols=(0), dtype='|S10')
    good_parameters = check_good_parameters(good_id, good_parameters, m_factor)
    good_parameters_0 = good_parameters.copy()

print 

emcee_plots = os.path.join(cli.full_path,'plots')
if (not os.path.isdir(emcee_plots)):
  os.makedirs(emcee_plots)

#median_mode_parameters = []

#print ' Number of bins: %d' %(k)
for i in range(0, nfit):
  emcee_fig_file = os.path.join(emcee_plots, 'chain_' + parameter_names_emcee[i] + '.png')
  print ' %s' %(emcee_fig_file),
  fig, (axChain, axHist) = plt.subplots(nrows=1, ncols=2, figsize=(12,12))

  #axHist.hist(flatchain_posterior_1[:,i], bins=k, range=(parameter_boundaries[i,0], parameter_boundaries[i,1]), orientation='horizontal', normed=False, histtype='bar', rwidth=0.9, color='gray', edgecolor='white')
  #if (k<100):
    #(counts, bins_val, patches) = axHist.hist(flatchain_posterior_1[:,i], bins=k, range=(parameter_boundaries[i,0], parameter_boundaries[i,1]), orientation='horizontal', normed=False, histtype='stepfilled', color='darkgrey', align='mid')
  #else:
    #(counts, bins_val, patches) = axHist.hist(flatchain_posterior_1[:,i], bins=k, range=(parameter_boundaries[i,0], parameter_boundaries[i,1]), orientation='horizontal', normed=False, histtype='step', color='darkgray')
  
  (counts, bins_val, patches) = axHist.hist(flatchain_posterior_1[:,i], bins=k, range=(flatchain_posterior_1[:,i].min(), flatchain_posterior_1[:,i].max()), orientation='horizontal', normed=True, stacked=True, histtype='stepfilled', color='darkgrey', edgecolor='lightgray', align='mid')

  xpdf = scipy_norm.pdf(flatchain_posterior_1[:,i], loc = flatchain_posterior_1[:,i].mean(), scale = flatchain_posterior_1[:,i].std())
  idx = np.argsort(flatchain_posterior_1[:,i])
  axHist.plot(xpdf[idx], flatchain_posterior_1[idx,i], color='black', marker='None', ls='-.', lw=1.5, label='pdf')

  axChain.plot(chains_T[:,:,i], '-', alpha=0.3)

  if (cli.good):
    axChain.axhline(good_parameters_0[i], marker='None', c='k',ls='-', lw=1., label='best guess')

  median  = median_parameters[i]
  lower = median_parameters[i]+median_confint[0,i]
  upper = median_parameters[i]+median_confint[3,i]

  # compute mean of higher peak/bin
  axChain.axhline(mode_parameters[i], color='forestgreen', ls='-.', lw=2.1, alpha=0.7, label='mode')
  axChain.axhline(mode_parameters[i]+mode_confint[0,i], color='forestgreen', ls='-.', lw=1.1, alpha=0.7, label='lower mode')
  axChain.axhline(mode_parameters[i]+mode_confint[3,i], color='forestgreen', ls='-.', lw=1.1, alpha=0.7, label='upper mode')
  
  #axChain.axhline(max_lnprob_parameters[i], marker='None', c='cyan',ls='-', lw=1.5, label='max lnprob fit', alpha=0.65)
  axChain.axhline(median, marker='None', c='b',ls='--', lw=2., label='median fit')
  axChain.axhline(lower, marker='None', c='coral',ls='--', lw=2., label='lower median')
  axChain.axhline(upper, marker='None', c='r',ls='--', lw=2., label='upper median')
  #axChain.axhline(median_mode, marker='None', c='y',ls='--', lw=1., label='median mode fit')
  
  axChain.ticklabel_format(useOffset=False)
  axChain.set_xlabel('$N_\mathrm{steps}$')
  axChain.set_ylabel(kel_legends[i])
  
  y_min, y_max = compute_limits(flatchain_posterior_1[:,i], 0.05)
  if(y_min == y_max):
    y_min = parameter_boundaries[i,0]
    y_max = parameter_boundaries[i,1]
  
  #axChain.set_ylim([parameter_boundaries[i,0], parameter_boundaries[i,1]])
  #axChain.set_ylim([flatchain_posterior_1[:,i].min(), flatchain_posterior_1[:,i].max()])
  axChain.set_ylim([y_min, y_max])
  axChain.set_title('Full chain %s:=[%.3f , %.3f]' %(kel_legends[i], parameter_boundaries[i,0], parameter_boundaries[i,1]))
  plt.draw()

  axHist.ticklabel_format(useOffset=False)
  #axHist.set_ylim([parameter_boundaries[i,0], parameter_boundaries[i,1]])
  #axHist.set_ylim([flatchain_posterior_1[:,i].min(), flatchain_posterior_1[:,i].max()])
  axHist.set_ylim([y_min, y_max])

  if (cli.good):
    axHist.axhline(good_parameters_0[i], marker='None', c='k',ls='-', lw=1., label='best guess')

  axHist.axhline(mode_parameters[i], color='forestgreen', ls='-.', lw=2.1, alpha=0.7, label='mode')
  axHist.axhline(mode_parameters[i]+mode_confint[0,i], color='forestgreen', ls='-.', lw=1.1, alpha=0.7, label='lower mode')
  axHist.axhline(mode_parameters[i]+mode_confint[3,i], color='forestgreen', ls='-.', lw=1.1, alpha=0.7, label='upper mode')
  axHist.axhline(median, marker='None', c='b',ls='--', lw=2., label='median fit')
  axHist.axhline(lower, marker='None', c='coral',ls='--', lw=2., label='lower median')
  axHist.axhline(upper, marker='None', c='r',ls='--', lw=2., label='upper median')
  
  axHist.set_title('Distribution of posterior chain')
  axHist.legend()
  plt.draw()

  fig.savefig(emcee_fig_file, bbox_inches='tight', dpi=150)
  print 'saved'
  if (cli.good):
    print ' guess = ',good_parameters_0[i],
  print 'median = ', median,             ' lower = ', lower,                                ' upper = ', upper
  print '  mode = ', mode_parameters[i], ' lower = ', mode_parameters[i]+mode_confint[0,i], ' upper = ', mode_parameters[i]+mode_confint[3,i]
  print

lnprob_burnin = lnprobability[:,nburnin:]
#if ('lnprobability' in data_names):
fig = plt.figure(figsize=(12,12))
plt.plot(lnprob_burnin.T, '-', alpha=0.8)
#plt.plot(lnprobability[:,:].T, '-', alpha=0.8)
plt.xlabel('$N_\mathrm{steps}$')
plt.ylabel('logprob')
min_lnp = np.min(lnprobability[:,:].T, axis=0).min()
max_lnp = np.max(lnprobability[:,:].T, axis=0).max()
#plt.ylim((min_lnp, max_lnp))
#plt.ylim((-100., max_lnp))
y_min, y_max = compute_limits(np.asarray([min_lnp, max_lnp]), 0.05)
plt.ylim((y_min, y_max))
plt.draw()
fig.savefig(os.path.join(emcee_plots, 'emcee_lnprobability.png'), bbox_inches='tight', dpi=150)
print ' %s saved' %(os.path.join(emcee_plots, 'emcee_lnprobability.png'))






