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
kel_legends, labels_list = anc.keplerian_legend(parameter_names_emcee, cli.m_type)

nfit, nwalkers, nruns, npost, nruns_sel = anc.get_emcee_parameters(chains, cli.temp_status, cli.npost, completed_steps)

anc.print_memory_usage(chains)

chains_T, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries, chains)

# create a flat array of the posterior: from (nruns_sel, nwalkers, nfit) -> (nruns_sel * nwalkers, nfit)
flatchain_posterior_0 = chains_T[:,:,:].reshape((-1, nfit))
flatchain_posterior_1 = flatchain_posterior_0

sel_min = np.arange(10000,completed_steps,5000)
sel_max = sel_min + 5000
if(sel_max[-1] > completed_steps): sel_max[-1] = completed_steps
n_sel = sel_min.shape[0]
print 'completed_steps = %d and n_sel = %d' %(completed_steps, n_sel)
m_full_posterior = flatchain_posterior_0[:,0]

fig = plt.figure(figsize=(12,12))
#fig, axHist = plt.subplots(nrows=sel_min.shape[0], ncols=1, figsize=(12,12))

for isel in range(0, sel_min.shape[0]):
  print 'isel = %d with selection: %d : %d' %(isel, sel_min[isel], sel_max[isel])
  m_sel = m_full_posterior[sel_min[isel]:sel_max[isel]]
  k = np.ceil(2. * m_sel.shape[0]**(1./3.)).astype(int)
  if(k>50): k=50
  print 'selected slide of the chain'
  axHist = plt.subplot2grid((n_sel, 1), (isel,0))
  (counts, bins_val, patches) = axHist.hist(m_sel, bins=k, range=(m_sel.min(), m_sel.max()), orientation='vertical', normed=True, stacked=True, histtype='stepfilled', color='darkgrey', edgecolor='lightgray', align='mid', label='%s [steps %d : %d]' %(kel_legends[0],sel_min[isel], sel_max[isel]))
  print 'plotted histogram'
  xpdf = scipy_norm.pdf(m_sel, loc = m_sel.mean(), scale = m_sel.std())
  idx = np.argsort(m_sel)
  axHist.plot(m_sel[idx], xpdf[idx], color='black', marker='None', ls='-.', lw=1.5, label='pdf')
  print 'computed and plotted probability distribution function pdf'
  #axHist.set_xlabel('%s [steps %d : %d]' %(parameter_names_emcee[0].strip(),sel_min[isel], sel_max[isel]))
  axHist.legend(loc='center left', fontsize=9, bbox_to_anchor=(1, 0.5))
  print 'set xlabel'
plt.draw()

print 'saving into file ...'
emcee_plots = os.path.join(cli.full_path,'plots')
if (not os.path.isdir(emcee_plots)):
  os.makedirs(emcee_plots)
emcee_fig_file = os.path.join(emcee_plots, 'sel_hist_' + parameter_names_emcee[0].strip() + '.png')
fig.savefig(emcee_fig_file, bbox_inches='tight', dpi=150)
print 'saved into %s' %(emcee_fig_file)
  





  
