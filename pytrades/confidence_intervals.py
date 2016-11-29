#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import sys
import argparse
import time
import os
import numpy as np # array
import h5py
import constants as cst # local constants module
import ancillary as anc
import pytrades_lib

#from matplotlib import use as mpluse
#mpluse("Agg")
##mpluse("Qt4Agg")
#import matplotlib.pyplot as plt
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('text', usetex=True)
#from matplotlib import rcParams
#rcParams['text.latex.unicode']=True

# read command line (cli) arguments

# main

def main():
  print
  print ' TRADES: EMCEE confidence intervals'
  print

  cli = anc.get_args()

  # init trades
  pytrades_lib.pytrades.initialize_trades(os.path.join(cli.full_path, ''), '')
  MR_star = pytrades_lib.pytrades.mr_star

  nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case_list = anc.get_fitted(cli.full_path)

  # read emcee data
  emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(cli.full_path, cli.temp_status)
  # get data from the hdf5 file
  names_par, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)
  # print Memory occupation of ...
  anc.print_memory_usage(chains)

  nfit, nwalkers, nruns, npost, nruns_sel = anc.get_emcee_parameters(chains, cli.temp_status, cli.npost, completed_steps)
  #chains_T, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, names_par, parameter_boundaries, chains)
  chains_T = np.zeros((nruns_sel, nwalkers, nfit))
  for ii in xrange(0,nfit):
    chains_T[:,:,ii] = chains[:,npost:nruns,ii].T # transpose after removing the burnin steps
  flatchain_posterior_0 = chains_T[:,:,:].reshape((-1, nfit))

  # computes mass conversion factor
  #m_factor = anc.mass_conversion_factor(cli.m_type)
  m_factor, mass_unit = anc.mass_type_factor(MR_star[0,0], cli.m_type, True)
  # set label and legend names
  #kel_legends, labels_list = anc.keplerian_legend(names_par, cli.m_type)

  top_header, header = anc.get_header(anc.percentile_val)

  

  def get_intervals(full_path, id_sim, parameter_names, parameters, flatchain_posterior, derived_type=None, full_output=False):
    
    out_folder = os.path.join(os.path.join(full_path, '%04d_sim' %(id_sim)), '')
    if (not os.path.isdir(out_folder)):
      os.makedirs(out_folder)
    out_file = os.path.join(out_folder, 'parameters_summary.txt')
    out = open(out_file, 'w')
    pytrades_lib.pytrades.path_change(out_folder)

    anc.print_both(' ORIGINAL PARAMETER VALUES -> %d' %(id_sim), out)
    fitness, lgllhd, check = pytrades_lib.pytrades.write_summary_files(id_sim, parameters)

    kel_file, kep_elem = anc.elements(out_folder, id_sim, lmf=0)

    sigma_par = anc.compute_intervals(flatchain_posterior, parameters, anc.percentile_val)
    units_par = anc.get_units(names_par, mass_unit)
    
    names_derived, der_posterior = anc.compute_derived_posterior(parameter_names, kep_elem, id_fit, case_list, cols_list, flatchain_posterior, conv_factor=m_factor)
    #der_posterior_T = der_posterior.T
    der_posterior_T = der_posterior
    

    if(str(derived_type).strip().lower() == 'median'):
      # MEDIAN PARAMETERS ID == 1050
      derived_par = np.percentile(der_posterior_T, 50., axis=0)
      par_type = 'MEDIAN'
    elif(str(derived_type).strip().lower() == 'mode'):
      # MODE-LIKE PARAMETERS -> id 3050
      k = np.ceil(2. * flatchain_posterior_0.shape[0]**(1./3.)).astype(int)
      if(k>50): k=50
      der_bin, derived_par = anc.get_mode_parameters(der_posterior_T, k)
      par_type = 'MODE'
    else:
      # ORIGINAL FITTING PARAMETERS ID == 0
      # or
      # MAX LNPROBABILITY -> id 2050
      names_derived, derived_par = anc.compute_derived_parameters(parameter_names, kep_elem, id_fit, case_list, cols_list, parameters, conv_factor=m_factor)
      derived_par, der_posterior_T = anc.adjust_derived_parameters(names_derived, derived_par, der_posterior_T)
      if(id_sim == 0):
        par_type = 'ORIGINAL FIT'
      else:
        par_type = 'MAX LNPROB'
    
    sigma_derived = anc.compute_intervals(der_posterior_T, derived_par, anc.percentile_val)
    units_der = anc.get_units(names_derived, mass_unit)
    
    anc.print_both('\n# SUMMARY: %s' %(par_type), out)
    anc.print_both('# FITTED PARAMETERS', out)
    anc.print_parameters(top_header, header, names_par, units_par, parameters, sigma_par, out)
    
    anc.print_both('# DERIVED PARAMETERS', out)
    anc.print_parameters(top_header, header, names_derived, units_der, derived_par, sigma_derived, out)
    out.close()
    if(full_output):
      return out_folder, names_derived, der_posterior_T
    else:
      return out_folder


  # ************************************
  ## ORIGINAL FITTING PARAMETERS ID == 0
  # save initial_fitting parameters into array
  original_fit_parameters = pytrades_lib.pytrades.fitting_parameters # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
  folder_0 = get_intervals(cli.full_path, 0, names_par, original_fit_parameters, flatchain_posterior_0, derived_type=None)
  # ************************************
  
  print
  print
  
  # ************************************
  ## MAX LNPROBABILITY AND PARAMETERS -> id 2050
  max_lnprob, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint = anc.get_maxlnprob_parameters(npost, nruns, lnprobability, chains_T, flatchain_posterior_0)
  folder_2050, names_derived, der_posterior = get_intervals(cli.full_path, 2050, names_par, max_lnprob_parameters, flatchain_posterior_0, derived_type=None, full_output=True)
  # write out the derived names and posterior into an hdf5 file
  der_post_file = os.path.join(cli.full_path, 'derived_posterior.hdf5')
  h5f = h5py.File(der_post_file, 'w')
  h5f.create_dataset('derived_names', data=names_derived, dtype='S10')
  h5f.create_dataset('derived_posterior', data=der_posterior, dtype=np.float64)
  h5f.close()
  # ************************************
  
  print
  print
  
  # ************************************
  ## MEDIAN PARAMETERS ID == 1050
  median_parameters, median_perc68, median_confint = anc.get_median_parameters(flatchain_posterior_0)
  folder_1050 = get_intervals(cli.full_path, 1050, names_par, median_parameters, flatchain_posterior_0, derived_type='median')
  # ************************************
  
  print
  print

  # ************************************
  ## MODE-LIKE PARAMETERS -> id 3050
  ## take the mean of 5 bin centered to the higher bin
  k = np.ceil(2. * flatchain_posterior_0.shape[0]**(1./3.)).astype(int)
  if(k>50): k=50
  mode_bin, mode_parameters = anc.get_mode_parameters(flatchain_posterior_0, k)
  folder_3050 = get_intervals(cli.full_path, 3050, names_par, mode_parameters, flatchain_posterior_0, derived_type='mode')
  # ************************************
  
  print
  print
  

  pytrades_lib.pytrades.deallocate_variables()

  return

if __name__ == "__main__":
  main()
