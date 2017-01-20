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

# main

def main():
  print
  print ' TRADES: EMCEE confidence intervals'
  print

  cli = anc.get_args()

  # init trades
  pytrades_lib.pytrades.initialize_trades(os.path.join(cli.full_path, ''), '')

  nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case_list = anc.get_fitted(cli.full_path)
  ndata = pytrades_lib.pytrades.ndata
  dof = pytrades_lib.pytrades.dof
  
  # read emcee data
  emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(cli.full_path, cli.temp_status)
  # get data from the hdf5 file
  names_par, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)
  # print Memory occupation of ...
  anc.print_memory_usage(chains)

  nfit, nwalkers, nruns, npost, nruns_sel = anc.get_emcee_parameters(chains, cli.temp_status, cli.npost, completed_steps)
  #chains_T, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, names_par, parameter_boundaries, chains)
  chains_T_full = np.zeros((nruns, nwalkers, nfit))
  for ii in xrange(0,nfit):
    chains_T_full[:,:,ii] = chains[:,:nruns,ii].T # transpose after removing the burnin steps
    
  chains_T, flatchain_posterior_0, lnprob_burnin, thin_steps = anc.thin_the_chains(cli.use_thin, npost, nruns, nruns_sel, autocor_time, chains_T_full, lnprobability, burnin_done=False)
  
  # computes mass conversion factor
  #m_factor = anc.mass_conversion_factor(cli.m_type)
  MR_star = pytrades_lib.pytrades.mr_star
  m_factor, mass_unit = anc.mass_type_factor(Ms=1.0, mtype=cli.m_type, mscale=False)
  np.random.seed(seed=cli.seed)
  Ms_gaussian = MR_star[0,0] + np.random.normal(0., 1., size=(np.shape(flatchain_posterior_0)[0]))*MR_star[0,1] # if exists an error on the mass, it creates a Normal distribution for the values and use it to re-scale mp/Ms to mp.
  m_factor_boot = m_factor * Ms_gaussian # given the factor from Msun to mass_unit it multiply it by the Normal Mstar.
  m_factor = m_factor * MR_star[0,0]

  # set label and legend names
  #kel_legends, labels_list = anc.keplerian_legend(names_par, cli.m_type)

  top_header, header = anc.get_header(anc.percentile_val)

  # GET INTERVALS
  def get_intervals(full_path, id_sim, parameter_names, parameters, flatchain_posterior, derived_type=None, full_output=False, idx_sample=None, summary_file_hdf5=None):
    
    out_folder = os.path.join(os.path.join(full_path, '%04d_sim' %(id_sim)), '')
    if (not os.path.isdir(out_folder)):
      os.makedirs(out_folder)
    out_file = os.path.join(out_folder, 'parameters_summary.txt')
    out = open(out_file, 'w')
    pytrades_lib.pytrades.path_change(out_folder)
    
    anc.print_both(' ', out)
    anc.print_both(' --------------------------------- ', out)
    anc.print_both(' PARAMETER VALUES -> %d' %(id_sim), out)
    fitness, lgllhd, check = pytrades_lib.pytrades.write_summary_files(id_sim, parameters)
    print 'CHECK = ',check
      
    kel_file, kep_elem = anc.elements(out_folder, id_sim, lmf=0)

    sigma_par = anc.compute_intervals(flatchain_posterior, parameters, anc.percentile_val)
    units_par = anc.get_units(names_par, mass_unit)
    
    if(not bool(check)):
      print 'WRTING WARNING FILE: %s' %(os.path.join(out_folder,'WARNING.txt'))
      warn_o = open(os.path.join(out_folder,'WARNING.txt'), 'w')
      warn_o.write('*******\nWARNING: FITTED PARAMETERS COULD NOT BE PHYSICAL!\nWARNING: BE VERY CAREFUL WITH THIS PARAMETER SET!\n*******')
      warn_o.close()
      if(full_output):
        return out_folder, None, None
      else:
        return out_folder
    
    
    names_derived, der_posterior = anc.compute_derived_posterior(parameter_names, kep_elem, id_fit, case_list, cols_list, flatchain_posterior, conv_factor=m_factor_boot)
    
    #der_posterior_T = der_posterior
    der_posterior_T = anc.derived_posterior_check(names_derived, der_posterior)

    if(str(derived_type).strip().lower() == 'median'):
      # MEDIAN PARAMETERS ID == 1050
      derived_par = np.percentile(der_posterior_T, 50., axis=0, interpolation='midpoint')
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
      elif(id_sim == 1051):
        par_type = 'MEDIAN PARAMETERS TO DERIVED'
      elif(id_sim == 2050):
        par_type = 'MAX LNPROB'
      elif(id_sim == 3051):
        par_type = 'MODE PARAMETERS TO DERIVED'
      elif(id_sim == 666):
        par_type = 'SELECTED SAMPLE'
        if(idx_sample is not None):
          par_type = '%s <-> idx = %d' %(par_type, idx_sample)
          derived_par = der_posterior_T[idx_sample, :]
      else:
        par_type = 'AD HOC'
    
    sigma_derived = anc.compute_intervals(der_posterior_T, derived_par, anc.percentile_val)
    units_der = anc.get_units(names_derived, mass_unit)
    
    if(s_h5f is not None):
      s_id_sim = '%04d' %(id_sim)
      #print s_id_sim
      s_h5f.create_dataset('parameters/%s/fitted/parameters' %(s_id_sim), data=parameters, dtype=np.float64)
      s_h5f.create_dataset('parameters/%s/fitted/names' %(s_id_sim), data=names_par, dtype='S10')
      s_h5f.create_dataset('parameters/%s/fitted/units' %(s_id_sim), data=units_par, dtype='S15')
      s_h5f.create_dataset('parameters/%s/fitted/sigma' %(s_id_sim), data=sigma_par.T, dtype=np.float64)
      s_h5f['parameters/%s/fitted/sigma' %(s_id_sim)].attrs['percentiles'] = anc.percentile_val
    #if(s_h5f is not None):
      #s_id_sim = '%04d' %(id_sim)
      s_h5f.create_dataset('parameters/%s/derived/parameters' %(s_id_sim), data=derived_par, dtype=np.float64)
      s_h5f.create_dataset('parameters/%s/derived/names' %(s_id_sim), data=names_derived, dtype='S10')
      s_h5f.create_dataset('parameters/%s/derived/units' %(s_id_sim), data=units_der, dtype='S15')
      s_h5f.create_dataset('parameters/%s/derived/sigma' %(s_id_sim), data=sigma_derived.T, dtype=np.float64)
      s_h5f['parameters/%s/derived/sigma' %(s_id_sim)].attrs['percentiles'] = anc.percentile_val
      s_h5f['parameters/%s' %(s_id_sim)].attrs['info'] = '%s ==> %s' %(s_id_sim, par_type)
      s_h5f['parameters/%s' %(s_id_sim)].attrs['fitness'] = fitness
      s_h5f['parameters/%s' %(s_id_sim)].attrs['lgllhd'] = lgllhd
      s_h5f['parameters/%s' %(s_id_sim)].attrs['check'] = check
      if(idx_sample is not None):
        s_h5f['parameters/%s' %(s_id_sim)].attrs['idx_sample'] = idx_sample
    
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
  
  ## CREATE A HDF5 FILE WITH CONFIDNCE INTERVALS AND ALL THE SUMMARY PARAMETERS
  summary_file = os.path.join(cli.full_path, 'summary_parameters.hdf5')
  s_h5f = h5py.File(summary_file, 'w')
  
  ## COMPUTE CONFIDENCE INTERVALS OF THE FITTED PARAMETER DISTRIBUTIONS
  ci_fitted = np.percentile(flatchain_posterior_0, anc.percentile_val[2:], axis=0, interpolation='midpoint') # (n_percentile-1 x nfit) ==> skip first item, the 68.27th
  units_par = anc.get_units(names_par, mass_unit)
  
  s_h5f.create_dataset('confidence_intervals/fitted/ci', data=ci_fitted.T, dtype=np.float64)
  s_h5f.create_dataset('confidence_intervals/fitted/names', data=names_par, dtype='S10')
  s_h5f.create_dataset('confidence_intervals/fitted/units', data=units_par, dtype='S15')
  s_h5f.create_dataset('confidence_intervals/fitted/percentiles', data=np.array(anc.percentile_val[2:]), dtype=np.float64)
  
  s_h5f['confidence_intervals/fitted'].attrs['nfit'] = nfit
  s_h5f['confidence_intervals/fitted'].attrs['ndata'] = ndata
  s_h5f['confidence_intervals/fitted'].attrs['dof'] = dof
  
  
  ## ORIGINAL FITTING PARAMETERS ID == 0
  # save initial_fitting parameters into array
  original_fit_parameters = pytrades_lib.pytrades.fitting_parameters # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
  folder_0 = get_intervals(cli.full_path, 0, names_par, original_fit_parameters, flatchain_posterior_0, derived_type=None, summary_file_hdf5=s_h5f)
  # ************************************
  
  print
  print
  
  # ************************************
  ## MAX LNPROBABILITY AND PARAMETERS -> id 2050
  max_lnprob, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint = anc.get_maxlnprob_parameters(lnprob_burnin, chains_T, flatchain_posterior_0)
  max_id1, max_id2 = anc.get_max_indices(lnprob_burnin)
  folder_2050, names_derived, der_posterior = get_intervals(cli.full_path, 2050, names_par, max_lnprob_parameters, flatchain_posterior_0, derived_type=None, full_output=True, summary_file_hdf5=s_h5f)
  units_der = anc.get_units(names_derived, mass_unit)
  # write out the derived names and posterior into an hdf5 file
  der_post_file = os.path.join(cli.full_path, 'derived_posterior.hdf5')
  h5f = h5py.File(der_post_file, 'w')
  h5f.create_dataset('derived_names', data=names_derived, dtype='S10')
  h5f.create_dataset('derived_posterior', data=der_posterior, dtype=np.float64)
  h5f.create_dataset('units_derived', data=units_der, dtype='S15')
  h5f.close()
  # ************************************
  
  ## COMPUTE CONFIDENCE INTERVALS OF THE DERIVED PARAMETER DISTRIBUTIONS
  ci_derived = np.percentile(der_posterior, anc.percentile_val[2:], axis=0, interpolation='midpoint') # (n_percentile-1 x nfit) ==> skip first item, the 68.27th
  
  s_h5f.create_dataset('confidence_intervals/derived/ci', data=ci_derived.T, dtype=np.float64)
  s_h5f.create_dataset('confidence_intervals/derived/names', data=names_derived, dtype='S10')
  s_h5f.create_dataset('confidence_intervals/derived/units', data=units_der, dtype='S15')
  s_h5f.create_dataset('confidence_intervals/derived/percentiles', data=np.array(anc.percentile_val[2:]), dtype=np.float64)
  
  print
  print
  
  # ************************************
  ## MEDIAN PARAMETERS ID == 1050
  median_parameters, median_perc68, median_confint = anc.get_median_parameters(flatchain_posterior_0)
  folder_1050 = get_intervals(cli.full_path, 1050, names_par, median_parameters, flatchain_posterior_0, derived_type='median', summary_file_hdf5=s_h5f)
  ## MEDIAN PARAMETERS ID == 1051
  folder_1051 = get_intervals(cli.full_path, 1051, names_par, median_parameters, flatchain_posterior_0, derived_type=None, summary_file_hdf5=s_h5f)
  # ************************************
  
  print
  print

  # ************************************
  ## MODE-LIKE PARAMETERS -> id 3050
  ## take the mean of 5 bin centered to the higher bin
  k = np.ceil(2. * flatchain_posterior_0.shape[0]**(1./3.)).astype(int)
  if(k>50): k=50
  mode_bin, mode_parameters = anc.get_mode_parameters(flatchain_posterior_0, k)
  folder_3050 = get_intervals(cli.full_path, 3050, names_par, mode_parameters, flatchain_posterior_0, derived_type='mode', summary_file_hdf5=s_h5f)
  ## MODE-LIKE PARAMETERS -> id 3051
  folder_3051 = get_intervals(cli.full_path, 3051, names_par, mode_parameters, flatchain_posterior_0, derived_type=None, summary_file_hdf5=s_h5f)
  # ************************************
  
  print
  print
  
  name_par, name_excluded = anc.get_sample_list(cli.sample_str, names_par)
  sample_parameters, idx_sample = anc.pick_sample_parameters(flatchain_posterior_0, names_par, name_par = name_par, name_excluded = name_excluded)
  if(sample_parameters is not None):
    folder_666 = get_intervals(cli.full_path, 666, names_par, sample_parameters, flatchain_posterior_0, idx_sample=idx_sample, summary_file_hdf5=s_h5f)
    s_h5f['parameters/%04d' %(666)].attrs['par_selection'] = name_par
    if(name_excluded is not None):
      s_h5f['parameters/%04d' %(666)].attrs['par_excluded'] = name_excluded
  else:
    print 'NONE SAMPLE PARAMETERS!!!'

  ## SELECT AD HOC PARAMETERS: K-19 median b,c, mode d
  #adhoc_par = median_parameters.copy()
  ##adhoc_par[10:] = mode_parameters[10:].copy()
  #adhoc_par[12] = mode_parameters[12].copy()
  #folder_777 = get_intervals(cli.full_path, 777, names_par, adhoc_par, flatchain_posterior_0, derived_type=None, summary_file_hdf5=s_h5f)
  
  
  s_h5f.close()

  print

  # print into file CONFIDENCE INTERVALS of fitted and derived parameters
  ci_file = os.path.join(cli.full_path, 'confidence_intervals.dat')
  oci = open(ci_file, 'w')
  anc.print_both('# CONFIDENCE INTERVALS',oci)
  
  anc.print_both('## FITTED PARAMETERS',oci)
  anc.print_confidence_intervals(anc.percentile_val[2:], conf_interv=ci_fitted, name_parameters=names_par, unit_parameters=units_par, output=oci)
  
  anc.print_both('## DERIVED PARAMETERS',oci)
  anc.print_confidence_intervals(anc.percentile_val[2:], conf_interv=ci_derived, name_parameters=names_derived, unit_parameters=units_der, output=oci)
  
  oci.close()

  pytrades_lib.pytrades.deallocate_variables()

  #print
  #print thin_steps
  #print np.shape(flatchain_posterior_0)
  #print np.shape(der_posterior)


  return

if __name__ == "__main__":
  main()
