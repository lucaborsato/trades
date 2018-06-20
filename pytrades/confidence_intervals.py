#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import sys
import argparse
import time
import os
import numpy as np # array
import h5py

# custom modules
script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), '../'))
sys.path.append(module_path)
import constants as cst # local constants module
import ancillary as anc
import pytrades_lib

# ==============================================================================

def save_ttra_and_rv_from_samples(samples_fit_par, names_par, NB, n_samples, smp_h5):
  
  n_samples, nfit = np.shape(samples_fit_par)
  Pephem = pytrades_lib.pytrades.pephem
  nt0_full, nrv_nmax = pytrades_lib.pytrades.get_max_nt0_nrv(Pephem,NB)
  
  for ismp in range(0, n_samples):
    smp_name = '%04d' %(ismp)
    gr = smp_h5.create_group(smp_name)
    gr.attrs['tepoch'] = pytrades_lib.pytrades.tepoch
    trades_fit_par = anc.sqrte_to_e_fitting(samples_fit_par[ismp], names_par)
    ttra_full, dur_full, id_ttra_full, stats_ttra, time_rv_nmax, rv_nmax, stats_rv = pytrades_lib.pytrades.fit_par_to_ttra_rv(trades_fit_par, nt0_full, nrv_nmax)
    
    stats_rv = np.array(stats_rv).astype(bool)
    time_rvx, rvx = time_rv_nmax[stats_rv], rv_nmax[stats_rv]
    id_rv = np.argsort(time_rvx)
    time_rv, rv = time_rvx[id_rv], rvx[id_rv]
    gr.create_dataset('time_rv_mod', data=time_rv,
                      dtype=np.float64, compression='gzip')
    gr.create_dataset('rv_mod', data=rv,
                      dtype=np.float64, compression='gzip')
    
    stats_ttra = np.array(stats_ttra).astype(bool)
    for ipl in range(2, NB+1):
      sel_tra = np.logical_and(np.logical_and(id_ttra_full == ipl, stats_ttra),
                               ttra_full > -9.e10
                               )
      ttra = np.array(ttra_full)[sel_tra]
      dur = np.array(dur_full)[sel_tra]
      idx_tra = np.argsort(ttra)
      ttra = ttra[idx_tra]
      dur = dur[idx_tra]
      gr.create_dataset('TTs_%d' %(ipl), data=ttra,
                        dtype=np.float64, compression='gzip'
                        )
      gr.create_dataset('T41s_%d' %(ipl), data=dur,
                        dtype=np.float64, compression='gzip'
                        )
  return

# ==============================================================================

# main

def main():
  print
  print ' TRADES: EMCEE confidence intervals'
  print

  cli = anc.get_args()

  # init trades
  pytrades_lib.pytrades.initialize_trades(os.path.join(cli.full_path, ''), '', 1)

  nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case_list = anc.get_fitted(cli.full_path)
  ndata = pytrades_lib.pytrades.ndata
  nfree = pytrades_lib.pytrades.nfree
  dof = pytrades_lib.pytrades.dof
  
  # read emcee data
  emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(cli.full_path, cli.temp_status)
  # get data from the hdf5 file
  names_par, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)
  # print Memory occupation of ...
  anc.print_memory_usage(chains)

  nfit, nwalkers, nruns, nburnin, nruns_sel = anc.get_emcee_parameters(chains, cli.temp_status, cli.nburnin, completed_steps)
  
  #chains_T, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, nburnin, nruns, nruns_sel, m_factor, names_par, parameter_boundaries, chains)
  chains_T_full = np.zeros((nruns, nwalkers, nfit))
  for ii in xrange(0,nfit):
    chains_T_full[:,:,ii] = chains[:,:nruns,ii].T # transpose
    
  chains_T, flatchain_posterior_0, lnprob_burnin, thin_steps = anc.thin_the_chains(cli.use_thin, nburnin, nruns, nruns_sel, autocor_time, chains_T_full, lnprobability, burnin_done=False)
  
  # computes mass conversion factor
  #m_factor = anc.mass_conversion_factor(cli.m_type)
  MR_star = pytrades_lib.pytrades.mr_star
  m_factor_0, mass_unit = anc.mass_type_factor(Ms=1.0, mtype=cli.m_type, mscale=False)
  np.random.seed(seed=cli.seed)
  Ms_gaussian = MR_star[0,0] + np.random.normal(0., 1., size=(np.shape(flatchain_posterior_0)[0]))*MR_star[0,1] # if exists an error on the mass, it creates a Normal distribution for the values and use it to re-scale mp/Ms to mp.
  m_factor_boot = m_factor_0 * Ms_gaussian # given the factor from Msun to mass_unit it multiply it by the Normal Mstar.
  m_factor = m_factor_0 * MR_star[0,0]

  # set label and legend names
  #kel_legends, labels_list = anc.keplerian_legend(names_par, cli.m_type)
  
  flatchain_posterior = flatchain_posterior_0.copy()
  for ifit in range(0,nfit):
    if('Ms' in names_par[ifit]):
      flatchain_posterior[:,ifit] = m_factor_0 * flatchain_posterior[:,ifit]
  posterior_file = os.path.join(cli.full_path, 'posterior.hdf5')
  p_h5f = h5py.File(posterior_file, 'w')
  p_h5f.create_dataset('posterior', data=flatchain_posterior, dtype=np.float64)
  p_h5f['posterior'].attrs['nfit'] = nfit
  p_h5f['posterior'].attrs['nposterior'] = np.shape(flatchain_posterior)[0]
  p_h5f.create_dataset('parameter_names', data=names_par, dtype='S10')
  p_h5f.close()
  anc.print_both(' Saved posterior file: %s' %(posterior_file))

  top_header, header = anc.get_header(anc.percentile_val)

# ==============================================================================
# ==============================================================================
# 2017-01-26 EMCEE NOW USED sqrt(e)cos(w), sqrt(e)sin(w)
# GET INTERVALS
# ==============================================================================
# ==============================================================================
  def get_intervals(full_path, id_sim, names_par_in, parameters_in, flatchain_posterior_in, derived_type=None, full_output=False, idx_sample=None, summary_file_hdf5=None):
    
    names_trades = anc.emcee_names_to_trades(names_par_in) # emcee to trades
    parameters_trades = anc.sqrte_to_e_fitting(parameters_in, names_par_in) # emcee to trades
    
    names_par = names_par_in # emcee kind
    parameters = parameters_in # emcee kind
    flatchain_posterior = flatchain_posterior_in # emcee kind
    
    loglhdx, checkx = pytrades_lib.pytrades.fortran_loglikelihood(np.array(parameters_trades,dtype=np.float64))
    loglhdx = loglhdx + ln_err_const
    
    out_folder = os.path.join(os.path.join(full_path, '%04d_sim' %(id_sim)), '')
    if (not os.path.isdir(out_folder)):
      os.makedirs(out_folder)
    out_file = os.path.join(out_folder, 'parameters_summary.txt')
    out = open(out_file, 'w')
    pytrades_lib.pytrades.path_change(out_folder)
    
    anc.print_both(' ', out)
    anc.print_both(' --------------------------------- ', out)
    anc.print_both(' PARAMETER VALUES -> %d' %(id_sim), out)
    fitness, lgllhd, check = pytrades_lib.pytrades.write_summary_files(id_sim, parameters_trades)
    
    kel_file, kep_elem = anc.elements(out_folder, id_sim, lmf=0)

    #sigma_par = anc.compute_intervals(flatchain_posterior, parameters, anc.percentile_val)
    sigma_par = anc.compute_sigma_hdi(flatchain_posterior, parameters) # uses HDI
    sigma_par = sigma_par.T
    units_par = anc.get_units(names_par, mass_unit)
    
    if(not bool(check)):
      print 'WRTING WARNING FILE: %s' %(os.path.join(out_folder,'WARNING.txt'))
      warn_o = open(os.path.join(out_folder,'WARNING.txt'), 'w')
      warn_o.write('*******\nWARNING: FITTED PARAMETERS COULD NOT BE PHYSICAL!\nWARNING: BE VERY CAREFUL WITH THIS PARAMETER SET!\n*******')
      warn_o.close()
    
    nbins = anc.get_auto_bins(flatchain_posterior_0)
    
    names_derived, der_posterior = anc.compute_derived_posterior(names_par, kep_elem, id_fit, case_list, cols_list, flatchain_posterior, conv_factor=m_factor_boot)
    
    #der_posterior_T = der_posterior
    der_posterior_T = anc.derived_posterior_check(names_derived, der_posterior)

    par_type = ''
    descr = ''
    if(str(derived_type).strip().lower() == 'median'):
      # MEDIAN PARAMETERS ID == 1050
      derived_par = np.percentile(der_posterior_T, 50., axis=0, interpolation='midpoint')
      par_type = 'MEDIAN:'
      descr = 'median of posterior and median of derived posterior'
    elif(str(derived_type).strip().lower() == 'mode'):
      # MODE-LIKE PARAMETERS -> id 3050
      #k = anc.get_bins(flatchain_posterior, rule='doane')
      
      der_bin, derived_par = anc.get_mode_parameters(der_posterior_T, nbins)
      par_type = 'MODE'
      descr = 'mode of posterior and mode of derived posterior'
    else:      
      # ORIGINAL FITTING PARAMETERS ID == 0
      # or
      # MAX LNPROBABILITY -> id 2050
      names_derived, derived_par = anc.compute_derived_parameters(names_par, kep_elem, id_fit, case_list, cols_list, parameters, conv_factor=m_factor)
      derived_par, der_posterior_T = anc.adjust_derived_parameters(names_derived, derived_par, der_posterior_T)
      if(id_sim == 0):
        par_type = 'ORIGINAL FIT:'
        descr = 'initial set of parameters'
      elif(id_sim == 1051):
        par_type = 'MEDIAN PARAMETERS TO DERIVED:'
        descr = 'median of posterior and converted to derived parameter'
      elif(id_sim == 2050):
        par_type = 'MAX LNPROB'
      elif(id_sim == 3051):
        par_type = 'MODE PARAMETERS TO DERIVED:'
        descr = 'mode of posterior and converted to derived parameter'
      elif(id_sim == 666):
        par_type = 'SELECTED SAMPLE WITHIN HDI'
        # ***COMMENTED 2017-02-02: TO CHECK IF REALLY NEEDED
        #if(idx_sample is not None):
          #par_type = '%s <-> idx = %d' %(par_type, idx_sample)
          #derived_par = der_posterior_T[idx_sample, :]
          #for ider in range(0,np.shape(derived_par)[0]):
            ##print ider, names_derived[ider], names_derived[ider][0], names_derived[ider][1]
            #if(names_derived[ider][0] == 'm' and names_derived[ider][1] != 'A'):
              ##print 'doing'
              #derived_par[ider] = der_posterior_T[idx_sample, ider]*m_factor/m_factor_boot[idx_sample]
      elif(id_sim == 667):
        par_type = 'SELECTED SAMPLE CLOSE TO MEDIAN LGLLHD WITHIN POSTERIOR HDI'
        descr = ""
      elif(id_sim == 668):
        par_type = 'MAX LGLLHD WITHIN POSTERIOR HDI:'
        descr = "Select posterior within HDI and take the parameter set with higher loglikelihood."
      else:
        par_type = 'AD HOC'
        descr = "from input file"
    
    par_type = '%s %s' %(par_type, descr)
    #sigma_derived = anc.compute_intervals(der_posterior_T, derived_par, anc.percentile_val)
    sigma_derived = anc.compute_sigma_hdi(der_posterior_T, derived_par)
    sigma_derived = sigma_derived.T
    
    units_der = anc.get_units(names_derived, mass_unit)
    
    if(s_h5f is not None):
      s_id_sim = '%04d' %(id_sim)
      s_h5f.create_dataset('parameters/%s/fitted/parameters' %(s_id_sim), data=parameters,
                           dtype=np.float64, compression='gzip')
      s_h5f.create_dataset('parameters/%s/fitted/names' %(s_id_sim), data=names_par,
                           dtype='S10', compression='gzip')
      s_h5f.create_dataset('parameters/%s/fitted/units' %(s_id_sim), data=units_par,
                           dtype='S15', compression='gzip')
      s_h5f.create_dataset('parameters/%s/fitted/sigma' %(s_id_sim), data=sigma_par.T,
                           dtype=np.float64, compression='gzip')
      s_h5f['parameters/%s/fitted/sigma' %(s_id_sim)].attrs['percentiles'] = anc.percentile_val

      s_h5f.create_dataset('parameters/%s/derived/parameters' %(s_id_sim), data=derived_par,
                           dtype=np.float64, compression='gzip')
      s_h5f.create_dataset('parameters/%s/derived/names' %(s_id_sim), data=names_derived,
                           dtype='S10', compression='gzip')
      s_h5f.create_dataset('parameters/%s/derived/units' %(s_id_sim), data=units_der,
                           dtype='S15', compression='gzip')
      s_h5f.create_dataset('parameters/%s/derived/sigma' %(s_id_sim), data=sigma_derived.T,
                           dtype=np.float64, compression='gzip')
      s_h5f['parameters/%s/derived/sigma' %(s_id_sim)].attrs['percentiles'] = anc.percentile_val

      s_h5f['parameters/%s' %(s_id_sim)].attrs['info'] = '%s ==> %s' %(s_id_sim, par_type)
      s_h5f['parameters/%s' %(s_id_sim)].attrs['fitness'] = fitness
      s_h5f['parameters/%s' %(s_id_sim)].attrs['lgllhd'] = lgllhd
      s_h5f['parameters/%s' %(s_id_sim)].attrs['check'] = check
      if(idx_sample is not None):
        s_h5f['parameters/%s' %(s_id_sim)].attrs['idx_sample'] = idx_sample
    
    #print '\nComputed sigma_par with shape ',np.shape(sigma_par)
    #print 'Computed sigma_derived with shape ',np.shape(sigma_derived)
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
# ==============================================================================
# ==============================================================================

# ==============================================================================
## CREATE A HDF5 FILE WITH CONFIDNCE INTERVALS AND ALL THE SUMMARY PARAMETERS
# ==============================================================================
  summary_file = os.path.join(cli.full_path, 'summary_parameters.hdf5')
  s_h5f = h5py.File(summary_file, 'w')
  
  ### COMPUTE CONFIDENCE INTERVALS OF THE FITTED PARAMETER DISTRIBUTIONS
  #ci_fitted = np.percentile(flatchain_posterior_0, anc.percentile_val[2:], axis=0, interpolation='midpoint') # (n_percentile-2 x nfit) ==> skip 1st and 2nd items, the 68.27th and 50th
# ==============================================================================
# HDI INSTEAD OF CREDIBLE INTERVALS
# ==============================================================================
  nbins = anc.get_auto_bins(flatchain_posterior_0)
  hdi_ci, mode_parameters = anc.compute_hdi_full(flatchain_posterior_0, mode_output=True)
  ci_fitted = hdi_ci.T
  print ' shape: hdi_ci = ', np.shape(hdi_ci), ' ci_fitted = ', np.shape(ci_fitted)
  # hdi_ci: nfit x nci
  # ci_fitted: nci x nfit
  # nci -> -1sigma(0) +1sigma(1) -2sigma(2) +2sigma(3) -3sigma(4) +3sigma(5)
  
  #sys.exit()
  
  units_par = anc.get_units(names_par, mass_unit)
  
  s_h5f.create_dataset('confidence_intervals/fitted/ci', data=ci_fitted.T, 
                       dtype=np.float64, compression='gzip')
  s_h5f.create_dataset('confidence_intervals/fitted/names', data=names_par,
                       dtype='S10', compression='gzip')
  s_h5f.create_dataset('confidence_intervals/fitted/units', data=units_par,
                       dtype='S15', compression='gzip')
  s_h5f.create_dataset('confidence_intervals/fitted/percentiles', data=np.array(anc.percentile_val[2:]), 
                       dtype=np.float64, compression='gzip') # now it not true...
  
  s_h5f['confidence_intervals/fitted'].attrs['nfit'] = nfit
  s_h5f['confidence_intervals/fitted'].attrs['nfree'] = nfree
  s_h5f['confidence_intervals/fitted'].attrs['ndata'] = ndata
  s_h5f['confidence_intervals/fitted'].attrs['dof'] = dof
# ==============================================================================
  
# ==============================================================================
## ORIGINAL FITTING PARAMETERS ID == 0
# ==============================================================================
  # save initial_fitting parameters into array
  original_fit_parameters = pytrades_lib.pytrades.fitting_parameters # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
  #folder_0 = get_intervals(cli.full_path, 0, names_par, original_fit_parameters, flatchain_posterior_0, derived_type=None, summary_file_hdf5=s_h5f)
  # WARNING: original_fit_parameters from TRADES have to be converted to emcee parameters:
  # (ecosw,esinw) => (sqrecosw, sqrtesinw)
  trades_names = anc.emcee_names_to_trades(names_par)
  original_parameters = anc.e_to_sqrte_fitting(original_fit_parameters, trades_names)
  folder_0 = get_intervals(cli.full_path, 0, names_par, original_parameters, flatchain_posterior_0, derived_type=None, summary_file_hdf5=s_h5f)
# ==============================================================================
  
  print
  print
  
# ==============================================================================
## MAX LNPROBABILITY AND PARAMETERS -> id 2050
# ==============================================================================
  max_lnprob, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint = anc.get_maxlnprob_parameters(lnprob_burnin, chains_T, flatchain_posterior_0)
  max_id1, max_id2 = anc.get_max_indices(lnprob_burnin)
  folder_2050, names_derived, der_posterior = get_intervals(cli.full_path, 2050, names_par, max_lnprob_parameters, flatchain_posterior_0, derived_type=None, full_output=True, summary_file_hdf5=s_h5f)
  units_der = anc.get_units(names_derived, mass_unit)
  # write out the derived names and posterior into an hdf5 file
  der_post_file = os.path.join(cli.full_path, 'derived_posterior.hdf5')
  h5f = h5py.File(der_post_file, 'w')
  h5f.create_dataset('derived_names', data=names_derived,
                     dtype='S10', compression='gzip')
  h5f.create_dataset('derived_posterior', data=der_posterior,
                     dtype=np.float64, compression='gzip')
  h5f.create_dataset('units_derived', data=units_der,
                     dtype='S15', compression='gzip')
  h5f.close()
# ==============================================================================
  
  ### COMPUTE CONFIDENCE INTERVALS OF THE DERIVED PARAMETER DISTRIBUTIONS
  #ci_derived = np.percentile(der_posterior, anc.percentile_val[2:], axis=0, interpolation='midpoint') # (n_percentile-1 x nder) ==> skip first item, the 68.27th
# ==============================================================================
# HDI INSTEAD OF CREDIBLE INTERVALS
# ==============================================================================
  #npost_der, nder = np.shape(der_posterior)
  #k_der = anc.get_auto_bins(der_posterior)
  hdi_ci_derived = anc.compute_hdi_full(der_posterior, mode_output=False)
  ci_derived = hdi_ci_derived.T
  
  s_h5f.create_dataset('confidence_intervals/derived/ci', 
                       data=ci_derived.T,
                       dtype=np.float64, compression='gzip')
  s_h5f.create_dataset('confidence_intervals/derived/names', 
                       data=names_derived, 
                       dtype='S10', compression='gzip')
  s_h5f.create_dataset('confidence_intervals/derived/units', 
                       data=units_der,
                       dtype='S15', compression='gzip')
  s_h5f.create_dataset('confidence_intervals/derived/percentiles',
                       data=np.array(anc.percentile_val[2:]),
                       dtype=np.float64, compression='gzip')
# ==============================================================================  
  print
  print
  
# ==============================================================================
## MEDIAN PARAMETERS ID == 1050
# ==============================================================================
  median_parameters, median_perc68, median_confint = anc.get_median_parameters(flatchain_posterior_0)
  folder_1050 = get_intervals(cli.full_path, 1050, names_par, median_parameters, flatchain_posterior_0, derived_type='median', summary_file_hdf5=s_h5f)
  ## MEDIAN PARAMETERS ID == 1051
  folder_1051 = get_intervals(cli.full_path, 1051, names_par, median_parameters,
                              flatchain_posterior_0, derived_type=None,
                              summary_file_hdf5=s_h5f
                              )
# ==============================================================================
  
  print
  print

# ==============================================================================
# select n_samples from the posterior within the CI
# ==============================================================================
  if(cli.n_samples > 0):
    anc.print_both(' Selecting %d samples from the posterior ...' %(cli.n_samples))
    sys.stdout.flush()
    samples_fit_par = anc.take_n_samples(flatchain_posterior_0, ci_fitted[0:2,:],
                                         n_samples=cli.n_samples
                                         )
    samples_fit_par[0,:] = median_parameters # first sample as the median of the posterior
    anc.print_both(' Running TRADES and computing the T0s and RVs ...')
    samples_file = os.path.join(cli.full_path, 'samples_ttra_rv.hdf5')
    anc.print_both(' Saving into %s' %(samples_file))
    smp_h5 = h5py.File(samples_file, 'w')
    save_ttra_and_rv_from_samples(samples_fit_par, names_par, NB, cli.n_samples, smp_h5)
    #tra_gr = smp_h5.create_group('T0')
    #for key in ttra_samples.keys():
      #tra_gr.create_dataset(key, data=ttra_samples[key], dtype=np.float64, compression='gzip')
    #rv_gr =  smp_h5.create_group('RV')
    #for key in rv_samples.keys():
      #rv_gr.create_dataset(key, data=rv_samples[key], dtype=np.float64, compression='gzip')
    #rv_gr['time_rv_mod'].attrs['tepoch'] = pytrades_lib.pytrades.tepoch
    smp_h5.close()
    anc.print_both(' ... done')
    sys.stdout.flush()
  #sys.exit()
# ==============================================================================

  print
  print

# ==============================================================================
## MODE-LIKE PARAMETERS -> id 3050
# ==============================================================================
  ## take the mean of 5 bin centered to the higher bin

  #anc.print_both('nbins = %d' %(nbins))
  #sys.stdout.flush()
  #mode_bin, mode_parameters = anc.get_mode_parameters(flatchain_posterior_0, nbins)
  # mode_parameters computed at the beginning with hdi
  
  if(np.any(np.isnan(mode_parameters))):
    print 'Some values are Nan, skip the mode parameters'
  else:
    folder_3050 = get_intervals(cli.full_path, 3050, names_par, mode_parameters,
                                flatchain_posterior_0, derived_type='mode', 
                                summary_file_hdf5=s_h5f
                                )
    ## MODE-LIKE PARAMETERS -> id 3051
    folder_3051 = get_intervals(cli.full_path, 3051, names_par, mode_parameters,
                                flatchain_posterior_0, derived_type=None,
                                 summary_file_hdf5=s_h5f
                                 )
# ==============================================================================
  
  print
  print

# ==============================================================================
# ONE SAMPLE PARAMETER SET --> 666
# ==============================================================================
  name_par, name_excluded = anc.get_sample_list(cli.sample_str, names_par)
  sample_parameters, idx_sample = anc.pick_sample_parameters(flatchain_posterior_0,
                                                             names_par, 
                                                             name_par = name_par,
                                                             name_excluded = name_excluded,
                                                             post_ci=ci_fitted[0:2,:]
                                                             )
  if(sample_parameters is not None):
    folder_666 = get_intervals(cli.full_path, 666, names_par, sample_parameters,
                               flatchain_posterior_0, idx_sample=idx_sample, 
                               summary_file_hdf5=s_h5f
                               )
    s_h5f['parameters/%04d' %(666)].attrs['par_selection'] = name_par
    if(name_excluded is not None):
      s_h5f['parameters/%04d' %(666)].attrs['par_excluded'] = name_excluded
  else:
    print 'NONE SAMPLE PARAMETERS!!!'
# ==============================================================================

# ==============================================================================
## SELECT AD HOC PARAMETERS:
# ==============================================================================
  #adhoc_par = median_parameters.copy()
  ##adhoc_par[10:] = mode_parameters[10:].copy()
  #adhoc_par[12] = mode_parameters[12].copy()
  #if(cli.overplot is not None):
  if( cli.adhoc is not None):
    print cli.overplot, cli.adhoc
    adhoc_names, adhoc_par_trades = anc.read_fitted_file(cli.adhoc)
    adhoc_par = anc.e_to_sqrte_fitting(adhoc_par_trades, adhoc_names)
    folder_777 = get_intervals(cli.full_path, 777, names_par, adhoc_par,
                               flatchain_posterior_0, derived_type=777,
                                summary_file_hdf5=s_h5f
                                )
# ==============================================================================

# ==============================================================================
# select the sample within post_ci and close to median lgllhd --> 667
# ==============================================================================
  sample2_parameters, sample2_lgllhd = anc.get_sample_by_sorted_lgllhd(
                                       flatchain_posterior_0,
                                       lnprob_burnin.T,
                                       #post_ci = ci_fitted[0:2,:]
                                       post_ci = ci_fitted.T
                                       )
  folder_667 = get_intervals(cli.full_path, 667, names_par, sample2_parameters, flatchain_posterior_0, derived_type=667, summary_file_hdf5=s_h5f)
# ==============================================================================
  
# ==============================================================================
# another kind of selection: parameter set withi HDI, then take the max(loglikelihood) --> 668
# ==============================================================================
  name_par, name_excluded = anc.get_sample_list(cli.sample_str, names_par)
  #sample3_parameters, sample3_lgllhd = anc.get_sample_by_par_and_lgllhd(flatchain_posterior_0,
                                       #lnprob_burnin.T, 
                                       #names_par, 
                                       #post_ci = ci_fitted[0:2,:], 
                                       #name_par= name_par)
  sample3_parameters, sample3_lgllhd = \
                             anc.select_maxlglhd_with_hdi(flatchain_posterior_0,
                                                          #ci_fitted[0:2,:],
                                                          ci_fitted.T,
                                                          lnprob_burnin.T
                                                          )
  folder_668 = get_intervals(cli.full_path, 668, names_par, sample3_parameters,
                             flatchain_posterior_0, derived_type=668,
                             summary_file_hdf5=s_h5f
                             )
# ==============================================================================

  s_h5f.close()

  print

# ==============================================================================
# print into file CONFIDENCE INTERVALS of fitted and derived parameters
# ==============================================================================
  ci_file = os.path.join(cli.full_path, 'confidence_intervals.dat')
  oci = open(ci_file, 'w')
  anc.print_both('\n# SUMMARY:\n# CONFIDENCE INTERVALS',oci)
  
  anc.print_both('## FITTED PARAMETERS',oci)
  #anc.print_confidence_intervals(anc.percentile_val[2:], conf_interv=ci_fitted, name_parameters=names_par, unit_parameters=units_par, output=oci)
  anc.print_hdi(conf_interv=ci_fitted, name_parameters=names_par,
                unit_parameters=units_par, output=oci)
  
  anc.print_both('## DERIVED PARAMETERS',oci)
  #anc.print_confidence_intervals(anc.percentile_val[2:], conf_interv=ci_derived, name_parameters=names_derived, unit_parameters=units_der, output=oci)
  anc.print_hdi(conf_interv=ci_derived, name_parameters=names_derived,
                unit_parameters=units_der, output=oci)
  
  oci.close()
# ==============================================================================

  pytrades_lib.pytrades.deallocate_variables()

  return

# ==============================================================================
# ==============================================================================

if __name__ == "__main__":
  main()
