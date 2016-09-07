#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import sys
import argparse
import os
import numpy as np # array
import h5py
#import trades_lib
#import random
import constants as cst # local constants module

# common variables needed to create labels and parameter names
kel_id = ['$M$', '$R$', '$P$', '$e$', '$\omega$', '$\\nu$', '$i$', '$\Omega$']
kel_units = ['$[M_{\oplus}]$', '$[R_\mathrm{Jup}]$', '[days]', '', '', '[deg]', '[deg]', '[deg]']

kel_id_2 = ['$M/M_\star$', '$R$', '$P$', '$e\cos\omega$', '$e\sin\omega$', '$\lambda$', '$i\cos\Omega$', '$i\sin\Omega$']
kel_units_2 = ['', '$[R_\mathrm{Jup}]$', '[days]', '', '', '[deg]', '', '']

kel_fmt = ['%.3f', '%.4f', '%.3f', '%.1f', '%.1f', '%.1f', '%.3f', '%.3f']
nelem = len(kel_id)
letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'l', 'm', 'n', 'o', 'p']

def set_bool_argument(arg):
  if (arg != False):
    if (arg.lower() in ['t', 'tr', 'tru', 'true']):
      arg = True
    else:
      arg = False
  return arg

# read command line (cli) arguments
def get_args():
  parser = argparse.ArgumentParser(description='TRADES+EMCEE PLOT')
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')
  parser.add_argument('-nb', '--nburn', '-np', '--npost', action='store', dest='npost', required=True, help='The number of posterior/burn in steps to discard at the beginning of each chain. It has to be > 0')
  parser.add_argument('-m', '--mtype', '--mass-type', action='store', dest='m_type', default='e', help='Mass type: j = Jupiter, e = Earth, s = Sun. Default is Earth = e.')
  parser.add_argument('-g', '--good-parameters', action='store', dest='good', default=False, help='If you want to use a previous solution, set it to True and it will search for a good_parameters.dat file with first column the name of the parameter and the value in the second. Mass parameters in Jupiter mass. Default is False')
  parser.add_argument('-t', '--temp-file', action='store', dest='temp_status', default=False, help='If you want to read temporary emcee_temp.hdf5 file, because simulation is not finished yet. Default is False')
  
  cli = parser.parse_args()
  cli.full_path = os.path.abspath(cli.full_path)
  cli.m_type = cli.m_type.lower()
  
  cli.good = set_bool_argument(cli.good)

  cli.temp_status = set_bool_argument(cli.temp_status)
  
  return cli

# given the mass flag letter it computes the proper mass conversion factor
def mass_conversion_factor(m_type):
  if (m_type in ['j', 'jup', 'jupiter']):
    m_factor = cst.Msjup
  elif (m_type in ['s', 'sun']):
    m_factor = 1.
  else:
    m_factor = cst.Msear
  return m_factor

def mass_factor_unit(m_type):
  if (m_type in ['j', 'jup', 'jupiter']):
    m_factor = cst.Msjup
    m_unit = 'M_Jup'
  elif (m_type in ['s', 'sun']):
    m_factor = 1.
    m_unit = 'M_Sun'
  else:
    m_factor = cst.Msear
    m_unit = 'M_Earth'
  return m_factor, m_unit


# prepare the labels of the Keplerian orbital elements for the legend
def keplerian_legend(parameter_names, m_type):
  if (m_type in ['j', 'jup', 'jupiter']):
    kel_units[0] = '$[M_\mathrm{Jup}]$'
  elif (m_type in ['s', 'sun']):
    kel_units[0] = '$[M_{\odot}]$'
  nfit = parameter_names.shape[0]
  kel_legends = np.zeros((nfit), dtype='|S256')
  labels_list = []
  for i in range(0, nfit):
    #print parameter_names[i], ' ==> ',
    parameter_names[i] = parameter_names[i].strip()
    #print parameter_names[i],
    
    if ('Ms' in parameter_names[i]):
      planet_id = int(parameter_names[i].split('m')[1].split('Ms')[0]) - 1
      kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id_2[0], letters[planet_id], kel_units_2[0])
      labels_list.append(r'$%s_%s$' %(kel_id_2[0].strip('$'), letters[planet_id]))
      
    elif (parameter_names[i][0] == 'm' and parameter_names[i][-1] != 's'):
      planet_id = int(parameter_names[i][1:]) - 1
      kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id[0], letters[planet_id], kel_units[0])
      labels_list.append(r'$%s_%s$' %(kel_id[0].strip('$'), letters[planet_id]))
    
    elif (parameter_names[i][0] == 'R'):
      planet_id = int(parameter_names[i][1:]) - 1
      kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id[1], letters[planet_id], kel_units[1])
      labels_list.append(r'$%s_%s$' %(kel_id_2[1].strip('$'), letters[planet_id]))
    
    elif (parameter_names[i][0] == 'P'):
      planet_id = int(parameter_names[i][1:]) - 1
      kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id[2], letters[planet_id], kel_units[2])
      labels_list.append(r'$%s_%s$' %(kel_id_2[2].strip('$'), letters[planet_id]))
      
    elif ('e' in parameter_names[i]):
      if (parameter_names[i][0:4] == 'ecos'):
        planet_id = int(parameter_names[i].split('w')[1].strip()) - 1
        kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id_2[3], letters[planet_id], kel_units_2[3])
        labels_list.append(r'$%s_%s$' %(kel_id_2[3].strip('$'), letters[planet_id]))
      elif (parameter_names[i][0:4] == 'esin'):
        planet_id = int(parameter_names[i].split('w')[1].strip()) - 1
        kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id_2[4], letters[planet_id], kel_units_2[4])
        labels_list.append(r'$%s_%s$' %(kel_id_2[4].strip('$'), letters[planet_id]))
      else:
        planet_id = int(parameter_names[i].split('e')[1].strip()) - 1
        kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id[3], letters[planet_id], kel_units[3])
        labels_list.append(r'$%s_%s$' %(kel_id[3].strip('$'), letters[planet_id]))
    
    elif (parameter_names[i][0] == 'w'):
      planet_id = int(parameter_names[i].split('w')[1].strip()) - 1
      kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id[4], letters[planet_id], kel_units[4])
      labels_list.append(r'$%s_%s$' %(kel_id[4].strip('$'), letters[planet_id]))
    
    elif (parameter_names[i][0:2] == 'mA'):
      planet_id = int(parameter_names[i][2:]) - 1
      kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id[5], letters[planet_id], kel_units[5])
      labels_list.append(r'$%s_%s$' %(kel_id[5].strip('$'), letters[planet_id]))
      
    elif ('lambda' in parameter_names[i]):
      planet_id = int(parameter_names[i].split('lambda')[1].strip()) - 1
      kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id_2[5], letters[planet_id], kel_units_2[5])
      labels_list.append(r'$%s_%s$' %(kel_id_2[5].strip('$'), letters[planet_id]))
      
    elif (parameter_names[i][0] == 'i'):
      if('icos' in parameter_names[i]):
        planet_id = int(parameter_names[i].split('lN')[1].strip()) - 1
        kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id_2[6], letters[planet_id], kel_units_2[6])
        labels_list.append(r'$%s_%s$' %(kel_id_2[6].strip('$'), letters[planet_id]))
      elif('isin' in parameter_names[i]):
        planet_id = int(parameter_names[i].split('lN')[1].strip()) - 1
        kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id_2[7], letters[planet_id], kel_units_2[7])
        labels_list.append(r'$%s_%s$' %(kel_id_2[7].strip('$'), letters[planet_id]))
      else:
        planet_id = int(parameter_names[i][1:]) - 1
        kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id[6], letters[planet_id], kel_units[6])
        labels_list.append(r'$%s_%s$' %(kel_id_2[6].strip('$'), letters[planet_id]))
      
    elif (parameter_names[i][0:2] == 'lN'):
      planet_id = int(parameter_names[i][2:]) - 1
      kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id[7], letters[planet_id], kel_units[7])
      labels_list.append(r'$%s_%s$' %(kel_id_2[7].strip('$'), letters[planet_id]))
      
      
    #print ' ==> label = ',labels_list[-1]
    
    
  kel_legends = [i.strip() for i in kel_legends]
  return kel_legends, labels_list

# only needed by check_good_parameters to convert Mjup to flagged mass
def check_good_parameters(good_id, good_parameters_in, m_factor, nfit):
  good_parameters_out = np.zeros(nfit) + good_parameters_in
  for i in range(0, nfit):
    if (good_id[i].strip()[0] == 'm' and good_id[i].strip()[1] != 'A'):
      good_parameters_out[i] = good_parameters_in[i] * cst.Mjups * m_factor
  return good_parameters_out

def read_check_good_parameters(full_path, good_status, m_factor, nfit):
  # read good parameters file
  good_parameters_0, good_parameters_1 = False, False
  if (good_status):
    good_file = os.path.join(full_path, 'good_parameters.dat')
    if (os.path.exists(good_file)):
      good_parameters_0 = np.genfromtxt(good_file, usecols=(1), dtype=np.float64)
      good_id = np.genfromtxt(good_file, usecols=(0), dtype='|S10')
      good_parameters_1 = check_good_parameters(good_id, good_parameters_0, m_factor, nfit)
      #for ii in range(0,nfit):
        #if (parameter_names_emcee[ii][:2] == 'mA'):
          #good_parameters_0[ii] = rescale_angle(good_parameters[ii])
  return good_parameters_0, good_parameters_1


# sometimes the angle distribution are bimodal because they are centered in the wrong position
# e.g.: -1° = 359° etc...
def rescale_angle(angle):

  # with arctan2
  cosa = np.cos(angle*np.pi/180.)
  sina = np.sin(angle*np.pi/180.)
  new_angle = np.arctan2(sina,cosa)*180./np.pi
  return new_angle

# renormalize the angular parameter mean Anomaly mA
def renormalize_parameters(parameters, parameter_names):
  new_parameters = parameters.copy()
  nfit = np.array(new_parameters).shape[0]
  for i in range(0, nfit):
    if (parameter_names[i][:2] == 'mA'):
      new_parameters[i] = (parameters[i] + 360.) % 360.
  return new_parameters

# get indices of max of an array
def get_max_indices(array_values):
  (idx1, idx2) = np.unravel_index(array_values.argmax(), array_values.shape)
  return idx1, idx2

# prepare file and best folder emcee
def get_emcee_file_and_best(emcee_folder,temp_status):
  if (temp_status):
    folder_best = 'best_temp'
    emcee_file = os.path.join(emcee_folder, 'emcee_temp.hdf5')
    emcee_best = os.path.join(emcee_folder, folder_best)
  else:
    folder_best = 'best'
    emcee_file = os.path.join(emcee_folder, 'emcee_summary.hdf5')
    emcee_best = os.path.join(emcee_folder, folder_best)
  return emcee_file, emcee_best, folder_best

def get_devol_file(emcee_folder):
  devol_file = os.path.join(emcee_folder, 'best_devol.hdf5')
  return devol_file

def get_percentile_angle(angle_posterior):
  cosa_posterior = np.cos(angle_posterior*np.pi/180.)
  sina_posterior = np.sin(angle_posterior*np.pi/180.)
  temp_cos = np.percentile(cosa_posterior, 50)
  temp_sin = np.percentile(sina_posterior, 50)
  median_angle = ((np.arctan2(temp_sin, temp_cos)*180./np.pi)+360.)%360.
  lower_angle = np.percentile(angle_posterior, 16)
  upper_angle = np.percentile(angle_posterior, 84)
  return median_angle, lower_angle, upper_angle


def get_data(emcee_file, temp_status):
  # read data from hdf5 file
  completed_steps = 0
  #print ' reading file', emcee_file
  f_read = h5py.File(emcee_file, 'r')
  data_names = [str(i) for i in f_read.keys()]
  parameter_names_emcee, parameter_boundaries = [], []
  chains = []
  acceptance_fraction, autocor_time, lnprobability = [], [], []

  if ('parameter_names' in data_names):   parameter_names_emcee = np.array(f_read['parameter_names'], dtype='S10')
  if ('boundaries' in data_names):        parameter_boundaries = np.array(f_read['boundaries'], dtype=np.float64)
  #if ('final_parameters' in data_names):  final_parameters = np.array(f_read['final_parameters'], dtype=np.float64)
  if ('chains' in data_names):            chains = np.array(f_read['chains'], dtype=np.float64) # shape (nwalkers, nrusn, nfit)
  if ('acceptance_fraction' in data_names): acceptance_fraction = np.array(f_read['acceptance_fraction'], dtype=np.float64)
  if ('autocor_time' in data_names): autocor_time = np.array(f_read['autocor_time'], dtype=np.float64)
  if ('lnprobability' in data_names): lnprobability = np.array(f_read['lnprobability'], dtype=np.float64)
  try:
    ln_err_const = f_read['lnprobability'].attrs['ln_err_const']
  except:
    ln_err_const = 0.0
  if (temp_status):  completed_steps = f_read['chains'].attrs['completed_steps']
  # close hdf5 file
  f_read.close()
  
  return parameter_names_emcee, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps

def get_emcee_parameters(chains, temp_status, npost_input, completed_steps):
  # determine nwalkers, nruns, npost (nburn), and nfit parameters
  nwalkers = chains.shape[0]
  nruns = chains.shape[1]
  if (temp_status):
    nruns = int(completed_steps)
  # select posterior chains, without burn in steps
  npost = 0
  if (npost_input < 0):
    #print ' WARNING: npost <= 0. It will be set to: npost = nruns * 10% = %d * 10% = %d' %(nruns, nruns*0.1)
    npost = int(nruns*0.1)
  else:
    npost = int(npost_input)
  nfit = chains.shape[2]
  nruns_sel = nruns - npost
  return nfit, nwalkers, nruns, npost, nruns_sel

def print_memory_usage(array_values):
  print ' MEMORY USAGE: array_values = %d bytes = %.2d MBytes = %.4f GBytes' %(array_values.nbytes, array_values.nbytes/(1024.**2), array_values.nbytes/(1024.**3))
  
def select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries_in, chains):
  # chain is transposed: needed to plot quicker chains for each walker: nruns vs value of parameter
  parameter_boundaries = parameter_boundaries_in.copy()
  chains_T = np.zeros((nruns_sel, nwalkers, nfit))
  for ii in xrange(0,nfit):
    chains_T[:,:,ii] = chains[:,npost:nruns,ii].T # transpose after removing the burnin steps
    if (parameter_names_emcee[ii][0] == 'm' and parameter_names_emcee[ii][1] != 'A'):
      chains_T[:,:,ii] = chains_T[:,:,ii] * m_factor
      parameter_boundaries[ii,:] = parameter_boundaries[ii,:] * m_factor

  return chains_T, parameter_boundaries

def prepare_plot_folder(full_path):
  plot_folder = prepare_emcee_plot_folder(full_path)
  return plot_folder

def prepare_emcee_plot_folder(full_path):
  emcee_plots = os.path.join(full_path, 'plots')
  if (not os.path.isdir(emcee_plots)):
    os.makedirs(emcee_plots)
  return emcee_plots
  
def computation_time(elapsed):
  elapsed_d = elapsed / 86400.
  elapsed_h = (elapsed_d - int(elapsed_d)) * 24.
  elapsed_m = (elapsed_h - int(elapsed_h)) * 60.
  elapsed_s = (elapsed_m - int(elapsed_m)) * 60.
  return int(elapsed_d), elapsed_h, elapsed_m, elapsed_s

def get_pso_data(pso_file):
  of_pso = h5py.File(pso_file, 'r')
  population = np.array(of_pso['population'], dtype=np.float64)
  population_fitness = np.array(of_pso['population_fitness'], dtype=np.float64)
  pso_parameters = np.array(of_pso['pso_parameters'], dtype=np.float64)
  pso_fitness = np.array(of_pso['pso_fitness'], dtype=np.float64)
  if ('pso_best_evolution' in of_pso.keys()):
    pso_best_evolution = np.array(of_pso['pso_best_evolution'], dtype=np.float64)
  else:
    pso_best_evolution = False
  if ('parameters_minmax' in of_pso.keys()):
    parameters_minmax = np.array(of_pso['parameters_minmax'], dtype=np.float64)
  else:
    parameters_minmax = False
  if ('parameter_names' in of_pso.keys()):
    parameter_names = np.array(of_pso['parameter_names'], dtype='S10')
  else:
    parameter_names = False
  of_pso.close()
  pop_shape = population.shape
  return population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution, parameters_minmax, parameter_names, pop_shape
  

def compute_limits(vec_a, delta=0.05):
  a_min = np.min(np.array(vec_a))
  a_max = np.max(np.array(vec_a))
  da = np.abs(a_max - a_min)
  lim_min = a_min - da*delta
  lim_max = a_max + da*delta
  return lim_min, lim_max
  


def get_sigmas(best_parameters, flatchain_posterior):
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  sigma_perc68 = np.percentile(np.abs(flatchain_posterior-best_parameters), 68.27, axis=0, interpolation='midpoint')
  # retrieve confidence intervals of the residual distribution
  sigma_confint = np.percentile(flatchain_posterior-best_parameters, sigmas_percentiles, axis=0, interpolation='midpoint')
  
  return sigma_perc68, sigma_confint

def get_maxlnprob_parameters(npost, nruns, lnprobability, chains_T, flatchain_posterior):
  
  lnprob_burnin = lnprobability[:,npost:nruns]
  maxlnprob_row, maxlnprob_col = get_max_indices(lnprob_burnin)
  maxlnprob = lnprob_burnin[maxlnprob_row, maxlnprob_col]
  # retrieve best parameters <-> best lnprobability
  maxlnprob_parameters = chains_T[maxlnprob_col, maxlnprob_row, :]
  # retrieve 1sigma as 68.27th percentile of the absolute residual distribution
  # retrieve confidence intervals of the residual distribution
  maxlnprob_perc68, maxlnprob_confint = get_sigmas(maxlnprob_parameters, flatchain_posterior)
  
  return maxlnprob, maxlnprob_parameters, maxlnprob_perc68, maxlnprob_confint
  
  
def get_median_parameters(flatchain_posterior):
  median_parameters = np.percentile(flatchain_posterior, 50., axis=0, interpolation='midpoint')
  median_perc68, median_confint = get_sigmas(median_parameters, flatchain_posterior)
  
  return median_parameters, median_perc68, median_confint
  
def get_parameters_median_fitness(nwalkers, npost, nruns, lnprobability, flatchain_posterior, ln_err_const):
  nruns_sel = nruns - npost
  lnprob_burnin = lnprobability[:,npost:nruns]
  flat_lnprob = lnprob_burnin.T.reshape((nwalkers*nruns_sel))
  flat_fitness = -2.*(flat_lnprob-ln_err_const)
  n_med = int(nwalkers*nruns_sel*0.5)
  id_med = np.argsort(flat_fitness)[n_med]
  # retrieve median fitness
  #median_fitness = np.percentile(flat_fitness, 50., interpolation='midpoint')
  median_fitness = flat_fitness[id_med]
  # retrieve parameters at id_med
  medfit_parameters = flatchain_posterior[id_med,:]
  medfit_perc68, medfit_confint = get_sigmas(medfit_parameters, flatchain_posterior)
  
  return median_fitness, medfit_parameters, medfit_perc68, medfit_confint
  

#def compute_max_mean(data_vec, k):
  #hist_counts, bin_edges = np.histogram(data_vec, bins=k)
  #max_bin = np.argmax(np.array(hist_counts))
  #max_mean = np.mean(data_vec[np.logical_and(data_vec>=bin_edges[max_bin], data_vec<=bin_edges[max_bin+1])])
  #return max_mean, max_bin
  
def compute_max_mean(data_vec, k):
  hist_counts, bin_edges = np.histogram(data_vec, bins=k)
  max_bin = np.argmax(np.array(hist_counts))
  if (max_bin == 0 or max_bin == k):
    ext_bin = 0
  elif (max_bin == 1 or max_bin == k-1):
    ext_bin = 1
  else:
    ext_bin = 2

  max_mean = np.mean(data_vec[np.logical_and(data_vec>=bin_edges[max_bin-ext_bin], data_vec<=bin_edges[max_bin+ext_bin])])
  
  return max_mean, max_bin

def get_mode_parameters(flatchain_posterior, k):
  nfit = flatchain_posterior.shape[1]
  mode_parameters = np.zeros((nfit))
  mode_bin = np.zeros((nfit)).astype(int)
  for i_fit in range(0, nfit):
    data_vec = flatchain_posterior[:,i_fit]
    mode_parameters[i_fit], mode_bin[i_fit] = compute_max_mean(data_vec, k)
  mode_perc68, mode_confint = get_sigmas(mode_parameters, flatchain_posterior)
  
  return mode_bin, mode_parameters, mode_perc68, mode_confint

def print_parameters(parameter_names, parameters, perc68, confint, par_type):
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  nsigmas = ['-1','-2','-3','+1','+2','+3']
  
  f_w = 23
  f_d = 16
  
  print
  line_68 = '+/-1sigma(68.27|res|)'.rjust(f_w)
  line_s = '%s' %(' '.join(['%2ssigma(%5.2fres)'.rjust(f_w) %(nsigmas[j],sigmas_percentiles[j]) for j in range(len(sigmas_percentiles))]))
  header = '%s %23s %s %s' %('#'.ljust(10), par_type, line_68, line_s)
  print header
  
  for i_fit in range(0, parameters.shape[0]):
    sigma_line = ''.join(['%23.16f' %(confint[j_perc, i_fit]) for j_perc in range(len(sigmas_percentiles))])
    line = '%10s %23.16f %23.16f %s' %(parameter_names[i_fit], parameters[i_fit], perc68[i_fit], sigma_line)
    print line
  print

  return

def print_parameters_logtxt(of_run, parameter_names, parameters, perc68, confint, par_type):
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  nsigmas = ['-1','-2','-3','+1','+2','+3']
  
  f_w = 23
  f_d = 16
  
  print
  of_run.write('\n')
  line_68 = '+/-1sigma(68.27|res|)'.rjust(f_w)
  line_s = '%s' %(' '.join(['%2ssigma(%5.2fres)'.rjust(f_w) %(nsigmas[j],sigmas_percentiles[j]) for j in range(len(sigmas_percentiles))]))
  header = '%s %23s %s %s' %('#'.ljust(10), par_type, line_68, line_s)
  print header
  of_run.write(header + '\n')
  
  for i_fit in range(0, parameters.shape[0]):
    sigma_line = ''.join(['%23.16f' %(confint[j_perc, i_fit]) for j_perc in range(len(sigmas_percentiles))])
    line = '%10s %23.16f %23.16f %s' %(parameter_names[i_fit], parameters[i_fit], perc68[i_fit], sigma_line)
    print line
    of_run.write(line + '\n')
  print
  of_run.write('\n')

  return


def print_parameters_logger(logger, parameter_names, parameters, perc68, confint, par_type):
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  nsigmas = ['-1','-2','-3','+1','+2','+3']
  
  f_w = 23
  f_d = 16
  
  logger.info('')
  line_68 = '+/-1sigma(68.27|res|)'.rjust(f_w)
  line_s = '%s' %(' '.join(['%2ssigma(%5.2fres)'.rjust(f_w) %(nsigmas[j],sigmas_percentiles[j]) for j in range(len(sigmas_percentiles))]))
  header = '%s %23s %s %s' %('#'.ljust(10), par_type, line_68, line_s)
  logger.info(header)
  
  for i_fit in range(0, parameters.shape[0]):
    sigma_line = ''.join(['%23.16f' %(confint[j_perc, i_fit]) for j_perc in range(len(sigmas_percentiles))])
    line = '%10s %23.16f %23.16f %s' %(parameter_names[i_fit], parameters[i_fit], perc68[i_fit], sigma_line)
    logger.info(line)
  logger.info('')

  return



def get_derived_posterior_parameters(parameter_names, chains_T, flatchain_posterior):
  nfit = flatchain_posterior.shape[1]
  derived_names = []
  derived_chains = []
  derived_posterior = []
  
 
  for i_fit in range(0, nfit):
    if('ecosw' in parameter_names[i_fit]):
      body_id = parameter_names[i_fit].split('ecosw')[1]
      derived_names.append('e%s' %(body_id))
      derived_names.append('w%s' %(body_id))
      
      temp_e, temp_w = 0., 0.
      
      temp_e = np.sqrt(chains_T[:,:,i_fit]**2 + chains_T[:,:,i_fit+1]**2)
      temp_w = np.arctan2(chains_T[:,:,i_fit+1], chains_T[:,:,i_fit])*180./np.pi
      derived_chains.append(temp_e)
      derived_chains.append(temp_w)
      
      temp_e, temp_w = 0., 0.
      
      temp_e = np.sqrt(flatchain_posterior[:,i_fit]**2 + flatchain_posterior[:,i_fit+1]**2)
      temp_w = np.arctan2(flatchain_posterior[:,i_fit+1], flatchain_posterior[:,i_fit])*180./np.pi
      derived_posterior.append(temp_e)
      derived_posterior.append(temp_w)

    if('icoslN' in parameter_names[i_fit]):
      body_id = parameter_names[i_fit].split('icoslN')[1]
      derived_names.append('i%s' %(body_id))
      derived_names.append('lN%s' %(body_id))
      
      temp_i, temp_lN = 0., 0.
      
      temp_i = np.sqrt(chains_T[:,:,i_fit]**2 + chains_T[:,:,i_fit+1]**2)
      temp_lN = np.arctan2(chains_T[:,:,i_fit+1], chains_T[:,:,i_fit])*180./np.pi
      derived_chains.append(temp_i)
      derived_chains.append(temp_lN)
      
      temp_i, temp_lN = 0., 0.
      
      temp_i = np.sqrt(flatchain_posterior[:,i_fit]**2 + flatchain_posterior[:,i_fit+1]**2)
      temp_lN = np.arctan2(flatchain_posterior[:,i_fit+1], flatchain_posterior[:,i_fit])*180./np.pi
      derived_posterior.append(temp_i)
      derived_posterior.append(temp_lN)
      
  nder = len(derived_names)
  derived_chains_T = np.zeros((chains_T.shape[0], chains_T.shape[1], nder))
  for i_der in range(0, nder):
    derived_chains_T[:,:,i_der] = np.array(derived_chains[i_der])
  
  return derived_names, derived_chains_T, np.array(derived_posterior).T
  
  
def get_proper_mass(m_type, parameter_names_emcee, full_path):

  m_factor = mass_conversion_factor(m_type)
  
  nfit = np.size(parameter_names_emcee)
  m_fit_logical = [False]*nfit
  
  for ifit in range(0, nfit):
    if('Ms' in parameter_names_emcee[ifit]): m_fit_logical[ifit] = True

  MR_star = np.ones((2,2))

  if(any(m_fit_logical)):
    bodies = os.path.join(full_path, 'bodies.lst')
    obf = open(bodies, 'r')
    line_star = obf.readline()
    obf.close()
    file_star = line_star.split()[0].strip()
    MR_star = np.genfromtxt(file_star)
    if(len(MR_star.shape) == 1): MR_star = np.column_stack((MR_star, np.zeros((2))))
  
  m_factor = m_factor * MR_star[0,0]
  
  return m_factor
  
  
