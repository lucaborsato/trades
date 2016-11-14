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
#from emcee import autocorr as acor
import acor

# common variables needed to create labels and parameter names
kel_id = ['$M$', '$R$', '$P$', '$e$', '$\omega$', '$\\nu$', '$i$', '$\Omega$']
kel_units = ['$[M_{\oplus}]$', '$[R_\mathrm{Jup}]$', '[days]', '', '', '[deg]', '[deg]', '[deg]']

kel_id_2 = ['$M/M_\star$', '$R$', '$P$', '$e\cos\omega$', '$e\sin\omega$', '$\lambda$', '$i\cos\Omega$', '$i\sin\Omega$']
kel_units_2 = ['', '$[R_\mathrm{Jup}]$', '[days]', '', '', '[deg]', '', '']

kel_fmt = ['%.3f', '%.4f', '%.3f', '%.1f', '%.1f', '%.1f', '%.3f', '%.3f']
nelem = len(kel_id)
letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'l', 'm', 'n', 'o', 'p']

deg2rad = np.pi/180.0
rad2deg = 180.0/np.pi

eps64bit = np.finfo(np.float64(1.0)).eps
eps32bit = np.finfo(np.float32(1.0)).eps

# sigma =          1      0    -1       1     -2       2     -3       3
# sigma pos        0      1     2       3      4       5      6       7
percentile_val = [68.27, 50.0, 15.865, 84.135, 2.265, 97.735, 0.135, 99.865]

def set_bool_argument(arg_in):
  if (str(arg_in).lower() in ['t', 'tr', 'tru', 'true', 'y', 'ye', 'yes', '1']):
    arg_out = True
  else:
    arg_out = False
  return arg_out

def set_int_argument(arg_in):
  try:
    arg_out = int(arg_in)
  except:
    arg_out = 0
  return arg_out
  
# read command line (cli) arguments
def get_args():
  parser = argparse.ArgumentParser(description='TRADES+EMCEE PLOT')
  
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')
  
  parser.add_argument('-nb', '--nburn', '-np', '--npost', action='store', dest='npost', required=True, help='The number of posterior/burn in steps to discard at the beginning of each chain. It has to be > 0')
  
  parser.add_argument('-m', '--mtype', '--mass-type', action='store', dest='m_type', default='e', help='Mass type: j = Jupiter, e = Earth, s = Sun. Default is Earth = e.')

  parser.add_argument('-t', '--temp-file', action='store', dest='temp_status', default=False, help='If you want to read temporary emcee_temp.hdf5 file, because simulation is not finished yet. Default is False')
  parser.add_argument('-b', '--boot-file', '--save-like-boot', action='store', dest='boot_id', default=0, help='If you want to save flatchain posterior as the posterior file, in order to be read by read_finalpar_v2.py. If set to 0 (default) it doesn\'t save anything. Bootstrap file output: boot_id_posterior_sim.dat')
  
  parser.add_argument('-c', '--cumulative', '--cumulative-histogram', action='store', dest='cumulative', default='False', help='If you want to plot the cumulative histogram instead of normal histogram. Set it to "True/Yes/1" Default is False. ')
  
  parser.add_argument('-s', '--steps', '--gelmanrubin-steps', action='store', dest='sel_steps', default=0, help='Set to positive integer. It will allow to compute GelmanRubin/Geweke statistics of only sel_steps, and not for each step. Just to debug and overall view of the chains. Default is 0.')
  
  cli = parser.parse_args()

  cli.full_path = os.path.abspath(cli.full_path)
  cli.m_type = cli.m_type.lower()
  cli.temp_status = set_bool_argument(cli.temp_status)
  cli.cumulative = set_bool_argument(cli.cumulative)
  cli.boot_id = set_int_argument(cli.boot_id)
  
  return cli

# given the mass flag letter it computes the proper mass conversion factor
def mass_type_factor(Ms=1, mtype='earth', mscale=True):
  if (mtype.lower() in ['s', 'su', 'sun']):
    conv_factor = np.float64(1.0)
    scale_factor = Ms
    mass_unit = 'M_Sun'
  elif (mtype.lower() in ['e', 'ea', 'ear', 'eart', 'earth']):
    conv_factor = cst.Msear
    scale_factor = Ms * cst.Msear
    mass_unit = 'M_Earth'
  elif (mtype.lower() in ['n', 'ne', 'nep', 'nept', 'neptu', 'neptun', 'neptune']):
    conv_factor = cst.Msnep
    scale_factor = Ms * cst.Msnep
    mass_unit = 'M_Nep'
  else:
    conv_factor = cst.Msjup
    scale_factor = Ms * cst.Msjup
    mass_unit = 'M_Jup'
  if(mscale):
    return scale_factor, mass_unit
  else:
    return conv_factor, mass_unit
  
def get_proper_mass(m_type, parameter_names_emcee, full_path):

  nfit = np.size(parameter_names_emcee)

  MR_star = np.ones((2,2))

  bodies = os.path.join(full_path, 'bodies.lst')
  obf = open(bodies, 'r')
  line_star = obf.readline()
  obf.close()
  file_star = line_star.split()[0].strip()
  MR_star = np.genfromtxt(file_star)
  if(len(MR_star.shape) == 1): MR_star = np.column_stack((MR_star, np.zeros((2))))
  
  m_factor, m_unit = mass_type_factor(MR_star[0,0], cli.m_type, True)
  
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
  cosa = np.cos(angle*deg2rad)
  sina = np.sin(angle*deg2rad)
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

def compute_autocor_time(chains):
  nr, nw, nfit = chains.shape

  autocor_time = np.zeros((nfit), dtype=np.float64)
  for ifit in range(0, nfit):
    x = chains[:,:,ifit]
    #tau, mean, sigma = acor.acor(x)
    temp_acor = np.mean(np.array([acor.acor(x[:,iw]) for iw in range(0, nw)]), axis=0)
    autocor_time[ifit] = temp_acor[0]
    
  return autocor_time

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
  
def select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries_in, chains, m_star=1.0):
  # chain is transposed: needed to plot quicker chains for each walker: nruns vs value of parameter
  parameter_boundaries = parameter_boundaries_in.copy()
  chains_T = np.zeros((nruns_sel, nwalkers, nfit))
  for ii in xrange(0,nfit):
    chains_T[:,:,ii] = chains[:,npost:nruns,ii].T # transpose after removing the burnin steps
    if (parameter_names_emcee[ii][0] == 'm' and parameter_names_emcee[ii][1] != 'A'):
      if('Ms' in parameter_names_emcee[ii]):
        m_conv = m_factor * m_star
      else:
        m_conv = np.float64(1.)
      chains_T[:,:,ii] = chains_T[:,:,ii] * m_conv
      parameter_boundaries[ii,:] = parameter_boundaries[ii,:] * m_conv

  return chains_T, parameter_boundaries

def posterior_back_to_msun(m_factor, parameter_names_emcee, flatchain_posterior_in):
  nfit = flatchain_posterior_in.shape[1]
  flatchain_posterior_out = flatchain_posterior_in.copy()
  for ifit in range(0,nfit):
    if(parameter_names_emcee[ifit][0] == 'm' and parameter_names_emcee[ifit][1] != 'A'):
      flatchain_posterior_out[:,ifit] = flatchain_posterior_out[:,ifit]/m_factor
  return flatchain_posterior_out

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

def get_mode_parameters_full(flatchain_posterior, k):
  nfit = flatchain_posterior.shape[1]
  mode_parameters = np.zeros((nfit))
  mode_bin = np.zeros((nfit)).astype(int)
  for i_fit in range(0, nfit):
    data_vec = flatchain_posterior[:,i_fit]
    mode_parameters[i_fit], mode_bin[i_fit] = compute_max_mean(data_vec, k)
  mode_perc68, mode_confint = get_sigmas(mode_parameters, flatchain_posterior)
  
  return mode_bin, mode_parameters, mode_perc68, mode_confint
  
def get_mode_parameters(flatchain_posterior, k):
  nfit = flatchain_posterior.shape[1]
  mode_parameters = np.zeros((nfit))
  mode_bin = np.zeros((nfit)).astype(int)
  for i_fit in range(0, nfit):
    data_vec = flatchain_posterior[:,i_fit]
    mode_parameters[i_fit], mode_bin[i_fit] = compute_max_mean(data_vec, k)
  
  return mode_bin, mode_parameters

def print_parameters_nolog(parameter_names, parameters, perc68, confint, par_type):
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
    
def save_posterior_like(out_folder, boot_id, parameter_names, flatchain_posterior):
  nfit = flatchain_posterior.shape[1]
  nboot = flatchain_posterior.shape[0]
  header_0000 = ' iboot %s' %(' '.join(parameter_names))
  fmt_dp = ' '.join( [ '%s' %('%23.16e') for ii in range(0,nfit) ] )
  fmt_full = '%s %s' %('%6d', fmt_dp)
  iboot_fake = np.arange(1,nboot+1,1)
  boot_fake = np.column_stack((iboot_fake, flatchain_posterior))
  boot_file = os.path.join(out_folder, '%d_posterior_sim.dat' %(boot_id))
  np.savetxt(boot_file, boot_fake, fmt=fmt_full, header=header_0000)
  return boot_file

def GelmanRubin_test_1(chains_T):
  n, M = np.shape(chains_T)
  theta_m = np.mean(chains_T, axis=0)
  
  B = np.sum((theta_m - np.mean(theta_m))**2)*n/(M-1)
  
  W = np.mean(np.var(chains_T, axis=0, ddof=1))
  
  Var = (W*(n-1) + B)/n
  Rc = np.sqrt(Var / W)
  return Rc
  
def GelmanRubin_test_2(chains_T):

  n, M = np.shape(chains_T)
  
  theta_m = [np.mean(chains_T[:,i_m]) for i_m in range(0, M)]
  theta = np.mean(theta_m)
  
  d_theta2 = (theta_m - theta)**2
  B_n = np.sum(d_theta2) / (M-1)
  
  arg_W = [np.sum((chains_T[:,i_m] - theta_m[i_m])**2) / (n-1) for i_m in range(0, M)]
  W = np.mean(arg_W)
    
  n_frac = (n-1)/n
  var_plus = n_frac*W + B_n
  Var = var_plus + (B_n/M)

  Rc = np.sqrt(Var / W)
  return Rc


def GelmanRubin_PyORBIT(chains_input):
  n_iters, n_chain = np.shape(chains_input)
  W = np.asarray(0., dtype=np.double)
  z_pc = np.sum(chains_input, axis=0) / n_iters  # eq 20
  for nn in xrange(0, n_chain):
    W += (np.sum(np.power(chains_input[:, nn] - z_pc[nn], 2))) / ((n_iters - 1) * n_chain)  # eq 21
  z_pp = np.sum(chains_input) / (n_chain * n_iters)
  B = np.sum(np.power(z_pc - z_pp, 2)) * (n_chain / (n_iters - 1))
  var = W * (n_chain - 1) / n_chain + B / n_chain
  
  return np.sqrt(var / W)

def GelmanRubin_pymc(x):
  """ Returns estimate of R for a set of traces.

  The Gelman-Rubin diagnostic tests for lack of convergence by comparing
  the variance between multiple chains to the variance within each chain.
  If convergence has been achieved, the between-chain and within-chain
  variances should be identical. To be most effective in detecting evidence
  for nonconvergence, each chain should have been initialized to starting
  values that are dispersed relative to the target distribution.

  Parameters
  ----------
  x : array-like
    An array containing the 2 or more traces of a stochastic parameter. That is, an array of dimension m x n x k, where m is the number of traces, n the number of samples, and k the dimension of the stochastic.

  Returns
  -------
  Rhat : float
    Return the potential scale reduction factor, :math:`\hat{R}`

  Notes
  -----

  The diagnostic is computed by:

    .. math:: \hat{R} = \frac{\hat{V}}{W}

  where :math:`W` is the within-chain variance and :math:`\hat{V}` is
  the posterior variance estimate for the pooled traces. This is the
  potential scale reduction factor, which converges to unity when each
  of the traces is a sample from the target posterior. Values greater
  than one indicate that one or more chains have not yet converged.

  References
  ----------
  Brooks and Gelman (1998)
  Gelman and Rubin (1992)"""

  if np.shape(x) < (2,):
    raise ValueError(
      'Gelman-Rubin diagnostic requires multiple chains of the same length.')

  try:
    m, n = np.shape(x)
  except ValueError:
    return [gelman_rubin(np.transpose(y)) for y in np.transpose(x)]

  # Calculate between-chain variance
  B_over_n = np.sum((np.mean(x, 1) - np.mean(x)) ** 2) / (m - 1)

  # Calculate within-chain variances
  W = np.sum([(x[i] - xbar) ** 2 for i, xbar in enumerate(np.mean(x,1))]) / (m * (n - 1))

  # (over) estimate of variance
  s2 = W * (n - 1) / n + B_over_n

  # Pooled posterior variance estimate
  V = s2 + B_over_n / m

  # Calculate PSRF
  R = V / W

  return R


# Geweke 1992 test 1
def geweke_test(chains_T, start_frac=0.05, n_sel_steps=5):
  n_steps, n_chains = chains_T.shape
  half = int(0.5*n_steps)
  
  int_b = chains_T[half:,:]
  mean_b = np.mean(int_b, axis=0)
  var_b = np.var(int_b, axis=0, ddof=1)
  
  start_step = int(start_frac * n_steps)
  sel_a = np.linspace(start=start_step, stop=half, num=n_sel_steps, endpoint=False, dtype=np.int)
  
  zscore = np.zeros((n_sel_steps, n_chains))
  for i_a in range(0, n_sel_steps):
    int_a = chains_T[sel_a[i_a]:half,:]
    mean_a = np.mean(int_a, axis=0)
    var_a = np.var(int_a, axis=0, ddof=1)
    z_temp = (mean_a - mean_b) / np.sqrt(var_a + var_b)
    zscore[i_a,:] = z_temp
  
  return sel_a, zscore

# ------------------------------
# function from read_finalpar_v2

def print_both(line, output=None):
  print line
  if(output is not None):
    output.write(line + '\n')
  return

def get_units(names, mass_unit):
  n_names = len(names)
  units_par = []
  
  for i in range(0, n_names):
    
    if(str(names[i])[0] == 'm'):
      if('Ms' in names[i]):
        units_par.append('[M_pl/M_star]')
      elif('mA' in names[i]):
        units_par.append('[deg]')
      else:
        units_par.append('[%s]' %(mass_unit))
    
    if('R' in str(names[i])):
      units_par.append('[R_sun]')
    
    if('P' in str(names[i])):
      units_par.append('[days]')
      
    if(str(names[i])[0] == 'e'): # the same for e, ecosw, esinw
      units_par.append(' ')
    
    if(str(names[i])[0] == 'w'):
      units_par.append('[deg]')
      
    if('lambda' in str(names[i])):
      units_par.append('[deg]')
      
    if(str(names[i])[0] == 'i'):
      if('cos' in names[i] or 'sin' in names[i]):
        units_par.append(' ')
      else:
        units_par.append('[deg]')
    
    if(str(names[i])[0:2] == 'lN'):
      units_par.append('[deg]')
    
  return units_par

def elements(fpath, idsim, lmf=0):
  #kel_file = os.path.join(fpath, str(idsim) + "_" + str(lmf) + "_initialElements.dat")
  kel_file = os.path.join(fpath, '%d_%d_initialElements.dat' %(idsim, lmf))
  try:
    kep_elem = np.genfromtxt(kel_file) # (NB-1) x (M_Msun R_Rsun P_d a_AU ecc w_deg mA_deg inc_deg lN_deg)
  except:
    print ' KEPLERIAN ELEMENTS FILE NOT FOUND %s' %(kel_file)
    sys.exit()
  return kel_file, kep_elem

def get_case(id_body, fit_body):
  fit_n = np.arange(1,11,1)
  id_fit = [fit_n[j] for j in range(len(fit_body)) if (fit_body[j]=='1')]
  id_all = [8*(id_body-1) + 2 + id_fit[j] for j in range(len(id_fit))]
  
  if (6 not in id_fit): # not fitting lambda
    case = [0]
  
  else: #(6 in id_fit) # fitting lambda
    if ( all( x in id_fit for x in [4, 5, 6, 7, 8]) ):  # fit ecosw, esinw, lambda, icoslN, isinlN
      case = [4]
    elif ( all( x in id_fit for x in [4, 5, 6]) ): # fit ecosw, esinw, lambda
      case = [2]
    elif ( all( x in id_fit for x in [6, 7, 8]) ): # fit lambda, icoslN, isinlN
      case = [3]
    else:
      case = [1] # fit lambda & w || lambda & lN || lamda & w & lN
  
  return id_fit, id_all, case

def get_fitted(full_path):
  of = open(os.path.join(full_path, 'bodies.lst'), 'r')
  lines = of.readlines()
  of.close()
  NB = len(lines)
  n_pl = NB - 1
  bodies_file = []
  fit_string = ''
  fit_list = []
  for line in lines:
    fname = line.strip().split()[0]
    bodies_file.append(fname)
    temp_string = line.strip().split(fname)[1].split('#')[0].strip()
    if('F' in temp_string):
      temp_string = temp_string.split('F')[0].strip()
    elif('f' in temp_string):
      temp_string = temp_string.split('f')[0].strip()
    elif('T' in temp_string):
      temp_string = temp_string.split('T')[0].strip()
    elif('t' in temp_string):
      temp_string = temp_string.split('t')[0].strip()

    fit_string = '%s %s' %(fit_string, temp_string)
    fit_list.append(temp_string.split())

  nfit = np.sum(np.array(fit_string.split(), dtype=np.int))
  
  case = []
  id_fit = []
  id_all = []
  nfit_list = []
  cols_list = []
  nfit_cnt = 0
  for j in range(1,NB+1):
    id_fit_j, id_all_j, case_j= get_case(j, fit_list[j-1])
    id_fit.append(id_fit_j)
    id_all.append(id_all_j)
    case.append(case_j)
    nfit_j = len(id_fit_j)
    nfit_list.append(nfit_j)
    cols_list.append([cc for cc in range(nfit_cnt, nfit_cnt+nfit_j)])
    nfit_cnt += nfit_j
  cols = [jj for jj in range(0,nfit)]

  
  return nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case


def compute_intervals(flatchain, parameters, percentiles):
  sigma_par = np.percentile(np.subtract(flatchain, parameters), percentiles, axis=0, interpolation='midpoint') # (n_percentile x nfit)
  sigma_par[0] = np.percentile(np.abs(np.subtract(flatchain, parameters)), percentiles[0], axis=0, interpolation='midpoint') # 68.27th
  sigma_par[1] = np.percentile(np.abs(np.subtract(flatchain, parameters)), percentiles[1], axis=0, interpolation='midpoint') # MAD
  #sigma_par = np.percentile(flatchain, parameters, percentiles, axis=0) # (n_percentile x nfit)
  #sigma_par[0] = np.percentile(np.abs(flatchain, parameters), percentiles[0], axis=0) # 68.27th
  #sigma_par[1] = np.percentile(np.abs(flatchain, parameters), percentiles[1], axis=0) # MAD
  
  return sigma_par

def get_good_distribution(posterior_scale, posterior_mod):
  par_scale = np.median(posterior_scale)
  p68_scale = np.percentile(np.abs(posterior_scale-par_scale), 68.27, interpolation='midpoint')
  
  par_mod = np.median(posterior_mod)
  p68_mod = np.percentile(np.abs(posterior_mod-par_mod), 68.27, interpolation='midpoint')
  
  if(np.abs(p68_scale - p68_mod) <= eps64bit):
    return par_mod, posterior_mod
  else:
    if(p68_mod < p68_scale):
      return par_mod, posterior_mod
    else:
      return par_scale, posterior_scale
  return par_mod, posterior_mod

def get_proper_posterior_correlated(posterior, col=None):
  
  if(col is not None):
    posterior_scale = np.arctan2(posterior[:,col+1], posterior[:,col]) * rad2deg
    posterior_mod = (posterior_scale+360.)%360.
  else:
    posterior_scale = np.arctan2(np.sin(posterior*deg2rad), np.cos(posterior*deg2rad)) * rad2deg
    posterior_mod = (posterior_scale+360.)%360.
  
  par_out, boot_out = get_good_distribution(posterior_scale, posterior_mod)
  return par_out, boot_out

def get_trigonometric_posterior(id_NB, id_l, col, idpar, posterior): # ecosw,esinw -> (e,w) && icoslN,isinlN -> (i,lN) || No lambda
  
  cosboot = posterior[:,col]
  sinboot = posterior[:,col+1]
  aboot = np.sqrt( (cosboot*cosboot) + (sinboot*sinboot) )
  temp, bboot = get_proper_posterior_correlated(posterior, col)
  if('e' in idpar[id_l]):
    sel_ezero = aboot <= eps64bit
    bboot[sel_ezero] = 90.
        
  if ('ecosw' in idpar[id_l]):
    id_aa = 'e'
    id_bb = 'w'
  else:
    id_aa = 'i'
    id_bb = 'lN'
        
  id_a0 = '%s%s' %(id_aa,id_NB+1)
  id_b0 = '%s%s' %(id_bb,id_NB+1)
    
  return id_a0, aboot, id_b0, bboot

def get_trigonometric_parameters(id_NB, id_l, col, idpar, parameters): # ecosw,esinw -> (e,w) && icoslN,isinlN -> (i,lN) || No lambda
  
  cosp = parameters[col]
  sinp = parameters[col+1]
  a_par = np.sqrt( (cosp*cosp) + (sinp*sinp) )
  b_par = np.arctan2(sinp, cosp)*rad2deg

  if('e' in idpar[id_l]):
    if(a_par <= eps64bit):
      b_par = 90.
        
  if ('ecosw' in idpar[id_l]):
    id_aa = 'e'
    id_bb = 'w'
  else:
    id_aa = 'i'
    id_bb = 'lN'
        
  id_a0 = '%s%s' %(id_aa,id_NB+1)
  id_b0 = '%s%s' %(id_bb,id_NB+1)
    
  return id_a0, a_par, id_b0, b_par


def derived_posterior_case0(idpar_NB, id_fit_NB, i_NB, cols, posterior): # not fitting lambda
  name_der = []
  der_posterior = []
  
  nlist = len(id_fit_NB)
  for i_l in range(nlist):
    if ('ecosw' in idpar_NB[i_l] or 'icoslN' in idpar_NB[i_l]):
      id_a0, aboot, id_b0, bboot = get_trigonometric_posterior(i_NB, i_l, cols[i_l], idpar_NB, posterior)
      name_der.append(id_a0)
      name_der.append(id_b0)
      der_posterior.append(aboot)
      der_posterior.append(bboot)
      
  return name_der, der_posterior

def derived_parameters_case0(idpar_NB, id_fit_NB, i_NB, cols, parameters): # not fitting lambda
  name_der = []
  der_par = []
  
  nlist = len(id_fit_NB)
  for i_l in range(nlist):
    if ('ecosw' in idpar_NB[i_l] or 'icoslN' in idpar_NB[i_l]):
      id_a0, a_par, id_b0, b_par = get_trigonometric_parameters(i_NB, i_l, cols[i_l], idpar_NB, parameters)
      name_der.append(id_a0)
      name_der.append(id_b0)
      der_par.append(a_par)
      der_par.append(b_par)
      
  return name_der, der_par


def derived_posterior_case1(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, posterior): #  fitting lambda || lambda & w || lambda & lN || lamda & w & lN
  name_der = []
  der_posterior = []
  
  nlist = len(id_fit_NB)
  for i_l in range(0,nlist):
    
    if(id_fit_NB[i_l] == 6):
      if(idpar_NB[i_l+2][0:2] != 'mA'):
        
        aboot = (posterior[:,cols[i_l]] - kep_elem[5] - kep_elem[8])%360.
        temp, aboot = get_proper_posterior_correlated(aboot)
          
        name_der.append('%s%s' %('mA',i_NB+1))
        der_posterior.append(aboot)
        
  return name_der, der_posterior

def derived_parameters_case1(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, parameters): #  fitting lambda || lambda & w || lambda & lN || lamda & w & lN
  name_der = []
  der_par = []
  
  nlist = len(id_fit_NB)
  for i_l in range(0,nlist):
    
    if(id_fit_NB[i_l] == 6):
      if(idpar_NB[i_l+2][0:2] != 'mA'):
        
        a_par = (parameters[cols[i_l]] - kep_elem[5] - kep_elem[8])%360.
        name_der.append('%s%s' %('mA',i_NB+1))
        der_par.append(a_par)
        
  return name_der, der_par


def derived_posterior_case2(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, posterior): # fit ecosw, esinw, lambda
  name_der = []
  der_posterior = []
  
  nlist = len(id_fit_NB)
  for i_l in range(0,nlist):
    
    if(id_fit_NB[i_l] == 4):
      id_e0, eboot, id_w0, wboot = get_trigonometric_posterior(i_NB, i_l, cols[i_l], idpar_NB, posterior)
      name_der.append(id_e0)
      name_der.append(id_w0)
      der_posterior.append(eboot)
      der_posterior.append(wboot)
      
      if(idpar_NB[i_l+2][0:2] != 'mA'):
        mAboot = (posterior[:,cols[i_l+2]] - wboot - kep_elem[8])%360. # mA = lambda - w - lN
        temp, mAboot = get_proper_posterior_correlated(mAboot)
        name_der.append('%s%s' %('mA',i_NB+1))
        der_posterior.append(mAboot)
        
  return name_der, der_posterior

def derived_parameters_case2(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, parameters): # fit ecosw, esinw, lambda
  name_der = []
  der_par = []
  
  nlist = len(id_fit_NB)
  for i_l in range(0,nlist):
    
    if(id_fit_NB[i_l] == 4):
      id_e0, e_par, id_w0, w_par = get_trigonometric_parameters(i_NB, i_l, cols[i_l], idpar_NB, parameters)
      name_der.append(id_e0)
      name_der.append(id_w0)
      der_par.append(e_par)
      der_par.append(w_par)
      
      if(idpar_NB[i_l+2][0:2] != 'mA'):
        
        mA_par = (parameters[cols[i_l+2]] - w_par - kep_elem[8])%360. # mA = lambda - w - lN
        name_der.append('%s%s' %('mA',i_NB+1))
        der_par.append(mA_par)
        
  return name_der, der_par


def derived_posterior_case3(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, posterior): # fit lambda, icoslN, isinlN
  name_der = []
  der_posterior = []
  
  nlist = len(id_fit_NB)
  for i_l in range(0,nlist):
    
    if(id_fit_NB[i_l] == 7):
      id_i0, iboot, id_lN0, lNboot = get_trigonometric_posterior(i_NB, i_l, cols[i_l], idpar_NB, posterior)
      name_der.append(id_i0)
      name_der.append(id_lN0)
      der_posterior.append(iboot)
      der_posterior.append(lNboot)
      
      if(idpar_NB[i_l+2][0:2] != 'mA'):
        mAboot = (posterior[:,cols[i_l-1]] - kep_elem[5] - lNboot)%360. # mA = lambda - w - lN
        temp, mAboot = get_proper_posterior_correlated(mAboot)
        name_der.append('%s%s' %('mA',i_NB+1))
        der_posterior.append(mAboot)
        
  return name_der, der_posterior

def derived_parameters_case3(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, parameters): # fit lambda, icoslN, isinlN
  name_der = []
  der_par = []
  
  nlist = len(id_fit_NB)
  for i_l in range(0,nlist):
    
    if(id_fit_NB[i_l] == 7):
      id_i0, i_par, id_lN0, lN_par = get_trigonometric_parameters(i_NB, i_l, cols[i_l], idpar_NB, parameters)
      name_der.append(id_i0)
      name_der.append(id_lN0)
      der_par.append(i_par)
      der_par.append(lN_par)
      
      if(idpar_NB[i_l+2][0:2] != 'mA'):
        mA_par = (parameters[cols[i_l-1]] - kep_elem[5] - lN_par)%360. # mA = lambda - w - lN
        name_der.append('%s%s' %('mA',i_NB+1))
        der_par.append(mA_par)
        
  return name_der, der_par


def derived_posterior_case4(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, posterior): # fit ecosw, esinw, lambda, icoslN, isinlN
  name_der = []
  der_posterior = []
  
  nlist = len(id_fit_NB)
  for i_l in range(0,nlist):
    
    if(id_fit_NB[i_l] == 4):
      id_e0, eboot, id_w0, wboot = get_trigonometric_posterior(i_NB, i_l, cols[i_l], idpar_NB, posterior)
      name_der.append(id_e0)
      name_der.append(id_w0)
      der_posterior.append(eboot)
      der_posterior.append(wboot)
      
      id_i0, iboot, id_lN0, lNboot = get_trigonometric_posterior(i_NB, i_l+3, cols[i_l+3], idpar_NB, posterior)
      name_der.append(id_i0)
      name_der.append(id_lN0)
      der_posterior.append(iboot)
      der_posterior.append(lNboot)
      
      if(idpar_NB[i_l+2][0:2] != 'mA'):
        mAboot = (posterior[:,cols[i_l+2]] - wboot - lNboot)%360.# mA = lambda - w - lN
        temp, mAboot = get_proper_posterior_correlated(mAboot)
        name_der.append('%s%s' %('mA',i_NB+1))
        der_posterior.append(mAboot)
        
  return name_der, der_posterior

def derived_parameters_case4(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, parameters): # fit ecosw, esinw, lambda, icoslN, isinlN
  name_der = []
  der_par = []
  
  nlist = len(id_fit_NB)
  for i_l in range(0,nlist):
    
    if(id_fit_NB[i_l] == 4):
      id_e0, e_par, id_w0, w_par = get_trigonometric_parameters(i_NB, i_l, cols[i_l], idpar_NB, parameters)
      name_der.append(id_e0)
      name_der.append(id_w0)
      der_par.append(e_par)
      der_par.append(w_par)
      
      id_i0, i_par, id_lN0, lN_par = get_trigonometric_parameters(i_NB, i_l+3, cols[i_l+3], idpar_NB, parameters)
      name_der.append(id_i0)
      name_der.append(id_lN0)
      der_par.append(i_par)
      der_par.append(lN_par)
      
      if(idpar_NB[i_l+2][0:2] != 'mA'):
        mA_par = (parameters[cols[i_l+2]] - w_par - lN_par)%360.# mA = lambda - w - lN
        name_der.append('%s%s' %('mA',i_NB+1))
        der_par.append(mA_par)
        
  return name_der, der_par

def compute_derived_posterior(idpar, kep_elem, id_fit, case_list, cols_list, posterior, conv_factor=1.):
  NB = len(case_list)
  nfit = len(id_fit)
  #print NB, nfit
  
  posterior_out = posterior.copy()
    
  id_derived = []
  der_posterior = []
  
  cc_fit = 0
  
  for i_NB in range(0, NB):
    nfit_NB = len(id_fit[i_NB])
    
    if(nfit_NB > 0):
      
      idpar_NB = idpar[cols_list[i_NB][0]:cols_list[i_NB][-1]+1] # names of the parameter for the body i_NB
      id_fit_NB = id_fit[i_NB] # integers that identify the proper type of the fitted parameters
      
      for i_fit in range(0, nfit_NB):
        
        if(id_fit[i_NB][i_fit] == 1): # convert Mp and Mp/Ms into mass unit specified by the user <-> mtype
          
          if('Ms' in idpar[cc_fit]):
            mboot = posterior[:,cc_fit]*conv_factor
            xid = 'm%d' %(i_NB+1)
            id_derived.append(xid)
            der_posterior.append(mboot)
        
        cc_fit += 1
      
      id_temp, der_temp = [], []
      
      if(case_list[i_NB][0] == 0):
        id_temp, der_temp = derived_posterior_case0(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], posterior)
      
      elif(case_list[i_NB][0] == 1):
        id_temp, der_temp = derived_posterior_case1(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], kep_elem[i_NB-1,:], posterior)
        
      elif(case_list[i_NB][0] == 2):
        id_temp, der_temp = derived_posterior_case2(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], kep_elem[i_NB-1,:], posterior)
        
      elif(case_list[i_NB][0] == 3):
        id_temp, der_temp = derived_posterior_case3(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], kep_elem[i_NB-1,:], posterior)
        
      elif(case_list[i_NB][0] == 4):
        id_temp, der_temp = derived_posterior_case4(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], kep_elem[i_NB-1,:], posterior)
        
      id_derived.append(id_temp)
      der_posterior.append(der_temp)
  
  n_der = 0
  n_xder = len(der_posterior)
  if(n_xder > 0):
    n_der_single  = []
    derived_post   = []
    names_derived = []
    
    for ii in range(0, n_xder):
      temp_names = np.array(id_derived[ii],dtype=str)
      ntemp = np.size(temp_names)
      n_der_single.append(ntemp)
      n_der += ntemp
      temp_der = np.array(der_posterior[ii],dtype=np.float64)
      
      if(ntemp > 0):
        
        if(ntemp == 1):
          derived_post.append(temp_der)
          names_derived.append(temp_names)
        else:
          for jj in range(0, ntemp):
            derived_post.append(temp_der[jj])
            names_derived.append(temp_names[jj])
  
  return np.array(names_derived, dtype=str), np.array(derived_post, dtype=np.float64)

def compute_derived_parameters(idpar, kep_elem, id_fit, case_list, cols_list, parameters, conv_factor=1.):
  NB = len(case_list)
  nfit = len(id_fit)
  #print NB, nfit
  
  id_derived = []
  der_par = []
  
  cc_fit = 0
  
  for i_NB in range(0, NB):
    nfit_NB = len(id_fit[i_NB])
    
    if(nfit_NB > 0):
      
      idpar_NB = idpar[cols_list[i_NB][0]:cols_list[i_NB][-1]+1] # names of the parameter for the body i_NB
      id_fit_NB = id_fit[i_NB] # integers that identify the proper type of the fitted parameters
      
      for i_fit in range(0, nfit_NB):
        
        if(id_fit[i_NB][i_fit] == 1): # convert Mp and Mp/Ms into mass unit specified by the user <-> mtype
          
          if('Ms' in idpar[cc_fit]):
            m_par = parameters[cc_fit]*conv_factor
            xid = 'm%d' %(i_NB+1)
            id_derived.append(xid)
            der_par.append(m_par)
        
        cc_fit += 1
      
      id_temp, der_temp = [], []
      
      if(case_list[i_NB][0] == 0):
        id_temp, der_temp = derived_parameters_case0(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], parameters)
      
      elif(case_list[i_NB][0] == 1):
        id_temp, der_temp = derived_parameters_case1(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], kep_elem[i_NB-1,:], parameters)
        
      elif(case_list[i_NB][0] == 2):
        id_temp, der_temp = derived_parameters_case2(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], kep_elem[i_NB-1,:], parameters)
        
      elif(case_list[i_NB][0] == 3):
        id_temp, der_temp = derived_parameters_case3(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], kep_elem[i_NB-1,:], parameters)
        
      elif(case_list[i_NB][0] == 4):
        id_temp, der_temp = derived_parameters_case4(idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], kep_elem[i_NB-1,:], parameters)
        
      id_derived.append(id_temp)
      der_par.append(der_temp)
  
  n_der = 0
  n_xder = len(der_par)
  if(n_xder > 0):
    n_der_single  = []
    derived_par   = []
    names_derived = []
    
    for ii in range(0, n_xder):
      temp_names = np.array(id_derived[ii],dtype=str)
      ntemp = np.size(temp_names)
      n_der_single.append(ntemp)
      n_der += ntemp
      temp_der = np.array(der_par[ii],dtype=np.float64)
            
      if(ntemp > 0):
        
        if(ntemp == 1):
          derived_par.append(temp_der)
          names_derived.append(temp_names)
        else:
          for jj in range(0, ntemp):
            derived_par.append(temp_der[jj])
            names_derived.append(temp_names[jj])

  return np.array(names_derived, dtype=str), np.array(derived_par, dtype=np.float64).T


def compare_par_post_angle(der_par, der_post):
  cospar = np.cos(der_par*deg2rad)
  sinpar = np.sin(der_par*deg2rad)
  
  par_scale = np.arctan2(sinpar, cospar)*rad2deg
  par_mod = (par_scale+360.)%360.
  
  cospos = np.cos(der_post*deg2rad)
  sinpos = np.sin(der_post*deg2rad)
  
  pos_scale = np.arctan2(sinpos, cospos)*rad2deg
  pos_mod = (pos_scale+360.)%360.

  delta_scale = np.abs(np.mean(pos_scale) - par_scale)
  delta_mod = np.abs(np.mean(pos_mod) - par_mod)
  
  if(np.abs(delta_scale - delta_mod) <= eps64bit):
    return par_mod, pos_mod
  else:
    if(delta_scale > delta_mod):
      return par_mod, pos_mod
    else:
      return par_scale, pos_scale
  return par_mod, pos_mod

def adjust_derived_parameters(derived_names, derived_par, derived_post):
  n_der = np.array(derived_par).shape[0]
  derived_par_adj = derived_par.copy()
  derived_post_adj = derived_post.copy()
  
  for i_d in range(0, n_der):
    if(derived_names[i_d][0] in ['w' ,'i'] or derived_names[i_d][0:2] in ['lN', 'mA']):
      derived_par_adj[i_d], derived_post_adj[:,i_d] = compare_par_post_angle(derived_par[i_d], derived_post[:,i_d])
      
  return derived_par_adj, derived_post_adj


def get_header(perc_val):

  top_header = '%1s %15s %20s %15s %15s %15s %15s %15s %15s %15s %15s %15s' %('#', '', '', '', '+/-1sigma', 'MAD', '-1sigma', '+1sigma', '-2sigma', '+2sigma', '-3sigma', '+3sigma')
  header = '%1s %15s %20s %15s' %('#', 'name', 'unit', 'parameter')
  perc_str = ' '.join(['%15s' %('%4.2f-th' %(perc_val[i_p])) for i_p in range(len(perc_val))])
  header = '%s %s' %(header, perc_str)
  
  return top_header, header

def print_parameters(top_header, header, name_parameters, unit_parameters, parameters, sigma_parameters, output=None):
  
  print_both(top_header, output)
  print_both(header, output)
  #n_par = parameters.shape[0]
  n_par = len(name_parameters)
  for i_p in range(0, n_par):
    sigma_line = ' '.join(['%15.10f' %(sigma_parameters[ii,i_p]) for ii in range(0,len(percentile_val))])
    line = '%17s %20s %15.10f %s' %(name_parameters[i_p], unit_parameters[i_p], parameters[i_p], sigma_line)
    print_both(line, output)
  
  return

# ------------------------------
