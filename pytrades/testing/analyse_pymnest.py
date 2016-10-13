#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
import h5py
import sys
from ancillary import *
from constants import Mjups, Msjup, Mears, Msear
from pytrades_lib import pytrades as trades
import logging
import warnings

from matplotlib import use as mpluse
mpluse("Agg")
#mpluse("Qt4Agg")
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)
import matplotlib.colors as colors

warnings.simplefilter('ignore', np.RankWarning)

def get_args():
  parser = argparse.ArgumentParser(description='TRADES+pyMultiNest')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')

  parser.add_argument('-l', '--ln-err', '--ln-err-const', action='store', dest='ln_flag', default=True, help='Computes or not constant to add to the lglikelihood: SUM( ln 2pi * sigma_obs^2 ). Default = True. Set to False if not use this value.')

  parser.add_argument('-pf', '--prefix', action='store', dest='prefix', default='trades_mnest_', help='Files prefix. Default = trades_mnest_')
  
  parser.add_argument('-m', '--mtype', '--mass-type', action='store', dest='m_type', default='e', help='Mass type: j = Jupiter, e = Earth, s = Sun. Default is Earth = e.')

  cli = parser.parse_args()
  
  #cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.full_path = os.path.abspath(cli.full_path)
  
  cli.ln_flag = set_bool_argument(cli.ln_flag)
  
  cli.m_type = cli.m_type.lower()
  
  return cli

def calculate_ecc_argp_phi(ecosw, esinw, mA):
  tolerance = np.finfo(float).eps
  ecc = np.asarray(np.sqrt(np.asarray([ecosw])**2 + np.asarray([esinw])**2), dtype=np.float64)
  ecc[ecc <= 0.] = tolerance
  ecc[ecc >= 1.] = 1. - tolerance
  argp = ((np.arctan2(np.asarray([esinw]), np.asarray([ecosw]))*180./np.pi)+360.)/360.
  argp[ecc <= tolerance] = 90.
  phi = (argp + mA + 360.)%360.
  return ecc, argp, phi

def compute_posterior_derived_parameters(parameter_names, posterior_parameters, NB):
  n_fit = np.asarray(parameter_names).shape[0]
  n_posterior = posterior_parameters.shape[0]
  n_derived = NB*3
  derived_parameters = np.zeros((n_posterior, n_derived))
  derived_names = []
  cnt_der = 0
  cnt_nb = 2
  for i in range(0, n_fit):
    if('ecosw' in parameter_names[i]):
      ecc, argp, phi = calculate_ecc_argp_phi(posterior_parameters[:,i], posterior_parameters[:,i+1], posterior_parameters[:,i+2])
      #print ecc.shape, argp.shape, phi.shape
      derived_parameters[:,cnt_der] = ecc
      derived_names.append('e%d' %(cnt_nb))
      derived_parameters[:,cnt_der+1] = argp
      derived_names.append('w%d' %(cnt_nb))
      derived_parameters[:,cnt_der+2] = phi
      derived_names.append('phi%d' %(cnt_nb))
      cnt_der = cnt_der + 3
      cnt_nb += 1

  return derived_parameters, derived_names

def convert_masses(parameter_names, parameters_minmax_in, posterior_parameters_in, m_factor):
  n_fit = parameter_names.shape[0]
  parameters_minmax_out = parameters_minmax_in.copy()
  posterior_parameters_out = posterior_parameters_in.copy()
  for i in range(0, n_fit):
    if (parameter_names[i][0] == 'm' and parameter_names[i][1] != 'A'):
      parameters_minmax_out[i,:] = parameters_minmax_in[i,:] * m_factor
      posterior_parameters_out[:,i] = posterior_parameters_in[:,i] * m_factor
  return parameters_minmax_out, posterior_parameters_out

def select_id_parameters(n_fit, n_derived, id_to_sel, posterior_parameters, derived_parameters):
    sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
    sel_parameters = np.zeros((n_fit,7))
    sel_parameters[:,0] = posterior_parameters[id_to_sel, :]
    for i in range(0, n_fit):
      sel_parameters[i,1:] = np.percentile(posterior_parameters[:,i] - sel_parameters[i,0], sigmas_percentiles, interpolation='lower')
    
    sel_derived = np.zeros((n_derived, 7))
    sel_derived[:,0] = derived_parameters[id_to_sel, :]
    for i in range(0, n_derived):
      sel_derived[i,1:] = np.percentile(derived_parameters[:,i] - sel_derived[i,0], sigmas_percentiles, interpolation='lower')
    
    return sel_parameters, sel_derived

def compute_mode_median(n_fit, n_derived, posterior_parameters, k, parameter_names, derived_parameters):
  mode_parameters = np.zeros((n_fit,7))
  mode_derived = np.zeros((n_derived, 7))
  median_parameters = np.zeros((n_fit,7))
  median_derived = np.zeros((n_derived, 7))
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  cnt_der = 0
  
  for i in range(0, n_fit):
    
    data = posterior_parameters[:,i]
    max_mean, max_bin = compute_max_mean(data, k)
    mode_parameters[i,0] = max_mean
    res = posterior_parameters[:,i] - max_mean
    sigma_mode = np.percentile(res, sigmas_percentiles, interpolation='lower')
    mode_parameters[i,1:] = sigma_mode
    median_parameters[i,0] = np.percentile(posterior_parameters[:,i], [50], interpolation='lower')
    median_parameters[i,1:] = np.percentile(posterior_parameters[:,i]-median_parameters[i,0], sigmas_percentiles, interpolation='lower')
    
    if('mA' in parameter_names[i]):
      #ecc, argp, phi = calculate_ecc_argp_phi(mode_parameters[i-2,0], mode_parameters[i-1,0], mode_parameters[i,0])

      ecc = derived_parameters[:,cnt_der]
      argp = derived_parameters[:,cnt_der+1]
      phi = derived_parameters[:,cnt_der+2]

      mode_ecc, max_bin = compute_max_mean(ecc, k)
      mode_argp, max_bin = compute_max_mean(argp, k)
      mode_phi, max_bin = compute_max_mean(phi, k)
      
      mode_derived[cnt_der,0] = mode_ecc
      mode_derived[cnt_der,1:] = np.percentile(ecc-mode_ecc, sigmas_percentiles, interpolation='lower')
      median_derived[cnt_der,0] = np.percentile(ecc, [50], interpolation='lower')
      median_derived[cnt_der,1:] = np.percentile(ecc-median_derived[cnt_der,0], sigmas_percentiles, interpolation='lower')
      
      mode_derived[cnt_der+1,0] = mode_argp
      mode_derived[cnt_der+1,1:] = np.percentile(argp-mode_argp, sigmas_percentiles,interpolation='lower')
      median_derived[cnt_der+1,0] = np.percentile(argp, [50], interpolation='lower')
      median_derived[cnt_der+1,1:] = np.percentile(argp-median_derived[cnt_der+1,0], sigmas_percentiles, interpolation='lower')
      
      mode_derived[cnt_der+2,0] = mode_phi
      mode_derived[cnt_der+2,1:] = np.percentile(phi-mode_phi, sigmas_percentiles,interpolation='lower')
      median_derived[cnt_der+2,0] = np.percentile(phi, [50], interpolation='lower')
      median_derived[cnt_der+2,1:] = np.percentile(phi-median_derived[cnt_der+2,0], sigmas_percentiles, interpolation='lower')
      
      cnt_der = cnt_der + 3
    
    #if ('esinw' in parameter_names[i]):
  return mode_parameters, mode_derived, median_parameters, median_derived

def print_values(n_fit, n_derived, parameter_names, values_parameters, values_derived, derived_names, kind, logger):
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  sigmas = ['-1','-2','-3','+1','+2','+3']
  line = '%10s %15s %s' %(' ', kind, ' '.join(['%ssigma(%5.2fres)'.rjust(19) %(sigmas[j],sigmas_percentiles[j]) for j in range(len(sigmas_percentiles))]))
  logger.info(line)
  for i in range(0, n_fit):
    sigma_line = ''.join(['%19.7f' %(values_parameters[i,j+1]) for j in range(len(sigmas_percentiles))])
    line = '%10s %15.7f %s' %(parameter_names[i], values_parameters[i,0], sigma_line)
    logger.info(line)
  logger.info('-'* n_line)
  for i in range(0, n_derived):
    sigma_line = ''.join(['%19.7f' %(values_derived[i,j+1]) for j in range(len(sigmas_percentiles))])
    line = '%10s %15.7f %s' %(derived_names[i], values_derived[i,0], sigma_line)
    logger.info(line)

def read_priors(priors_file):
  # Masses must be provided in M_earth!
  priors_type = []
  priors = []
  of = open(priors_file, 'r')
  lines = of.readlines()
  of.close()
  for line in lines:
    fields = line.strip().split()
    if(fields[0] != '#'):
      priors_type.append([fields[0], fields[1]])
      if(fields[1].lower() == 'g'):
        factor = 1.
        if(fields[0][0] == 'm' and fields[0][1] != 'A'): factor = Mears
        priors.append(np.asarray([fields[2], fields[3]], dtype=np.float64)*factor)
      else:
        priors.append(None)
  return priors, priors_type

# LNPRIOR TO BE ADDED TO LOGLHD
# it can use all the variables defined before this point!
def lnprior(fitting_parameters, ndim, fitting_priors, fitting_priors_type, derived_priors, derived_priors_type):
  lnprior_value = 0.
  i_der = 0
  for i in range(0, ndim):
    ln_temp = 0.
    # calculate the LogLikelihood<->prior of fitting parameter
    if(fitting_priors_type[i][1].lower() == 'g'):
      ln_temp = -0.5*(((fitting_parameters[i]-fitting_priors[i][0])/fitting_priors[i][1])**2)
      lnprior_value = lnprior_value + ln_temp
    # calculate the LogLikelihood<->prior of derived parameter
    if('mA' in parameter_names[i]):
      ln_temp = 0.
      ecc = np.sqrt(fitting_parameters[i-2]**2 + fitting_parameters[i-1]**2)
      if (ecc <= 0.):
        ecc = np.finfo(float).eps
      elif (ecc > 1.):
        ecc = 1.- np.finfo(float).eps
      # ecc prior
      if(derived_priors_type[i_der][1].lower() == 'g'):
        ln_temp = -0.5*(((ecc-derived_priors[i_der][0])/derived_priors[i_der][1])**2)
        lnprior_value = lnprior_value + ln_temp
      # phi prior
      if(derived_priors_type[i_der+1][1].lower() == 'g'):
        if(ecc <= np.finfo(float).eps):
          argp = 90.
        else:
          argp = ((np.arctan2(fitting_parameters[i-1], fitting_parameters[i-2])*180./np.pi)+360.)%360.
        phi = (argp + fitting_parameters[i] + 360.)%360.
        ln_temp = 0.
        ln_temp = -0.5*(((phi-derived_priors[i_der+1][0])/derived_priors[i_der+1][1])**2)
        lnprior_value = lnprior_value + ln_temp
        i_der = i_der + 2
  return lnprior_value


# ---
# initialize logger
logger = logging.getLogger("Main_log")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(message)s")

# READ COMMAND LINE ARGUMENTS
cli = get_args()
# RENAME 
working_path = cli.full_path

plot_folder = prepare_plot_folder(working_path)
log_file = os.path.join(plot_folder, 'pymnest_triangle_log.txt')

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
m_factor = mass_conversion_factor(cli.m_type)

trades_hdf5 = h5py.File(os.path.join(working_path, 'system_summary.hdf5'), 'r')
#parameter_names = trades_hdf5['parameter_names'][:]
#parameters_minmax = trades_hdf5['parameters_minmax'][:,:]
#ln_err_const = trades_hdf5['ln_err_const'][:]
parameter_names = np.array(trades_hdf5['parameter_names'])
parameters_minmax = np.array(trades_hdf5['parameters_minmax'])
ln_err_const = np.array(trades_hdf5['ln_err_const'])
n_fit = parameter_names.shape[0]
#for i in range(0, n_fit):
  #print '%15s %12.7f %12.7f' %(parameter_names[i], parameters_minmax[i,0], parameters_minmax[i,1])

pymnest_file = os.path.join(working_path, '%s.txt' %(cli.prefix))

posterior_pymnest = np.genfromtxt(pymnest_file)
if(len(posterior_pymnest.shape) == 1):
  logger.info('n_posterior = 1 stop script')
  sys.exit()
posterior_parameters = posterior_pymnest[:,2:]
n_posterior = posterior_parameters.shape[0]
if(n_posterior < 7):
  logger.info('n_posterior < 7 stop script')
  sys.exit()


fitting_prior_file = os.path.join(os.path.dirname(working_path), 'fitting_priors.dat')
fitting_priors, fitting_priors_type = read_priors(fitting_prior_file)

derived_prior_file = os.path.join(os.path.dirname(working_path), 'derived_priors.dat')
derived_priors, derived_priors_type = read_priors(derived_prior_file)
n_der_priors = len(derived_priors)

NB=3
derived_parameters, derived_names = compute_posterior_derived_parameters(parameter_names, posterior_parameters, NB)
n_derived = NB * 3

#n_ecc = count_eccentricities(parameter_names)
#ecc, ln_ecc, ln_ecc_value = compute_ecc_ln_ecc(parameter_names, posterior_parameters, n_ecc)

lnprior_values = np.zeros((n_posterior))
for i in range(0, n_posterior):
  lnprior_values[i] = lnprior(posterior_parameters[i,:], n_fit, fitting_priors, fitting_priors_type, derived_priors, derived_priors_type)

# Msun to Meart
parameters_minmax, posterior_parameters = convert_masses(parameter_names, parameters_minmax, posterior_parameters, m_factor)

#fitness = posterior_pymnest[:,1] - (-2.*ln_err_const) - (-2.*ln_ecc_value)
lnlkhd = -0.5*posterior_pymnest[:,1]
fitness = -2.* ( lnlkhd - ln_err_const - lnprior_values )
fitness_min_id = np.argmin(fitness)
lnlkhd_max_id = np.argmax(lnlkhd)

k = np.ceil(2. * n_posterior**(1./3.)).astype(int)
if(k>40): k=40

mode_parameters, mode_derived, median_parameters, median_derived = compute_mode_median(n_fit, n_derived, posterior_parameters, k, parameter_names, derived_parameters)

#min_fitness_parameters[:,0] = posterior_parameters[fitness_min_id,:]
#min_fitness_derived = derived_parameters[fitness_min_id,:]
min_fitness_parameters, min_fitness_derived = select_id_parameters(n_fit, n_derived, fitness_min_id, posterior_parameters, derived_parameters)

#maxlnlklhd_parameters = posterior_parameters[lnlkhd_max_id,:]
#maxlnlklhd_derived = derived_parameters[lnlkhd_max_id,:]
maxlnlklhd_parameters, maxlnlklhd_derived = select_id_parameters(n_fit, n_derived, lnlkhd_max_id, posterior_parameters, derived_parameters)

logger.info('')
logger.info('ANALYSIS OF TRADES-PYMNEST POSTERIOR FILE')
logger.info('')

# PLOT CORRELATIONS

# set label and legend names
kel_legends, labels_list = keplerian_legend(parameter_names, cli.m_type)

label_separation=-0.75
label_pad = 16
label_size = 8
ticklabel_size = 5


def set_xaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_label, ticks_formatter):
  ax.get_xaxis().set_visible(True)
  ax.xaxis.set_tick_params(labelsize=ticklabel_size)
  ax.xaxis.set_label_coords(0.5, label_separation)
  ax.ticklabel_format(style='plain', axis='both', useOffset=False)
  plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
  ax.xaxis.labelpad = label_pad
  ax.set_xlabel(kel_label, fontsize=label_size)
  tick_step = (ticks_formatter[1] - ticks_formatter[0]) / ticks_formatter[2]
  ax.xaxis.set_ticks(np.arange(ticks_formatter[0], ticks_formatter[1], tick_step))

def set_yaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_label, ticks_formatter):
  ax.get_yaxis().set_visible(True)
  ax.yaxis.set_tick_params(labelsize=ticklabel_size)
  ax.yaxis.set_label_coords(label_separation,0.5)
  ax.ticklabel_format(style='plain', axis='both', useOffset=False)
  ax.yaxis.labelpad = label_pad
  ax.set_ylabel(kel_label, fontsize=label_size)
  tick_step = (ticks_formatter[1] - ticks_formatter[0]) / ticks_formatter[2]
  ax.yaxis.set_ticks(np.arange(ticks_formatter[0], ticks_formatter[1], tick_step))


fig = plt.figure(figsize=(12,12))
fig.subplots_adjust(hspace=0.05, wspace=0.05)

for ii in range(0, n_fit, 1):
  x_data = posterior_parameters[:, ii]
  x_med = np.median(x_data)
  x_16 = np.percentile(x_data, 16)
  x_84 = np.percentile(x_data, 84)
  x_min, x_max = compute_limits(x_data, 0.05)
  
  # compute mean of higher peak/bin
  x_max_mean, x_max_bin = compute_max_mean(x_data, k)

  for jj in range(n_fit-1, -1, -1):
    y_data = posterior_parameters[:, jj]
    y_med = np.median(y_data)
    y_min, y_max = compute_limits(y_data, 0.05)
    
    # compute mean of higher peak/bin
    y_max_mean, y_max_bin = compute_max_mean(y_data, k)
    
    if(jj > ii): # correlation plot
      logger.info('%s vs %s' %(labels_list[ii], labels_list[jj]))
      #print ' (%.4f) %.4f <= X <= %.4f (%.4f)' %(parameters_minmax[ii,0], x_data.min(), x_data.max(), parameters_minmax[ii,1])
      #print ' (%.4f) %.4f <= Y <= %.4f (%.4f)' %(parameters_minmax[jj,0], y_data.min(), y_data.max(), parameters_minmax[jj,1])

      ax = plt.subplot2grid((n_fit+1, n_fit), (jj,ii))
      
      hist2d_counts, xedges, yedges, image2d = ax.hist2d(x_data, y_data, bins=k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], cmap=plt.get_cmap('Greys'), normed=True)
      
      # plot mean_mode
      ax.axvline(x_max_mean, color='forestgreen', ls='-.', lw=0.8, alpha=0.7)
      ax.axhline(y_max_mean, color='forestgreen', ls='-.', lw=0.8, alpha=0.7)
      
      new_k = k
      hist2d_counts_2, xedges_2, yedges_2 = np.histogram2d(x_data, y_data, bins=new_k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], normed=True)
      x_bins = [0.5*(xedges_2[i]+xedges_2[i+1]) for i in range(0, new_k)]
      y_bins = [0.5*(yedges_2[i]+yedges_2[i+1]) for i in range(0, new_k)]
      ax.contour(x_bins, y_bins, hist2d_counts_2.T, 3, colors=('forestgreen', 'royalblue', 'red'), linestyle='solid', linewidths=(0.5, 0.5, 0.5))
      
      ax.axvline(x_med, color='dimgrey', ls='-', lw=0.5)
      ax.axhline(y_med, color='dimgrey', ls='-', lw=0.5)
      
      ax.get_xaxis().set_visible(False)
      ax.get_yaxis().set_visible(False)
      if(jj == n_fit-1):
        set_xaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_legends[ii], [xedges[0], xedges[-1], 3])
      if(ii == 0): 
        set_yaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_legends[jj], [yedges[0], yedges[-1], 5])
      
      ax.set_ylim([y_min, y_max])
      ax.set_xlim([x_min, x_max])
      plt.draw()
    
    elif(jj == ii): # distribution plot
      logger.info('%s histogram' %(labels_list[ii]))
      #print ' (%.4f) %.4f <= X <= %.4f (%.4f)' %(parameters_minmax[ii,0], x_data.min(), x_data.max(), parameters_minmax[ii,1])
      
      ax = plt.subplot2grid((n_fit+1, n_fit), (ii,ii))
      if (ii == n_fit-1):
        hist_orientation='horizontal'
      else:
        hist_orientation='vertical'
        
      hist_counts, edges, patces = ax.hist(x_data, bins=k, range=[x_data.min(), x_data.max()], histtype='stepfilled', color='darkgrey', align='mid', orientation=hist_orientation, normed=True)
      
      if (ii == n_fit-1):
        ax.set_ylim([y_min, y_max])
        ax.axhline(x_med, color='dimgrey', ls='-', lw=0.5)
        ax.axhline(x_16, color='dimgrey', ls='--', lw=0.5)
        ax.axhline(x_84, color='dimgrey', ls='--', lw=0.5)
        # plot mean_mode
        ax.axhline(x_max_mean, color='forestgreen', ls='-.', lw=0.8, alpha=0.7)
      else:
        ax.set_xlim([x_min, x_max])
        ax.axvline(x_med, color='dimgrey', ls='-', lw=0.5)
        ax.axvline(x_16, color='dimgrey', ls='--', lw=0.5)
        ax.axvline(x_84, color='dimgrey', ls='--', lw=0.5)
        # plot mean_mode
        ax.axvline(x_max_mean, color='forestgreen', ls='-.', lw=0.8, alpha=0.7)
      
      ax.get_xaxis().set_visible(False)
      ax.get_yaxis().set_visible(False)
      ax.set_title(kel_legends[ii], fontsize=label_size)

      plt.draw()

logger.info('saving plot')
correlation_fig_file = os.path.join(plot_folder, 'pymnest_triangle.png')
fig.savefig(correlation_fig_file, bbox_inches='tight', dpi=300)
plt.close(fig)
logger.info('')

n_line = 120

logger.info('')
logger.info('posterior sample: %d elements' %(posterior_pymnest.shape[0]))
logger.info('Number of bins: k = %d' %(k))
logger.info('%15s = %19.7f %3s = %19.7f [id = %5d]'%('lnlkhd min', lnlkhd.min(), 'max',lnlkhd.max(), lnlkhd_max_id))
logger.info('%15s = %19.7f %3s = %19.7f [id = %5d]'%('fitness max', fitness.max(), 'min', fitness.min(), fitness_min_id))
logger.info('-'* n_line)
#logger.info('%15s = %20s | %20s' %('id','par:min(fitness)', 'par:max(lnlkhd)'))
#for i in range(0, n_fit):
  #logger.info('%15s = %20.12f | %20.12f' %(parameter_names[i], posterior_parameters[fitness_min_id,i], posterior_parameters[lnlkhd_max_id,i]))
print_values(n_fit, n_derived, parameter_names, min_fitness_parameters, min_fitness_derived, derived_names, 'min fitness', logger)
logger.info('-'* n_line)
logger.info('%15s = %20.12f' %('ln_err_const', ln_err_const[0]))
logger.info('%15s = %20.12f' %('ln_priors', lnprior_values[fitness_min_id]))
logger.info('-'* n_line)
logger.info('%15s = %20.12f' %('fitness', fitness[fitness_min_id]))
logger.info('%15s = %20.12f' %('lnlkhd', lnlkhd[fitness_min_id]))
logger.info('-'* n_line)
print_values(n_fit, n_derived, parameter_names, maxlnlklhd_parameters, maxlnlklhd_derived, derived_names, 'max lklkhd', logger)
logger.info('-'* n_line)
logger.info('%15s = %20.12f' %('ln_err_const', ln_err_const[0]))
logger.info('%15s = %20.12f' %('ln_priors', lnprior_values[lnlkhd_max_id]))
logger.info('-'* n_line)
logger.info('%15s = %20.12f' %('fitness', fitness[lnlkhd_max_id]))
logger.info('%15s = %20.12f' %('lnlkhd', lnlkhd[lnlkhd_max_id]))
logger.info('-'* n_line)
logger.info('-'* n_line)
print_values(n_fit, n_derived, parameter_names, mode_parameters, mode_derived, derived_names, 'mode', logger)
logger.info('-'* n_line)
logger.info('-'* n_line)
print_values(n_fit, n_derived, parameter_names, median_parameters, median_derived, derived_names, 'median', logger)
logger.info('-'* n_line)
logger.info('-'* n_line)

# INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
trades.initialize_trades(os.path.join(os.path.dirname(working_path), ''), os.path.basename(working_path))
trades.path_change(os.path.join(working_path, ''))

# Mearth to Msun: needed by trades
parameters_minmax, posterior_parameters = convert_masses(parameter_names, parameters_minmax, posterior_parameters, 1./m_factor)

mode_parameters, mode_derived, median_parameters, median_derived = compute_mode_median(n_fit, n_derived, posterior_parameters, k, parameter_names, derived_parameters)
min_fitness_parameters = posterior_parameters[fitness_min_id,:]
maxlnlklhd_parameters = posterior_parameters[lnlkhd_max_id,:]

logger.info('')

# functions needed to redirecting fortran stdout to file and read to print with python
def save_temp(working_path):
  temp_out_file_name = os.path.join(working_path, 'temp.txt')
  outfile = os.open(temp_out_file_name, os.O_RDWR|os.O_CREAT)
  # save the current file descriptor
  save = os.dup(1)
  # put outfile on 1
  os.dup2(outfile, 1)
  # end magic
  return save, outfile, temp_out_file_name

def print_and_remove_temp(save, outfile, temp_out_file_name, logger):
  # restore the standard output file descriptor
  os.dup2(save, 1)
  # close the output file
  os.close(outfile)
  # read it
  of = open(temp_out_file_name, 'r')
  for line in of.readlines():
    logger.info(line.strip())
  of.close()
  # delete it
  os.remove(temp_out_file_name)
  return

logger.info('')
logger.info('**FITNESS MIN PARAMETERS FIT**')
save, outfile, temp_out_file_name = save_temp(working_path)
fitness_1001, lgllhd_1001, check_1001 = trades.write_summary_files(1001, min_fitness_parameters)
print_and_remove_temp(save, outfile, temp_out_file_name, logger)

logger.info('')
logger.info('**LNLKHD MAX PARAMETERS FIT**')
save, outfile, temp_out_file_name = save_temp(working_path)
fitness_1999, lgllhd_1999, check_1999 = trades.write_summary_files(1999, maxlnlklhd_parameters)
print_and_remove_temp(save, outfile, temp_out_file_name, logger)

logger.info('')
logger.info('**MEAN OF MAX BIN PARAMETERS FIT**')
save, outfile, temp_out_file_name = save_temp(working_path)
fitness_1666, lgllhd_1666, check_1666 = trades.write_summary_files(1666, mode_parameters[:,0])
print_and_remove_temp(save, outfile, temp_out_file_name, logger)

logger.info('')
logger.info('**MEDIAN PARAMETERS FIT**')
save, outfile, temp_out_file_name = save_temp(working_path)
fitness_1050, lgllhd_1050, check_1050 = trades.write_summary_files(1050, median_parameters[:,0])
print_and_remove_temp(save, outfile, temp_out_file_name, logger)
save, outfile, temp_out_file_name = save_temp(working_path)

trades.deallocate_variables()

logger.info('')
logger.info('FINISHED')
logger.info('')
