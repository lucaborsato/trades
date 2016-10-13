#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os
import sys
import numpy as np # array
import ancillary as anc

from matplotlib import use as mpluse
mpluse("Agg")
#mpluse("Qt4Agg")
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)
import matplotlib.colors as colors
#from matplotlib import rcParams
#rcParams['text.latex.unicode']=True
from matplotlib.ticker import FormatStrFormatter

    
label_separation=-0.25
label_pad = 16
label_size = 10
ticklabel_size = 5


def set_xaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_label, ticks_formatter, tick_fmt='%.4f'):
  ax.get_xaxis().set_visible(True)
  ax.xaxis.set_tick_params(labelsize=ticklabel_size)
  ax.xaxis.set_label_coords(0.5, label_separation)
  #ax.ticklabel_format(style='plain', axis='both', useOffset=False)
  plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
  ax.xaxis.labelpad = label_pad
  ax.set_xlabel(kel_label, fontsize=label_size)
  tick_step = (ticks_formatter[1] - ticks_formatter[0]) / ticks_formatter[2]
  ax.xaxis.set_ticks(np.arange(ticks_formatter[0], ticks_formatter[1], tick_step))
  tick_formatter = FormatStrFormatter(tick_fmt)
  ax.xaxis.set_major_formatter(tick_formatter)
  return
  
def set_yaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_label, ticks_formatter, tick_fmt='%.4f'):
  ax.get_yaxis().set_visible(True)
  ax.yaxis.set_tick_params(labelsize=ticklabel_size)
  ax.yaxis.set_label_coords(label_separation,0.5)
  #ax.ticklabel_format(style='plain', axis='both', useOffset=False)
  ax.yaxis.labelpad = label_pad
  ax.set_ylabel(kel_label, fontsize=label_size)
  tick_step = (ticks_formatter[1] - ticks_formatter[0]) / ticks_formatter[2]
  ax.yaxis.set_ticks(np.arange(ticks_formatter[0], ticks_formatter[1], tick_step))
  tick_formatter = FormatStrFormatter(tick_fmt)
  ax.yaxis.set_major_formatter(tick_formatter)
  return
    

def plot_correlation_K19_ilN_3_4(plot_id, full_path, nder, sel_derived_flatposterior, kel_legends, divide_plot):

  median_parameters, median_perc68, median_confint = anc.get_median_parameters(sel_derived_flatposterior)
  median = {'parameters': median_parameters, 'sigma': median_perc68, 'intervals': median_confint}

  k = np.ceil(2. * sel_derived_flatposterior.shape[0]**(1./3.)).astype(int)
  if(k>11): k=11
  mode_bin, mode_parameters, mode_perc68, mode_confint = anc.get_mode_parameters(sel_derived_flatposterior, k)
  mode = {'parameters': mode_parameters, 'sigma': mode_perc68, 'intervals': mode_confint}
  
  fig = plt.figure(figsize=(12,12))
  fig.subplots_adjust(hspace=0.05, wspace=0.05)

  
  for ider in range(0, nder, 1):
    x_data = sel_derived_flatposterior[:, ider]
    x_med = median_parameters[ider]
    x_min, x_max = anc.compute_limits(x_data, 0.05)
    x_max_mean = mode_parameters[ider]
    divide_x = divide_plot[ider]
    
    for jder in range(nder-1, -1, -1):
      y_data = sel_derived_flatposterior[:, jder]
      y_med = median_parameters[jder]
      y_min, y_max = anc.compute_limits(y_data, 0.05)
      y_max_mean = mode_parameters[jder]
      divide_y = divide_plot[jder]
      
      if(jder > ider):
        ax = plt.subplot2grid((nder+1, nder), (jder,ider))
      
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
        
        # to remove in case
        if(divide_x is not None): ax.axvline(divide_x, color='red', ls='-', lw=1)
        if(divide_y is not None): ax.axhline(divide_y, color='red', ls='-', lw=1)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        if(jder == nder-1):
          set_xaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_legends[ider], [xedges[0], xedges[-1], 3])
        if(ider == 0): 
          set_yaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_legends[jder], [yedges[0], yedges[-1], 5])
        
        ax.set_ylim([y_min, y_max])
        ax.set_xlim([x_min, x_max])
        plt.draw()
      
      elif(jder == ider): # distribution plot
        ax = plt.subplot2grid((nder+1, nder), (ider,ider))
        if (ider == nder-1):
          hist_orientation='horizontal'
        else:
          hist_orientation='vertical'
          
        hist_counts, edges, patces = ax.hist(x_data, bins=k, range=[x_data.min(), x_data.max()], histtype='stepfilled', color='darkgrey', align='mid', orientation=hist_orientation, normed=True)
        
        if (ider == nder-1):
          ax.set_ylim([y_min, y_max])
          ax.axhline(x_med, color='dimgrey', ls='-', lw=0.5)
          # plot mean_mode
          ax.axhline(x_max_mean, color='forestgreen', ls='-.', lw=0.8, alpha=0.7)
          
          # to remove in case
          if(divide_x is not None): ax.axhline(divide_x, color='red', ls='-', lw=1)
          
        else:
          ax.set_xlim([x_min, x_max])
          ax.axvline(x_med, color='dimgrey', ls='-', lw=0.5)
          # plot mean_mode
          ax.axvline(x_max_mean, color='forestgreen', ls='-.', lw=0.8, alpha=0.7)
        
          # to remove in case
          if(divide_y is not None): ax.axvline(divide_y, color='red', ls='-', lw=1)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title(kel_legends[ider], fontsize=label_size)
        if(divide_x is not None): ax.set_title('%s [red line at %.2f]' %(kel_legends[ider], divide_x), fontsize=label_size)

        plt.draw()
  
  emcee_plots = anc.prepare_emcee_plot_folder(full_path)
  emcee_fig_file = os.path.join(emcee_plots, '%s_derived_emcee_triangle.png' %(plot_id))
  fig.savefig(emcee_fig_file, bbox_inches='tight', dpi=300)
  #fig.savefig(emcee_fig_file, dpi=300)
  emcee_fig_file = os.path.join(emcee_plots, '%s_derived_emcee_triangle.pdf' %(plot_id))
  fig.savefig(emcee_fig_file, bbox_inches='tight', dpi=300)
  #fig.savefig(emcee_fig_file, dpi=300)
  
  return median, mode

def derived_correlation_K19_ilN_3_4():
  
  # take arguments by command line
  cli = anc.get_args()
  
  emcee_folder = cli.full_path
  trades_folder = os.path.join(os.path.dirname(cli.full_path), '')
  # and best folder
  emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(emcee_folder, cli.temp_status)
  
  # read emcee data
  parameter_names_emcee, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)
  
  nfit, nwalkers, nruns, npost, nruns_sel = anc.get_emcee_parameters(chains, cli.temp_status, cli.npost, completed_steps)
  
  m_factor = anc.get_proper_mass(cli.m_type, parameter_names_emcee, trades_folder)

  chains_T, parameter_boundaries_2 = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries, chains)

  flatchain_posterior_0 = chains_T[:,:,:].reshape((-1, nfit))
  
  derived_names, derived_chains_T, derived_flatposterior = anc.get_derived_posterior_parameters(parameter_names_emcee, chains_T, flatchain_posterior_0)
  nder = len(derived_names)
  
  kel_legends = ['$i_\\mathrm{c}$', '$\Omega_\\mathrm{c}$', '$i_\\mathrm{d}$', '$\Omega_\\mathrm{d}$']
  
  inc3_lim = 88.724
  inc4_lim = 89.245
  
  divide_plot = [None, 181., None, None]  
  print_divide = ' '.join([str(divide_plot[ider]) for ider in range(0, nder)])
  print 'divide plot at: %s' %(print_divide)
  
  #sel_inc = np.logical_and(derived_flatposterior[:,0]<inc3_lim, derived_flatposterior[:,2]<inc4_lim)
  sel_inc = derived_flatposterior[:,0]<inc3_lim

  sel_derived_flatposterior = derived_flatposterior[sel_inc, :]
  
  print 'PLOTTING CORRELATION PLOT WITH ANGLES DEFINED BETWEEN -180 AND 180 DEGREES ... ',
  median_scale, mode_scale = plot_correlation_K19_ilN_3_4('scale', cli.full_path, nder, sel_derived_flatposterior, kel_legends, divide_plot)
  print 'DONE'
  
  sel_derived_flatposterior_module = sel_derived_flatposterior.copy()
  sel_derived_flatposterior_module[:,1] = (sel_derived_flatposterior_module[:,1]+360.)%360.
  sel_derived_flatposterior_module[:,3] = (sel_derived_flatposterior_module[:,3]+360.)%360.
  
  print 'PLOTTING CORRELATION PLOT WITH ANGLES DEFINED BETWEEN 0 AND 360 DEGREES ... ',
  median_module, mode_module = plot_correlation_K19_ilN_3_4('module', cli.full_path, nder, sel_derived_flatposterior_module, kel_legends, divide_plot)
  print 'DONE'
  
  sel_lN3_1 = sel_derived_flatposterior_module[:,1] <= divide_plot[1]
  print 'PLOTTING CORRELATION PLOT WITH ANGLES DEFINED BETWEEN 0 AND 360 DEGREES AND OMEGA < %.3f... ' %(divide_plot[1]),
  median_lN3_1, mode_lN3_1 = plot_correlation_K19_ilN_3_4('lN3_1', cli.full_path, nder, sel_derived_flatposterior_module[sel_lN3_1,:], kel_legends, [None]*4)
  print 'DONE'
  anc.print_parameters(derived_names, mode_lN3_1['parameters'], mode_lN3_1['sigma'], mode_lN3_1['intervals'], 'lN3 <= 181 deg')
  
  
  print
  sel_lN3_2 = sel_derived_flatposterior_module[:,1] > divide_plot[1]
  print 'PLOTTING CORRELATION PLOT WITH ANGLES DEFINED BETWEEN 0 AND 360 DEGREES AND OMEGA > %.3f... ' %(divide_plot[1]),
  median_lN3_2, mode_lN3_2 = plot_correlation_K19_ilN_3_4('lN3_2', cli.full_path, nder, sel_derived_flatposterior_module[sel_lN3_2,:], kel_legends, [None]*4)
  print 'DONE'
  anc.print_parameters(derived_names, mode_lN3_2['parameters'], mode_lN3_2['sigma'], mode_lN3_2['intervals'], 'lN3 > 181 deg')
  
  
  return


def derived_correlation_K19_ilN_3_4_000():
  # take arguments by command line
  cli = anc.get_args()
  
  emcee_folder = cli.full_path
  trades_folder = os.path.join(os.path.dirname(cli.full_path), '')
  # and best folder
  emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(emcee_folder, cli.temp_status)
  
  # read emcee data
  parameter_names_emcee, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)
  
  nfit, nwalkers, nruns, npost, nruns_sel = anc.get_emcee_parameters(chains, cli.temp_status, cli.npost, completed_steps)
  
  m_factor = anc.get_proper_mass(cli.m_type, parameter_names_emcee, trades_folder)

  chains_T, parameter_boundaries_2 = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries, chains)

  flatchain_posterior_0 = chains_T[:,:,:].reshape((-1, nfit))
  
  derived_names, derived_chains_T, derived_flatposterior = anc.get_derived_posterior_parameters(parameter_names_emcee, chains_T, flatchain_posterior_0)
  nder = len(derived_names)
  
  kel_legends = ['$i_\\mathrm{c}$', '$\Omega_\\mathrm{c}$', '$i_\\mathrm{d}$', '$\Omega_\\mathrm{d}$']

  divide_plot = [90., None, None, None]  
  print_divide = ' '.join([str(divide_plot[ider]) for ider in range(0, nder)])
  print 'divide plot at: %s' %(print_divide)
  
  print 'PLOTTING CORRELATION PLOT WITH ANGLES DEFINED BETWEEN -180 AND 180 DEGREES ... ',
  median_scale, mode_scale = plot_correlation_K19_ilN_3_4('scale', cli.full_path, nder, derived_flatposterior, kel_legends, divide_plot)
  print 'DONE'
  
  derived_flatposterior_module = derived_flatposterior.copy()
  derived_flatposterior_module[:,1] = (derived_flatposterior_module[:,1]+360.)%360.
  derived_flatposterior_module[:,3] = (derived_flatposterior_module[:,3]+360.)%360.
  
  print 'PLOTTING CORRELATION PLOT WITH ANGLES DEFINED BETWEEN 0 AND 360 DEGREES ... ',
  median_module, mode_module = plot_correlation_K19_ilN_3_4('module', cli.full_path, nder, derived_flatposterior_module, kel_legends, divide_plot)
  print 'DONE'
  
  anc.print_parameters(derived_names, median_module['parameters'], median_module['sigma'], median_module['intervals'], 'median_module')
  
  anc.print_parameters(derived_names, mode_module['parameters'], mode_module['sigma'], mode_module['intervals'], 'mode_module')
  
  sel_inc2 = derived_flatposterior_module[:,0] < divide_plot[0]
  print 'PLOTTING CORRELATION PLOT WITH ANGLES DEFINED BETWEEN 0 AND 360 DEGREES AND INC_C < ', divide_plot[0],' ... ',
  median_module_sel, mode_module_sel = plot_correlation_K19_ilN_3_4('module_sel_lt90', cli.full_path, nder, derived_flatposterior_module[sel_inc2, :], kel_legends, divide_plot)
  print 'DONE'
  anc.print_parameters(derived_names, median_module_sel['parameters'], median_module_sel['sigma'], median_module_sel['intervals'], 'median_module')
  anc.print_parameters(derived_names, mode_module_sel['parameters'], mode_module_sel['sigma'], mode_module_sel['intervals'], 'mode_module_sel')
  
  sel_inc2 = derived_flatposterior_module[:,0] > divide_plot[0]
  print 'PLOTTING CORRELATION PLOT WITH ANGLES DEFINED BETWEEN 0 AND 360 DEGREES AND INC_C >', divide_plot[0],' ... ',
  median_module_sel, mode_module_sel = plot_correlation_K19_ilN_3_4('module_sel_gt90', cli.full_path, nder, derived_flatposterior_module[sel_inc2, :], kel_legends, divide_plot)
  print 'DONE'
  anc.print_parameters(derived_names, median_module_sel['parameters'], median_module_sel['sigma'], median_module_sel['intervals'], 'median_module')
  anc.print_parameters(derived_names, mode_module_sel['parameters'], mode_module_sel['sigma'], mode_module_sel['intervals'], 'mode_module_sel')
  
  return

# 
# Here run the proper function
#

#derived_correlation_K19_ilN_3_4()
derived_correlation_K19_ilN_3_4_000()

  
  
