#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import sys
import argparse
import os
import numpy as np # array
import ancillary as anc # not so good...but fast
import random
import logging
import warnings
import h5py
from scipy.stats import norm as scipy_norm

import matplotlib.cm as cm
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

warnings.simplefilter('ignore', np.RankWarning)

def main():
  # ---
  # initialize logger

  logger = logging.getLogger("Main_log")
  logger.setLevel(logging.DEBUG)
  formatter = logging.Formatter("%(asctime)s - %(message)s")

  # global variables

  label_separation=-0.90
  label_pad = 16
  label_size = 7
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

  print 
  print ' ================== '
  print ' CORRELATION PLOTS'
  print ' ================== '
  print

  # read cli arguments
  cli = anc.get_args()

  #plot_folder = prepare_plot_folder(working_path)
  emcee_plots = anc.prepare_emcee_plot_folder(cli.full_path)
  log_file = os.path.join(emcee_plots, 'emcee_triangle_log.txt')
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
  #m_factor = anc.mass_conversion_factor(cli.m_type)
  m_factor, m_unit = anc.mass_type_factor(1., cli.m_type, False)

  # set emcee and trades folder
  emcee_folder = cli.full_path
  trades_folder = os.path.join(os.path.dirname(cli.full_path), '')
  # and best folder
  emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(emcee_folder, cli.temp_status)

  # get data from the hdf5 file
  parameter_names_emcee, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)
  # print Memory occupation of ...
  anc.print_memory_usage(chains)

  nfit, nwalkers, nruns, npost, nruns_sel = anc.get_emcee_parameters(chains, cli.temp_status, cli.npost, completed_steps)
  logger.info('nfit(%d), nwalkers(%d), nruns(%d), npost(%d), nruns_sel(%d)' %(nfit, nwalkers, nruns, npost, nruns_sel))

  # test label_separation
  #if (nfit <= 3): label_separation = -0.1
  if(nfit > 2):
    label_separation = -0.1 - ( 0.075 * (nfit-2) )
  else:
    label_separation = -0.15

  # set label and legend names
  kel_plot_labels = anc.keplerian_legend(parameter_names_emcee, cli.m_type)

  chains_T_full, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries, chains)

  chains_T, flatchain_posterior_0, lnprob_burnin, thin_steps = anc.thin_the_chains(cli.use_thin, npost, nruns, nruns_sel, autocor_time, chains_T_full, lnprobability, burnin_done=False)

  if(cli.boot_id > 0):
    flatchain_posterior_msun = anc.posterior_back_to_msun(m_factor,parameter_names_emcee,flatchain_posterior_0)
    boot_file = anc.save_bootstrap_like(emcee_folder, cli.boot_id, parameter_names_emcee, flatchain_posterior_msun)
    logger.info('saved bootstrap like file: %s' %(boot_file))
    del flatchain_posterior_msun

  ## TO BE SUBSTITUTED WITH PROPER READ OF summary_parameters.hdf5 FILE
  ## GET MAX LNPROBABILITY AND PARAMETERS -> id 40XX
  #max_lnprob, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint = anc.get_maxlnprob_parameters(lnprob_burnin, chains_T, flatchain_posterior_0)

  ## std way: median of the posterior parameter distribution -> id 10XX
  #median_parameters, median_perc68, median_confint = anc.get_median_parameters(flatchain_posterior_0)

  ## MODE-LIKE PARAMETERS -> id 30XX
  ## take the mean of 5 bin centered to the higher bin
  #k = np.ceil(2. * flatchain_posterior_0.shape[0]**(1./3.)).astype(int)
  ##if(k>50): k=50
  #if(k>20): k=20
  k = anc.get_bins(flatchain_posterior_0, rule='doane')
  
  #mode_bin, mode_parameters, mode_perc68, mode_confint = anc.get_mode_parameters_full(flatchain_posterior_0, k)

  #name_par, name_excluded = anc.get_sample_list(cli.sample_str, parameter_names_emcee)
  #sample_parameters, idx_sample = anc.pick_sample_parameters(flatchain_posterior_0, parameter_names_emcee, name_par = name_par, name_excluded = name_excluded)
  
  if(cli.overplot is not None):
    
    if(cli.adhoc is not None):
      overp_names, read_par = anc.read_fitted_file(cli.adhoc)
      cli.overplot = 777
    else:
      ## OPEN summary_parameters.hdf5 FILE
      s_h5f = h5py.File(os.path.join(cli.full_path, 'summary_parameters.hdf5'), 'r')
      # take only the selected sample
      s_overplot = '%04d' %(cli.overplot)
      read_par = s_h5f['parameters/%s/fitted/parameters' %(s_overplot)][...]
      s_h5f.close()
      
    # fitted parameters has always Mp/Ms in Msun/Mstar, so it is needed to rescale it properly
    overp_par = read_par.copy()
    for ii in range(0, nfit):
      if('Ms' in parameter_names_emcee[ii]):
      #if('Ms' in overp_names[ii]):
        overp_par[ii] = overp_par[ii]*m_factor

  #fig, ax = plt.subplots(nrows = nfit-1, ncols=nfit, figsize=(12,12))
  fig = plt.figure(figsize=(12,12))
  fig.subplots_adjust(hspace=0.05, wspace=0.05)

  for ix in range(0, nfit, 1):
    x_data = flatchain_posterior_0[:, ix]
    ##x_med = median_parameters[ix]
    
    x_min, x_max = anc.compute_limits(x_data, 0.05)
    if(x_min == x_max):
      x_min = parameter_boundaries[ix,0]
      x_max = parameter_boundaries[ix,1]

    #x_max_mean = mode_parameters[ix]

    for iy in range(nfit-1, -1, -1):
      y_data = flatchain_posterior_0[:, iy]
      y_min, y_max = anc.compute_limits(y_data, 0.05)
      if(y_min == y_max):
        y_min = parameter_boundaries[iy,0]
        y_max = parameter_boundaries[iy,1]
      
      #y_max_mean = mode_parameters[iy]
      
      if(iy > ix): # correlation plot
        logger.info('%s vs %s' %(parameter_names_emcee[ix], parameter_names_emcee[iy]))

        ax = plt.subplot2grid((nfit+1, nfit), (iy,ix))
        
        #hist2d_counts, xedges, yedges, image2d = ax.hist2d(x_data, y_data, bins=k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], cmap=cm.gray_r, normed=True)
        hist2d_counts, xedges, yedges, image2d = ax.hist2d(x_data, y_data, bins=k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], cmap=cm.gray_r, normed=False)
        
       #new_k = int(k/3)
        new_k = k
        hist2d_counts_2, xedges_2, yedges_2 = np.histogram2d(x_data, y_data, bins=new_k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], normed=True)
        
        x_bins = [0.5*(xedges_2[i]+xedges_2[i+1]) for i in range(0, new_k)]
        y_bins = [0.5*(yedges_2[i]+yedges_2[i+1]) for i in range(0, new_k)]

        #ax.contour(x_bins, y_bins, hist2d_counts_2.T, 3, cmap=cm.gray, linestyle='solid', linewidths=(0.7, 0.7, 0.7))
        nl = 5
        levels = [1.-np.exp(-0.5*ii) for ii in range(0,nl)] # 2D sigmas: 0sigma, 1sigma, 2sigma, 3sigma, ..
        #ax.contour(x_bins, y_bins, hist2d_counts_2.T, levels, cmap=cm.gray, linestyle='solid', linewidths=(0.7, 0.7, 0.7))
        #ax.contour(x_bins, y_bins, hist2d_counts_2.T, levels, cmap=cm.viridis, linestyle='solid', linewidths=1.)
        #ax.contour(x_bins, y_bins, hist2d_counts_2.T, cmap=cm.viridis, linestyle='solid', linewidths=0.7)
        
        #ax.contour(x_bins, y_bins, hist2d_counts_2.T, levels, cmap=cm.viridis, linestyle='solid', linewidths=0.7, norme=True)
        ax.contour(x_bins, y_bins, hist2d_counts_2.T, nl, cmap=cm.viridis, linestyle='solid', linewidths=0.5, norme=True)
        
        
        if(cli.overplot is not None):
          # plot selected overplot sample
          ax.axvline(overp_par[ix], color='blue', ls='--', lw=1.1, alpha=0.7)
          ax.axhline(overp_par[iy], color='blue', ls='--', lw=1.1, alpha=0.7)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        if(iy == nfit-1):
          set_xaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_plot_labels[ix], [xedges[0], xedges[-1], 3])
        if(ix == 0): 
          set_yaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_plot_labels[iy], [yedges[0], yedges[-1], 5])
        
        ax.set_ylim([y_min, y_max])
        ax.set_xlim([x_min, x_max])
        plt.draw()
      
      elif(iy == ix): # distribution plot
        logger.info('%s histogram' %(parameter_names_emcee[ix]))

        ax = plt.subplot2grid((nfit+1, nfit), (ix,ix))
        if (ix == nfit-1):
          hist_orientation='horizontal'
        else:
          hist_orientation='vertical'
          
        
        idx = np.argsort(x_data)
        
        if(not cli.cumulative):
          # HISTOGRAM and pdf
          hist_counts, edges, patces = ax.hist(x_data, bins=k, range=[x_data.min(), x_data.max()], histtype='stepfilled', color='darkgrey', edgecolor='lightgray', align='mid', orientation=hist_orientation, normed=True, stacked=True)
          
          ##x_pdf = scipy_norm.pdf(x_data[idx], loc=x_data.mean(), scale=x_data.std())
          #x_pdf = scipy_norm.pdf(x_data, loc=x_data.mean(), scale=x_data.std())
          #if(ix == nfit-1):
            #ax.plot(x_pdf[idx], x_data[idx], marker='None', color='black', ls='-', lw=1.1, alpha=0.99)
          #else:
            #ax.plot(x_data[idx], x_pdf[idx] , marker='None', color='black', ls='-', lw=1.1, alpha=0.99)
        
        else:
          # CUMULATIVE HISTOGRAM and cdf
          hist_counts, edges, patces = ax.hist(x_data, bins=k, range=[x_data.min(), x_data.max()], histtype='stepfilled', color='darkgrey', edgecolor='lightgray', align='mid', orientation=hist_orientation, normed=True, stacked=True, cumulative=True)
          
          ##x_cdf = scipy_norm.cdf(x_data[idx], loc=x_data.mean(), scale=x_data.std())
          #x_cdf = scipy_norm.cdf(x_data, loc=x_data.mean(), scale=x_data.std())
          #if(ix == nfit-1):
            #ax.plot(x_cdf[idx], x_data[idx], marker='None', color='black', ls='-', lw=1.1, alpha=0.99)
          #else:
            #ax.plot(x_data[idx], x_cdf[idx] , marker='None', color='black', ls='-', lw=1.1, alpha=0.99)
        
        if (ix == nfit-1):
          ax.set_ylim([y_min, y_max])
          if(cli.overplot is not None):
            # plot selected overplot sample
            ax.axhline(overp_par[ix], color='blue', ls='--', lw=1.1, alpha=0.7)
        else:
          ax.set_xlim([x_min, x_max])
          if(cli.overplot is not None):
            # plot selected overplot sample
            ax.axvline(overp_par[ix], color='blue', ls='--', lw=1.1, alpha=0.7)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title(kel_plot_labels[ix], fontsize=label_size)

        plt.draw()


  logger.info('saving plot')
  emcee_fig_file = os.path.join(emcee_plots, 'emcee_triangle.png')
  fig.savefig(emcee_fig_file, bbox_inches='tight', dpi=300)
  logger.info('png done')
  emcee_fig_file = os.path.join(emcee_plots, 'emcee_triangle.pdf')
  fig.savefig(emcee_fig_file, bbox_inches='tight', dpi=300)
  logger.info('pdf done')
  plt.close(fig)

  logger.info('')

  return

if __name__ == "__main__":
  main()
