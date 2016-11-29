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

from scipy.stats import gaussian_kde
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

  # set label and legend names
  kel_plot_labels = anc.keplerian_legend(parameter_names_emcee, cli.m_type)

  chains_T, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee, parameter_boundaries, chains)

  flatchain_posterior_0 = chains_T[:,:,:].reshape((-1, nfit))

  if(cli.boot_id > 0):
    flatchain_posterior_msun = anc.posterior_back_to_msun(m_factor,parameter_names_emcee,flatchain_posterior_0)
    boot_file = anc.save_bootstrap_like(emcee_folder, cli.boot_id, parameter_names_emcee, flatchain_posterior_msun)
    logger.info('saved bootstrap like file: %s' %(boot_file))
    del flatchain_posterior_msun

  # GET MAX LNPROBABILITY AND PARAMETERS -> id 40XX
  max_lnprob, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint = anc.get_maxlnprob_parameters(npost, nruns, lnprobability, chains_T, flatchain_posterior_0)

  # std way: median of the posterior parameter distribution -> id 10XX
  median_parameters, median_perc68, median_confint = anc.get_median_parameters(flatchain_posterior_0)

  # MODE-LIKE PARAMETERS -> id 30XX
  # take the mean of 5 bin centered to the higher bin
  k = np.ceil(2. * flatchain_posterior_0.shape[0]**(1./3.)).astype(int)
  #if(k>11): k=11
  if(k>50): k=50
  mode_bin, mode_parameters, mode_perc68, mode_confint = anc.get_mode_parameters_full(flatchain_posterior_0, k)

  #fig, ax = plt.subplots(nrows = nfit-1, ncols=nfit, figsize=(12,12))
  fig = plt.figure(figsize=(12,12))
  fig.subplots_adjust(hspace=0.05, wspace=0.05)

  for ii in range(0, nfit, 1):
    x_data = flatchain_posterior_0[:, ii]
    #x_med = np.median(x_data)
    #x_16 = np.percentile(x_data, 16)
    #x_84 = np.percentile(x_data, 84)
    x_med = median_parameters[ii]
    x_min, x_max = anc.compute_limits(x_data, 0.05)
    if(x_min == x_max):
      x_min = parameter_boundaries[ii,0]
      x_max = parameter_boundaries[ii,1]

    
    #x_max_mean, x_max_bin = compute_max_mean(x_data, k)
    x_max_mean = mode_parameters[ii]

    for jj in range(nfit-1, -1, -1):
      y_data = flatchain_posterior_0[:, jj]
      #y_med = np.median(y_data)
      y_med = median_parameters[jj]
      y_min, y_max = anc.compute_limits(y_data, 0.05)
      if(y_min == y_max):
        y_min = parameter_boundaries[jj,0]
        y_max = parameter_boundaries[jj,1]
      
      #y_max_mean, y_max_bin = compute_max_mean(y_data, k)
      y_max_mean = mode_parameters[jj]
      
      if(jj > ii): # correlation plot
        logger.info('%s vs %s' %(parameter_names_emcee[ii], parameter_names_emcee[jj]))
        logger.info('(%.4f) %.4f <= X <= %.4f (%.4f)' %(parameter_boundaries[ii,0], x_data.min(), x_data.max(), parameter_boundaries[ii,1]))
        logger.info('(%.4f) %.4f <= Y <= %.4f (%.4f)' %(parameter_boundaries[jj,0], y_data.min(), y_data.max(), parameter_boundaries[jj,1]))

        ax = plt.subplot2grid((nfit+1, nfit), (jj,ii))
        
        #hist2d_counts, xedges, yedges, image2d = ax.hist2d(x_data, y_data, bins=k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], cmap=plt.get_cmap('Greys'), normed=True)
        hist2d_counts, xedges, yedges, image2d = ax.hist2d(x_data, y_data, bins=k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], cmap=cm.gray_r, normed=True)
        
        #hist2d_counts, xedges, yedges, image2d = ax.hist2d(x_data, y_data, bins=k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], cmap=plt.get_cmap('Greys'))
              
        new_k = int(k/3)
        hist2d_counts_2, xedges_2, yedges_2 = np.histogram2d(x_data, y_data, bins=new_k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], normed=True)
        
        #hist2d_counts_2, xedges_2, yedges_2 = np.histogram2d(x_data, y_data, bins=new_k, range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]])
        
        ## test gaussian_kde scatter plot !!VERY VERY SLOw!!
        #print 'applying gaussian kde ...',
        #sys.stdout.flush()
        #xy = np.vstack([x_data,y_data])
        #z = gaussian_kde(xy)(xy)
        #idx = np.argsort(z)
        #ax.scatter(x_data[idx], y_data[idx], c=z[idx], s=2, edgecolor='', cmap=plt.get_cmap('Greys'))
        #print 'done'
        #sys.stdout.flush()
        #ax.scatter(x_data, y_data, color='black', s=2, edgecolor='', alpha=0.33) # SLOW!!
        
        x_bins = [0.5*(xedges_2[i]+xedges_2[i+1]) for i in range(0, new_k)]
        y_bins = [0.5*(yedges_2[i]+yedges_2[i+1]) for i in range(0, new_k)]
        #ax.contour(x_bins, y_bins, hist2d_counts_2.T, 3, colors=('forestgreen', 'royalblue', 'red'), linestyle='solid', linewidths=(0.5, 0.5, 0.5))
        #ax.contour(x_bins, y_bins, hist2d_counts_2.T, 3, colors=('black', 'forestgreen', 'lightgray'), linestyle='solid', linewidths=(0.5, 0.5, 0.5))
        ax.contour(x_bins, y_bins, hist2d_counts_2.T, 3, cmap=cm.gray, linestyle='solid', linewidths=(0.7, 0.7, 0.7))
        
        # plot mean_mode
        ax.axvline(x_max_mean, color='red', ls='-', lw=0.9, alpha=0.7)
        ax.axhline(y_max_mean, color='red', ls='-', lw=0.9, alpha=0.7)
        # plot median
        ax.axvline(x_med, color='blue', ls='--', lw=1.1, alpha=0.7)
        ax.axhline(y_med, color='blue', ls='--', lw=1.1, alpha=0.7)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        if(jj == nfit-1):
          set_xaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_plot_labels[ii], [xedges[0], xedges[-1], 3])
        if(ii == 0): 
          set_yaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_plot_labels[jj], [yedges[0], yedges[-1], 5])
        
        ax.set_ylim([y_min, y_max])
        ax.set_xlim([x_min, x_max])
        plt.draw()
      
      elif(jj == ii): # distribution plot
        logger.info('%s histogram' %(parameter_names_emcee[ii]))
        logger.info('(%.4f) %.4f <= X <= %.4f (%.4f)' %(parameter_boundaries[ii,0], x_data.min(), x_data.max(), parameter_boundaries[ii,1]))

        ax = plt.subplot2grid((nfit+1, nfit), (ii,ii))
        if (ii == nfit-1):
          hist_orientation='horizontal'
        else:
          hist_orientation='vertical'
          
        
        idx = np.argsort(x_data)
        
        if(not cli.cumulative):
          # HISTOGRAM and pdf
          hist_counts, edges, patces = ax.hist(x_data, bins=k, range=[x_data.min(), x_data.max()], histtype='stepfilled', color='darkgrey', edgecolor='lightgray', align='mid', orientation=hist_orientation, normed=True, stacked=True)
          
          #x_pdf = scipy_norm.pdf(x_data[idx], loc=x_data.mean(), scale=x_data.std())
          x_pdf = scipy_norm.pdf(x_data, loc=x_data.mean(), scale=x_data.std())
          if(ii == nfit-1):
            ax.plot(x_pdf[idx], x_data[idx], marker='None', color='black', ls='-', lw=1.1, alpha=0.99)
          else:
            ax.plot(x_data[idx], x_pdf[idx] , marker='None', color='black', ls='-', lw=1.1, alpha=0.99)
        
        else:
          # CUMULATIVE HISTOGRAM and cdf
          hist_counts, edges, patces = ax.hist(x_data, bins=k, range=[x_data.min(), x_data.max()], histtype='stepfilled', color='darkgrey', edgecolor='lightgray', align='mid', orientation=hist_orientation, normed=True, stacked=True, cumulative=True)
          
          
          #x_cdf = scipy_norm.cdf(x_data[idx], loc=x_data.mean(), scale=x_data.std())
          x_cdf = scipy_norm.cdf(x_data, loc=x_data.mean(), scale=x_data.std())
          if(ii == nfit-1):
            ax.plot(x_cdf[idx], x_data[idx], marker='None', color='black', ls='-', lw=1.1, alpha=0.99)
          else:
            ax.plot(x_data[idx], x_cdf[idx] , marker='None', color='black', ls='-', lw=1.1, alpha=0.99)
        
        
        if (ii == nfit-1):
          ax.set_ylim([y_min, y_max])
          # plot mean_mode
          ax.axhline(x_max_mean, color='red', ls='-', lw=0.9, alpha=0.7)
          # plot median
          ax.axhline(x_med, color='blue', ls='--', lw=1.1, alpha=0.7)
        else:
          ax.set_xlim([x_min, x_max])
          # plot mean_mode
          ax.axvline(x_max_mean, color='red', ls='-', lw=0.9, alpha=0.7)
          # plot median
          ax.axvline(x_med, color='blue', ls='--', lw=1.1, alpha=0.7)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title(kel_plot_labels[ii], fontsize=label_size)

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

  anc.print_parameters_logger(logger, parameter_names_emcee, median_parameters, median_perc68, median_confint, 'median')
  anc.print_parameters_logger(logger, parameter_names_emcee, mode_parameters, mode_perc68, mode_confint, 'mode')

  return

if __name__ == "__main__":
  main()
