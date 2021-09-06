#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import argparse
import os
import sys
import numpy as np # array
import h5py
import ancillary as anc
from scipy.stats import norm as scipy_norm

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# matplotlib rc params
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.figsize'] = [6, 6]
plt.rcParams["figure.facecolor"] = 'white'
plt.rcParams["savefig.facecolor"] = 'white'
plt.rcParams["figure.dpi"]  = 200
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"]   = 14
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter

def set_xaxis(ax, label_size, label_separation, label_pad, ticklabel_size, kel_label, ticks_formatter, tick_fmt='%.4f'):
  ax.get_xaxis().set_visible(True)
  ax.xaxis.set_tick_params(labelsize=ticklabel_size)
  ax.xaxis.set_label_coords(0.5, label_separation)
  #ax.ticklabel_format(style='plain', axis='both', useOffset=False)
  plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
  #ax.xaxis.labelpad = label_pad
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
  #ax.yaxis.labelpad = label_pad
  ax.set_ylabel(kel_label, fontsize=label_size)
  tick_step = (ticks_formatter[1] - ticks_formatter[0]) / ticks_formatter[2]
  ax.yaxis.set_ticks(np.arange(ticks_formatter[0], ticks_formatter[1], tick_step))
  tick_formatter = FormatStrFormatter(tick_fmt)
  ax.yaxis.set_major_formatter(tick_formatter)
  return
    

## PLOT DERIVED POSTERIOR CORRELATION PLOT

def main():
  cli = anc.get_args()
  # read derived posterior file
  derived_file = os.path.join(cli.full_path, 'derived_posterior.hdf5')
  h5f = h5py.File(derived_file, 'r')
  # derived_names = np.array(h5f['derived_names'], dtype='S10')
  derived_names = anc.decode_list(h5f['derived_names'][...])
  derived_posterior_in = np.array(h5f['derived_posterior'], dtype=np.float64)
  h5f.close()

  n_der = np.shape(derived_names)[0]
  # n_flatchain = np.shape(derived_posterior_in)[0]

  derived_posterior = anc.derived_posterior_check(derived_names, derived_posterior_in)

      
  label_separation=-0.90 # if uses this, comment ax.xyaxis.labelpad = label_pad
  label_pad = 12 # it uses this, comment ax.xyaxis.set_label_coords()...
  label_size = 8
  ticklabel_size = 4

  if(n_der > 2):
    #label_separation = -0.1 - ( 0.075 * (n_der-2) )
    label_separation = -0.15 - ( 0.125 * (n_der-2) )
  #else:
    #label_separation = -0.15

  #label_size = label_size - 1 * int(n_der / 10.)
  #label_size = label_size - 1 * int(n_der / 5.)
  label_size = label_size - 1 * int(n_der / 2.5)

  labels_list = anc.derived_labels(derived_names, cli.m_type)

  k = anc.get_auto_bins(derived_posterior)
  
  if(cli.overplot is not None):
    ## OPEN summary_parameters.hdf5 FILE
    s_h5f = h5py.File(os.path.join(cli.full_path, 'summary_parameters.hdf5'), 'r')
    # take only the selected sample
    s_overplot = '%04d' %(cli.overplot)
    #overp_der = s_h5f['parameters/%s/derived/parameters' %(s_overplot)][...]
    read_der = s_h5f['parameters/{0:s}/derived/parameters'.format(s_overplot)][...]
    s_h5f.close()
  
    overp_der = anc.derived_parameters_check(derived_names, read_der, derived_posterior)
  
  #fig = plt.figure(figsize=(12,12))
  fig = plt.figure(figsize=(6,6))
  fig.subplots_adjust(hspace=0.05, wspace=0.05)
  
  for ix in range(0, n_der):
    x_data = derived_posterior[:,ix]
    x_min, x_max = anc.compute_limits(x_data, 0.05)
    if(x_min == x_max):
      x_min = x_min - 1.
      x_max = x_max + 1.
   
    for iy in range(0, n_der):
      y_data = derived_posterior[:,iy]
      y_min, y_max = anc.compute_limits(y_data, 0.05)
      if(y_min == y_max):
        y_min = y_min - 1.
        y_max = y_max + 1.
        
      if(iy > ix): # correlation plot
        anc.print_both('correlation %s vs %s' %(derived_names[ix], derived_names[iy]) )
        ax = plt.subplot2grid((n_der+1, n_der), (iy,ix))
        
        hist2d_counts, xedges, yedges, image2d = ax.hist2d(\
          x_data, y_data, bins=k,
          range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]], 
          cmap=cm.gray_r,
          density=False
          )
        
        #new_k = int(k/3)
        new_k = k
        hist2d_counts_2, xedges_2, yedges_2 = np.histogram2d(\
          x_data, y_data, bins=new_k,
          range=[[x_data.min(), x_data.max()],[y_data.min(), y_data.max()]],
          density=True
          )
        
        x_bins = [0.5*(xedges_2[i]+xedges_2[i+1]) for i in range(0, new_k)]
        y_bins = [0.5*(yedges_2[i]+yedges_2[i+1]) for i in range(0, new_k)]
        
        nl = 5
        levels = [1.-np.exp(-0.5*ii) for ii in range(0,nl)] # 2D sigmas: 0sigma, 1sigma, 2sigma, 3sigma, ..
        ax.contour(x_bins, y_bins, hist2d_counts_2.T, 
                   nl, cmap=cm.viridis,
                   linestyles='solid', linewidths=0.5,
                   )
        
        if(cli.overplot is not None):
          # plot selected overplot sample
          # check angle and plot %360 and %-360...
          if('w' in derived_names[ix] or 
             'lN' in derived_names[ix] or 
             'mA' in derived_names[ix]):
            ax.axvline(overp_der[ix]%360., color='C0', ls='--', lw=1.1, alpha=0.7)
            ax.axvline(overp_der[ix]%-360., color='C0', ls='--', lw=1.1, alpha=0.7)
          else:
            ax.axvline(overp_der[ix], color='C0', ls='--', lw=1.1, alpha=0.7)
          if('w' in derived_names[iy] or 
             'lN' in derived_names[iy] or 
             'mA' in derived_names[iy]):
            ax.axhline(overp_der[iy]%360., color='C0', ls='--', lw=1.1, alpha=0.7)
            ax.axhline(overp_der[iy]%-360., color='C0', ls='--', lw=1.1, alpha=0.7)
          else:
            ax.axhline(overp_der[iy], color='C0', ls='--', lw=1.1, alpha=0.7)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        if(iy == n_der-1):
          set_xaxis(ax, label_size, label_separation, label_pad, ticklabel_size, labels_list[ix], [xedges[0], xedges[-1], 4])
        if(ix == 0): 
          set_yaxis(ax, label_size, label_separation, label_pad, ticklabel_size, labels_list[iy], [yedges[0], yedges[-1], 5])
        
        ax.set_ylim([y_min, y_max])
        ax.set_xlim([x_min, x_max])
        plt.draw()
  
      elif(iy == ix): # distribution plot
        anc.print_both('%s histogram' %(derived_names[ix]))
        ax = plt.subplot2grid((n_der+1, n_der), (ix,ix))
        if (ix == n_der-1):
          hist_orientation='horizontal'
        else:
          hist_orientation='vertical'
        
        idx = np.argsort(x_data)
        
        if(not cli.cumulative):
          # HISTOGRAM
          hist_counts, edges, patces = ax.hist(x_data, bins=k,
                                               range=[x_data.min(), x_data.max()], 
                                               histtype='stepfilled', 
                                               color='darkgrey', 
                                               #edgecolor='lightgray',
                                               edgecolor='None',
                                               align='mid', 
                                               orientation=hist_orientation, 
                                               density=True,
                                               stacked=True
                                               )
          
        else:
          # CUMULATIVE HISTOGRAM
          hist_counts, edges, patces = ax.hist(x_data, bins=k,
                                               range=[x_data.min(), x_data.max()],
                                               histtype='stepfilled', 
                                               color='darkgrey', 
                                               #edgecolor='lightgray',
                                               edgecolor='None',
                                               align='mid', 
                                               orientation=hist_orientation, 
                                               density=True,
                                               stacked=True,
                                               cumulative=True
                                               )
        
        #print parameter_names_emcee[ix], overp_der[ix]
        if (ix == n_der-1):
          if(cli.overplot is not None):
            # check angle and plot %360 and %-360...
            if('w' in derived_names[ix] or 
               'lN' in derived_names[ix] or 
               'mA' in derived_names[ix]):
              ax.axhline(overp_der[ix]%360., color='C0', ls='--', lw=1.1, alpha=0.7)
              ax.axhline(overp_der[ix]%-360., color='C0', ls='--', lw=1.1, alpha=0.7)
            else:
              # plot selected overplot sample
              ax.axhline(overp_der[ix], color='C0', ls='--', lw=1.1, alpha=0.7)
          ax.set_ylim([y_min, y_max])
        else:
          if(cli.overplot is not None):
            if('w' in derived_names[ix] or 
               'lN' in derived_names[ix] or 
               'mA' in derived_names[ix]):
              ax.axvline(overp_der[ix]%360., color='C0', ls='--', lw=1.1, alpha=0.7)
              ax.axvline(overp_der[ix]%-360., color='C0', ls='--', lw=1.1, alpha=0.7)
            else:
              # plot selected overplot sample
              ax.axvline(overp_der[ix], color='C0', ls='--', lw=1.1, alpha=0.7)
          ax.set_xlim([x_min, x_max])
        if(cli.overplot is not None):
          print(derived_names[ix], ' overplot val = ', overp_der[ix], ' min = ', x_data.min(), ' max = ', x_data.max())
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title(labels_list[ix], fontsize=label_size)
        plt.draw()
  
  plot_folder = os.path.join(cli.full_path, 'plots')
  if (not os.path.isdir(plot_folder)):
      os.makedirs(plot_folder)
  correlation_file = os.path.join(plot_folder, 'derived_triangle.png')
  fig.savefig(correlation_file, bbox_inches='tight', dpi=300)
  anc.print_both('png done')
  correlation_file = os.path.join(plot_folder, 'derived_triangle.pdf')
  fig.savefig(correlation_file, bbox_inches='tight', dpi=96)    
  anc.print_both('pdf done')
  plt.close(fig)
      
  return

# 
# Here run the proper function
#
if __name__ == "__main__":
  main()

  
  
