#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import sys
import argparse
import time
import os
import numpy as np # array
import h5py
import constants as cst # local constants module
import ancillary as anc
import corner

import matplotlib.pyplot as plt
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

def read_posterior(folder_path):
  
  posterior_file = os.path.join(folder_path, 'posterior.hdf5')
  p_h5f = h5py.File(posterior_file, 'r')
  posterior = p_h5f['posterior'][...]
  parameter_names = p_h5f['parameter_names'][...]
  p_h5f.close()
  nposterior, nfit = np.shape(posterior)
  
  return parameter_names, posterior, nposterior, nfit

def main():

  cli = anc.get_args()
  parameter_names, posterior, nposterior, nfit = read_posterior(cli.full_path)
  # set label and legend names
  Kep_labels = anc.keplerian_legend(parameter_names, cli.m_type)
  plot_labels = ["%s" %(Kep_labels[ii]) for ii in range(0, nfit)]
  #print Kep_labels
  #print plot_labels
  plot_folder = os.path.join(cli.full_path, 'plots')
  if (not os.path.isdir(plot_folder)):
      os.makedirs(plot_folder)
  corner_file = os.path.join(plot_folder, 'corner_plot_fitted')
  
  #k = np.ceil(2. * nposterior**(1./3.)).astype(int)
  #if(k>20): k=20
  k = anc.get_bins(posterior, rule='doane')
  
  nl = 4
  levels = [1.-np.exp(-0.5*ii) for ii in range(0,nl)] # 2D sigmas: 0sigma, 1sigma, 2sigma, 3sigma, ..
  
  #fig = corner.corner(posterior, bins=k,
                      #labels = parameter_names,
                      #quantiles = [0.16, 0.5, 0.84],
                      #show_titles = True,
                      #title_kwargs = {'fontsize': 10},
                      #levels = (1.-np.exp(-0.5), 1.-np.exp(-1.), 1.-np.exp(-1.5)) # 2D sigmas: 1sigma, 2sigma, 3sigma
                      #)
  
  fig = corner.corner(posterior,
                      labels = plot_labels,
                      quantiles = [0.16, 0.5, 0.84],
                      show_titles = True,
                      levels = levels,
                      )

  
  fig.savefig('%s.png' %(corner_file), dpi=300)
  print('Saved plot %s.png' %(corner_file))
  fig.savefig('%s.pdf' %(corner_file), dpi=150)
  print('Saved plot %s.pdf' %(corner_file))
  #fig.savefig('%s.png' %(corner_file), bbox_inches='tight', dpi=300)
  #fig.savefig('%s.pdf' %(corner_file), bbox_inches='tight', dpi=150)
  
  return


if __name__ == "__main__":
  main()
  
  
