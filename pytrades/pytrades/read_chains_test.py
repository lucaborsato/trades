#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import os
import sys
import numpy as np # array
import h5py


def main():

  hdf5_file = '/data2/borsato/TRADESsimulations/Kepler-19/EMCEE/2016-12-01_K19_emcee_001/2016-12-01_emcee_001/emcee_summary.hdf5'
  
  f_h5 = h5py.File(hdf5_file, 'r')
  full_chains = f_h5['chains'][...]
  f_h5.close()

  nwalkers, nsteps, nfit = np.shape(full_chains)

  nburnin = 125000
  thin = 200
  nuse = nsteps-nburnin
  npost = nwalkers*(nuse)
  sel_thin = np.arange(0, nuse+thin, thin).astype(int)
  n_thin = np.shape(sel_thin)[0]
  if(sel_thin[-1] >= nuse): sel_thin[-1] = nuse-1
  npost_thin = nwalkers * n_thin
  
  print 'nfit         = ', nfit
  print 'nwalkers     = ', nwalkers
  print 'nsteps       = ', nsteps
  print 'nburnin      = ', nburnin
  print 'nuse         = ', nuse
  print 'npost        = ', npost
  print 'thin         = ', thin
  print 'n_thin       = ', n_thin
  print 'npost_thin   = ', npost_thin
  print 'sel_thin[0]  = ',sel_thin[0]
  print 'sel_thin[-1] = ',sel_thin[-1]
  
  

  # full chains of Pb
  Pb_full = full_chains[:,:,1].reshape((nwalkers*nsteps))
  # select only steps after the burnin for each walker
  Pb_chains_posterior = full_chains[:,nburnin:,1]
  # create the flat posterior of Pb
  Pb_posterior = Pb_chains_posterior[:,:].reshape((npost))
  # select only the posterior thinned for each walker
  Pb_chains_thinned = Pb_chains_posterior[:,sel_thin]
  # create the flat thinned posterior of Pb
  Pb_thinned = Pb_chains_thinned[:,:].reshape(npost_thin)
  
  print np.percentile(Pb_thinned, 16.), np.median(Pb_thinned), np.percentile(Pb_thinned, 84.)
  
  np.savetxt('/data2/borsato/TRADESsimulations/Kepler-19/EMCEE/2016-12-01_K19_emcee_001/2016-12-01_emcee_001/Pb_posterior_thinned.dat', Pb_thinned, fmt='%23.16e', header='P_b_posterior_thinned')
  print 'SAVE Pb_thinned'
  
  np.savetxt('/data2/borsato/TRADESsimulations/Kepler-19/EMCEE/2016-12-01_K19_emcee_001/2016-12-01_emcee_001/Pb_sorted_numpy.dat', np.sort(Pb_thinned), fmt='%23.16e')
  print 'SAVE Pb_thinned sorted'
    
  np.savetxt('/data2/borsato/TRADESsimulations/Kepler-19/EMCEE/2016-12-01_K19_emcee_001/2016-12-01_emcee_001/Pb_posterior.dat', Pb_posterior, fmt='%23.16e', header='P_b_posterior')
  print 'SAVE Pb_posterior'
  
  np.savetxt('/data2/borsato/TRADESsimulations/Kepler-19/EMCEE/2016-12-01_K19_emcee_001/2016-12-01_emcee_001/Pb_full.dat', Pb_full, fmt='%23.16e', header='P_b_full')
  print 'SAVE Pb_full'


  return


if __name__ == "__main__":
  main()

