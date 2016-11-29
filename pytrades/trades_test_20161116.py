#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
#import argparse
import os
import numpy as np # array
from pytrades_lib import pytrades
#import emcee
#import h5py
#import sys
#import time
#import shutil
from constants import Mjups
import ancillary as anc
import glob

def main():
  # integration parameters
  t_start = np.float64(5999.)
  t_epoch = np.float64(6050.)
  t_int = np.float64(201.)
  step_in = np.float64(1.e-3)

  print 'set integration parameters'

  # number of bodies in the system
  n_body = 2

  print 'set n_body'

  # read rv data
  rv_path = os.path.abspath('/home/borsato/Dropbox/Transfer/LMa_TestTTVFast/2016-11-08_data/rv_data_1p.dat')
  rv_data = np.genfromtxt(rv_path, usecols=(0,1,2)) # t_rv, rv_obs, erv_obs
  n_rv = rv_data.shape[0]

  t_rv = rv_data[:,0]
  rv_obs  = rv_data[:,1]
  erv_obs = rv_data[:,2]

  print 'read and set rv data'

  # read t0 data
  t0b_path = os.path.abspath('/home/borsato/Dropbox/Transfer/LMa_TestTTVFast/2016-11-08_data/obsT0.dat')
  t0b_data = np.genfromtxt(t0b_path)
  n_t0b = t0b_data.shape[0]

  n_t0 = np.zeros((n_body)).astype(int)
  n_t0[1] = n_t0b

  n_max_t0 = np.max(n_t0)

  t0_num = np.zeros((n_max_t0,n_body)).astype(int)
  t0_obs = np.zeros((n_max_t0,n_body))
  et0_obs = np.zeros((n_max_t0,n_body))

  t0_num[0:n_t0[1],1] = t0b_data[:,0].astype(int)
  t0_obs[0:n_t0[1],1] = t0b_data[:,1]
  et0_obs[0:n_t0[1],1] = t0b_data[:,2]

  print 'read and set T0 data'
  print 'n_t0 = ', n_t0
  print 'n_max_t0 = ',n_max_t0

  pytrades.args_init(t_start,t_epoch,t_int,n_body,n_t0,t0_num,t0_obs,et0_obs)

  print 'initialised args for trades'

  transit_flag = np.array([False, True])

  M_msun   = np.array([1.0, 7.9792038535161035E-04], dtype=np.float64) # Mb=265.66443948Me=0.8357008336179664Mj
  R_rsun   = np.array([1.0, 1.0271839080459770E-01], dtype=np.float64) # Rb=11.208916409849234Re=1Rj
  P_day    = np.array([0, 33.52], dtype=np.float64)
  ecc      = np.array([0, 0.35598141891760265], dtype=np.float64)
  argp_deg = np.array([0, 334.10664447355475], dtype=np.float64)
  mA_deg   = np.array([0, 249.76203184284145], dtype=np.float64)
  inc_deg  = np.array([0, 90.], dtype=np.float64)
  lN_deg   = np.array([0, 180.], dtype=np.float64)

  print 'set orbital elements'

  rv_sim, t0_sim = pytrades.kelements_to_data(t_start,t_epoch,step_in,t_int,
                     M_msun,R_rsun,P_day,ecc,argp_deg,mA_deg,inc_deg,lN_deg,
                     t_rv,transit_flag,n_t0,t0_num)

  print
  print '## RV'
  print '# time_rv rv_obs rv_sim'
  for irv in range(0, n_rv):
    print '%18.6f %12.6f %12.6f' %(t_rv[irv], rv_obs[irv], rv_sim[irv])
  print 

  print '## T0'
  for inb in range(1,n_body):
    print '# BODY %d (planet %d)' %(inb+1, inb)
    print '# t0_num t0_obs t0_sim'
    for itt in range(0,n_t0[inb]):
      print '%5d %20.6f %16.6f' %(t0_num[itt,inb], t0_obs[itt,inb], t0_sim[itt,inb])

  return

if __name__ == "__main__":
  main()
  
