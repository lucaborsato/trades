#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
import pytrades_lib
#import emcee
import h5py
import sys
import time
import shutil
#from emcee.utils import MPIPool
from constants import Mjups, Msjup, Mears, Msear
from ancillary import *
#from pyde.de import DiffEvol as devol
import pymultinest
import logging
import warnings

#
def get_args():
  parser = argparse.ArgumentParser(description='TRADES+pyMultiNest')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--sub-folder', '--sb', action='store', dest='sub_folder', default='emcee_run', help='Sub-folder name, without full path. Default = emcee_run')
  
  # COMPUTE OR NOT CONSTANT TO ADD TO THE LGLIKELIHOOD: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 )
  parser.add_argument('-l', '--ln-err', '--ln-err-const', action='store', dest='ln_flag', default=True, help='Computes or not constant to add to the lglikelihood: SUM( ln 2pi * sigma_obs^2 ). Default = True. Set to False if not use this value.')

  parser.add_argument('--npop', action='store', dest='n_pop', default=None, help='Number of populations/configurations to use with pyMultiNest. Set to values less than nfit is not allowed. Default = nfit * 10.')

  parser.add_argument('--seed', action='store', dest='seed_value', default=-1, help='Seed for random numbers. Set it to -1 for automatic determination. Default = -1')

  parser.add_argument('--resume', action='store', dest='resume_flag', default=False, help='Flag for resume or not. Default = False')

  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  
  cli.ln_flag = set_bool_argument(cli.ln_flag)
  cli.resume_flag =  set_bool_argument(cli.resume_flag)
    
  return cli

#
# INITIALISE FOLDER AND LOG FILE
#
def init_folder(working_path, sub_folder):
  working_folder = os.path.join(working_path, sub_folder)
  if (not os.path.isdir(working_folder)):
      os.makedirs(working_folder)
  run_log = os.path.join(working_folder, "trades_run.log")
  of_run = open(run_log, 'w')
  of_run.write("# pyTRADES LOG FILE\n")
  of_run.write("# working_path = %s\n" %(working_path))
  of_run.write("# working_folder = %s\n" %(working_folder))
  of_run.write("# run_log = %s\n" %(run_log))
  
  return working_folder, run_log, of_run

#def compute_ln_err_const(ndata, dof, e_RVo, e_T0o, ln_flag):
  #if (ln_flag):
    #eRV = e_RVo[e_RVo > 0.]
    #eT0 = e_T0o[e_T0o > 0.]
    
    #ln_e_RVo = np.abs(np.sum(np.log(eRV*eRV)))
    #ln_e_T0o = np.abs(np.sum(np.log(eT0*eT0)))

    ##ln_err_const = -(0.5/dof) * (ndata*np.log(2.*np.pi) + ln_e_RVo + ln_e_T0o)
    ##ln_err_const = np.log(2.*np.pi) + ln_e_RVo + ln_e_T0o
    #ln_err_const = - 0.5 * dof * np.log(2.*np.pi) - 0.5 * ( ln_e_RVo + ln_e_T0o)
  #else:
    #ln_err_const = 0.
  #return ln_err_const


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

def convert_to_int(val_in):
  try:
    val_out = int(val_in)
  except:
    val_out = -1
  return val_out
    
# some global constants
Mear2jup = Mears*Msjup
Mjup2ear = Mjups*Msear

def main():
  # MAIN -- TRADES + pyMultiNest
  # ---
  # initialize logger
  logger = logging.getLogger("Main_log")
  logger.setLevel(logging.DEBUG)
  formatter = logging.Formatter("%(asctime)s - %(message)s")

  # READ COMMAND LINE ARGUMENTS
  cli = get_args()

  # STARTING TIME
  start = time.time()

  # RENAME 
  working_path = cli.full_path
  #nthreads=cli.nthreads

  log_file = os.path.join(working_path, '%s_log.txt' %(os.path.dirname(cli.sub_folder)))

  flog = logging.FileHandler(log_file, 'w')
  flog.setLevel(logging.DEBUG)
  flog.setFormatter(formatter)
  logger.addHandler(flog)
  # log screen
  slog = logging.StreamHandler()
  slog.setLevel(logging.DEBUG)
  slog.setFormatter(formatter)
  logger.addHandler(slog)


  fitting_priors, fitting_priors_type = read_priors(os.path.join(working_path, 'fitting_priors.dat'))
  derived_priors, derived_priors_type = read_priors(os.path.join(working_path, 'derived_priors.dat'))
  n_der_priors = len(derived_priors)


  # INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
  pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder, 1)

  # RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
  fitting_parameters = pytrades_lib.pytrades.fitting_parameters # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
  parameters_minmax = pytrades_lib.pytrades.parameters_minmax # PARAMETER BOUNDARIES
  delta_parameters = np.abs(parameters_minmax[:,1] - parameters_minmax[:,0]) # DELTA BETWEEN MAX AND MIN OF BOUNDARIES

  n_bodies = pytrades_lib.pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
  n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
  ndata = pytrades_lib.pytrades.ndata # TOTAL NUMBER OF DATA AVAILABLE
  npar  = pytrades_lib.pytrades.npar # NUMBER OF TOTAL PARAMATERS ~n_planets X 6
  nfit  = pytrades_lib.pytrades.nfit # NUMBER OF PARAMETERS TO FIT
  nfree  = pytrades_lib.pytrades.nfree # NUMBER OF FREE PARAMETERS (ie nrvset)
  dof   = pytrades_lib.pytrades.dof # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
  global inv_dof
  #inv_dof = np.float64(1.0 / dof)
  inv_dof = pytrades_lib.pytrades.inv_dof

  str_len = pytrades_lib.pytrades.str_len
  temp_names = pytrades_lib.pytrades.get_parameter_names(nfit,str_len)
  trades_names = anc.convert_fortran_charray2python_strararray(temp_names)
  parameter_names = anc.trades_names_to_emcee(trades_names)

  # RADIAL VELOCITIES SET
  n_rv = pytrades_lib.pytrades.nrv
  n_set_rv = pytrades_lib.pytrades.nrvset

  # TRANSITS SET
  n_t0 = pytrades_lib.pytrades.nt0
  n_t0_sum = pytrades_lib.pytrades.ntts
  n_set_t0 = 0
  for i in range(0, n_bodies-1):
    if (n_t0[i] > 0): n_set_t0 += 1


  # compute global constant for the loglhd
  global ln_err_const

  #try:
    #e_RVo = np.asarray(pytrades_lib.pytrades.ervobs[:], dtype=np.float64) # fortran variable RV in python will be rv!!!
  #except:
    #e_RVo = np.asarray([0.], dtype=np.float64)
  #try:
    #e_T0o = np.asarray(pytrades_lib.pytrades.et0obs[:,:], dtype=np.float64).reshape((-1))
  #except:
    #e_T0o = np.asarray([0.], dtype=np.float64)
  #ln_err_const = anc.compute_ln_err_const(ndata, dof, e_RVo, e_T0o, cli.ln_flag)
  ln_err_const = pytrades_lib.pytrades.ln_err_const
  
  # READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
  #reshaped_names = pytrades_lib.pytrades.parameter_names.reshape((10,nfit), order='F').T
  #parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]

  # INITIALISE SCRIPT FOLDER/LOG FILE
  working_folder, run_log, of_run = init_folder(working_path, cli.sub_folder)

  logger.info('')
  logger.info('==================== ')
  logger.info('pyTRADES-pyMultiNest')
  logger.info('==================== ')
  logger.info('')
  logger.info('WORKING PATH = %s' %(working_path))
  logger.info('dof = ndata(%d) - nfit(%d) = %d' %(ndata, nfit, dof))
  logger.info('Total N_RV = %d for %d set(s)' %(n_rv, n_set_rv))
  logger.info('Total N_T0 = %d for %d out of %d planet(s)' %(n_t0_sum, n_set_t0, n_planets))
  logger.info('%s = %.7f' %('log constant error = ', ln_err_const))

  #of_run.write(' dof = ndata(%d) - nfit(%d) = %d\n' %(ndata, nfit, dof))
  #of_run.write(' Total N_RV = %d for %d set(s)\n' %(n_rv, n_set_rv))
  #of_run.write(' Total N_T0 = %d for %d out of %d planet(s)\n' %(n_t0_sum, n_set_t0, n_planets))
  #of_run.write(' %s = %.7f\n' %('log constant error = ', ln_err_const))

  # save parameter_names and boundaries to be read by a script
  trades_hdf5 = h5py.File(os.path.join(working_folder, 'system_summary.hdf5'), 'w')
  trades_hdf5.create_dataset('parameter_names', data=parameter_names, dtype='S10')
  trades_hdf5.create_dataset('parameters_minmax', data=parameters_minmax, dtype=np.float64)
  trades_hdf5.create_dataset('ln_err_const', data=np.asarray([ln_err_const], dtype=np.float64), dtype=np.float64)
  trades_hdf5.close()

  # MULTINEST HERE
  #output_mnest = os.path.join(working_folder, 'trades_mnest_')
  #output_mnest = os.path.join(cli.sub_folder, 'trades_mnest_')
  output_mnest = 'trades_mnest_'
  os.chdir(working_folder)
  log_zero_value = -0.5*1.e8
  seed_value = convert_to_int(cli.seed_value)
  #seed_value = 392919
  #n_pop = nfit+int(nfit*0.5)
  n_pop = convert_to_int(cli.n_pop)
  if( n_pop < nfit): n_pop = nfit*10
  n_update = 5 # by argument
  resume_flag = cli.resume_flag
  multi_node_flag = True
  mpi_onoff = False
  logger.info('pyMultiNest parameters:')
  logger.info('folder = %s' %(output_mnest))
  logger.info('seed_value = %d , n_pop = %d , n_update = %d' %(seed_value, n_pop, n_update))
  logger.info('resume_flag = %r , multi_node_flag = %r, mpi_onoff = %r' %(resume_flag, multi_node_flag, mpi_onoff))
  #of_run.write(' pyMultiNest parameters:\n')
  #of_run.write(' folder = %s\n seed_value = %d , n_pop = %d , n_update = %d\n resume_flag = %r , multi_node_flag = %r, mpi_onoff = %r\n' %(output_mnest, seed_value, n_pop, n_update, resume_flag, multi_node_flag, mpi_onoff))

  #
  # RESCALE PARAMETERS FUNCTION NEEDED BY LNLIKE
  #
  def trades_rescale(fitting_parameters, ndim, nparams):
    for i in range(0,ndim):
      fitting_parameters[i] = parameters_minmax[i,0] + fitting_parameters[i]*delta_parameters[i]
    return fitting_parameters

  # LNPRIOR TO BE ADDED TO LOGLHD
  # it can use all the variables defined before this point!
  def lnprior(fitting_parameters, ndim):
    lnprior_value = 0.
    i_der = 0
    for i in range(0, ndim):
      #print i,parameter_names[i], fitting_priors_type[i]
      ln_temp = 0.
      
      # calculate the LogLikelihood<->prior of fitting parameter
      if(fitting_priors_type[i][1].lower() == 'g'):
        ln_temp = -0.5*(((fitting_parameters[i]-fitting_priors[i][0])/fitting_priors[i][1])**2)
        lnprior_value = lnprior_value + ln_temp
        #print '%18.12f %18.12f (%18.12f) => ln = %18.12f' %(fitting_parameters[i], fitting_priors[i][0], fitting_priors[i][1], ln_temp)
      
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
          #print derived_priors_type[i_der]
          #print '%18.12f %18.12f (%18.12f) => ln = %18.12f' %(ecc, derived_priors[i_der][0], derived_priors[i_der][1], ln_temp)
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
          #print derived_priors_type[i_der+1]
          #print '%18.12f (argp[%18.12f]+mA[%18.12f]) %18.12f (%18.12f) => ln = %18.12f' %(phi, argp, fitting_parameters[i], derived_priors[i_der+1][0], derived_priors[i_der+1][1], ln_temp)
          i_der = i_der + 2
      
    return lnprior_value

  # LNLIKEHOOD FUNCTION NEEDED BY MULTINEST
  def lnlike(fitting_parameters, ndim, nparams):
    loglhd = 0.
    check = 1
    trades_parameters = np.asarray([fitting_parameters[i] for i in range(0,ndim)], dtype=np.float64)
    loglhd, check = pytrades_lib.pytrades.fortran_loglikelihood(trades_parameters)
    if ( check == 0 ):
      loglhd = -0.5e10
    else:
      lnprior_value = lnprior(fitting_parameters, ndim)
      #lnprior_value = 0.0
      loglhd = loglhd + ln_err_const + lnprior_value # ln_err_const (global variable) & lnprior_value
    
    return loglhd


  # run MultiNest
  pymultinest.run(LogLikelihood = lnlike, Prior = trades_rescale, n_dims = nfit, n_params = nfit, outputfiles_basename=output_mnest, multimodal = multi_node_flag, log_zero=log_zero_value, seed = seed_value, n_live_points = n_pop, n_iter_before_update = n_update, resume = resume_flag, verbose = True, init_MPI = mpi_onoff)

  elapsed = time.time() - start
  elapsed_d, elapsed_h, elapsed_m, elapsed_s = computation_time(elapsed)

  logger.info('')
  logger.info('pyTRADES: pyMultiNest FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s))
  logger.info('')

  #of_run.write(' pyTRADES: pyMultiNest FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye\n' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s))
  #of_run.close()
  pytrades_lib.pytrades.deallocate_variables()

  return

if __name__ == "__main__":
  main()
