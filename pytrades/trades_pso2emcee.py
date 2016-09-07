#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
import pytrades_lib
import emcee
import h5py
import sys
import time
import shutil
#from emcee.utils import MPIPool
from constants import Mjups
from ancillary import *
#from pyde.de import DiffEvol as devol

#from matplotlib import use as mpluse
##mpluse("Agg")
#mpluse("Qt4Agg")
#import matplotlib.pyplot as plt
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('text', usetex=True)

#
#
#
def get_args():
  parser = argparse.ArgumentParser(description='TRADES+PSO+EMCEE')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--sub-folder', '--sb', action='store', dest='sub_folder', default='emcee_run', help='Sub-folder name, without full path. Default = emcee_run')
  
  # HOW TO RUN PSO OR SKIP
  parser.add_argument('-r', '--pso', '--pso-type', action='store', dest='pso_type', default='run', help='Define PSO run type: "skip" = do not run PSO, only emcee with random walkers; "run" = run PSO normally; "exists" = do not run PSO, but read previous PSO (pso_run.hdf5) file and start emcee from its population')
  
  # NUMBER OF CPU TO USE WITH EMCE !! WARNING !! TO USE PARALLEL IN PSO BEFORE THIS SCRIPT TYPE:
  # export OMP_NUM_THREADS=CPUNUMBER
  parser.add_argument('-c', '--cpu', '--nthreads', action='store', dest='nthreads', default=1, help='Number of threads to use. default nthreads = 1.')
  
  # NUMBER OF WALKERS TO USE WITH EMCEE
  parser.add_argument('-nw', '--nwalkers', '-np', '--npop', action='store', dest='nwalkers', default=1, help='Number of walkers (or number of chains) to use. default nwalkers = nfit*2')
  
  # NUMBER OF STEPS/RUNS TO DO FOR EACH WALKER OF EMCEE
  parser.add_argument('-nr', '--nruns', '-ns', '--nsteps', action='store', dest='nruns', default=10000, help='Number of runs/steps to use for each chain. default nruns = 10000.')
  
  # NUMBER OF STEPS/RUNS TO SAVE TEMPORARY EMCEE SIMULATION
  parser.add_argument('--isave', '--iter-save', '--iterations-save', action='store', dest='nsave', default='False', help='Number of iterations to do for save temporary chain. default each 0.1 of nruns.')
  
  # NUMBER OF BURNIN TO SKIP FOR PRELIMINARY ANALYSIS
  parser.add_argument('-nb', '--nburn', '--npost', action='store', dest='npost', default=1000, help='Number of burn in, or number of posterior to discard at the beginning of the each chain. default npost = 1000.')
  
  # COMPUTE OR NOT CONSTANT TO ADD TO THE LGLIKELIHOOD: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 )
  parser.add_argument('-l', '--ln-err', '--ln-err-const', action='store', dest='ln_flag', default=True, help='Computes or not constant to add to the lglikelihood: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 ). Default = True. Set to False if not use this value.')
  
  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  
  if (cli.pso_type not in ['skip', 'run', 'exists']):
    cli.pso_type = 'run'
  
  cli.nthreads = int(cli.nthreads)
  
  cli.nwalkers = int(cli.nwalkers)
  cli.nruns = int(cli.nruns)
  cli.npost = int(cli.npost)
  
  cli.ln_flag = set_bool_argument(cli.ln_flag)
    
  return cli


# 
# LOGPROBABILITY FUNCTION NEEDED BY EMCEE
#
def lnprob(fitting_parameters):
  loglhd = 0.
  check = 1
  loglhd, check = pytrades_lib.pytrades.fortran_loglikelihood(np.asarray(fitting_parameters,dtype=np.float64))
  #print loglhd, ln_err_const
  loglhd = loglhd + ln_err_const # ln_err_const: global variable
  #print loglhd
  if ( check == 0 ):
    loglhd = -np.inf
  return loglhd

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

def compute_ln_err_const(ndata, dof, e_RVo, e_T0o, ln_flag=False):
  if (ln_flag):
    eRV = e_RVo[e_RVo > 0.]
    eT0 = e_T0o[e_T0o > 0.]
    
    ln_e_RVo = np.sum(np.log(eRV*eRV))
    ln_e_T0o = np.sum(np.log(eT0*eT0))

    #ln_err_const = -(0.5*inv_dof) * (ndata*np.log(2.*np.pi) + ln_e_RVo + ln_e_T0o)
    ln_err_const = - 0.5 * dof * np.log(2.*np.pi) - 0.5 * ( ln_e_RVo + ln_e_T0o)
  else:
    ln_err_const = 0.
  return ln_err_const


#
# PSO SIMULATION TO EMCEE p0
#
def pso_to_emcee(nfit, nwalkers, pso_population, pso_population_fitness, pso_parameters, pso_fitness, pso_best_evolution):
  np_pso = pso_population_fitness.shape[0]
  ni_pso = pso_population_fitness.shape[1]
  np_ni = np_pso * ni_pso
  #print np_ni
  pso_fitness_p0 = np.zeros(nwalkers)
  p0 = np.zeros((nwalkers, nfit))
  
  flat_pop_fitness = pso_population_fitness.reshape((np_ni))
  #print flat_pop_fitness.shape
  flat_pop = pso_population.reshape((nfit, np_ni))
  #print flat_pop.shape
  sort_idx = np.argsort(flat_pop_fitness)
  
  if (str(pso_best_evolution) != 'False'):
    nw_best = int(nwalkers*0.25)+1
    #print nw_best
    sel_best = np.random.choice(ni_pso, nw_best)
    for i in range(1,nw_best+1):
      #print i, sel_best[i-1]
      p0[i,:] = pso_best_evolution[:nfit,sel_best[i-1]]
      pso_fitness_p0[i] = pso_best_evolution[-1,sel_best[i-1]]
      if (pso_fitness_p0[i] >= 1.e10):
        one_sel = np.random.choice(int(ni_pso*0.5), 1)
        p0[i,:] = pso_best_evolution[:nfit, one_sel[0]]
        pso_fitness_p0[i] = pso_best_evolution[-1, one_sel[0]]
    
    nw_pop = nwalkers - nw_best - 1
    #print nw_pop
    for i in range(0, nw_pop):
      ii=i+nw_best+1
      p0[ii,:] = flat_pop[:, sort_idx[i]]
      pso_fitness_p0[ii] = flat_pop_fitness[sort_idx[i]]

  else:
    for i in range(0,nwalkers-1):
      p0[i+1,:] = flat_pop[:, sort_idx[i]]
      pso_fitness_p0[i+1] = flat_pop_fitness[sort_idx[i]]
  
  #print pso_fitness, pso_parameters[1]
  p0[0,:] = np.asarray(pso_parameters, dtype=np.float64)
  pso_fitness_p0[0] = pso_fitness
  
  return p0, pso_fitness_p0

def print_p0(nwalkers, p0):
  lglhd_p0 = np.zeros(nwalkers)
  check_p0 = np.zeros(nwalkers).astype(int)
  for i in range(0, nwalkers):
    lglhd_temp, check_p0[i] = pytrades_lib.pytrades.fortran_loglikelihood(p0[i,:])
    lglhd_p0[i] = lnprob(p0[i])
    #line = ' p0[%4d] -> fitness(PSO) = %25.4f -> lnprob = %25.4f (lglhd = %25.4f, check = %5r, fitness = %25.4f)\n' %(i, pso_fitness_p0[i], lglhd_p0[i], lglhd_temp, check_p0[i].astype(bool), -2.*(lglhd_temp))
    line = ' p0[%4d] -> fitness(PSO) = %25.4f -> lnprob = %25.4f (lglhd = %25.4f, check = %5r, fitness = %25.4f)\n' %(i, pso_fitness_p0[i], lglhd_p0[i], lglhd_temp, check_p0[i].astype(bool), -2.*(lglhd_temp)*inv_dof) # inv_dof: global variable
    of_run.write(line)
    print line.strip('\n')


def get_emcee_arguments(cli,nfit):
  #
  # NUMBER OF WALKERS
  if (cli.nwalkers < nfit*2):
    nwalkers = nfit * 2
  else:
    nwalkers = cli.nwalkers
  if (nwalkers % 2) != 0:
    #print ' Provided odd nwalkers = %d ==> using nwalkers+1 = %d' %(nwalkers, nwalkers+1)
    nwalkers += 1

  # NUMBER OF STEPS/RUNS FOR EACH WALKER
  if (cli.nruns < 1):
    nruns = 10000
  else:
    nruns = cli.nruns

  # NUMBER OF SAVE STEPS
  if (cli.nsave != 'False'):
    if (int(cli.nsave) > 0 and int(cli.nsave) < nruns):
      nsave = int(cli.nsave)
    elif (int(cli.nsave) <= 0):
      nsave = False
    else:
      nsave = nruns/10
  else:
    nsave = False

  # NUMBER OF BURNIN/POSTERIOR TO DISCARD
  if (cli.npost < 0):
    npost = 1000
  else:
    npost = cli.npost
  #print nwalkers, nruns, npost
  return nwalkers, nruns, nsave, npost

# MAIN -- TRADES + EMCEE

# READ COMMAND LINE ARGUMENTS
cli = get_args()

# STARTING TIME
start = time.time()

# RENAME 
working_path = cli.full_path
nthreads=cli.nthreads

# INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder)

# RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
fitting_parameters = pytrades_lib.pytrades.fitting_parameters # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
parameters_minmax = pytrades_lib.pytrades.parameters_minmax # PARAMETER BOUNDARIES
delta_parameters = np.abs(parameters_minmax[:,1] - parameters_minmax[:,0]) # DELTA BETWEEN MAX AND MIN OF BOUNDARIES


#global n_bodies, n_planets, ndata, npar, nfit, dof, inv_dof
n_bodies = pytrades_lib.pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
ndata = pytrades_lib.pytrades.ndata # TOTAL NUMBER OF DATA AVAILABLE
npar  = pytrades_lib.pytrades.npar # NUMBER OF TOTAL PARAMATERS ~n_planets X 6
nfit  = pytrades_lib.pytrades.nfit # NUMBER OF PARAMETERS TO FIT
dof   = pytrades_lib.pytrades.dof # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
global inv_dof
inv_dof = np.float64(1.0 / dof)

# RADIAL VELOCITIES SET
n_rv = pytrades_lib.pytrades.nrv
n_set_rv = pytrades_lib.pytrades.nrvset

# TRANSITS SET
n_t0 = pytrades_lib.pytrades.nt0
n_t0_sum = np.sum(n_t0)
n_set_t0 = 0
for i in range(0, n_bodies):
  if (np.sum(n_t0[i]) > 0): n_set_t0 += 1

# compute global constant for the loglhd
global ln_err_const

try:
  e_RVo = np.asarray(pytrades_lib.pytrades.ervobs[:], dtype=np.float64) # fortran variable RV in python will be rv!!!
except:
  e_RVo = np.asarray([0.], dtype=np.float64)
try:
  e_T0o = np.asarray(pytrades_lib.pytrades.et0obs[:,:], dtype=np.float64).reshape((-1))
except:
  e_T0o = np.asarray([0.], dtype=np.float64)
ln_err_const = compute_ln_err_const(ndata, dof, e_RVo, e_T0o, cli.ln_flag)

# SET EMCEE PARAMETERS:
nwalkers, nruns, nsave, npost = get_emcee_arguments(cli,nfit)

# READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
reshaped_names = pytrades_lib.pytrades.parameter_names.reshape((10,nfit), order='F').T
parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]

# INITIALISE SCRIPT FOLDER/LOG FILE
working_folder, run_log, of_run = init_folder(working_path, cli.sub_folder)

print 
print ' ======== '
print ' pyTRADES'
print ' ======== '
print
print ' WORKING PATH = %s' %(working_path)
print ' NUMBER OF THREADS = %d' %(nthreads)
print ' dof = ndata(%d) - nfit(%d) = %d' %(ndata, nfit, dof)
print ' Total N_RV = %d for %d set(s)' %(n_rv, n_set_rv)
print ' Total N_T0 = %d for %d out of %d planet(s)' %(n_t0_sum, n_set_t0, n_planets)
print ' %s = %.7f' %('log constant error = ', ln_err_const)
print ' %s = %.7f' %('IN FORTRAN log constant error = ', pytrades_lib.pytrades.ln_err_const)
of_run.write(' NUMBER OF THREADS = %d\n' %(nthreads))
of_run.write(' dof = ndata(%d) - nfit(%d) = %d\n' %(ndata, nfit, dof))
of_run.write(' Total N_RV = %d for %d set(s)\n' %(n_rv, n_set_rv))
of_run.write(' Total N_T0 = %d for %d out of %d planet(s)\n' %(n_t0_sum, n_set_t0, n_planets))
of_run.write(' %s = %.7f\n' %('log constant error = ', ln_err_const))
of_run.write(' %s = %.7f\n' %('IN FORTRAN log constant error = ', pytrades_lib.pytrades.ln_err_const))

# INITIALISE PSO ARGUMENTS FROM pso.opt FILE
pytrades_lib.pytrades.init_pso(1,working_path) # read PSO options
# PSO VARIABLES
np_pso = pytrades_lib.pytrades.np_pso
nit_pso = pytrades_lib.pytrades.nit_pso
n_global = pytrades_lib.pytrades.n_global
#n_global = 1
print ' PSO n_global = %d npop = %d ngen = %d' %(n_global, np_pso, nit_pso)
of_run.write("# PSO\n n_global = %d\n npop = %d\n ngen = %d\n" %(n_global, np_pso, nit_pso))

# RUN PSO+EMCEE n_global TIMES
for iter_global in range(0,n_global):

  # CREATES PROPER WORKING PATH AND NAME
  i_global = iter_global + 1
  pso_path = os.path.join(os.path.join(working_folder, '%04d_pso' %(i_global)), '')
  pytrades_lib.pytrades.path_change(pso_path)
  
  print '\n\n GLOBAL RUN %04d INTO PATH: %s' %(i_global, pso_path)
  of_run.write('\n\n GLOBAL RUN %04d INTO PATH: %s\n' %(i_global, pso_path))

  if (cli.pso_type == 'run'):
    # RUN PSO
    print ' RUN PSO'

    pso_start = time.time()
    if(not os.path.exists(pso_path)): os.makedirs(pso_path)

    # CALL RUN_PSO SUBROUTINE FROM TRADES_LIB: RUNS PSO AND COMPUTES THE BEST SOLUTION, SAVING ALL THE POPULATION EVOLUTION
    pso_parameters = np.zeros((nfit)) + fitting_parameters
    pso_fitness = 0.
    pso_parameters, pso_fitness = pytrades_lib.pytrades.pyrun_pso(nfit,i_global)
    print ' completed run_pso'
    
    print ' fortran system_parameters'
    print pytrades_lib.pytrades.system_parameters
    
    pso_best_evolution = np.asarray(pytrades_lib.pytrades.pso_best_evolution[...], dtype=np.float64)
    print ' pso_best_evolution retrieved'
    
    print ' last pso_best_evolution'
    best_parameters = np.asarray(pso_best_evolution[:nfit,-1],dtype=np.float64)
    print best_parameters
    best_fitness = pso_best_evolution[-1,-1].astype(np.float64)
    print ' fitness = ',best_fitness
    
    
    print ' Best pso_parameters'
    of_run.write(' Best pso_parameters\n')
    for ii in range(0,nfit):
      print ' %23s %23.16E' %(parameter_names[ii],pso_parameters[ii])
      of_run.write(' %23s %23.16E\n' %(parameter_names[ii],pso_parameters[ii]))
    #print pso_parameters
    print ' Best PSO fitness = %23.16E' %(pso_fitness)
    of_run.write(" Best PSO fitness = %23.16E\n" %(pso_fitness))
    
    # test fortran_loglikelihood and lnprob
    print '---'
    lglhd_fortran, check_fortran = pytrades_lib.pytrades.fortran_loglikelihood(np.asarray(pso_best_evolution[:nfit,-1], dtype=np.float64))
    print 'lglhd_fortran = %23.16E check_fortran = %r' %(lglhd_fortran, check_fortran)
    lglhd_lnprob = lnprob(np.asarray(pso_best_evolution[:nfit,-1], dtype=np.float64))
    print 'lglhd_lnprob = %23.16E' %(lglhd_lnprob)
    print '---'
    lglhd_fortran, check_fortran = pytrades_lib.pytrades.fortran_loglikelihood(best_parameters)
    print 'lglhd_fortran = %23.16E check_fortran = %r' %(lglhd_fortran, check_fortran)
    lglhd_lnprob = lnprob(best_parameters)
    print 'lglhd_lnprob = %23.16E' %(lglhd_lnprob)
    print '---'
    
    
    # SAVE PSO SIMULATION IN pso_run.hdf5 FILE
    print ' Creating pso hdf5 file: %s' %(os.path.join(pso_path, 'pso_run.hdf5'))
    pso_hdf5 = h5py.File(os.path.join(pso_path, 'pso_run.hdf5'), 'w')
    pso_hdf5.create_dataset('population', data=pytrades_lib.pytrades.population, dtype=np.float64)
    pso_hdf5.create_dataset('population_fitness', data=pytrades_lib.pytrades.population_fitness, dtype=np.float64)
    pso_hdf5.create_dataset('pso_parameters', data=pso_parameters, dtype=np.float64)
    pso_hdf5.create_dataset('pso_fitness', data=np.array(pso_fitness), dtype=np.float64)
    pso_hdf5.create_dataset('pso_best_evolution', data=pso_best_evolution, dtype=np.float64)
    pso_hdf5.create_dataset('parameters_minmax', data=parameters_minmax, dtype=np.float64)
    pso_hdf5.create_dataset('parameter_names', data=parameter_names, dtype='S10')
    pso_hdf5.close()

    population = np.asarray(pytrades_lib.pytrades.population, dtype=np.float64)
    population_fitness = np.asarray(pytrades_lib.pytrades.population_fitness, dtype=np.float64)
    p0, pso_fitness_p0 = pso_to_emcee(nfit, nwalkers, population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution)
    
    print_p0(nwalkers, p0)
    
    print
    fitness_iter, lgllhd_iter, check_iter = pytrades_lib.pytrades.write_summary_files(i_global, pso_parameters)
    elapsed = time.time() - pso_start
    elapsed_d, elapsed_h, elapsed_m, elapsed_s = computation_time(elapsed)
    print
    print ' PSO FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s)
    of_run.write(' PSO FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye\n' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s))
    
  
  elif (cli.pso_type == 'exists'):
    # READ PREVIOUS PSO_RUN.HDF5 FILE AND INITIALISE POPULATION FOR EMCEE
    print ' READ PREVIOUS PSO_RUN.HDF5 FILE AND INITIALISE POPULATION FOR EMCEE'
    
    population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution, pso_parameters_minmax, pso_parameter_names, pop_shape = get_pso_data(os.path.join(pso_path, 'pso_run.hdf5'))
    
    fitness_iter, lgllhd_iter, check_iter = pytrades_lib.pytrades.write_summary_files(i_global, pso_parameters)
    
    print ' read pso_run.hdf5 file with '
    print ' best pso_fitness = %.7f' %(pso_fitness)
    
    p0, pso_fitness_p0 = pso_to_emcee(nfit, nwalkers, population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution)
        
    print_p0(nwalkers, p0)
    
  elif (cli.pso_type == 'skip'):
    # DO NOT RUN PSO, ONLY EMCEE
    print ' DO NOT RUN PSO, ONLY EMCEE'
    
    p0 = [parameters_minmax[:,0] + np.random.random(nfit)*delta_parameters for i in range(0, nwalkers)]

  #sys.exit('MERDACANEPORCOBIO')

  print ' emcee chain: nwalkers = %d nruns = %d' %(nwalkers, nruns)
  of_run.write("# EMCEE chain:\n nwalkers = %d\n nruns = %d\n" %(nwalkers, nruns))
  print ' sampler ... '
  sampler = emcee.EnsembleSampler(nwalkers, nfit, lnprob, threads=nthreads)
  print ' ready to go'
  print ' with nsave = %r' %(nsave)
  sys.stdout.flush()
  
  #sys.exit()
  
  if (nsave != False):
    # save temporary sampling during emcee every nruns*10%
    if(os.path.exists(os.path.join(pso_path, 'emcee_temp.hdf5')) and os.path.isfile(os.path.join(pso_path, 'emcee_temp.hdf5'))):
      os.remove(os.path.join(pso_path, 'emcee_temp.hdf5'))
    of_temp = h5py.File(os.path.join(pso_path, 'emcee_temp.hdf5'), 'a')
    of_temp.create_dataset('parameter_names', data=parameter_names, dtype='S10')
    of_temp.create_dataset('boundaries', data=parameters_minmax, dtype=np.float64)
    temp_dset = of_temp.create_dataset('chains', (nwalkers, nruns, nfit), dtype=np.float64)
    temp_lnprob = of_temp.create_dataset('lnprobability', (nwalkers, nruns), dtype=np.float64)
    temp_lnprob.attrs['ln_err_const'] = ln_err_const
    of_temp.close()
    pos = p0
    nchains = int(nruns/nsave)
    state=None
    print ' Running emcee with temporary saving'
    sys.stdout.flush()
    for i in range(0, nchains):
      print
      print ' iter: ',i+1,
      aaa = i*nsave
      bbb = aaa+nsave
      pos, prob, state = sampler.run_mcmc(pos, N=nsave, rstate0=state)
      print 'completed %d steps of %d' %(bbb, nruns)
      of_temp = h5py.File(os.path.join(pso_path, 'emcee_temp.hdf5'), 'a')
      temp_dset = of_temp['chains'] #[:,:,:]
      temp_dset[:,aaa:bbb,:] = sampler.chain[:, aaa:bbb, :]
      #of_temp['chains'].attrs['completed_steps'] = bbb
      temp_dset.attrs['completed_steps'] = bbb
      temp_lnprob = of_temp['lnprobability'] #[:,:]
      temp_lnprob[:, aaa:bbb] = sampler.lnprobability[:, aaa:bbb]
      shape_lnprob = sampler.lnprobability.shape
      #print 'aaa = %6d bbb = %6d -> sampler.lnprobability.shape = (%6d , %6d)' %(aaa, bbb, shape_lnprob[0], shape_lnprob[1])
      of_temp.close()
      sys.stdout.flush()
    print
    print '...done with saving temporary total shape = ', sampler.chain.shape
    print
    sys.stdout.flush()

  # RUN EMCEE AND RESET AFTER REMOVE BURN-IN
  #pos, prob, state = sampler.run_mcmc(p0, npost)
  #sampler.reset()
  #sampler.run_mcmc(pos, nruns, rstate0=state)
  else:
    # GOOD COMPLETE SINGLE RUNNING OF EMCEE, WITHOUT REMOVING THE BURN-IN
    print ' Running full emcee ...',
    sys.stdout.flush()
    sampler.run_mcmc(p0, nruns)
    print 'done'
    print
    sys.stdout.flush()

  flatchains = sampler.chain[:, :, :].reshape((nwalkers*nruns, nfit)) # full chain values
  acceptance_fraction = sampler.acceptance_fraction
  mean_acceptance_fraction = np.mean(acceptance_fraction)
  autocor_time = sampler.acor
  lnprobability = sampler.lnprobability

  # save chains with original shape as hdf5 file
  f_hdf5 = h5py.File(os.path.join(pso_path, 'emcee_summary.hdf5'), 'w')
  f_hdf5.create_dataset('chains', data=sampler.chain, dtype=np.float64)
  f_hdf5.create_dataset('parameter_names', data=parameter_names, dtype='S10')
  f_hdf5.create_dataset('boundaries', data=parameters_minmax, dtype=np.float64)
  #f_hdf5.create_dataset('final_parameters', data=final_parameters, dtype=np.float64)
  #f_hdf5.create_dataset('max_lnprob_parameters', data=max_lnprob_parameters, dtype=np.float64)
  #f_hdf5.create_dataset('max_lnprob', data=max_lnprob, dtype=np.float64)
  f_hdf5.create_dataset('acceptance_fraction', data=acceptance_fraction, dtype=np.float64)
  f_hdf5.create_dataset('autocor_time', data=autocor_time, dtype=np.float64)
  f_hdf5.create_dataset('lnprobability', data=lnprobability, dtype=np.float64)
  f_hdf5['lnprobability'].attrs['ln_err_const'] = ln_err_const
  f_hdf5.close()

  nruns_sel = nruns - npost
  # chain is transposed: needed to plot quicker chains for each walker: nruns vs value of parameter
  chains_T, parameter_boundaries = select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, 1., parameter_names, parameters_minmax, sampler.chain)
   
  flatchain_posterior_0 = chains_T[:,:,:].reshape((nruns_sel*nwalkers, nfit))

  of_run.write(" Mean_acceptance_fraction should be between [0.25-0.5] = %.6f\n" %(mean_acceptance_fraction))
  print(" Mean acceptance fraction should be between [0.25-0.5]: {0:.6f}".format(mean_acceptance_fraction))

  print

  derived_names, derived_chains_T, derived_posterior = get_derived_posterior_parameters(parameter_names, chains_T, flatchain_posterior_0)
  nder = len(derived_names)

  # GET MAX LNPROBABILITY AND PARAMETERS -> id 40XX
  max_lnprob, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint = get_maxlnprob_parameters(npost, nruns, lnprobability, chains_T, flatchain_posterior_0)
  max_lnprob_der, max_lnprob_parameters_der, max_lnprob_perc68_der, max_lnprob_confint_der = get_maxlnprob_parameters(npost, nruns, lnprobability, derived_chains_T, derived_posterior)
  
  # MEDIAN PARAMETERS
  
  # std way: median of the posterior parameter distribution -> id 10XX
  median_parameters, median_perc68, median_confint = get_median_parameters(flatchain_posterior_0)
  median_parameters_der, median_perc68_der, median_confint_der = get_median_parameters(derived_posterior)
  
  # parameters linked to the median of the fitness =  - 2 * (lglhd + ln_err_const) -> id 20XX
  median_fitness, medfit_parameters, medfit_perc68, medfit_confint = get_parameters_median_fitness(nwalkers, npost, nruns, lnprobability, flatchain_posterior_0, ln_err_const)
  median_fitness_der, medfit_parameters_der, medfit_perc68_der, medfit_confint_der = get_parameters_median_fitness(nwalkers, npost, nruns, lnprobability, derived_posterior, ln_err_const)
  
  # MODE-LIKE PARAMETERS -> id 30XX
  # take the mean of 5 bin centered to the higher bin
  k = 11
  mode_bin, mode_parameters, mode_perc68, mode_confint = get_mode_parameters(flatchain_posterior_0, k)
  mode_bin_der, mode_parameters_der, mode_perc68_der, mode_confint_der = get_mode_parameters(derived_posterior, k)


  print ' MAX LNPROBABILITY PARAMETER VALUES -> 4000'
  of_run.write(' MAX LNPROBABILITY PARAMETER VALUES -> 4000\n')
  fitness_4000, lgllhd_4000, check_4000 = pytrades_lib.pytrades.write_summary_files(4000, max_lnprob_parameters)
  
  print ' MEDIAN PARAMETER VALUES -> 1050'
  of_run.write(' MEDIAN PARAMETER VALUES -> 1050\n')
  fitness_1050, lgllhd_1050, check_1050 = pytrades_lib.pytrades.write_summary_files(1050, median_parameters)
  
  print ' MEDFIT PARAMETER VALUES -> 2050'
  of_run.write(' MEDFIT PARAMETER VALUES -> 2050\n')
  fitness_2050, lgllhd_2050, check_2050 = pytrades_lib.pytrades.write_summary_files(2050, medfit_parameters)

  print ' MODE PARAMETER VALUES -> 3050'
  of_run.write(' MODE PARAMETER VALUES -> 3050\n')
  fitness_3050, lgllhd_3050, check_3050 = pytrades_lib.pytrades.write_summary_files(3050, mode_parameters)

  
  print ' MAX LNPROBABILITY PARAMETER VALUES -> 4000'
  of_run.write(' MAX LNPROBABILITY PARAMETER VALUES -> 4000\n')
  print_parameters_logtxt(of_run, parameter_names, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint, 'maxlnpr')
  print_parameters_logtxt(of_run, derived_names, max_lnprob_parameters_der, max_lnprob_perc68_der, max_lnprob_confint_der, 'maxlnpr_der')
  print 'FITNESS (4000) = %23.16f LOGLIKELIHOOD = %23.16f\n' %(fitness_4000, lgllhd_4000+ln_err_const)
  of_run.write('FITNESS (4000) = %23.16f LOGLIKELIHOOD = %23.16f\n\n' %(fitness_4000, lgllhd_4000+ln_err_const))
  sys.stdout.flush()

  print ' MEDIAN PARAMETER VALUES -> 1050'
  of_run.write(' MEDIAN PARAMETER VALUES -> 1050\n')
  print_parameters_logtxt(of_run, parameter_names, median_parameters, median_perc68, median_confint, 'median')
  print_parameters_logtxt(of_run, derived_names, median_parameters_der, median_perc68_der, median_confint_der, 'median_der')
  print 'FITNESS (1050) = %23.16f LOGLIKELIHOOD = %23.16f\n' %(fitness_1050, lgllhd_1050+ln_err_const)
  of_run.write('FITNESS (1050) = %23.16f LOGLIKELIHOOD = %23.16f\n\n' %(fitness_1050, lgllhd_1050+ln_err_const))
  sys.stdout.flush()
  
  print ' MEDFIT PARAMETER VALUES -> 2050'
  of_run.write(' MEDFIT PARAMETER VALUES -> 2050\n')
  print_parameters_logtxt(of_run, parameter_names, medfit_parameters, medfit_perc68, medfit_confint, 'medfit')
  print_parameters_logtxt(of_run, derived_names, medfit_parameters_der, medfit_perc68_der, medfit_confint_der, 'medfit_der')
  print 'FITNESS (2050) = %23.16f LOGLIKELIHOOD = %23.16f\n' %(fitness_2050, lgllhd_2050+ln_err_const)
  of_run.write('FITNESS (2050) = %23.16f LOGLIKELIHOOD = %23.16f\n\n' %(fitness_2050, lgllhd_2050+ln_err_const))
  sys.stdout.flush()
  
    
  print ' MODE PARAMETER VALUES -> 3050'
  of_run.write(' MODE PARAMETER VALUES -> 3050\n')
  print_parameters_logtxt(of_run, parameter_names, mode_parameters, mode_perc68, mode_confint, 'mode')
  print_parameters_logtxt(of_run, derived_names, mode_parameters_der, mode_perc68_der, mode_confint_der, 'mode_der')
  print 'FITNESS (3050) = %23.16f LOGLIKELIHOOD = %23.16f\n' %(fitness_3050, lgllhd_3050+ln_err_const)
  of_run.write('FITNESS (3050) = %23.16f LOGLIKELIHOOD = %23.16f\n\n' %(fitness_3050, lgllhd_3050+ln_err_const))
  sys.stdout.flush()
  
  del chains_T
  del flatchain_posterior_0

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = computation_time(elapsed)

print
print ' pyTRADES: PSO+EMCEE FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s)
of_run.write(' pyTRADES: PSO+EMCEE FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye\n' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s))
of_run.close()
print 
pytrades_lib.pytrades.deallocate_variables()




