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
import ancillary as anc
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
  parser = argparse.ArgumentParser(description='TRADES+EMCEE')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--sub-folder', '--sb', action='store', dest='sub_folder', default='emcee_run', help='Sub-folder name, without full path. Default = emcee_run')
  
  # NUMBER OF CPU TO USE WITH EMCE !!
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
  
  parser.add_argument('-ds', '--d-sigma', '--delta-sigma', action='store', dest='delta_sigma', default=1.e-4, type=np.float64, help='Value of the sigma to compute initial walkers from initial solution. Default=1.e-4')
  
  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  
  cli.nthreads = int(cli.nthreads)
  
  cli.nwalkers = int(cli.nwalkers)
  cli.nruns = int(cli.nruns)
  cli.npost = int(cli.npost)
  
  cli.ln_flag = anc.set_bool_argument(cli.ln_flag)
  
  try:
    cli.delta_sigma = np.float64(cli.delta_sigma)
  except:
    cli.delta_sigma = np.float64(1.e-4)
    
  return cli


# 
# LOGPROBABILITY FUNCTION NEEDED BY EMCEE
#
def lnprob(fitting_parameters):
  loglhd = 0.
  check = 1
  loglhd, check = pytrades_lib.pytrades.fortran_loglikelihood(np.array(fitting_parameters,dtype=np.float64))
  #print loglhd, ln_err_const
  loglhd = loglhd + ln_err_const # ln_err_const: global variable
  #print loglhd
  if ( check == 0 ):
    loglhd = -np.inf
  return loglhd

def print_both(line, log_file=None):
  print line
  if(log_file is not None):
    log_file.write(line + '\n')
  return


#
# INITIALISE FOLDER AND LOG FILE
#
def init_folder(working_path, sub_folder):
  working_folder = os.path.join(working_path, sub_folder)
  if (not os.path.isdir(working_folder)):
      os.makedirs(working_folder)
  bodies_file = os.path.join(working_path, 'bodies.lst')
  shutil.copy(bodies_file, os.path.join(working_folder,''))
  obd = open(bodies_file, 'r')
  for line in obd.readlines():
    shutil.copy(os.path.join(working_path, line.strip().split()[0]), os.path.join(working_folder,''))
  obd.close()
  run_log = os.path.join(working_folder, "trades_run.log")
  of_run = open(run_log, 'w')
  print_both("# pyTRADES LOG FILE", of_run)
  print_both("# working_path = %s" %(working_path), of_run)
  print_both("# working_folder = %s" %(working_folder), of_run)
  print_both("# run_log = %s" %(run_log), of_run)
  
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


def compute_initial_walkers(nfit, nwalkers, fitting_parameters, delta_parameters, delta_sigma, of_run):
  # initial walkers as input fitting_parameters + N(loc=0.,sigma=1.,size=nwalkers)*delta_sigma
  #p0 = [parameters_minmax[:,0] + np.random.random(nfit)*delta_parameters for i in range(0, nwalkers)]
  print_both(' Inititializing walkers with delta_sigma = %.16f' %(cli.delta_sigma),of_run)
  p0 = []
  i_p0 = 0
  print_both(' good p0:', of_run)
  while True:
    test_p0 = fitting_parameters + np.random.normal(loc=0., scale=1., size=nfit)*cli.delta_sigma
    #test_p0 = fitting_parameters + fitting_parameters*np.random.normal(loc=0., scale=1., size=nfit)*cli.delta_sigma
    #test_p0 = fitting_parameters + delta_parameters*np.random.normal(loc=0., scale=1., size=nfit)*cli.delta_sigma
    test_lg = lnprob(test_p0)
    if(not np.isinf(test_lg)):
      i_p0 +=1
      p0.append(test_p0)
      print i_p0,
      if(i_p0 == nwalkers): break
  print
  print_both(' done initial walkers.', of_run)
  return p0


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
  e_RVo = np.array(pytrades_lib.pytrades.ervobs[:], dtype=np.float64) # fortran variable RV in python will be rv!!!
except:
  e_RVo = np.array([0.], dtype=np.float64)
try:
  e_T0o = np.array(pytrades_lib.pytrades.et0obs[:,:], dtype=np.float64).reshape((-1))
except:
  e_T0o = np.array([0.], dtype=np.float64)
ln_err_const = compute_ln_err_const(ndata, dof, e_RVo, e_T0o, cli.ln_flag)

# SET EMCEE PARAMETERS:
nwalkers, nruns, nsave, npost = get_emcee_arguments(cli,nfit)

# READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
reshaped_names = pytrades_lib.pytrades.parameter_names.reshape((10,nfit), order='F').T
parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]

# INITIALISE SCRIPT FOLDER/LOG FILE
working_folder, run_log, of_run = init_folder(working_path, cli.sub_folder)

print_both('',of_run)
print_both(' ======== ',of_run)
print_both(' pyTRADES' ,of_run)
print_both(' ======== ',of_run)
print_both('',of_run)
print_both(' WORKING PATH = %s' %(working_path),of_run)
print_both(' NUMBER OF THREADS = %d' %(nthreads),of_run)
print_both(' dof = ndata(%d) - nfit(%d) = %d' %(ndata, nfit, dof),of_run)
print_both(' Total N_RV = %d for %d set(s)' %(n_rv, n_set_rv),of_run)
print_both(' Total N_T0 = %d for %d out of %d planet(s)' %(n_t0_sum, n_set_t0, n_planets),of_run)
print_both(' %s = %.7f' %('log constant error = ', ln_err_const),of_run)
print_both(' %s = %.7f' %('IN FORTRAN log constant error = ', pytrades_lib.pytrades.ln_err_const),of_run)

# save initial_fitting parameters into array
original_fit_parameters = fitting_parameters.copy()

p0 = compute_initial_walkers(nfit, nwalkers, fitting_parameters, delta_parameters, cli.delta_sigma, of_run)

print_both(' emcee chain: nwalkers = %d nruns = %d' %(nwalkers, nruns), of_run)
print_both(' sampler ... ',of_run)
sampler = emcee.EnsembleSampler(nwalkers, nfit, lnprob, threads=nthreads)
print_both(' ready to go', of_run)
print_both(' with nsave = %r' %(nsave), of_run)
sys.stdout.flush()

#sys.exit()

if (nsave != False):
  # save temporary sampling during emcee every nruns*10%
  if(os.path.exists(os.path.join(working_folder, 'emcee_temp.hdf5')) and os.path.isfile(os.path.join(working_folder, 'emcee_temp.hdf5'))):
    os.remove(os.path.join(working_folder, 'emcee_temp.hdf5'))
  of_temp = h5py.File(os.path.join(working_folder, 'emcee_temp.hdf5'), 'a')
  of_temp.create_dataset('parameter_names', data=parameter_names, dtype='S10')
  of_temp.create_dataset('boundaries', data=parameters_minmax, dtype=np.float64)
  temp_dset = of_temp.create_dataset('chains', (nwalkers, nruns, nfit), dtype=np.float64)
  temp_lnprob = of_temp.create_dataset('lnprobability', (nwalkers, nruns), dtype=np.float64)
  temp_lnprob.attrs['ln_err_const'] = ln_err_const
  of_temp.close()
  pos = p0
  nchains = int(nruns/nsave)
  state=None
  print_both(' Running emcee with temporary saving', of_run)
  sys.stdout.flush()
  for i in range(0, nchains):
    print_both('', of_run)
    print_both(' iter: %6d ' %(i+1), of_run)
    aaa = i*nsave
    bbb = aaa+nsave
    pos, prob, state = sampler.run_mcmc(pos, N=nsave, rstate0=state)
    print_both('completed %d steps of %d' %(bbb, nruns), of_run)
    of_temp = h5py.File(os.path.join(working_folder, 'emcee_temp.hdf5'), 'a')
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
  print_both('', of_run)
  print_both('...done with saving temporary total shape = %s' %(str(sampler.chain.shape)), of_run)
  print_both('', of_run)
  sys.stdout.flush()

# RUN EMCEE AND RESET AFTER REMOVE BURN-IN
#pos, prob, state = sampler.run_mcmc(p0, npost)
#sampler.reset()
#sampler.run_mcmc(pos, nruns, rstate0=state)
else:
  # GOOD COMPLETE SINGLE RUNNING OF EMCEE, WITHOUT REMOVING THE BURN-IN
  print_both(' Running full emcee ...', of_run)
  sys.stdout.flush()
  sampler.run_mcmc(p0, nruns)
  print_both('done', of_run)
  print_both('', of_run)
  sys.stdout.flush()

flatchains = sampler.chain[:, :, :].reshape((nwalkers*nruns, nfit)) # full chain values
acceptance_fraction = sampler.acceptance_fraction
mean_acceptance_fraction = np.mean(acceptance_fraction)
autocor_time = sampler.acor
lnprobability = sampler.lnprobability

# save chains with original shape as hdf5 file
f_hdf5 = h5py.File(os.path.join(working_folder, 'emcee_summary.hdf5'), 'w')
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
chains_T, parameter_boundaries = anc.select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, 1., parameter_names, parameters_minmax, sampler.chain)
 
flatchain_posterior_0 = chains_T[:,:,:].reshape((nruns_sel*nwalkers, nfit))

print_both(" Mean_acceptance_fraction should be between [0.25-0.5] = %.6f" %(mean_acceptance_fraction), of_run)
print_both('', of_run)

derived_names, derived_chains_T, derived_posterior = anc.get_derived_posterior_parameters(parameter_names, chains_T, flatchain_posterior_0)
nder = len(derived_names)

# GET MAX LNPROBABILITY AND PARAMETERS -> id 40XX
max_lnprob, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint = anc.get_maxlnprob_parameters(npost, nruns, lnprobability, chains_T, flatchain_posterior_0)
max_lnprob_der, max_lnprob_parameters_der, max_lnprob_perc68_der, max_lnprob_confint_der = anc.get_maxlnprob_parameters(npost, nruns, lnprobability, derived_chains_T, derived_posterior)

# MEDIAN PARAMETERS

# std way: median of the posterior parameter distribution -> id 10XX
median_parameters, median_perc68, median_confint = anc.get_median_parameters(flatchain_posterior_0)
median_parameters_der, median_perc68_der, median_confint_der = anc.get_median_parameters(derived_posterior)

# parameters linked to the median of the fitness =  - 2 * (lglhd + ln_err_const) -> id 20XX
median_fitness, medfit_parameters, medfit_perc68, medfit_confint = anc.get_parameters_median_fitness(nwalkers, npost, nruns, lnprobability, flatchain_posterior_0, ln_err_const)
median_fitness_der, medfit_parameters_der, medfit_perc68_der, medfit_confint_der = anc.get_parameters_median_fitness(nwalkers, npost, nruns, lnprobability, derived_posterior, ln_err_const)

# MODE-LIKE PARAMETERS -> id 30XX
# take the mean of 5 bin centered to the higher bin
k = 11
mode_bin, mode_parameters, mode_perc68, mode_confint = anc.get_mode_parameters(flatchain_posterior_0, k)
mode_bin_der, mode_parameters_der, mode_perc68_der, mode_confint_der = anc.get_mode_parameters(derived_posterior, k)


print_both(' MAX LNPROBABILITY PARAMETER VALUES -> 4000', of_run)
fitness_4000, lgllhd_4000, check_4000 = pytrades_lib.pytrades.write_summary_files(4000, max_lnprob_parameters)

print_both(' MEDIAN PARAMETER VALUES -> 1050', of_run)
fitness_1050, lgllhd_1050, check_1050 = pytrades_lib.pytrades.write_summary_files(1050, median_parameters)

print_both(' MEDFIT PARAMETER VALUES -> 2050', of_run)
fitness_2050, lgllhd_2050, check_2050 = pytrades_lib.pytrades.write_summary_files(2050, medfit_parameters)

print_both(' MODE PARAMETER VALUES -> 3050', of_run)
fitness_3050, lgllhd_3050, check_3050 = pytrades_lib.pytrades.write_summary_files(3050, mode_parameters)


print_both(' MAX LNPROBABILITY PARAMETER VALUES -> 4000', of_run)
anc.print_parameters_logtxt(of_run, parameter_names, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint, 'maxlnpr')
anc.print_parameters_logtxt(of_run, derived_names, max_lnprob_parameters_der, max_lnprob_perc68_der, max_lnprob_confint_der, 'maxlnpr_der')
print_both('FITNESS (4000) = %23.16f LOGLIKELIHOOD = %23.16f\n' %(fitness_4000, lgllhd_4000+ln_err_const), of_run)
sys.stdout.flush()

print_both(' MEDIAN PARAMETER VALUES -> 1050', of_run)
anc.print_parameters_logtxt(of_run, parameter_names, median_parameters, median_perc68, median_confint, 'median')
anc.print_parameters_logtxt(of_run, derived_names, median_parameters_der, median_perc68_der, median_confint_der, 'median_der')
print_both('FITNESS (1050) = %23.16f LOGLIKELIHOOD = %23.16f' %(fitness_1050, lgllhd_1050+ln_err_const), of_run)
sys.stdout.flush()

print_both(' MEDFIT PARAMETER VALUES -> 2050', of_run)
anc.print_parameters_logtxt(of_run, parameter_names, medfit_parameters, medfit_perc68, medfit_confint, 'medfit')
anc.print_parameters_logtxt(of_run, derived_names, medfit_parameters_der, medfit_perc68_der, medfit_confint_der, 'medfit_der')
print_both('FITNESS (2050) = %23.16f LOGLIKELIHOOD = %23.16f\n' %(fitness_2050, lgllhd_2050+ln_err_const), of_run)
sys.stdout.flush()

  
print_both(' MODE PARAMETER VALUES -> 3050', of_run)
anc.print_parameters_logtxt(of_run, parameter_names, mode_parameters, mode_perc68, mode_confint, 'mode')
anc.print_parameters_logtxt(of_run, derived_names, mode_parameters_der, mode_perc68_der, mode_confint_der, 'mode_der')
print_both('FITNESS (3050) = %23.16f LOGLIKELIHOOD = %23.16f\n' %(fitness_3050, lgllhd_3050+ln_err_const), of_run)
sys.stdout.flush()

#print_both(' COMPUTE SIGMAS FROM INITIAL FITTING PARAMETERS')
#fit_68res, fit_confint = anc.get_sigmas(original_fit_parameters, flatchain_posterior_0)
print_both(' ORIGINAL PARAMETER VALUES -> 0000', of_run)
fitness_0000, lgllhd_0000, check_0000 = pytrades_lib.pytrades.write_summary_files(0, original_fit_parameters)
# save flatchain_posterior_0 as bootstrap simulation!
#header_0000 = ' iboot %s' %(' '.join(parameter_names))
#fmt_dp = ' '.join( [ '%s' %('%23.16e') for ii in range(0,nfit) ] )
#fmt_full = '%s %s' %('%6d', fmt_dp)
#iboot_fake = np.arange(1,flatchain_posterior_0.shape[0]+1,1)
#boot_fake = np.column_stack((iboot_fake, flatchain_posterior_0))
#boot_file = os.path.join(working_folder, '0_bootstrap_sim.dat')
#np.savetxt(boot_file, boot_fake, fmt=fmt_full, header=header_0000)
boot_file = anc.save_bootstrap_like(working_folder, 0, parameter_names, flatchain_posterior_0)
print_both(' saved fake bootstrap file: %s' %(boot_file), of_run)

del chains_T
del flatchain_posterior_0

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)

print_both('', of_run)
print_both(' pyTRADES: EMCEE FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s), of_run)
print_both('', of_run)
of_run.close()
pytrades_lib.pytrades.deallocate_variables()




