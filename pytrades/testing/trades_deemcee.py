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
from emcee.utils import MPIPool
from constants import Mjups
from ancillary import *
#from scipy.optimize import differential_evolution as devol
from pyde.de import DiffEvol as devol
from matplotlib import use as mpluse
#mpluse("Agg")
mpluse("Qt4Agg")
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)

#def set_bool_argument(arg):
  #if (arg != False):
    #if (arg.lower() in ['t', 'tr', 'tru', 'true']):
      #arg = True
    #else:
      #arg = False
  #return arg

def get_args():
  parser = argparse.ArgumentParser(description='TRADES+EMCEE')
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')
  parser.add_argument('-s', '--sub-folder', '--sb', action='store', dest='sub_folder', default='emcee_run', help='Sub-folder name, without full path. Default = emcee_run')
  parser.add_argument('-m', '--mpi', action='store', dest='mpi_set', default=False, help='MPI: use MPI with True, otherwise set it to False. Defautl is False.')
  parser.add_argument('-c', '--cpu', '--nthreads', action='store', dest='nthreads', default=1, help='Number of threads to use. default nthreads = 1.')
  parser.add_argument('-nw', '--nwalkers', '-np', '--npop', action='store', dest='nwalkers', default=1, help='Number of walkers (or number of chains) to use. default nwalkers = nfit*2')
  parser.add_argument('-nr', '--nruns', '-ns', '--nsteps', action='store', dest='nruns', default=10000, help='Number of runs/steps to use for each chain. default nruns = 10000.')
  parser.add_argument('--isave', '--iter-save', '--iterations-save', action='store', dest='nsave', default='False', help='Number of iterations to do for save temporary chain. default each 0.1 of nruns.')
  parser.add_argument('-nb', '--nburn', '--npost', action='store', dest='npost', default=1000, help='Number of burn in, or number of posterior to discard at the beginning of the each chain. default npost = 1000.')
  parser.add_argument('-g', '--good-parameters', action='store', dest='good', default=False, help='If you want to use a previous solution, set it to True and it will search for a good_parameters.dat file with first column the name of the parameter and the value in the second. Mass paramenters in Jupiter mass. Default is False')
  parser.add_argument('-d', '--de', '--devol', '--differential-evolution', action='store', dest='devol', default=False, help='If you want to use the Differential Evolution algorithm to create the initial set of walkers for emcee. Default: False')
  parser.add_argument('-ng', '--ngen', action='store', dest='ngen', default=1, help='Number of generation for the Differential Evolution. Minimum value, the Default, sets to nwalkers * 100.')
  
  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  cli.nthreads = int(cli.nthreads)
  
  #if (cli.mpi_set != False):
    #if (cli.mpi_set.lower() == in ['t', 'tr', 'tru', 'true']):
      #cli.mpi_set = True
    #else:
      #cli.mpi_set = False
  cli.mpi_set = set_bool_argument(cli.mpi_set)
    
  cli.nwalkers = int(cli.nwalkers)
  cli.nruns = int(cli.nruns)
  cli.npost = int(cli.npost)
  #if (cli.good != False):
    #if (cli.good.lower() in ['t', 'tr', 'tru', 'true']):
      #cli.good = True
    #else:
      #cli.good = False
  cli.good = set_bool_argument(cli.good)
  
  #if (cli.devol != False):
    #if (cli.devol.lower() in ['t', 'tr', 'tru', 'true']):
      #cli.devol = True
    #else:
      #cli.devol = False
  cli.devol = set_bool_argument(cli.devol)
  if (cli.devol): cli.ngen = int(cli.ngen)
    
  return cli

def lnprob(fitting_parameters):
  loglhd = 0.
  loglhd, check = pytrades_lib.pytrades.fortran_loglikelihood(fitting_parameters)
  # add stability check
  if (not check or check == 0): loglhd = -np.inf
  return loglhd

def fitness_de(fitting_parameters):
  fitness = 0.
  lglhd, check = pytrades_lib.pytrades.fortran_loglikelihood(fitting_parameters)
  # add stability check
  #if (not check or check == 0): loglhd = -np.inf
  fitness = -2.*lglhd
  return fitness


def init_emcee_folder(working_path, sub_folder):
  emcee_folder = os.path.join(working_path, sub_folder)
  if (not os.path.isdir(emcee_folder)):
      os.makedirs(emcee_folder)
  emcee_log = os.path.join(emcee_folder, "emcee_run.log")
  of_emcee = open(emcee_log, 'w')
  of_emcee.write("# TRADES-EMCEE LOG FILE\n")
  of_emcee.write("# working_path = %s\n" %(working_path))
  of_emcee.write("# emcee_folder = %s\n" %(emcee_folder))
  of_emcee.write("# emcee_log = %s\n" %(emcee_log))
  
  return emcee_folder, emcee_log, of_emcee

#def check_good_parameters(good_id, good_parameters_in):
  #nfit = np.asarray(good_parameters_in).shape[0]
  #good_parameters_out = np.zeros(np.asarray(good_parameters_in).shape[0]) + good_parameters_in
  #for i in range(0, nfit):
    #if (good_id[i].strip()[0] == 'm' and good_id[i].strip()[1] != 'A'):
      #good_parameters_out[i] = good_parameters_in[i] * Mjups
  #return good_parameters_out
    



# MAIN -- TRADES + EMCEE

cli = get_args()

start = time.time()

working_path = cli.full_path
nthreads=cli.nthreads

pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder)

fitting_parameters = pytrades_lib.pytrades.par
parameters_minmax = pytrades_lib.pytrades.parameters_minmax
delta_parameters = np.abs(parameters_minmax[:,1] - parameters_minmax[:,0])

ndata = pytrades_lib.pytrades.ndata
npar  = pytrades_lib.pytrades.npar
nfit  = pytrades_lib.pytrades.nfit
dof   = pytrades_lib.pytrades.dof

n_rv = pytrades_lib.pytrades.nrv
n_set_rv = pytrades_lib.pytrades.nrvset

n_t0 = pytrades_lib.pytrades.nt0
n_t0_sum = np.sum(n_t0)
n_bodies = n_t0.shape[0]
n_set_t0 = 0
for i in range(0, n_bodies):
  if (np.sum(n_t0[i]) > 0): n_set_t0 += 1

# set emcee parameters
if (cli.nwalkers < nfit*2):
  nwalkers = nfit * 2
else:
  nwalkers = cli.nwalkers
if (nwalkers % 2) != 0:
  #print ' Provided odd nwalkers = %d ==> using nwalkers+1 = %d' %(nwalkers, nwalkers+1)
  nwalkers += 1
# differential evolution ngen
if (cli.devol):
  ngen = nwalkers*100
  if (cli.ngen > ngen):
    ngen = cli.ngen
#ngen=10

if (cli.nruns < 1):
  nruns = 10000
else:
  nruns = cli.nruns
  
if (cli.nsave != 'False'):
  if (int(cli.nsave) > 0 and int(cli.nsave) < nruns):
    nsave = int(cli.nsave)
  elif (int(cli.nsave) <= 0):
    nsave = False
  else:
    nsave = nruns/10
else:
  nsave = False

if (cli.npost < 0):
  npost = 1000
else:
  npost = cli.npost
#print nwalkers, nruns, npost

reshaped_names = pytrades_lib.pytrades.parameter_names.reshape((10,nfit), order='F').T
parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]

#print 'MPI SET = ', cli.mpi_set
# initialize MPI - IT MUST BE HERE OR IT RETURNS SEGMENTATION FAULT ...
if (cli.mpi_set):
  pool_mpi = MPIPool(loadbalance=True)
  if not pool_mpi.is_master():
    pool_mpi.wait()
    sys.exit(0)
emcee_folder, emcee_log, of_emcee = init_emcee_folder(working_path, cli.sub_folder)
if (cli.devol): devol_file = get_devol_file(emcee_folder)

print 
print ' ============ '
print ' TRADES+EMCEE'
print ' ============ '
print
print ' WORKING PATH = %s' %(working_path)
print ' NUMBER OF THREADS = %d' %(nthreads)
print ' dof = ndata(%d) - nfit(%d) = %d' %(ndata, nfit, dof)
print ' Total N_RV = %d for %d set(s)' %(n_rv, n_set_rv)
print ' Total N_T0 = %d for %d out of %d planet(s)' %(n_t0_sum, n_set_t0, n_bodies-1)
print ' Differential Evolution npop = %d ngen = %d' %(nwalkers, ngen)
print ' emcee chain: nwalkers = %d nruns = %d' %(nwalkers, nruns)
of_emcee.write(' NUMBER OF THREADS = %d\n' %(nthreads))
of_emcee.write(' dof = ndata(%d) - nfit(%d) = %d\n' %(ndata, nfit, dof))
of_emcee.write(' Total N_RV = %d for %d set(s)\n' %(n_rv, n_set_rv))
of_emcee.write(' Total N_T0 = %d for %d out of %d planet(s)\n' %(n_t0_sum, n_set_t0, n_bodies))
if (cli.devol): of_emcee.write("# Differential Evolution\n npop = %d ngen = %d\n" %(nwalkers, ngen))
of_emcee.write("# emcee chain:\nnwalkers = %d\nnruns = %d\n" %(nwalkers, nruns))

if (not cli.good):
  
  if (cli.devol):
    if (os.path.isfile(devol_file)): os.remove(devol_file)
    of_devol = h5py.File(devol_file, 'a')
    population = of_devol.create_dataset('population', (ngen, nwalkers, nfit), dtype=np.float64)
    of_devol['population'].dims[0].label = 'ngen'
    of_devol['population'].dims[1].label = 'nwalkers'
    of_devol['population'].dims[2].label = 'nfit'
    best_population = of_devol.create_dataset('best_population', (ngen, nfit), dtype=np.float64)
    best_fitness_de = of_devol.create_dataset('best_fitness', (ngen,), dtype=np.float64)
    t_start_de = time.time()
    best_de = np.zeros((ngen, nfit))
    llhd_de = np.zeros((ngen))
    print ' SETTING p0 WITH DIFFERENTIAL EVOLUTION (ngen = %d)' %(ngen)
    #diff_evol = devol(lnprob, parameters_minmax, nwalkers, maximize=True)
    diff_evol = devol(fitness_de, parameters_minmax, nwalkers)
    #diff_evol.optimize(ngen)
    for i_de, res_de in enumerate(diff_evol(ngen)):
      best_de[i_de, :] = diff_evol.minimum_location.copy()
      llhd_de[i_de] = diff_evol.minimum_value.copy()
      population = of_devol['population']
      population[i_de, :, :] = diff_evol.population.copy()
      best_population = of_devol['best_population']
      best_population[i_de, :] = best_de[i_de, :]
      best_fitness_de = of_devol['best_fitness']
      best_fitness_de[i_de] = llhd_de[i_de]
      print ' gen = %6d : best fitness = %.3f <-> parameters: %s' %(i_de, llhd_de[i_de], ' '.join(map(str, best_de[i_de,:])))
    time_de = time.time() - t_start_de
    of_devol.create_dataset('timing_de_seconds', data=np.array([time_de]), dtype=np.float64)
    of_devol.close()
    t_de_d, t_de_h, t_de_m, t_de_s = computation_time(time_de)
    print ' FINISHED DIFFERENTIAL EVOLUTION in %2d day %02d hour %02d min %5.2f sec' %(int(t_de_d), int(t_de_h), int(t_de_m), t_de_s)
    p0 = diff_evol.population
  else:
    print ' SETTING p0 WITH UNIFORM/RANDOM INITIALIZATION'
    p0 = [parameters_minmax[:,0] + np.random.random(nfit)*delta_parameters for i in range(0, nwalkers)]
    
else:
  
  good_file = os.path.join(working_path, 'good_parameters.dat')
  if (os.path.exists(good_file)):
    good_parameters = np.genfromtxt(good_file, usecols=(1), dtype=np.float64)
    good_id = np.genfromtxt(good_file, usecols=(0), dtype='|S10')
    good_parameters = check_good_parameters(good_id, good_parameters)

    good_copy_file = os.path.join(os.path.join(working_path, cli.sub_folder),'good_parameters.dat')
    shutil.copy(good_file, good_copy_file)
    good_lglhd = lnprob(good_parameters)

    p0 = [good_parameters + np.random.randn(nfit)*np.abs(delta_parameters*0.05) for i in range(0, nwalkers)]
    for i in range(0,nwalkers):
      par_temp = p0[i]
      for j in range(0,nfit):
        if (par_temp[j] < parameters_minmax[j,0] or par_temp[j] > parameters_minmax[j,1]):
          if (parameter_names[j][:2] == 'mA'):
            p0[i][j] = p0[i][j]%360.
          else:
            p0[i][j] = parameters_minmax[j,0] + np.random.random(1)*delta_parameters[j]
  
  else:
  
    print ' ERROR: good_parameters.dat not found in working folder, create it or run with --good False'
    sys.exit()

if (cli.good):
  print ' Good: lglhd = %.6f ==> Fitness = %.6f' %(good_lglhd, -2.*good_lglhd)
  #print ' Good: lglhd = %.6f ==> Fitness = %.6f' %(good_lglhd, 1./good_lglhd)
print ' Running emcee with nfit = %d nwalkers = %d nruns = %d npost = %d ...' %(nfit, nwalkers, nruns, npost)
if (cli.mpi_set):
  sampler = emcee.EnsembleSampler(nwalkers, nfit, lnprob, threads=nthreads, pool=pool_mpi)
else:
  sampler = emcee.EnsembleSampler(nwalkers, nfit, lnprob, threads=nthreads)

# first loglikelihood print
#for i in range(0,nwalkers):
  #lnprob_p0 = lnprob(p0[i])
  #print 'walker = %d ==> lnprob = %20.4f' %(i, lnprob_p0)

if (nsave != False):
  # save temporary sampling during emcee every nruns*10%
  if(os.path.exists(os.path.join(emcee_folder, 'emcee_temp.hdf5')) and os.path.isfile(os.path.join(emcee_folder, 'emcee_temp.hdf5'))):
    os.remove(os.path.join(emcee_folder, 'emcee_temp.hdf5'))
  of_temp = h5py.File(os.path.join(emcee_folder, 'emcee_temp.hdf5'), 'a')
  of_temp.create_dataset('parameter_names', data=parameter_names, dtype='S10')
  of_temp.create_dataset('boundaries', data=parameters_minmax, dtype=np.float64)
  temp_dset = of_temp.create_dataset('chains', (nwalkers, nruns, nfit), dtype=np.float64)
  temp_lnprob = of_temp.create_dataset('lnprobability', (nwalkers, nruns), dtype=np.float64)
  of_temp.close()
  pos = p0
  nchains = int(nruns/nsave)
  state=None
  for i in range(0, nchains):
    print
    print ' iter: ',i+1,
    aaa = i*nsave
    bbb = aaa+nsave
    pos, prob, state = sampler.run_mcmc(pos, N=nsave, rstate0=state)
    print 'completed %d steps of %d' %(bbb, nruns)
    of_temp = h5py.File(os.path.join(emcee_folder, 'emcee_temp.hdf5'), 'a')
    temp_dset = of_temp['chains'] #[:,:,:]
    temp_dset[:,aaa:bbb,:] = sampler.chain[:, aaa:bbb, :]
    #of_temp['chains'].attrs['completed_steps'] = bbb
    temp_dset.attrs['completed_steps'] = bbb
    temp_lnprob = of_temp['lnprobability'] #[:,:]
    temp_lnprob[:, aaa:bbb] = sampler.lnprobability[:, aaa:bbb]
    of_temp.close()
  print '...done with saving temporary total shape = ', sampler.chain.shape

# RUN EMCEE AND RESET AFTER REMOVE BURN-IN
#pos, prob, state = sampler.run_mcmc(p0, npost)
#sampler.reset()
#sampler.run_mcmc(pos, nruns, rstate0=state)
else:
  # GOOD COMPLETE SINGLE RUNNING OF EMCEE, WITHOUT REMOVING THE BURN-IN
  sampler.run_mcmc(p0, nruns)
  print 'done'

flatchains = sampler.chain[:, :, :].reshape((-1, nfit)) # full chain values
samples = sampler.chain[:, npost:, :].reshape((-1, nfit)) # discard first npost iterations of the chains

acceptance_fraction = sampler.acceptance_fraction
mean_acceptance_fraction = np.mean(acceptance_fraction)
autocor_time = sampler.acor
lnprobability = sampler.lnprobability

of_emcee.write(" Mean_acceptance_fraction [0.25-0.5] = %.6f\n" %(mean_acceptance_fraction))
print(" Mean acceptance fraction should be between [0.25-0.5]: {0:.3f}".format(mean_acceptance_fraction))

final_parameters = np.zeros((nfit,3))
final_parameters_median = np.percentile(samples, 50, axis=0)
best_fitness = lnprob(final_parameters_median)
print ' median_loglikelihood = ', best_fitness
print ' median_fitness = ', -2.*best_fitness
#print ' best_fitness = ', 1./best_fitness
final_parameters[:,1] = final_parameters_median
final_parameters[:,0] = np.percentile(samples, 16, axis=0)
final_parameters[:,2] = np.percentile(samples, 84, axis=0)

nruns_sel = nruns - npost
# get max lnprobability and associated parameters
# chain is transposed: needed to plot quicker chains for each walker: nruns vs value of parameter
chains_T, parameter_boundaries = select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, 1., parameter_names_emcee, parameter_boundaries, chains)
#chains_T = np.zeros((nruns-npost, nwalkers, nfit))
#for ii in xrange(0,nfit):
  #chains_T[:,:,ii] = sampler.chain[:,npost:nruns,ii].T
 
flatchain_posterior_0 = chains_T[npost:,:,:].reshape((-1, nfit))
#transposed_lnprob = lnprobability[:,npost:].T
#lnprob_shape = transposed_lnprob.shape
#flat_lnprob = transposed_lnprob[:,:].reshape(lnprob_shape[0]*lnprob_shape[1])
#max_lnprob_id = np.argmax(flat_lnprob)
#max_lnprob_parameters = flatchain_posterior_0[max_lnprob_id,:]
#print ' shape of max_lnprob_parameters = ', max_lnprob_parameters.shape
#max_lnprob = flat_lnprob[max_lnprob_id]
#max_lnprob_one_sigma = np.percentile(flatchain_posterior_0-max_lnprob_parameters, 68.3, axis=0)

# GET MAX LNPROBABILITY AND PARAMETERS
# original lnprobability dimensions: (nwalkers, nruns)
lnprob_burnin = lnprobability[:,npost:nruns]
max_lnprob_row, max_lnprob_col = get_max_indices(lnprob_burnin)
max_lnprob = lnprob_burnin[max_lnprob_row, max_lnprob_col]
print ' max_lnprob = %.12f => fitness = %.12f' %(max_lnprob, -2.*max_lnprob)
max_lnprob_parameters = chains_T[max_lnprob_col, max_lnprob_row, :]
max_lnprob_one_sigma = np.percentile(np.abs(flatchain_posterior_0-max_lnprob_parameters), 68.3, axis=0)
max_lnprob_parameters = renormalize_parameters(max_lnprob_parameters, parameter_names_emcee)

# save chains with original shape as hdf5 file
f_hdf5 = h5py.File(os.path.join(emcee_folder, 'emcee_summary.hdf5'), 'w')
f_hdf5.create_dataset('chains', data=sampler.chain, dtype=np.float64)
f_hdf5.create_dataset('parameter_names', data=parameter_names, dtype='S10')
f_hdf5.create_dataset('boundaries', data=parameters_minmax, dtype=np.float64)
f_hdf5.create_dataset('final_parameters', data=final_parameters, dtype=np.float64)
f_hdf5.create_dataset('max_lnprob_parameters', data=max_lnprob_parameters, dtype=np.float64)
f_hdf5.create_dataset('max_lnprob', data=max_lnprob, dtype=np.float64)
f_hdf5.create_dataset('acceptance_fraction', data=acceptance_fraction, dtype=np.float64)
f_hdf5.create_dataset('autocor_time', data=autocor_time, dtype=np.float64)
f_hdf5.create_dataset('lnprobability', data=lnprobability, dtype=np.float64)
f_hdf5.close()

of_emcee.write("# parameters boundaries:\n# parameter_name \t min \t max\n")
for i in range(0, nfit):
  line_log = '%10s %19.14f %19.14f\n' %(parameter_names[i], parameters_minmax[i,0], parameters_minmax[i,1])
  of_emcee.write(line_log)
of_emcee.write('# final summary\n')
of_emcee.write('max_loglikelihood = %.6f\nmin_fitness = %.6f\n'%(max_lnprob, -2.*max_lnprob))
of_emcee.write('median_loglikelihood = %.6f\nmedian_fitness = %.6f\n'%(best_fitness, -2.*best_fitness))
#of_emcee.write('best_loglikelihood = %.6f\nbest_fitness = %.6f\n'%(best_fitness, 1./best_fitness))
of_emcee.write("# final parameters [m -> mass in Msun, P -> period in days, R -> radius in Rsun, w -> arg. per. in degree, mA -> mean anomaly in degree, i -> inclination in degree, lN -> long. of node in degree]:\n")
of_emcee.write("# parameter_name \t value \t neg_err \t pos_err \t autocor_time\n")
for i in range(0, nfit):
  line_log = '%10s %20.14f %19.14f %19.14f %14.8f\n' %(parameter_names[i], final_parameters[i,1], final_parameters[i,1]-final_parameters[i,0], final_parameters[i,2]-final_parameters[i,1], autocor_time[i])
  of_emcee.write(line_log)

del chains_T
del flatchain_posterior_0


# create final files:
# 0_0_xxx files for max lnprobability parameter values
print
print ' MAX LNPROBABILITY PARAMETER VALUES'
fitness_0, lgllhd_0, check_0 = pytrades_lib.pytrades.write_summary_files(0, max_lnprob_parameters)

# 16_0_xxx files for lower parameter values
print
print ' LOWER PARAMETER VALUES'
fitness_16, lgllhd_16, check_16 = pytrades_lib.pytrades.write_summary_files(16, final_parameters[:,0])

# 50_0_xxx files for best parameter values
print
print ' MEDIAN PARAMETER VALUES'
fitness_50, lgllhd_50, check_50 = pytrades_lib.pytrades.write_summary_files(50, final_parameters[:,1])

# 84_0_xxx files for upper parameter values
print
print ' UPPER PARAMETER VALUES'
fitness_84, lgllhd_84, check_84 = pytrades_lib.pytrades.write_summary_files(84, final_parameters[:,2])

del chains_T
del flatchain_posterior_0

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = computation_time(elapsed)

print
print ' TRADES+EMCEE FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s)
of_emcee.write(' TRADES+EMCEE FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye\n' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s))
of_emcee.close()
print 
pytrades_lib.pytrades.deallocate_variables()
if (cli.mpi_set):
  pool_mpi.close()



