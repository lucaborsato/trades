#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
import pytrades_lib
import sys
import time
import shutil
import glob
#from emcee.utils import MPIPool
from constants import Mjups,Rjups
import ancillary as anc

import multiprocessing as mp

deg2rad = np.pi/180.0
rad2deg = 180.0/np.pi

#
def get_args():
  parser = argparse.ArgumentParser(description='TRADES ad hoc grid')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--sub-folder', '--sb', action='store', dest='sub_folder', default='emcee_run', help='Sub-folder name, without full path. Default = emcee_run')
  
  parser.add_argument('-c', '--cpu', '--threads', action='store', dest='nthreads', default=1, type=int, help='Number of threads/cpu to use. Default is 1.')
  
  parser.add_argument('-n', '--ns', '--n-sims', action='store', dest='nsims', default=10, type=int, help='Number of simulations to generate for each solution. Default is 10.')
  
  
  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  
  if(cli.nthreads <= 0 ):
    sys.exit(' INSERTED nthreads <= 0. Stop.')
    
  if(cli.nsims <= 0 ):
    sys.exit(' INSERTED nsims <= 0. Stop.')
    
  
  return cli

def print_both(line, of_log=None):
  print line
  if(of_log is not None):
    of_log.write(line + '\n')
  return

#
# INITIALISE FOLDER AND LOG FILE
#
def init_folder(working_path, sub_folder):
  working_folder = os.path.join(working_path, sub_folder)
  if (not os.path.isdir(working_folder)):
      os.makedirs(working_folder)
  run_log = os.path.join(working_folder, "trades_run.log")
  of_run = open(run_log, 'w')
  
  print_both("# pyTRADES LOG FILE", of_run)
  print_both("# working_path = %s" %(working_path), of_run)
  print_both("# working_folder = %s" %(working_folder), of_run)
  print_both("# run_log = %s" %(run_log), of_run)
  
  return working_folder, run_log, of_run

def copy_src_files(working_path, working_folder, of_run):
  files = glob.glob(os.path.join(working_path, '*.*'))
  for ifile in files:
    print_both(' COPYNG: %s TO %s' %(ifile, working_folder), of_run)
    shutil.copy(ifile, os.path.join(working_folder, ''))
  return

# function called by thread/process that perform trades simulation given the queue
def worker(input_q, full_id, full_fitness, full_check):
  while input_q.qsize() != 0:
    
    id_procs = mp.current_process().name
    
    set_parameters = input_q.get()
    i_sim = set_parameters['id_sim']
    fit_parameters = set_parameters['fitting_parameters']
    nfit = fit_parameters.shape[0]
    #fitness_run, check_run = pytrades_lib.pytrades.fortran_fitness_short(fit_parameters)
    
    all_parameters = set_parameters['all_parameters']
    npar = all_parameters.shape[0]
    
    #print 'all: %s' %(' '.join(['%14.8f' %(all_parameters[ii]) for ii in range(0,npar)]))
    #print 'fit: %s' %(' '.join(['%14.8f' %(fit_parameters[ii]) for ii in range(0,nfit)]))
    fitness_run, check_run = pytrades_lib.pytrades.fortran_fitness_long(all_parameters, fit_parameters)
    
    full_id[i_sim-1] = i_sim
    full_fitness[i_sim-1] = fitness_run
    full_check[i_sim-1] = check_run
    
    #print_both('id: %5d all: %s' %(i_sim, ' '.join(['%20.13f' %(all_parameters[ii]) for ii in range(0,npar)])), of_run)
    #print_both('id: %5d fit: %s' %(i_sim, ' '.join(['%20.13f' %(fit_parameters[ii]) for ii in range(0,nfit)])), of_run)
    
    print_both(' Thread == %s : Done simulation number %6d ==> check = %r, fitness = %14.8f\n' %(id_procs, i_sim, bool(check_run), fitness_run))

# Define the function for checking the process status
def status(proc):
  #print proc.is_alive(), type(proc.is_alive())
  if (proc.is_alive()==True):
    return 'alive'
  elif (proc.is_alive()==False):
    return 'dead'
  else:
    return proc.is_alive()
  
def read_solutions(working_path, of_run):
  file_sols = os.path.join(working_path, 'solutions_B14_H10.dat')
  if(os.path.exists(file_sols)):
    
    o_sol = open(file_sols, 'r')
    lines = o_sol.readlines()
    o_sol.close()
    names, B14, ner_B14, per_B14, H10, ner_H10, per_H10 = [], [], [], [], [], [], []
    for line in lines:
      if(len(line) > 0):
        if(line.strip()[0] != '#'):
          line_split = line.strip().split()
          names.append(line_split[0])
          if(line_split[0][0] == 'm' and line_split[0][1] != 'A'):
            factor = Mjups
          elif(line_split[0][0] == 'R'):
            factor = Rjups
          else:
            factor=np.float64(1.0)
           
          B14.append(    np.float64(line_split[1])*factor)
          ner_B14.append(np.float64(line_split[2])*factor)
          per_B14.append(np.float64(line_split[3])*factor)
          H10.append(    np.float64(line_split[4])*factor)
          ner_H10.append(np.float64(line_split[5])*factor)
          per_H10.append(np.float64(line_split[6])*factor)
          
    print_both(' READ FILE %s' %(file_sols), of_run)
    
    return names, np.array(B14, dtype=np.float64), np.array(ner_B14, dtype=np.float64), np.array(per_B14, dtype=np.float64), np.array(H10, dtype=np.float64), np.array(ner_H10, dtype=np.float64), np.array(per_H10, dtype=np.float64)

  else:
    print_both(' FILE %s NOT EXIST - STOP' %(file_sols), of_run)
    sys.exit()
  
  return

def get_normal(val,sigma):
  if(sigma <= 0.):
    new_val = val
  else:
    new_val = np.random.normal(val, sigma)
  return new_val

def generate_new_parameters(nfit, system_parameters, MR_star, sol_par, ner_sol, per_sol):
  
  npar = system_parameters.shape[0]
  all_parameters = np.array(system_parameters, dtype=np.float64).copy()
  n_sol = sol_par.shape[0]
  # generate new all_parameters
  all_parameters[0] = get_normal(MR_star[0,0], MR_star[0,1]) # Mstar
  all_parameters[1] = get_normal(MR_star[1,0], MR_star[1,1]) # Rstar
  for ii in range(0,n_sol):
    jj=ii+2
    sel_err = np.random.randint(2)
    #print sel_err, sol_par[ii], ner_sol[ii], per_sol[ii],
    if(sel_err == 0): # use the negative error
      all_parameters[jj] = sol_par[ii] - np.abs(get_normal(0., ner_sol[ii]))
    else: # use the positive error
      all_parameters[jj] = sol_par[ii] + np.abs(get_normal(0., per_sol[ii]))
    #print ' ==> ',all_parameters[jj]
  fit_parameters = pytrades_lib.pytrades.init_fit_parameters(all_parameters, n_fit=nfit)
  
  return all_parameters, fit_parameters

def main():
  # MAIN -- TRADES + GRID

  # READ COMMAND LINE ARGUMENTS
  cli = get_args()

  # STARTING TIME
  start_time = time.asctime()
  start = time.time()

  # RENAME 
  working_path = cli.full_path

  # INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
  pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder)

  # RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
  system_parameters = pytrades_lib.pytrades.system_parameters
  fitting_parameters = pytrades_lib.pytrades.fitting_parameters # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
  parameters_minmax = pytrades_lib.pytrades.parameters_minmax # PARAMETER BOUNDARIES
  MR_star = pytrades_lib.pytrades.mr_star
  MR_star_noerr = np.zeros((MR_star.shape))
  MR_star_noerr[:,0] = MR_star[:,0]

  #global n_bodies, n_planets, ndata, npar, nfit, dof, inv_dof
  n_bodies = pytrades_lib.pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
  n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
  ndata = pytrades_lib.pytrades.ndata # TOTAL NUMBER OF DATA AVAILABLE
  npar  = pytrades_lib.pytrades.npar # NUMBER OF TOTAL PARAMATERS ~n_planets X 6
  nfit  = pytrades_lib.pytrades.nfit # NUMBER OF PARAMETERS TO FIT
  dof   = pytrades_lib.pytrades.dof # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
  global inv_dof
  inv_dof = np.float64(1.0 / dof)

  # READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
  reshaped_names = pytrades_lib.pytrades.parameter_names.reshape((10,nfit), order='F').T
  parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]

  # INITIALISE SCRIPT FOLDER/LOG FILE
  working_folder, run_log, of_run = init_folder(working_path, cli.sub_folder)
  copy_src_files(working_path, working_folder, of_run)

  print_both(' STARTING DATE: %s' %(start_time), of_run)

  #all_min = pytrades_lib.pytrades.par_min
  #all_max = pytrades_lib.pytrades.par_max

  nthreads = int(cli.nthreads)

  print_both(' NUMBER OF THREADS = %d' %(nthreads), of_run)
  print_both(' dof = ndata(%d) - nfit(%d) = %d' %(ndata, nfit, dof), of_run)

  sol_names, B14, ner_B14, per_B14, H10, ner_H10, per_H10 = read_solutions(working_path, of_run)
  
  # create base parameter set for B14 and H10 without errors to write output files
  base_all_B14, base_fit_B14 = generate_new_parameters(nfit,system_parameters, MR_star_noerr, B14, np.zeros(len(sol_names)), np.zeros(len(sol_names)))
  fitness_B14, lgllhd_B14, check_B14 = pytrades_lib.pytrades.write_summary_files_long(0, base_all_B14, base_fit_B14)
  
  base_all_H10, base_fit_H10 = generate_new_parameters(nfit,system_parameters, MR_star_noerr, H10, np.zeros(len(sol_names)), np.zeros(len(sol_names)))
  fitness_H10, lgllhd_H10, check_H10 = pytrades_lib.pytrades.write_summary_files_long(1, base_all_H10, base_fit_H10)
  
  #print ' system: %s' %(' '.join(['%14.8f' %(system_parameters[ii]) for ii in range(0,npar)]))
  #print 'fit: %s' %(' '.join(['%14.8f' %(fit_parameters[ii]) for ii in range(0,nfit)]))

  print_both(''                                                                              , of_run)
  print_both(' B14: %s' %(' '.join(['%14.8f' %(base_all_B14[ii]) for ii in range(0,npar)])), of_run)
  print_both(' H10: %s' %(' '.join(['%14.8f' %(base_all_H10[ii]) for ii in range(0,npar)])), of_run)
  print_both(''                                                                              , of_run)

  
  # Initialize queues
  #input_queue  = mp.Queue()
  B14_queue = mp.Queue()
  H10_queue = mp.Queue()
  
  #set_parameters = {'fitting_parameters':grid_temp, 'id_grid':i_grid}
  #input_queue.put(set_parameters)
  for i_sim in range(0, cli.nsims):
    all_B14, fit_B14 = generate_new_parameters(nfit,system_parameters, MR_star, B14, ner_B14, per_B14)
    set_parameters_B14 = {'id_sim':i_sim+1, 'all_parameters':all_B14, 'fitting_parameters':fit_B14}
    B14_queue.put(set_parameters_B14)
    all_H10, fit_H10 = generate_new_parameters(nfit,system_parameters, MR_star, H10, ner_H10, per_H10)
    set_parameters_H10 = {'id_sim':i_sim+1, 'all_parameters':all_H10, 'fitting_parameters':fit_H10}
    H10_queue.put(set_parameters_H10)
    #print_both(' sim: %6d => all_B14: %s' %(i_sim, ' '.join(['%14.8f' %(all_B14[ii])           for ii in range(0,npar)])), of_run)
    #print_both(' sim: %6d => all_H10: %s' %(i_sim' '.join(['%14.8f' %(all_H10[ii])           for ii in range(0,npar)])), of_run)

  print_both(' INITIALISED QUEUE FOR B14 AND H10 SOLUTIONS', of_run)
  
  full_id_B14 = mp.Array('i',range(cli.nsims))
  full_fitness_B14  = mp.Array('d',range(cli.nsims))
  full_check_B14 = mp.Array('i',range(cli.nsims))
  
  full_id_H10 = mp.Array('i',range(cli.nsims))
  full_fitness_H10  = mp.Array('d',range(cli.nsims))
  full_check_H10 = mp.Array('i',range(cli.nsims))

  # B14 SIMULATIONS
  print_both(' RUN SIMULATIONS FOR B14\n', of_run)
  threads = []
  for i_th in range(0, nthreads):
    pp = mp.Process(target=worker, args=(B14_queue, full_id_B14, full_fitness_B14, full_check_B14))
    pp.start()
    threads.append(pp)

  for pp in threads:
    pp.join()
    
  # Check processes status
  for pp in threads:
    print "Process ",pp, " @ ", pp.pid, " is ", status(pp)
  # Wait processes to finish
  while len(threads) != 0:
    for pp in threads:
      if status(pp) == "dead":
        threads.pop(threads.index(pp))
    # loose 10 seconds before checking again
    time.sleep(5)
  
  print_both(' DONE SIMULATIONS FOR B14', of_run)

  np.savetxt(os.path.join(working_folder, 'simulations_B14.dat'), np.column_stack((full_id_B14, full_fitness_B14, full_check_B14)), fmt='%6d %23.16e %2d', header='id_sim fitness check\n%6d %23.16e %2d' %(0,fitness_B14,check_B14))
  print_both(' SAVED SIMULATIONS INTO FILE %s\n\n' %(os.path.join(working_folder, 'simulations_B14.dat')), of_run)
  
  
  # B14 SIMULATIONS
  print_both(' RUN SIMULATIONS FOR BH10\n', of_run)
  threads = []
  for i_th in range(0, nthreads):
    pp = mp.Process(target=worker, args=(H10_queue, full_id_H10, full_fitness_H10, full_check_H10))
    pp.start()
    threads.append(pp)

  for pp in threads:
    pp.join()
    
  # Check processes status
  for pp in threads:
    print "Process ",pp, " @ ", pp.pid, " is ", status(pp)
  # Wait processes to finish
  while len(threads) != 0:
    for pp in threads:
      if status(pp) == "dead":
        threads.pop(threads.index(pp))
    # loose 10 seconds before checking again
    time.sleep(5)
  
  print_both(' DONE SIMULATIONS FOR H10', of_run)
  
  np.savetxt(os.path.join(working_folder, 'simulations_H10.dat'), np.column_stack((full_id_H10, full_fitness_H10, full_check_H10)), fmt='%6d %23.16e %2d', header='id_sim fitness check\n%6d %23.16e %2d' %(0,fitness_H10,check_H10))
  print_both(' SAVED SIMULATIONS INTO FILE %s\n\n' %(os.path.join(working_folder, 'simulations_H10.dat')), of_run)
  

  pytrades_lib.pytrades.deallocate_variables()
  print_both(' THE END', of_run)
  of_run.close()
  
  return

if __name__ == "__main__":
  main()
