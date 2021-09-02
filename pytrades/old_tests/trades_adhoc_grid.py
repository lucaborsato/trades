#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
import pytrades_lib
import sys
import time
import shutil
#from emcee.utils import MPIPool
from constants import Mjups, deg2rad, rad2deg
import ancillary as anc

# deg2rad = np.pi/180.0
# rad2deg = 180.0/np.pi

#
def get_args():
  parser = argparse.ArgumentParser(description='TRADES ad hoc grid')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--sub-folder', '--sb', action='store', dest='sub_folder', default='emcee_run', help='Sub-folder name, without full path. Default = emcee_run')
  
  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')

  return cli

def print_both(line, of_log=None):
  print(line)
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

def generate_inc_same_size(n_single, all_min, all_max):
  n_ic = n_single
  n_id = n_single
  n_grid = n_ic*n_id
  s_ic, f_ic = 16, 10
  s_id, f_id = 24, 17
  inc_c = np.linspace(start=all_min[s_ic], stop=all_max[s_ic], num=n_ic, endpoint=True)
  inc_d = np.linspace(start=all_min[s_id], stop=all_max[s_id], num=n_id, endpoint=True)
  return inc_c, inc_d, n_ic, n_id, n_grid

def generate_inc_delta_inc(d_ic, d_id, all_min, all_max):
  s_ic, f_ic = 16, 10
  s_id, f_id = 24, 17
  inc_c = np.arange(start=all_min[s_ic], stop=all_max[s_ic]+0.5*d_ic, step=d_ic)
  inc_d = np.arange(start=all_min[s_id], stop=all_max[s_id]+0.5*d_id, step=d_id)
  n_ic = inc_c.shape[0]
  n_id = inc_d.shape[0]
  n_grid = n_ic * n_id
  return inc_c, inc_d, n_ic, n_id, n_grid

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
  pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder, 1)

  # RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
  system_parameters = pytrades_lib.pytrades.system_parameters
  fitting_parameters = pytrades_lib.pytrades.fitting_parameters # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
  parameters_minmax = pytrades_lib.pytrades.parameters_minmax # PARAMETER BOUNDARIES

  #global n_bodies, n_planets, ndata, npar, nfit, dof, inv_dof
  n_bodies = pytrades_lib.pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
  n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
  ndata = pytrades_lib.pytrades.ndata # TOTAL NUMBER OF DATA AVAILABLE
  npar  = pytrades_lib.pytrades.npar # NUMBER OF TOTAL PARAMATERS ~n_planets X 6
  nfit  = pytrades_lib.pytrades.nfit # NUMBER OF PARAMETERS TO FIT
  dof   = pytrades_lib.pytrades.dof # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
  global inv_dof, inv_dof_cor
  inv_dof = np.float64(1.0 / dof)
  inv_dof_cor = np.float64(1.0 / (dof+4.0))



  # READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
  reshaped_names = pytrades_lib.pytrades.parameter_names.reshape((10,nfit), order='F').T
  parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]

  # INITIALISE SCRIPT FOLDER/LOG FILE
  working_folder, run_log, of_run = init_folder(working_path, cli.sub_folder)

  print_both(' STARTING DATE: %s' %(start_time), of_run)
  print_both(' RUN THE INITIAL GUESS fitting_parameters WITH ID = 0', of_run)
  fitness_run0, lgllhd_run0, check_run0 = pytrades_lib.pytrades.write_summary_files(0, fitting_parameters)
  fitness_run0 = fitness_run0*inv_dof_cor/inv_dof
  print_both('fitness = %.8f lgllhd = %.8f check = %r' %(fitness_run0, lgllhd_run0, bool(check_run0)), of_run)
  print_both(' '.join([str(fitting_parameters[i]) for i in range(len(fitting_parameters))]), of_run)
  print_both(' GO ON WITH THE GRID', of_run)

  # WE WANT TO GENERATE A GRID ON INC_C AND INC_D BETWEEN THE BOUNDARIES
  s_ic, f_ic = 16, 10
  s_id, f_id = 24, 17
  all_min = pytrades_lib.pytrades.par_min
  all_max = pytrades_lib.pytrades.par_max

  #n_single=45
  #inc_c, inc_d, n_ic, n_id, n_grid = generate_inc_same_size(n_single, all_min, all_max)

  d_ic, d_id = 0.25, 0.5
  inc_c, inc_d, n_ic, n_id, n_grid = generate_inc_delta_inc(d_ic, d_id, all_min, all_max)

  # prepare lists to write summary
  full_inc_c, full_inc_d, full_fitness, full_check = [], [], [], []

  #print_both(' NUMBER OF THREADS = %d' %(nthreads), of_run)
  print_both(' dof = ndata(%d) - nfit(%d) = %d' %(ndata, nfit, dof), of_run)
  print_both(' n_grid = %d = n_ic x n_id = %d x %d' %(n_grid, n_ic, n_id), of_run)
  print_both(' INITIALIZING GRID', of_run)

  i_grid = 0
  #grid_fit = np.zeros((n_grid, nfit))
  for iic in range(0,n_ic):
    for iid in range(0,n_id):
      i_grid += 1
      print_both(' Setting grid parameters at id = %d (%d , %d): (%14.8f , %14.8f) ' %(i_grid, iic+1, iid+1, inc_c[iic], inc_d[iid]), of_run)
      
      grid_temp = np.zeros((nfit)) + fitting_parameters
      
      #grid_temp[f_ic]   = inc_c[iic]*np.cos(system_parameters[s_ic+1]*deg2rad)
      ##print_both('icoslN_c = %.8f' %(inc_c[iic]*np.cos(system_parameters[s_ic+1]*deg2rad)),of_run)
      
      #grid_temp[f_ic+1] = inc_c[iic]*np.sin(system_parameters[s_ic+1]*deg2rad)
      ##print_both('isinlN_c = %.8f' %(inc_c[iic]*np.sin(system_parameters[s_ic+1]*deg2rad)),of_run)
      
      #grid_temp[f_id]   = inc_d[iid]*np.cos(system_parameters[s_id+1]*deg2rad)
      ##print_both('icoslN_d = %.8f' %(inc_c[iid]*np.cos(system_parameters[s_id+1]*deg2rad)),of_run)
      
      #grid_temp[f_id+1] = inc_d[iid]*np.sin(system_parameters[s_id+1]*deg2rad)
      ##print_both('isinlN_d = %.8f' %(inc_c[iid]*np.sin(system_parameters[s_id+1]*deg2rad)),of_run)
      
      system_temp = system_parameters.copy()
      system_temp[s_ic] = inc_c[iic]
      system_temp[s_id] = inc_d[iid]
      
      set_parameters = {'fitting_parameters':grid_temp, 'id_grid':i_grid}
      
      # GOOD IN SERIAL
      #fitness_run, lgllhd_run, check_run = pytrades_lib.pytrades.write_summary_files(i_grid, grid_temp)
      #lgllhd_run, check_run = pytrades_lib.pytrades.fortran_loglikelihood(grid_temp)
      ##fitness_run = -2.0*lgllhd_run*inv_dof
      #fitness_run = -2.0*lgllhd_run*inv_dof_cor
      
      fitness_run, check_run = pytrades_lib.pytrades.fortran_fitness_long(system_temp, grid_temp)
      
      # save final lists
      full_inc_c.append(inc_c[iic])
      full_inc_d.append(inc_d[iid])
      full_fitness.append(fitness_run*inv_dof)
      full_check.append(check_run)
      
      print_both(' Done simulation number %d ==> check = %r, fitness = %.8f\n' %(i_grid, bool(check_run), fitness_run), of_run)

  print_both(' GRID COMPLETED', of_run)

  of_fin = open(os.path.join(working_folder, 'final_summary.dat'), 'w')
  of_fin.write('# iteration inc_c inc_d fitness check[0==False, 1==True]\n')
  of_fin.write('%4d %17.13f %17.13f %23.16e %2d\n' %(0, system_parameters[s_ic], system_parameters[s_id], fitness_run0, check_run0))

  for ii in range(0, n_grid):
    line = '%4d %17.13f %17.13f %23.16e %2d\n' %(ii+1, full_inc_c[ii], full_inc_d[ii], full_fitness[ii], full_check[ii])
    of_fin.write(line)
  of_fin.close()

  print_both(' WRITTEN SUMMARY FILE: %s' %(os.path.join(working_folder, 'final_summary.dat')), of_run)


  pytrades_lib.pytrades.deallocate_variables()
  print_both(' THE END', of_run)
  of_run.close()
  
  return

if __name__ == "__main__":
  main()
