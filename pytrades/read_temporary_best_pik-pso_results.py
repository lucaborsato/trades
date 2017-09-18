#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os
#import time
import glob
import sys
#import datetime
import numpy as np # array

# custom modules
script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), '../'))
sys.path.append(module_path)
import constants as cst
import ancillary as anc
from pytrades_lib import pytrades

# given folder it read best pso or pik results of each global simulations

# ==============================================================================
# read arguments
def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-p', '--p', '-path-folder', '--path-folder',
                      action='store', dest='fpath', required=True,
                      help='Folder path')
  parser.add_argument('-i', '--i', '-integrates', '--integrates',
                      action='store', dest='integrates', default='False',
                      help='Perform orbit integration of the temporary best solution. Set True/Yes or False/No. Default False.')
  parser.add_argument('-m', '--m', '-mass-type', '--mass-type',
                      action='store', dest='mass_type', default='e',
                      help='Needed only if --integrates is True/Yes. Provide mas type of further output: e/E/earth/Earth j/J/jup/Jup n/N/nep/Nep . Default Earth.')
  
  cli = parser.parse_args()
  cli.fpath = os.path.abspath(cli.fpath)
  if(str(cli.integrates).lower()[0] in ['t', 'y']):
    cli.integrates = True
  else:
    cli.integrates = False
    
  if(cli.mass_type.lower() in 'j ju jup jupi jupit jupite jupiter'.split()):
    cli.mass_type = 'j'
  elif(cli.mass_type.lower() in 'n ne nep nept neptu neptun neptune'.split()):
    cli.mass_type = 'n'
  else:
    cli.mass_type = 'e'
    
  return cli

# ==============================================================================
# check all the global Simulations and create the list
def create_global_list(full_path):
  temp_files = '*_best*par.out'
  best_list = sorted(glob.glob(os.path.join(full_path, temp_files)))
  return best_list
  
# ==============================================================================
# read last line of the file, keep first and the last column [iteration - fitness]
def read_last_results(best_file):
  of = open(best_file, 'r')
  line_empty = ['EMPTY', 'EMPTY']
  for line in of:
    pass
  try:
    line = line.strip().split()
  except:
    line = line_empty
  iteration = line[0]
  fitness = line[-1]
  file_name = os.path.basename(best_file)
  of.close()
  return file_name, iteration, fitness

# ==============================================================================
def read_temp_results(best_file):

  line_empty = ['EMPTY', 'EMPTY']
  file_name = os.path.basename(best_file)

  of = open(best_file, 'r')
  lines = of.readlines()
  of.close()
  nlines = len(lines)
  
  header = lines[0]
  par_names = header.split('part')[1].split('inv_fitness')[0].split()
  nfit = len(par_names)

  if (nlines <= 1):
    iteration = 'EMPTY'
    fitness   = 'EMPTY'
    par_best = 'EMPTY'
  else:
    last_row = lines[-1].split()
    iteration = last_row[0]
    fitness = last_row[-1]
    par_best = [np.float64(last_row[ifit]) for ifit in range(2,2+nfit)]
    
  return file_name, iteration, fitness, par_names, par_best

# ==============================================================================
def do_simulation(best_file, iteration_str, par_names, par_best, mass_type):

  # prepare new out_folder
  iteration = int(iteration_str)
  pso_folder = os.path.dirname(best_file)
  out_folder = os.path.join(pso_folder, '%s_best' %(iteration_str))
  out_folder = os.path.join(out_folder, '')
  if (not os.path.isdir(out_folder)):
      os.makedirs(out_folder)
  # init trades
  pytrades.initialize_trades(os.path.join(pso_folder, ''), '', 1)
  # change path to new folder to save trades outputs
  pytrades.path_change(out_folder)
  # convert properly the parameters
  par_trades = anc.sqrte_to_e_fitting(par_best, par_names)
  # run trades and save outputs
  fitness, lgllhd, check = pytrades.write_summary_files(iteration, par_trades)
  # warning if issues exist
  if(not bool(check)):
    print 'WRTING WARNING FILE: %s' %(os.path.join(out_folder,'WARNING.txt'))
    warn_o = open(os.path.join(out_folder,'WARNING.txt'), 'w')
    warn_o.write('*******\nWARNING: FITTED PARAMETERS COULD NOT BE PHYSICAL!\nWARNING: BE VERY CAREFUL WITH THIS PARAMETER SET!\n*******')
    warn_o.close()
  
  # read elements file and write it in selected mass type
  kel_file, kep_elem = anc.elements(out_folder, iteration, lmf=0)
  if(mass_type == 'j'):
    mass_out = 'M_Mjup'
    mass_conv = cst.Msjup
  elif(mass_type == 'n'):
    mass_out = 'M_Mnep'
    mass_conv = cst.Msnep
  else:
    mass_out = 'M_Mear'
    mass_conv = cst.Msear
  kel_file_mass = kel_file.replace('.dat', '%s.dat' %(mass_out))
  kep_header = '%s R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg' %(mass_out)
  kep_elem[:,0] = kep_elem[:,0] * mass_conv
  np.savetxt(kel_file_mass, kep_elem, fmt='%23.16e', header=kep_header)
  print
  
  return

# ==============================================================================
def main():
  # MAIN

  cli = get_args()
  
  best_list =  create_global_list(cli.fpath)
  
  for best_file in best_list:
    
    #file_name, iteration, fitness = read_last_results(best_file)
    file_name, iteration, fitness, par_names, par_best = read_temp_results(best_file)
    
    if(cli.integrates):
      do_simulation(best_file, iteration, par_names, par_best, cli.mass_type)
    
    line = ' File %s => iteration = %5s fitness = %s' %(file_name, iteration, fitness)
    if(fitness != 'EMPTY'):
      evalfit = np.asarray(fitness, dtype=np.float64)
      if(evalfit > 3. and evalfit <= 10.): # line = '%s ! (<10.)' %(line)
        line = '%s ! (<10.)' %(line)
      if(evalfit > 2. and evalfit <= 3.): # line = '%s !! (<3.)' %(line)
        line = '%s !! (<3.)' %(line)
      if(evalfit > 1. and evalfit <= 2.): # line = '%s !!! (<2.)' %(line)
        line = '%s !!! (<2.)' %(line)
      if(evalfit <= 1.): # line = '%s !!!! (<1.)' %(line)
        line = '%s !!!! (<1.)' %(line)
    print line
  
  return


# ==============================================================================
if __name__ == "__main__":
  main()
  
  
  
