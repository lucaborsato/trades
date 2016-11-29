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

# given folder it read best pso or pik results of each global simulations

# read arguments
def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-p', action='store', dest='fpath', required=True, help='Folder path')
  #parser.add_argument('-s', action='store', dest='idsim', required=True, help='Simulation ID number')
  
  cli = parser.parse_args()
  
  return os.path.abspath(cli.fpath)

# check all the global Simulations and create the list
def create_global_list(full_path):
  temp_files = '*_best*par.out'
  best_list = sorted(glob.glob(os.path.join(full_path, temp_files)))
  return best_list
  
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
  return file_name, iteration, fitness

def main():
  # MAIN

  full_path = get_args()
  best_list =  create_global_list(full_path)
  for best_file in best_list:
    file_name, iteration, fitness = read_last_results(best_file)
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

if __name__ == "__main__":
  main()
  
  
  
