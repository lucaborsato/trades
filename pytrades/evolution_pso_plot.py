#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import sys
import argparse
import os
import numpy as np # array
import h5py
import constants as cst # local constants module
import ancillary as anc
from matplotlib import use as mpluse
mpluse("Agg")
#mpluse("Qt4Agg")
import matplotlib.pyplot as plt
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

#from matplotlib import rcParams
#rcParams['text.latex.unicode']=True


# read command line (cli) arguments
def get_args():
  parser = argparse.ArgumentParser(description='TRADES+EMCEE PLOT')
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files.')
  parser.add_argument('-m', '--mtype', '--mass-type', action='store', dest='m_type', default='e', help='Mass type: j = Jupiter, e = Earth, s = Sun. Default is Earth = e.')
  
  cli = parser.parse_args()
  cli.full_path = os.path.abspath(cli.full_path)
  cli.m_type = cli.m_type.lower()
  
  return cli

def main():

  print() 
  print(' ================== ')
  print(' PSO PLOTS')
  print(' ================== ')
  print()

  # read cli arguments
  cli = get_args()
  # computes mass conversion factor
  #m_factor = mass_conversion_factor(cli.m_type)
  m_factor, m_unit = mass_type_factor(1., cli.mtype, False)

  # set pso_file
  pso_file = os.path.join(cli.full_path, 'pso_run.hdf5')
  population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution, parameters_minmax, parameter_names, pop_shape = anc.get_pso_data(pso_file)
  nfit = pop_shape[0]
  npop = pop_shape[1]
  niter = pop_shape[2]

  iteration = np.arange(0,niter)+1
  if(isinstance(parameters_minmax,type(population_fitness))):
    parameters_minmax_bck = parameters_minmax.copy()

  # set label and legend names
  kel_legends, labels_list = anc.keplerian_legend(parameter_names, cli.m_type)

  anc.print_memory_usage(population)
  anc.print_memory_usage(population_fitness)

  pso_plots = os.path.join(cli.full_path,'plots')
  if (not os.path.isdir(pso_plots)):
    os.makedirs(pso_plots)


  # parameter_names and parameters_minmax in pso_run.hdf5
  if(isinstance(parameter_names, type(population)) and isinstance(parameters_minmax,type(population_fitness))):
    for ii in range(0, nfit):
      print('parameter: %s' %(parameter_names[ii]))
      if (parameter_names[ii][0] == 'm' and parameter_names[ii][1] != 'A'):
        population[ii,:,:] = population[ii,:,:]*m_factor
        y_min = parameters_minmax[ii,0]*m_factor
        y_max = parameters_minmax[ii,1]*m_factor
      else:
        y_min = parameters_minmax[ii,0]
        y_max = parameters_minmax[ii,1]
        
      print('boundaries: [%.6f, %.6f]' %(y_min, y_max))  
      print('    minmax: [%.6f, %.6f]' %(np.min(population[ii,:,:]), np.max(population[ii,:,:])))
      pso_fig_file = os.path.join(pso_plots, 'evolution_%s.png' %(parameter_names[ii]))
      print(' %s' %(pso_fig_file), end=' ')
      fig = plt.figure(figsize=(12,12))
      for jj in range(0, npop):
        #print jj,
        plt.plot(iteration, population[ii,jj,:], marker='o', mfc='gray', mec='none', ls='', ms=4)
      plt.ylim(y_min, y_max)
    
      if(isinstance(pso_best_evolution, type(population_fitness))):
        if (parameter_names[ii][0] == 'm' and parameter_names[ii][1] != 'A'):
          plt.plot(iteration, pso_best_evolution[ii,:]*m_factor, marker='o', mfc='black', mec='white', mew=0.25, ls='-', ms=5)
          #print pso_best_evolution[-1,0],pso_best_evolution[-1,-1]
        else:
          plt.plot(iteration, pso_best_evolution[ii,:], marker='o', mfc='black', mec='white', mew=0.25, ls='-', ms=5)
    
      plt.xlabel('$N_\mathrm{iteration}$')
      plt.ylabel(kel_legends[ii])
      plt.draw()
      fig.savefig(pso_fig_file, bbox_inches='tight', dpi=150)
      print(' done')
  #elif ():
  
  return

if __name__ == "__main__":
  main()
