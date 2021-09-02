#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
import sys
# import h5py

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# matplotlib rc params
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.figsize'] = [6, 6]
plt.rcParams["figure.facecolor"] = 'white'
plt.rcParams["savefig.facecolor"] = 'white'
plt.rcParams["figure.dpi"]  = 200
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"]   = 14


# local constants module
# import constants as cst
import ancillary as anc

# read command line (cli) arguments
def get_args():
  parser = argparse.ArgumentParser(description='PSO PLOT')
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files.')
  parser.add_argument('-m', '--mtype', '--mass-type', action='store', dest='m_type', default='e', help='Mass type: j = Jupiter, e = Earth, s = Sun. Default is Earth = e.')
  
  cli = parser.parse_args()
  cli.full_path = os.path.abspath(cli.full_path)
  cli.m_type = str(cli.m_type).lower()
  
  return cli

# ======================================================================

def main():

  print() 
  print(' ================== ')
  print(' PSO PLOTS')
  print(' ================== ')
  print()

  # read cli arguments
  cli = get_args()

  # set pso_file
  pso_file = os.path.join(cli.full_path, 'pso_run.hdf5')
  population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution, parameters_minmax, parameter_names, pop_shape = anc.get_pso_data(pso_file)
  nfit = pop_shape[0]
  npop = pop_shape[1]
  niter = pop_shape[2]

  nwalkers = npop

  mass_star = 1.0
  m_factor = 1.0
  m_factor, m_unit = anc.mass_type_factor(
    Ms=mass_star, 
    mtype=cli.m_type,
    mscale=False
  )

  iteration = np.arange(0,niter)+1
  if(isinstance(parameters_minmax,type(population_fitness))):
    parameters_minmax_bck = parameters_minmax.copy()

  # set label and legend names
  kel_legends = anc.keplerian_legend(parameter_names, cli.m_type)

  anc.print_memory_usage(population)
  anc.print_memory_usage(population_fitness)
  fit_min = np.min(population_fitness)
  fit_max = fit_min + 100
  print('population fitness: min = {} (set max = min + 100 = {}'.format(fit_min, fit_max))

  pso_plots = os.path.join(cli.full_path,'plots')
  if (not os.path.isdir(pso_plots)):
    os.makedirs(pso_plots)

  # parameter_names and parameters_minmax in pso_run.hdf5
  for ii in range(0, nfit):
    print('parameter: {:s}'.format(parameter_names[ii]))
    if (parameter_names[ii][0] == 'm' and parameter_names[ii][1] != 'A'):
      pop_plt = population[ii,:,:]*m_factor
      y_min = parameters_minmax[ii,0]*m_factor
      y_max = parameters_minmax[ii,1]*m_factor
    else:
      pop_plt = population[ii,:,:]
      y_min = parameters_minmax[ii,0]
      y_max = parameters_minmax[ii,1]
      
    print('boundaries: [{:.6f}, {:.6f}]'.format(y_min, y_max))  
    print('    minmax: [{:.6f}, {:.6f}]'.format(np.min(population[ii,:,:]), np.max(population[ii,:,:])))
    pso_fig_file = os.path.join(pso_plots, 'pso_evolution_{:s}.png'.format(parameter_names[ii]))
    print(' {:s}'.format(pso_fig_file), end=' ')

    fig = plt.figure(figsize=(5,5))
    for jj in range(0, npop):
      plt.plot(iteration, pop_plt[jj,:],
        marker='o',
        mfc=(0.8, 0.8, 0.8, 0.8),
        mec='None',
        ls='',
        ms=2,
        zorder=6
      )
      # # SCATTER plot
      # sc = plt.scatter(iteration, pop_plt[jj,:],
      #   s=1,
      #   c=population_fitness[jj,:],
      #   marker='.',
      #   vmin=fit_min, vmax=fit_max,
      #   alpha=0.45
      # )

    plt.ylim(y_min, y_max)
  
    if (parameter_names[ii][0] == 'm' and parameter_names[ii][1] != 'A'):
      par_plt = pso_best_evolution[ii,:]*m_factor
    else:
      par_plt = pso_best_evolution[ii,:]

    # plt.plot(iteration, par_plt, 
    #   color='C0',
    #   marker='o', 
    #   mfc='C0', 
    #   mec='white', 
    #   mew=0.1, 
    #   ls='-', 
    #   ms=3,
    #   alpha=0.45
    # )

    # sel_min = np.argmin(population_fitness, axis=0)
    par_fit = np.min(population_fitness, axis=0)
    sc = plt.scatter(iteration, par_plt,
      s=2,
      c=par_fit,
      marker='o',
      vmin=fit_min, vmax=fit_max,
      alpha=1.0,
      zorder=7
    )

    plt.colorbar(sc)

    plt.xlabel(r'$N_\mathrm{iteration}$')
    plt.ylabel(kel_legends[ii])
    plt.draw()
    fig.savefig(pso_fig_file, bbox_inches='tight')
    print(' done')
  
  return

if __name__ == "__main__":
  main()
