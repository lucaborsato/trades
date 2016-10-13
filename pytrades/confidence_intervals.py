#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import sys
import argparse
import time
import os
import numpy as np # array
import h5py
import constants as cst # local constants module
from ancillary import *
from matplotlib import use as mpluse
mpluse("Agg")
#mpluse("Qt4Agg")
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)
#from matplotlib import rcParams
#rcParams['text.latex.unicode']=True

# read command line (cli) arguments
def get_args():
  parser = argparse.ArgumentParser(description='pytrades confidence intervals')
 
  #parser.add_argument('-s', '--sol', '--solution', action='store', dest='solution_file', default='False', help='The file path (absolute or relative) with solution file from trades: #_#_final#par*.dat')
  
  #parser.add_argument('-p', '--pso', '--pso-run-file', action='store', dest='pso_file', default='False', help='The file path (absolute or relative) with pso_run.hdf5 file with simulation data.')
 
  parser.add_argument('-e', '--emcee', '--emcee-file', action='store', dest='emcee_file', required=True, help='The file path (absolute or relative) of the trades-emcee hdf5 file: emcee_summary.hdf5 or emcee_temp.hdf5')
  
  parser.add_argument('-nb', '--nburn', '-np', '--nburnin', action='store', dest='nburnin', required=True, help='The number of posterior/burn in steps to discard at the beginning of each chain. It has to be > 0')
  
  parser.add_argument('-mi', '--mtype-in', '--mass-type-input', action='store', dest='m_type_in', default='j', help='Mass type: j = Jupiter, e = Earth, s = Sun. Default is Jupiter = j.')
  
  parser.add_argument('-mo', '--mtype-out', '--mass-type-output', action='store', dest='m_type_out', default='j', help='Mass type: j = Jupiter, e = Earth, s = Sun. Default is Jupiter = j.')
  
  #parser.add_argument('-o', '--output', '--output-file', action='store', dest='output_file', required=True, help='Output file with summary of the solution and confidence intervals with the specified mass type (default mass type is Jupiter-mass)')
  
  #parser.add_argument('-dof', '--dof', '--degree-of-freedom', action='store', dest='dof', required=True, help='Degrees of freedom = dof.')
  
  cli = parser.parse_args()

  #cli.solution_file = os.path.abspath(cli.solution_file)

  cli.emcee_file = os.path.abspath(cli.emcee_file)
  
  cli.m_type_in = cli.m_type_in.lower()
  cli.m_type_out = cli.m_type_out.lower()

  #cli.output_file = os.path.abspath(cli.output_file)

  return cli


# main

print
print ' TRADES: PSO-EMCEE confidence intervals'
print

start = time.time()

cli = get_args()

output_file = os.path.join(os.path.dirname(cli.emcee_file), 'emcee_confidence_intervals.txt')
of_out = open(output_file, 'w')

#dof = int(cli.dof)
m_factor_out, m_unit_out = mass_factor_unit(cli.m_type_out) # sun to mass out

# READ EMCEE FILE

if('temp' in os.path.basename(cli.emcee_file)):
  temp_status = True
else:
  temp_status = False
print ' temp_status = ',temp_status

print
parameter_names_emcee, parameter_boundaries_in, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = get_data(cli.emcee_file, temp_status)

print
print ' chains size:'
print_memory_usage(chains)

nfit, nwalkers, nruns, nburnin, nruns_sel = get_emcee_parameters(chains, temp_status, cli.nburnin, completed_steps)

chains_T, parameter_boundaries_out = select_transpose_convert_chains(nfit, nwalkers, nburnin, nruns, nruns_sel, m_factor_out, parameter_names_emcee, parameter_boundaries_in, chains)
print ' chains_T size:'
print_memory_usage(chains_T)

flatchain_posterior_0 = chains_T[:,:,:].reshape((nruns_sel*nwalkers, nfit))

derived_names, derived_chains_T, derived_posterior = get_derived_posterior_parameters(parameter_names_emcee, chains_T, flatchain_posterior_0)
nder = len(derived_names)

print

big_header = '# %s || nruns = %d nburnin = %d nwalkers = %d nfit = %d nder = %d' %(output_file, nruns, nburnin, nwalkers, nfit, nder)
print big_header
of_out.write(big_header + '\n\n')


# GET MAX LNPROBABILITY AND PARAMETERS -> id 40XX
# original lnprobability dimensions: (nwalkers, nruns)
max_lnprob, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint = get_maxlnprob_parameters(nburnin, nruns, lnprobability, chains_T, flatchain_posterior_0)
max_lnprob_der, max_lnprob_parameters_der, max_lnprob_perc68_der, max_lnprob_confint_der = get_maxlnprob_parameters(nburnin, nruns, lnprobability, derived_chains_T, derived_posterior)
print '# MAX LNPROBABILITY PARAMETER VALUES -> 4000'
of_out.write('# MAX LNPROBABILITY PARAMETER VALUES -> 4000\n')
print_parameters_logtxt(of_out, parameter_names_emcee, max_lnprob_parameters, max_lnprob_perc68, max_lnprob_confint, 'maxlnpr')
print_parameters_logtxt(of_out, derived_names, max_lnprob_parameters_der, max_lnprob_perc68_der, max_lnprob_confint_der, 'maxlnpr_der')
print '# ' + '-'*220
of_out.write('# ' + '-'*220 + '\n')

# MEDIAN PARAMETERS
# std way: median of the posterior parameter distribution -> id 10XX
median_parameters, median_perc68, median_confint = get_median_parameters(flatchain_posterior_0)
median_parameters_der, median_perc68_der, median_confint_der = get_median_parameters(derived_posterior)
print '# MEDIAN PARAMETER VALUES -> 1050'
of_out.write('# MEDIAN PARAMETER VALUES -> 1050\n')
print_parameters_logtxt(of_out, parameter_names_emcee, median_parameters, median_perc68, median_confint, 'median')
print_parameters_logtxt(of_out, derived_names, median_parameters_der, median_perc68_der, median_confint_der, 'median_der')
print '# ' + '-'*220
of_out.write('# ' + '-'*220 + '\n')

# parameters linked to the median of the fitness =  - 2 * (lglhd + ln_err_const) -> id 20XX
median_fitness, medfit_parameters, medfit_perc68, medfit_confint = get_parameters_median_fitness(nwalkers, nburnin, nruns, lnprobability, flatchain_posterior_0, ln_err_const)
median_fitness_der, medfit_parameters_der, medfit_perc68_der, medfit_confint_der = get_parameters_median_fitness(nwalkers, nburnin, nruns, lnprobability, derived_posterior, ln_err_const)
print '# MEDFIT PARAMETER VALUES -> 2050'
of_out.write('# MEDFIT PARAMETER VALUES -> 2050\n')
print_parameters_logtxt(of_out, parameter_names_emcee, medfit_parameters, medfit_perc68, medfit_confint, 'medfit')
print_parameters_logtxt(of_out, derived_names, medfit_parameters_der, medfit_perc68_der, medfit_confint_der, 'medfit_der')
print '# ' + '-'*220
of_out.write('# ' + '-'*220 + '\n')

# MODE-LIKE PARAMETERS -> id 30XX
# take the mean of 5 bin centered to the higher bin
k = 40
mode_bin, mode_parameters, mode_perc68, mode_confint = get_mode_parameters(flatchain_posterior_0, k)
mode_bin_der, mode_parameters_der, mode_perc68_der, mode_confint_der = get_mode_parameters(derived_posterior, k)
print '# MODE PARAMETER VALUES -> 3050'
of_out.write('# MODE PARAMETER VALUES -> 3050\n')
print_parameters_logtxt(of_out, parameter_names_emcee, mode_parameters, mode_perc68, mode_confint, 'mode')
print_parameters_logtxt(of_out, derived_names, mode_parameters_der, mode_perc68_der, mode_confint_der, 'mode_der')
print '# ' + '-'*220
of_out.write('# ' + '-'*220 + '\n')

of_out.close()

elapsed_d, elapsed_h, elapsed_m, elapsed_s = computation_time(time.time()-start)
