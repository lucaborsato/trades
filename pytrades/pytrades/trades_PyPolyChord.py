#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
from pytrades_lib import pytrades
import sys
import time

# PyPolyChord
import PyPolyChord as PC
import PyPolyChord.settings as PC_settings
#from PyPolyChord.priors import UniformPrior
import PyPolyChord.priors as PC_priors

# custom modules
script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), './'))
sys.path.append(module_path)
from constants import Mjups
import ancillary as anc

#from matplotlib import use as mpluse
##mpluse("Agg")
#mpluse("Qt4Agg")
#import matplotlib.pyplot as plt
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('text', usetex=True)

# ==============================================================================

# 1) INITIALISE TRADES
# 2) SET PARAMETER NAMES IN 
#    TRADES (ecosw,esinw) AND FITTING (sqrt(e)cosw,sqrt(e)sinw)
# 3) DEFINE LOGLIKELIHOOD FOR POLYCHORD SUCH THAT:
#    - loglikelihood(fitting_parameters, derived_parameters): theta,phi in PyPolyChord
#    - trades_names and fitting_names have to be 'visible' from loglikelihood
#    - within loglikelihood change fitting_parameters to trades_parameters
#    - call normal fitness loglikelihood such as in emcee
# 4) RUN POLYCHORD AS IN THE ORIGINAL SCRIPT run_PyPolyChord.py

# ==============================================================================

# ==============================================================================

def get_args():
  parser = argparse.ArgumentParser(description='TRADES+PYPOLYCHORD')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', '--p', '-path',  
                      action='store', dest='full_path', 
                      required=True, 
                      help='The path (absolute or relative) with simulation files for TRADES.')

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--s', '-sb', '--sb', '-sub-folder', '--sub-folder',
                      action='store', dest='sub_folder', 
                      required=True,
                      help='Sub-folder name, without full path.')
  
  parser.add_argument('-m', '--m', '-m-type', '--m-type',
                      action='store', dest='m_type', 
                      default='e',
                      help='Mass type (needed for the plot only). Values: "e ea ear eart earth" or "j ju jup jupi jupit jupite jupiter" or "n ne nep nept neptu neptun neptune".')
  
  parser.add_argument('--trades-fin', '--trades-final-previous', '--trades-final', 
                      action='store', dest='trades_previous', 
                      default='None', 
                      help='Define file from a previous TRADES simulation. File name and structure should be of type X_Y_finalNpar.dat. The parameters from this file will be the new original parameters')
  
  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  
  return cli

# ==============================================================================
#
# INITIALISE FOLDER AND LOG FILE
#
def init_folder(working_path, sub_folder):
  working_folder = os.path.join(working_path, sub_folder)
  if (not os.path.isdir(working_folder)):
      os.makedirs(working_folder)
  # copy files
  anc.copy_simulation_files(working_path, working_folder)
  
  run_log = os.path.join(working_folder, "trades_run.log")
  of_run = open(run_log, 'w')
  anc.print_both("# pyTRADES LOG FILE", of_run)
  anc.print_both("# working_path = %s" %(working_path), of_run)
  anc.print_both("# working_folder = %s" %(working_folder), of_run)
  anc.print_both("# run_log = %s" %(run_log), of_run)
  
  return working_folder, run_log, of_run

# ==============================================================================

def main():
  
  # READ COMMAND LINE ARGUMENTS
  cli = get_args()


  # STARTING TIME
  start = time.localtime()
  pc_output_dir = '%d-%02d-%02dT%02dh%02dm%02ds_' %(start.tm_year, 
                                                    start.tm_mon, 
                                                    start.tm_mday,
                                                    start.tm_hour, 
                                                    start.tm_min, 
                                                    start.tm_sec)
  pc_output_files = 'trades_pc'

  # RENAME 
  working_path = cli.full_path
  nthreads=1

  # INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
  pytrades.initialize_trades(working_path, cli.sub_folder, nthreads)

  # RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
  
  #global n_bodies, n_planets, ndata, npar, nfit, dof, inv_dof
  n_bodies = pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
  n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
  ndata = pytrades.ndata # TOTAL NUMBER OF DATA AVAILABLE
  npar  = pytrades.npar # NUMBER OF TOTAL PARAMATERS ~n_planets X 6
  nfit  = pytrades.nfit # NUMBER OF PARAMETERS TO FIT
  nfree  = pytrades.nfree # NUMBER OF FREE PARAMETERS (ie nrvset)
  dof   = pytrades.dof # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
  global inv_dof
  inv_dof = np.float64(1.0 / dof)
  
  # READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
  str_len = pytrades.str_len
  temp_names = pytrades.get_parameter_names(nfit,str_len)
  trades_names = anc.convert_fortran_charray2python_strararray(temp_names)
  fitting_names = anc.trades_names_to_emcee(trades_names)
  
  # save initial_fitting parameters into array
  original_fit_parameters = trades_parameters.copy()
  fitting_parameters = anc.e_to_sqrte_fitting(trades_parameters, trades_names)
  
  trades_minmax = pytrades.parameters_minmax # PARAMETER BOUNDARIES
  parameters_minmax = anc.e_to_sqrte_boundaries(trades_minmax, trades_names)

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
    # fortran variable RV in python will be rv!!!
    e_RVo = np.array(pytrades_lib.pytrades.ervobs[:], dtype=np.float64)
  except:
    e_RVo = np.array([0.], dtype=np.float64)
  try:
    e_T0o = np.array(pytrades_lib.pytrades.et0obs[:,:], dtype=np.float64).reshape((-1))
  except:
    e_T0o = np.array([0.], dtype=np.float64)
  ln_err_const = anc.compute_ln_err_const(dof, e_RVo, e_T0o, True)

  # INITIALISE SCRIPT FOLDER/LOG FILE
  working_folder, run_log, of_run = init_folder(working_path, cli.sub_folder)
  anc.print_both('',of_run)
  anc.print_both(' ======== ',of_run)
  anc.print_both(' pyTRADES' ,of_run)
  anc.print_both(' ======== ',of_run)
  anc.print_both('',of_run)
  anc.print_both(' WORKING PATH = %s' %(working_path),of_run)
  anc.print_both(' dof = ndata(%d) - nfit(%d) - nfree(%d) = %d' %(ndata, nfit, nfree, dof),of_run)
  anc.print_both(' Total N_RV = %d for %d set(s)' %(n_rv, n_set_rv),of_run)
  anc.print_both(' Total N_T0 = %d for %d out of %d planet(s)' %(n_t0_sum, n_set_t0, n_planets),of_run)
  anc.print_both(' %s = %.7f' %('log constant error = ', ln_err_const),of_run)
  
  # SET PYPOLYCHORD
  # needed to define number of derived parameters for PyPolyChord
  nder = 0
  
  # define the loglikelihood function for PyPolyChord
  def likelihood(fitting_par):
    
    # derived parameters
    derived_par = [0.0] * nder
    # convert fitting_par to trades_par
    trades_par = anc.sqrte_to_e_fitting(fitting_par, fitting_names)
    loglhd = 0.
    check = 1
    loglhd, check = pytrades.fortran_loglikelihood(np.array(trades_par, dtype=np.float64))
    #print loglhd, ln_err_const
    loglhd = loglhd + ln_err_const # ln_err_const: global variable
    
    return loglhd, derived_par

  # define the prior for the fitting parameters
  def prior(hypercube):
    """ Uniform prior from [-1,1]^D. """

    fitting_par = [0.0] * nfit
    for i, x in enumerate(hypercube):
        fitting_par[i] = PC_priors.UniformPrior(parameters_minmax[i,0], parameters_minmax[i,1])(x)

    return fitting_par

  # set PyPolyChord: the pc_settings define how to run PC, e.g. nlive, precision_criterio, etc.
  pc_settings = PC_settings.PolyChordSettings(nfit, nder)
  pc_settings.base_dir = cli.pc_output_dir
  pc_settings.file_root = cli.pc_output_files
  pc_settings.do_clustering = True
  # Possible PyPolyChord settings:
  #Keyword arguments
  #-----------------
  #nlive: int
      #(Default: nDims*25)
      #The number of live points.
      #Increasing nlive increases the accuracy of posteriors and evidences,
      #and proportionally increases runtime ~ O(nlive).

  #num_repeats : int
      #(Default: nDims*5)
      #The number of slice slice-sampling steps to generate a new point.
      #Increasing num_repeats increases the reliability of the algorithm.
      #Typically
      #* for reliable evidences need num_repeats ~ O(5*nDims).
      #* for reliable posteriors need num_repeats ~ O(nDims)

  #nprior : int
      #(Default: nlive)
      #The number of prior samples to draw before starting compression.

  #do_clustering : boolean
      #(Default: True)
      #Whether or not to use clustering at run time.

  #feedback : {0,1,2,3}
      #(Default: 1)
      #How much command line feedback to give

  #precision_criterion : float
      #(Default: 0.001)
      #Termination criterion. Nested sampling terminates when the evidence
      #contained in the live points is precision_criterion fraction of the
      #total evidence.

  #max_ndead : int
      #(Default: -1)
      #Alternative termination criterion. Stop after max_ndead iterations.
      #Set negative to ignore (default).

  #boost_posterior : float
      #(Default: 0.0)
      #Increase the number of posterior samples produced.  This can be set
      #arbitrarily high, but you won't be able to boost by more than
      #num_repeats
      #Warning: in high dimensions PolyChord produces _a lot_ of posterior
      #samples. You probably don't need to change this

  #posteriors : boolean
      #(Default: True)
      #Produce (weighted) posterior samples. Stored in <root>.txt.

  #equals : boolean
      #(Default: True)
      #Produce (equally weighted) posterior samples. Stored in
      #<root>_equal_weights.txt

  #cluster_posteriors : boolean
      #(Default: True)
      #Produce posterior files for each cluster?
      #Does nothing if do_clustering=False.

  #write_resume : boolean
      #(Default: True)
      #Create a resume file.

  #read_resume : boolean
      #(Default: True)
      #Read from resume file.

  #write_stats : boolean
      #(Default: True)
      #Write an evidence statistics file.

  #write_live : boolean
      #(Default: True)
      #Write a live points file.

  #write_dead : boolean
      #(Default: True)
      #Write a dead points file.

  #write_dead : boolean
      #(Default: True)
      #Write a prior points file.

  #update_files : int
      #(Default: nlive)
      #How often to update the files in <base_dir>.

  #base_dir : string
      #(Default: 'chains')
      #Where to store output files.

  #file_root : string
      #(Default: 'test')
      #Root name of the files produced.

  #grade_frac : List[float]
      #(Default: 1)
      #The amount of time to spend in each speed.

  #grade_dims : List[int]
      #(Default: 1)
      #The number of parameters within each speed.
  

  # RUN POLYCHORD
  pc_run = PC.run_polychord(likelihood, nfit, nder, pc_settings, prior)
  
  # set label and legend names
  kel_plot_labels = anc.keplerian_legend(fitting_names, cli.m_type)
  pc_paramnames = [('%s' %(fitting_names[i]), r'%s' %(kel_plot_labels[i])) for i in range(nfit)]
  #pc_paramnames += [('r*', 'r')]
  pc_run.make_paramnames_files(pc_paramnames)
  
  if(cli.pc_plot):
    import getdist.plots
    import matplotlib.pyplot as plt
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)
    posterior = pc_run.posterior
    g = getdist.plots.getSubplotPlotter()
    g.triangle_plot(posterior, filled=True)
    plt.show()
  
  return

# ==============================================================================
# ==============================================================================

if __name__ == "__main__":
  main()





