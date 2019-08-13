#!/usr/bin/env python
# -*- coding: utf-8 -*-

  # no more "zero" integer division bugs!:P
import sys
import argparse
import os
import numpy as np  # array
import h5py
# import trades_lib
# import random
import constants as cst  # local constants module

# from emcee import autocorr as acor
import emcee
# import acor

import scipy.optimize as sciopt
import scipy.odr as sodr
import sklearn.linear_model as sklm
import scipy.stats as scits
import glob
import shutil
import statsmodels.api as sm

# common variables needed to create labels and parameter names
kel_fmt = ['%.3f', '%.4f', '%.3f', '%.1f', '%.1f', '%.1f', '%.3f', '%.3f']

letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'l', 'm', 'n', 'o', 'p']

#deg2rad = np.pi / 180.0
#rad2deg = 180.0 / np.pi
deg2rad = cst.deg2rad
rad2deg = cst.rad2deg

eps64bit = np.finfo(np.float64(1.0)).eps
eps32bit = np.finfo(np.float32(1.0)).eps

# sigma =          1      0    -1       1     -2       2     -3       3
# sigma pos        0      1     2       3      4       5      6       7
percentile_val = [68.27, 50.0, 15.865, 84.135, 2.265, 97.735, 0.135, 99.865]


# ==============================================================================

def print_both(line, output=None):
  
  print(line)
  if (output is not None):
    output.write(line + '\n')

  return

def decode_list(alist):

  if (not 'str' in str(type(alist[0]))):
    blist = [a.decode('utf-8') for a in alist]
  else:
    blist = alist

  return blist

def encode_list(alist):

  if (not 'bytes' in str(type(alist[0]))):
    blist = [a.encode('utf-8') for a in alist]
  else:
    blist = alist

  return blist

# ==============================================================================
# INITIALISE FOLDER AND LOG FILE
# ==============================================================================
def copy_simulation_files(dir_src, dir_dest):
  # copy all the simulation files from the source directory to destination directory
  ## arg.in
  arg_file = os.path.join(dir_src, 'arg.in')
  shutil.copy(arg_file, os.path.join(dir_dest, ''))
  ## bodies.lst
  bodies_file = os.path.join(dir_src, 'bodies.lst')
  shutil.copy(bodies_file, os.path.join(dir_dest, ''))
  ## bodies file (inside bodies.lst)
  obd = open(bodies_file, 'r')
  for line in obd.readlines():
    shutil.copy(os.path.join(dir_src, line.strip().split()[0]), os.path.join(dir_dest, ''))
  obd.close()
  ## TT files
  t0files = glob.glob(os.path.join(dir_src, 'NB*_observations.dat'))
  for t0f in t0files:
    shutil.copy(t0f, os.path.join(dir_dest, ''))
  ## RV file
  if (os.path.exists(os.path.join(dir_src, 'obsRV.dat'))):
    shutil.copy(os.path.join(dir_src, 'obsRV.dat'), os.path.join(dir_dest, ''))

  ## lm.opt pso.opt pik.opt
  opt_files = 'lm.opt pso.opt pik.opt'.split()
  for opt in opt_files:
    if (os.path.exists(os.path.join(dir_src, opt))):
      shutil.copy(os.path.join(dir_src, opt), os.path.join(dir_dest, ''))

  return

# ==============================================================================


def init_folder(working_path, sub_folder):
  working_folder = os.path.join(working_path, sub_folder)
  if (not os.path.isdir(working_folder)):
      os.makedirs(working_folder)

  copy_simulation_files(working_path, working_folder)

  run_log = os.path.join(working_folder, "trades_run.log")
  of_run = open(run_log, 'w')
  print_both("# pyTRADES LOG FILE", of_run)
  print_both("# working_path = %s" %(working_path), of_run)
  print_both("# working_folder = %s" %(working_folder), of_run)
  print_both("# run_log = %s" %(run_log), of_run)

  return working_folder, run_log, of_run

# ==============================================================================

def set_bool_argument(arg_in):
  if (str(arg_in).lower() in ['t', 'tr', 'tru', 'true', 'y', 'ye', 'yes', '1']):
    arg_out = True
  else:
    arg_out = False
  return arg_out


# ==============================================================================

def set_int_argument(arg_in, default=0):
  try:
    arg_out = int(arg_in)
  except:
    arg_out = default
  return arg_out


# ==============================================================================

def set_overplot(arg_in):
  if (str(arg_in).lower() == 'none'):
    arg_out = None
  else:
    arg_out = os.path.abspath(arg_in)
    if (not os.path.isfile(arg_out)):
      try:
        arg_out = int(arg_in)
      except:
        arg_out = None
      if (arg_out not in [0, 666, 667, 668, 777, 1050, 1051, 2050, 3050, 3051]):
        arg_out = None
  return arg_out


# ==============================================================================

def set_adhoc_file(arg_in):
  if (str(arg_in).lower() == 'none'):
    arg_out = None
  else:
    arg_out = os.path.abspath(arg_in)
    if (not os.path.isfile(arg_out)):  arg_out = None

  return arg_out


# ==============================================================================

def set_int_or_none(arg_in):
  try:
    arg_out = int(arg_in)
  except:
    arg_out = None
  return arg_out


# ==============================================================================

# read command line (cli) arguments
def get_args():
  parser = argparse.ArgumentParser(description='TRADES+EMCEE PLOT')

  parser.add_argument('-p', '--path', 
                      action='store', dest='full_path', required=True,
                      help='The path (absolute or relative) with simulation files for TRADES.'
                      )

  parser.add_argument('-nb', '--nb', '--nburn', '-nburn', '--nburnin', '-nburnin',
                      action='store', dest='nburnin', required=True,
                      help='The number of posterior/burn in steps to discard at the beginning of each chain. It has to be > 0'
                      )

  parser.add_argument('-m', '--mtype', '--mass-type', 
                      action='store', dest='m_type', default='e',
                      help='Mass type: j = Jupiter, e = Earth, s = Sun. Default is Earth = e.'
                      )

  parser.add_argument('-t', '--temp-file', 
                      action='store', dest='temp_status', default=False,
                      help='If you want to read temporary emcee_temp.hdf5 file, because simulation is not finished yet. Default is False'
                      )
  
  parser.add_argument('-b', '--boot-file', '--save-like-boot',
                      action='store', dest='boot_id', default=0,
                      help='If you want to save flatchain posterior as the posterior file, in order to be read by read_finalpar_v2.py. If set to 0 (default) it doesn\'t save anything. Bootstrap file output: boot_id_posterior_sim.dat'
                      )

  parser.add_argument('-c', '--cumulative', '--cumulative-histogram',
                      action='store', dest='cumulative', default='False',
                      help='If you want to plot the cumulative histogram instead of normal histogram. Set it to "True/Yes/1" Default is False. '
                      )

  parser.add_argument('-s', '--steps', '--gelmanrubin-steps',
                      action='store', dest='sel_steps', default=0,
                      help='Set to positive integer. It will allow to compute GelmanRubin/Geweke statistics of only sel_steps, and not for each step. Just to debug and overall view of the chains. Default is 0.'
                      )

  parser.add_argument('-u', '--thin', '--use-thin',
                      action='store', dest='use_thin', default=False,
                      help='Set if you want to use or not the thinning parameter computed from the autocorrelation time. Default False'
                      )

  parser.add_argument('--sample',
                      action='store', dest='sample_str', default='None',
                      help='Set a list of parameter names to select a parameter sample within the Credible Intervals, first will be the "pivoting" parameters, then the parameters to not take into account. Default = None, it means first parameters fitted as "pivot", all other have to be within the Credible Intervals.'
                      )

  parser.add_argument('--seed', 
                      action='store', dest='seed', default='None',
                      help='Set the seed for random number generator. Default is None.'
                      )

  parser.add_argument('--overplot',
                      action='store', dest='overplot', default='None',
                      help='Define which parameter set to be overplotted on the correlation plot.\nInput the proper number:\n0 (starting parameters),\n666 (pick sample),\n777 (ad hoc sample),\n1050 (median, derived as median of the derived posterior),\n1051 (median, derived computed from median parameters), 2050 (max lnprob),\n3050 (mode, derived as mode of the derived posterior),\n3051 (mode, derived computed from mode parameters)')

  parser.add_argument('--ad-hoc', action='store', dest='adhoc', default='None',
                      help='Define adhoc file with fitted parameters to be overplotted with the id 777')

  parser.add_argument('--script-folder', action='store', dest='pyscript', default='./',
                      help='Folder of the python scripts... Default is "./"')

  parser.add_argument('-n-samples', '--n-samples', action='store', dest='n_samples', default=0,
                      help='Number of sample parameter within CI to generate T0s and RVs to be overplotted in the O-C plots. Defautl is 0')

  cli = parser.parse_args()

  cli.full_path = os.path.abspath(cli.full_path)
  cli.nburnin = set_int_or_none(cli.nburnin)
  cli.m_type = cli.m_type.lower()
  cli.temp_status = set_bool_argument(cli.temp_status)
  cli.cumulative = set_bool_argument(cli.cumulative)
  cli.boot_id = set_int_argument(cli.boot_id)
  cli.use_thin = set_int_argument(cli.use_thin, default=False)
  #cli.use_thin = set_bool_argument(cli.use_thin)
  cli.seed = set_int_or_none(cli.seed)
  cli.overplot = set_overplot(cli.overplot)
  cli.adhoc = set_adhoc_file(cli.adhoc)
  cli.pyscript = os.path.abspath(cli.pyscript)
  cli.n_samples = set_int_argument(cli.n_samples, default=0)

  return cli


# ==============================================================================

def compute_ln_err_const(dof, e_RVo, e_T0o, ln_flag=False):
  
  if (ln_flag):
    eRV = e_RVo[e_RVo > 0.]
    eT0 = e_T0o[e_T0o > 0.]
    
    ln_e_RVo = np.sum(np.log(eRV*eRV))
    ln_e_T0o = np.sum(np.log(eT0*eT0))

    ln_err_const = - 0.5 * dof * np.log(2.*np.pi) - 0.5 * ( ln_e_RVo + ln_e_T0o)
  else:
    ln_err_const = 0.
  
  return ln_err_const

# ==============================================================================

# given the mass flag letter it computes the proper mass conversion factor
def mass_type_factor(Ms=1.0, mtype='earth', mscale=True):
  if (mtype.lower() in ['s', 'su', 'sun']):
    conv_factor = np.float64(1.0)
    scale_factor = Ms
    mass_unit = 'M_Sun'
  elif (mtype.lower() in ['e', 'ea', 'ear', 'eart', 'earth']):
    conv_factor = cst.Msear
    scale_factor = Ms * cst.Msear
    mass_unit = 'M_Earth'
  elif (mtype.lower() in ['n', 'ne', 'nep', 'nept', 'neptu', 'neptun', 'neptune']):
    conv_factor = cst.Msnep
    scale_factor = Ms * cst.Msnep
    mass_unit = 'M_Nep'
  else:
    conv_factor = cst.Msjup
    scale_factor = Ms * cst.Msjup
    mass_unit = 'M_Jup'
  if (mscale):
    return scale_factor, mass_unit
  else:
    return conv_factor, mass_unit


# ==============================================================================

def get_proper_mass(m_type, parameter_names_emcee, full_path):

  # nfit = np.size(parameter_names_emcee)

  MR_star = np.ones((2, 2))

  bodies = os.path.join(full_path, 'bodies.lst')
  obf = open(bodies, 'r')
  line_star = obf.readline()
  obf.close()
  file_star = line_star.split()[0].strip()
  MR_star = np.genfromtxt(file_star)
  if (len(MR_star.shape) == 1): MR_star = np.column_stack((MR_star, np.zeros((2))))

  m_factor, m_unit = mass_type_factor(MR_star[0, 0], m_type, True)

  return m_factor, m_unit


# ==============================================================================

# prepare the labels of the Keplerian orbital elements for the legend
def keplerian_legend(parameter_names, m_type):
  _, m_unit = mass_type_factor(1.0, m_type, mscale=False)
  # if (m_type in ['j', 'jup', 'jupiter']):
  # unit_mass = '$[M_\mathrm{Jup}]$'
  # elif (m_type in ['s', 'sun']):
  # unit_mass = '$[M_{\odot}]$'
  if (m_unit == r'M_Jup'):
    unit_mass   = r'M_\mathrm{Jup}'
    unit_radius = r'R_\mathrm{Jup}'
  elif (m_unit == r'M_Earth'):
    unit_mass   = r'M_\oplus'
    unit_radius = r'R_\oplus'
  elif (m_unit == r'M_Sun'):
    unit_mass   = r'M_\odot'
    unit_radius = r'R_\odot'
  elif (m_unit == r'M_Nep'):
    unit_mass   = r'M_\mathrm{Nep}'
    unit_radius = r'R_\mathrm{Nep}'
  else:
    unit_mass = ' '

  nfit = np.shape(parameter_names)[0]
  kel_legends = np.zeros((nfit), dtype='|S256')

  for i in range(0, nfit):
    parameter_names[i] = parameter_names[i].strip()

    if ('Ms' in parameter_names[i]):
      planet_id = int(parameter_names[i].split('m')[1].split('Ms')[0]) - 1
      # kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id_2[0], letters[planet_id], kel_units_2[0])
      kel_legends[i] = r'$\displaystyle\frac{M_\mathrm{%s}}{M_\star} \left[\frac{%s}{M_\odot}\right]$' % (
      letters[planet_id], unit_mass)

    elif (parameter_names[i][0] == 'm' and parameter_names[i][-1] != 's'):
      planet_id = int(parameter_names[i][1:]) - 1
      kel_legends[i] = r'$M_\mathrm{%s} [%s]$' % (letters[planet_id], unit_mass)

    elif (parameter_names[i][0] == 'R'):
      planet_id = int(parameter_names[i][1:]) - 1
      kel_legends[i] = r'$R_\mathrm{%s}$ [%s]$' % (letters[planet_id], unit_radius)

    elif (parameter_names[i][0] == 'P'):
      planet_id = int(parameter_names[i][1:]) - 1
      kel_legends[i] = r'$P_\mathrm{%s}$ [days]' % (letters[planet_id])

    elif ('e' in parameter_names[i]):
      if (parameter_names[i][0:4] == 'ecos'):
        planet_id = int(parameter_names[i].split('w')[1].strip()) - 1
        kel_legends[i] = r'$e\cos\omega_\mathrm{%s}$' % (letters[planet_id])
      elif (parameter_names[i][0:6] == 'sqrtec'):
        planet_id = int(parameter_names[i].split('w')[1].strip()) - 1
        kel_legends[i] = r'$\sqrt{e}\cos\omega_\mathrm{%s}$' % (letters[planet_id])
      elif (parameter_names[i][0:4] == 'esin'):
        planet_id = int(parameter_names[i].split('w')[1].strip()) - 1
        kel_legends[i] = r'$e\sin\omega_\mathrm{%s}$' % (letters[planet_id])
      elif (parameter_names[i][0:6] == 'sqrtes'):
        planet_id = int(parameter_names[i].split('w')[1].strip()) - 1
        kel_legends[i] = r'$\sqrt{e}\sin\omega_\mathrm{%s}$' % (letters[planet_id])
      else:
        planet_id = int(parameter_names[i].split('e')[1].strip()) - 1
        kel_legends[i] = r'$e_\mathrm{%s}$' % (letters[planet_id])

    elif (parameter_names[i][0] == 'w'):
      planet_id = int(parameter_names[i].split('w')[1].strip()) - 1
      kel_legends[i] = r'$\omega_\mathrm{%s} [^\circ]$' % (letters[planet_id])

    elif (parameter_names[i][0:2] == 'mA'):
      planet_id = int(parameter_names[i][2:]) - 1
      kel_legends[i] = r'$\mathcal{M}_\mathrm{%s} [^\circ]$' % (letters[planet_id])

    elif ('lambda' in parameter_names[i]):
      planet_id = int(parameter_names[i].split('lambda')[1].strip()) - 1
      kel_legends[i] = r'$\lambda_\mathrm{%s} [^\circ]$' % (letters[planet_id])

    elif (parameter_names[i][0] == 'i'):
      if ('icos' in parameter_names[i]):
        planet_id = int(parameter_names[i].split('lN')[1].strip()) - 1
        kel_legends[i] = r'$i\cos\Omega_\mathrm{%s}$' % (letters[planet_id])
      elif ('isin' in parameter_names[i]):
        planet_id = int(parameter_names[i].split('lN')[1].strip()) - 1
        kel_legends[i] = r'$\i\sin\Omega_\mathrm{%s}$' % (letters[planet_id])
      else:
        planet_id = int(parameter_names[i][1:]) - 1
        kel_legends[i] = r'$i_\mathrm{%s} [^\circ]$' % (letters[planet_id])

    elif (parameter_names[i][0:2] == 'lN'):
      planet_id = int(parameter_names[i][2:]) - 1
      kel_legends[i] = r'$\Omega_\mathrm{%s} [^\circ]$' % (letters[planet_id])

  kel_legends = [r'{0}'.format(kl.replace('[', '(').replace(']', ')').strip())
                 for kl in decode_list(kel_legends)]

  return kel_legends


# ==============================================================================

def derived_labels(derived_names, m_type):
  _, m_unit = mass_type_factor(1., m_type, mscale=False)
  if (m_unit == 'M_Jup'):
    unit_mass = r'M_\mathrm{Jup}'
  elif (m_unit == 'M_Earth'):
    unit_mass = r'M_\oplus'
  elif (m_unit == 'M_Sun'):
    unit_mass = r'M_\odot'
  elif (m_unit == 'M_Nep'):
    unit_mass = r'M_\mathrm{Nep}'
  else:
    unit_mass = ' '

  labels_list = []
  n_der = np.shape(derived_names)[0]
  for ider in range(0, n_der):

    # mass
    if (derived_names[ider][0] == 'm' and derived_names[ider][1] != 'A'):
      planet_id = int(derived_names[ider].split('m')[1]) - 1
      labels_list.append(r'$M_\mathrm{%s} [%s]$' % (letters[planet_id], unit_mass))

    # eccentricity
    elif (derived_names[ider][0] == 'e'):
      planet_id = int(derived_names[ider].split('e')[1]) - 1
      labels_list.append(r'$e_\mathrm{%s}$' % (letters[planet_id]))

    # argument of pericentre
    elif (derived_names[ider][0] == 'w'):
      planet_id = int(derived_names[ider].split('w')[1]) - 1
      # labels_list.append(r'$\omega_\mathrm{%s} [\mathrm{deg}]$' %(letters[planet_id]))
      labels_list.append(r'$\omega_\mathrm{%s} [^\circ]$' % (letters[planet_id]))

    # mean anomaly
    elif (derived_names[ider][0:2] == 'mA'):
      planet_id = int(derived_names[ider].split('mA')[1]) - 1
      # labels_list.append(r'$\mathcal{M}_\mathrm{%s} [\mathrm{deg}]$' %(letters[planet_id]))
      labels_list.append(r'$\mathcal{M}_\mathrm{%s} [^\circ]$' % (letters[planet_id]))

    # inclination
    elif (derived_names[ider][0] == 'i'):
      planet_id = int(derived_names[ider].split('i')[1]) - 1
      # labels_list.append(r'$i_\mathrm{%s} [\mathrm{deg}]$' %(letters[planet_id]))
      labels_list.append(r'$i_\mathrm{%s} [^\circ]$' % (letters[planet_id]))

    # longitude of node
    elif (derived_names[ider][0:2] == 'lN'):
      planet_id = int(derived_names[ider].split('lN')[1]) - 1
      # labels_list.append(r'$\Omega_\mathrm{%s} [\mathrm{deg}]$' %(letters[planet_id]))
      labels_list.append(r'$\Omega_\mathrm{%s} [^\circ]$' % (letters[planet_id]))

  labels_list = [r'{0}'.format(l.replace('[', '(').replace(']', ')'))
                 for l in decode_list(labels_list)]

  return labels_list


# ==============================================================================

# only needed by check_good_parameters to convert Mjup to flagged mass
def check_good_parameters(good_id, good_parameters_in, m_factor, nfit):
  
  good_parameters_out = np.zeros(nfit) + good_parameters_in
  for i in range(0, nfit):
    if (good_id[i].strip()[0] == 'm' and good_id[i].strip()[1] != 'A'):
      good_parameters_out[i] = good_parameters_in[i] * cst.Mjups * m_factor
      
  return good_parameters_out


# ==============================================================================

def read_check_good_parameters(full_path, good_status, m_factor, nfit):
  
  # read good parameters file
  good_parameters_0, good_parameters_1 = False, False
  if (good_status):
    good_file = os.path.join(full_path, 'good_parameters.dat')
    if (os.path.exists(good_file)):
      good_parameters_0 = np.genfromtxt(good_file, usecols=(1), dtype=np.float64)
      good_id = np.genfromtxt(good_file, usecols=(0), dtype='|S10')
      good_parameters_1 = check_good_parameters(good_id, good_parameters_0, m_factor, nfit)
      # for ii in range(0,nfit):
      # if (parameter_names_emcee[ii][:2] == 'mA'):
      # good_parameters_0[ii] = rescale_angle(good_parameters[ii])
      
  return good_parameters_0, good_parameters_1


# ==============================================================================

# sometimes the angle distribution are bimodal because they are centered in the wrong position
# e.g.: -1° = 359° etc...
def rescale_angle(angle):
  
  # with arctan2
  cosa = np.cos(angle * deg2rad)
  sina = np.sin(angle * deg2rad)
  new_angle = np.arctan2(sina, cosa) * 180. / np.pi

  return new_angle


# ==============================================================================

# renormalize the angular parameter mean Anomaly mA
def renormalize_parameters(parameters, parameter_names):
  
  new_parameters = parameters.copy()
  nfit = np.array(new_parameters).shape[0]
  for i in range(0, nfit):
    if (parameter_names[i][:2] == 'mA'):
      new_parameters[i] = (parameters[i] + 360.) % 360.
      
  return new_parameters


# ==============================================================================

# get indices of max of an array
def get_max_indices(array_values):
  
  # print(' ++ from index = ', np.argmax(array_values))
  # idx = np.unravel_index(np.argmax(array_values), np.shape(array_values))
  lmax = np.max(array_values)
  print(' ++ lnL_max = {}'.format(lmax))
  xx = np.where(array_values == lmax)
  print(' ++ from index = ', xx)
  idx = [xx[0][0], xx[1][0]]
  print(' ++ to indexes = ', idx)

  return idx[0], idx[1]

# ==============================================================================

# prepare file and best folder emcee
def get_emcee_file_and_best(emcee_folder, temp_status):
  
  if (temp_status):
    folder_best = 'best_temp'
    emcee_file = os.path.join(emcee_folder, 'emcee_temp.hdf5')
    if (not os.path.exists(emcee_file)):
      emcee_file = os.path.join(emcee_folder, 'emcee_summary.hdf5')
    emcee_best = os.path.join(emcee_folder, folder_best)
  else:
    folder_best = 'best'
    emcee_file = os.path.join(emcee_folder, 'emcee_summary.hdf5')
    emcee_best = os.path.join(emcee_folder, folder_best)
    
  return emcee_file, emcee_best, folder_best


# ==============================================================================

def get_devol_file(emcee_folder):
  
  devol_file = os.path.join(emcee_folder, 'best_devol.hdf5')
  
  return devol_file


# ==============================================================================

def get_percentile_angle(angle_posterior):
  
  cosa_posterior = np.cos(angle_posterior * np.pi / 180.)
  sina_posterior = np.sin(angle_posterior * np.pi / 180.)
  temp_cos = np.percentile(cosa_posterior, 50.)
  temp_sin = np.percentile(sina_posterior, 50.)
  median_angle = (np.arctan2(temp_sin, temp_cos) * 180. / np.pi) % 360.
  lower_angle = np.percentile(angle_posterior, 16.)
  upper_angle = np.percentile(angle_posterior, 84.)
  
  return median_angle, lower_angle, upper_angle


# ==============================================================================

def get_data(emcee_file, temp_status):
  
  ## read data from hdf5 file
  completed_steps = 0
  # print ' reading file', emcee_file
  f_read = h5py.File(emcee_file, 'r')
  data_names = [str(i) for i in list(f_read.keys())]
  parameter_names_emcee, parameter_boundaries = [], []
  chains = []
  acceptance_fraction, autocor_time, lnprobability = [], [], []

  if ('parameter_names' in data_names):
    # parameter_names_emcee = np.array(f_read['parameter_names'], dtype='S15').astype(str)
    names_temp = f_read['parameter_names'][...]
    parameter_names_emcee = decode_list(names_temp)
    # parameter_names_emcee = names_temp

  if ('boundaries' in data_names):
    parameter_boundaries = np.array(f_read['boundaries'][...],
                                    dtype=np.float64)
  # if ('final_parameters' in data_names):  final_parameters = np.array(f_read['final_parameters'], dtype=np.float64)
  if ('chains' in data_names):
    chains = np.array(f_read['chains'],
                      dtype=np.float64)  # shape (nwalkers, nrusn, nfit)
  if ('acceptance_fraction' in data_names):
    acceptance_fraction = np.array(f_read['acceptance_fraction'],
                                   dtype=np.float64)
  if ('autocor_time' in data_names): 
    autocor_time = np.array(f_read['autocor_time'],
                            dtype=np.float64)
  if ('lnprobability' in data_names): 
    lnprobability = np.array(f_read['lnprobability'],
                             dtype=np.float64)
  try:
    ln_err_const = f_read['lnprobability'].attrs['ln_err_const']
  except:
    ln_err_const = 0.0
  if (temp_status):
    completed_steps = f_read['chains'].attrs['completed_steps']
  else:
    completed_steps = np.shape(chains)[1]
  # close hdf5 file
  f_read.close()

  return parameter_names_emcee, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps


# ==============================================================================

def get_last_emcee_iteration(emcee_file, nwalkers):
  
  f_read = h5py.File(emcee_file, 'r')
  chains = f_read['chains'][...]
  nw_c, nr_c, _ = np.shape(chains)
  done = True

  try:
    last_step = f_read['chains'].attrs['completed_steps']
    # print 'completed_steps'
  except:
    # print 'all steps'
    last_step = nr_c
  f_read.close()
  # print '*****', last_step,'*****'

  if (nwalkers == nw_c):
    # last_p0 = chains[:,last_step-1,:].copy()
    last_p0 = [chains[iw, last_step - 1, :].copy() for iw in range(0, nw_c)]
    done = True
  else:
    last_p0 = None
    done = False
  # elif(nwalkers > nw_c):

  # else:

  del chains

  return last_p0, nw_c, done


# ==============================================================================

def compute_autocor_time(chains, walkers_transposed=True):
  
  if (walkers_transposed):

    _, nw, nfit = np.shape(chains)
    autocor_time = np.zeros((nfit), dtype=np.float64)
    for ifit in range(0, nfit):
      x = chains[:, :, ifit]
      # tau, mean, sigma = acor.acor(x)
      #temp_acor = np.mean(np.array([acor.acor(x[:, iw]) for iw in range(0, nw)]), axis=0)
      temp_acor = np.mean(np.array([emcee.autocorr.integrated_time(x[:, iw], c=1) for iw in range(0, nw)]), axis=0)
      autocor_time[ifit] = temp_acor[0]

  else:

    nw, _, nfit = chains.shape
    autocor_time = np.zeros((nfit), dtype=np.float64)
    for ifit in range(0, nfit):
      x = chains[:, :, ifit]
      #temp_acor = np.mean(np.array([acor.acor(x[iw, :]) for iw in range(0, nw)]), axis=1)
      temp_acor = np.mean(np.array([emcee.autocorr.integrated_time(x[iw, :], c=1) for iw in range(0, nw)]), axis=1)
      autocor_time[ifit] = temp_acor[0]

  return autocor_time


# ==============================================================================
def compute_acor_time(sampler, steps_done=None):
  
  if (steps_done is None):
    actime_const = 100
  else:
    if (steps_done < 100):
      actime_const = steps_done
    else:
      actime_const = 100

  _, _, nfit = np.shape(sampler.chain)
  try:
    acor_time = sampler.acor
  except:
    acor_time = np.zeros((nfit), dtype=np.float64) + actime_const

  return acor_time


# ==============================================================================

def get_emcee_parameters(chains, temp_status, nburnin_in, completed_steps):
  
  # determine nwalkers, nruns, nburnin, and nfit parameters
  # nwalkers = chains.shape[0]
  # nruns = chains.shape[1]
  nwalkers, nruns, nfit = np.shape(chains)
  if (temp_status):
    nruns = int(completed_steps)
  # select posterior chains, without burn in steps
  nburnin = 0
  
  if (nburnin_in < 0):
    # print ' WARNING: nburnin <= 0. It will be set to: nburnin = nruns * 10% = %d * 10% = %d' %(nruns, nruns*0.1)
    nburnin = 0
  elif(nburnin_in >= nruns or nburnin_in is None):
    nburnin = np.rint(0.5*nruns).astype(int)
  else:
    nburnin = int(nburnin_in)
  nruns_sel = nruns - nburnin

  return nfit, nwalkers, nruns, nburnin, nruns_sel


# ==============================================================================

def print_memory_usage(array_values):
  
  print(' MEMORY USAGE: array_values = %d bytes = %.2d MBytes = %.4f GBytes' % (
  array_values.nbytes, array_values.nbytes / (1024. ** 2), array_values.nbytes / (1024. ** 3)))
  
  return


# ==============================================================================

def select_transpose_convert_chains(nfit, nwalkers, npost, nruns, nruns_sel, m_factor, parameter_names_emcee,
                                    parameter_boundaries_in, chains):
  
  # chain is transposed: needed to plot quicker chains for each walker: nruns vs value of parameter
  parameter_boundaries = parameter_boundaries_in.copy()
  chains_T = np.zeros((nruns, nwalkers, nfit))

  for ii in range(0, nfit):
    chains_T[:, :, ii] = chains[:, :nruns, ii].T  # transpose
    if (parameter_names_emcee[ii][0] == 'm' and parameter_names_emcee[ii][1] != 'A'):

      if ('Ms' in parameter_names_emcee[ii]):
        m_conv = m_factor
      else:
        m_conv = np.float64(1.)

      chains_T[:, :, ii] = chains_T[:, :, ii] * m_conv
      parameter_boundaries[ii, :] = parameter_boundaries[ii, :] * m_conv

  return chains_T, parameter_boundaries


# ==============================================================================

def posterior_back_to_msun(m_factor, parameter_names_emcee, flatchain_posterior_in):
  
  nfit = flatchain_posterior_in.shape[1]
  flatchain_posterior_out = flatchain_posterior_in.copy()
  for ifit in range(0, nfit):
    if (parameter_names_emcee[ifit][0] == 'm' and parameter_names_emcee[ifit][1] != 'A'):
      flatchain_posterior_out[:, ifit] = flatchain_posterior_out[:, ifit] / m_factor

  return flatchain_posterior_out


# ==============================================================================

def prepare_plot_folder(full_path):
  
  plot_folder = prepare_emcee_plot_folder(full_path)

  return plot_folder


# ==============================================================================

def prepare_emcee_plot_folder(full_path):
  
  emcee_plots = os.path.join(full_path, 'plots')
  if (not os.path.isdir(emcee_plots)):
    if (not os.path.exists(emcee_plots)):
      os.makedirs(emcee_plots)

  return emcee_plots


# ==============================================================================

def computation_time(elapsed):
  
  elapsed_d = elapsed / 86400.
  elapsed_h = (elapsed_d - int(elapsed_d)) * 24.
  elapsed_m = (elapsed_h - int(elapsed_h)) * 60.
  elapsed_s = (elapsed_m - int(elapsed_m)) * 60.

  return int(elapsed_d), elapsed_h, elapsed_m, elapsed_s


# ==============================================================================

def get_pso_data(pso_file):
  
  of_pso = h5py.File(pso_file, 'r')
  population = np.array(of_pso['population'], dtype=np.float64)
  population_fitness = np.array(of_pso['population_fitness'], dtype=np.float64)
  pso_parameters = np.array(of_pso['pso_parameters'], dtype=np.float64)
  pso_fitness = np.array(of_pso['pso_fitness'], dtype=np.float64)
  if ('pso_best_evolution' in list(of_pso.keys())):
    pso_best_evolution = np.array(of_pso['pso_best_evolution'], dtype=np.float64)
  else:
    pso_best_evolution = False
  if ('parameters_minmax' in list(of_pso.keys())):
    parameters_minmax = np.array(of_pso['parameters_minmax'], dtype=np.float64)
  else:
    parameters_minmax = False
  if ('parameter_names' in list(of_pso.keys())):
    parameter_names = np.array(of_pso['parameter_names'], dtype='S10')
  else:
    parameter_names = False
  of_pso.close()
  pop_shape = population.shape

  return population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution, parameters_minmax, parameter_names, pop_shape


# ==============================================================================

def compute_limits(vec_a, delta=0.05):
  
  a_min = np.min(np.array(vec_a))
  a_max = np.max(np.array(vec_a))
  da = np.abs(a_max - a_min)
  lim_min = a_min - da * delta
  lim_max = a_max + da * delta

  return lim_min, lim_max


# ==============================================================================

def thin_the_chains(use_thin, nburnin, nruns, nruns_sel, autocor_time,
                    chains_T_full, lnprobability, burnin_done=False,
                    full_chains_thinned=False):
  
  _, _, nfit = np.shape(chains_T_full)
  print(' nburnin = ', nburnin)
  print(' nruns = ', nruns)
  print(' nruns_sel = ', nruns_sel)
  # print ' np.shape(lnprobability) = ',np.shape(lnprobability)
  nr, nc = np.shape(lnprobability)

  if (not burnin_done):
    chains_T_posterior = chains_T_full[nburnin:nruns, :, :].copy()
    # lnprob_posterior = lnprobability[:, nburnin:nruns].copy()
    if(nr > nc):
      lnprob_posterior = lnprobability[nburnin:nruns, :].copy()
    else:
      lnprob_posterior = lnprobability.T[nburnin:nruns, :].copy()
  else:
    chains_T_posterior = chains_T_full[:, :nruns, :].copy()
    # lnprob_posterior = lnprobability[:, :nruns].copy()
    if(nr > nc):
      lnprob_posterior = lnprobability[:nruns, :].copy()
    else:
      lnprob_posterior = lnprobability.T[:nruns, :].copy()

  # print ' np.shape(lnprob_posterior) = ',np.shape(lnprob_posterior)

  if(not use_thin or use_thin <= 0): # use_thin is False or <= 0
    
    chains_T = chains_T_posterior
    flatchain_posterior_0 = chains_T[:, :, :].copy().reshape((-1, nfit))
    print(' posterior not thinned shape = ', np.shape(flatchain_posterior_0))
    lnprob_burnin = lnprob_posterior[:, :].copy()
    print(' lnprob_burnin not thinned shape = ', np.shape(lnprob_burnin))
    thin_steps = 0

  else: # use_thin is True or > 0
    if(type(use_thin) is bool and use_thin): # use_thin is True
      try:
        n_acor = autocor_time.shape[0]
      except:
        n_acor = len(autocor_time)
      if (n_acor == 0):
        # autocor_time = compute_autocor_time(flatchain_posterior_0)
        autocor_time = compute_autocor_time(chains_T_posterior)
      thin_steps = np.rint(np.mean(np.array(autocor_time, dtype=np.float64))).astype(int)
      print(' computed thin_steps = ', thin_steps)
      if (thin_steps > 1000):
        thin_steps = 1000
        print(' set thin_steps = 1000')
        
    else: # use_thin > 0
      thin_steps = use_thin

    sel_thin_steps = np.arange(0, nruns_sel + thin_steps, thin_steps).astype(int)
    if (sel_thin_steps[-1] >= nruns_sel): sel_thin_steps[-1] = nruns_sel - 1
    n_thin = np.shape(sel_thin_steps)[0]

    print(' n_thin = ', n_thin)

    chains_T = chains_T_posterior[sel_thin_steps, :, :]
    # create a flat array of the posterior: from (nruns_sel, nwalkers, nfit) -> (n_thin * nwalkers, nfit)
    flatchain_posterior_0 = chains_T[:, :, :].copy().reshape((-1, nfit))
    print(' posterior thinned shape = ', np.shape(flatchain_posterior_0))
    lnprob_burnin = lnprob_posterior[:, sel_thin_steps].copy()
    print(' lnprob_burnin thinned shape = ', np.shape(lnprob_burnin))

    if (full_chains_thinned):
      tempsel = np.arange(0, nruns + thin_steps, thin_steps).astype(int)
      if (tempsel[-1] >= nruns): tempsel[-1] = nruns - 1
      chains_T_full_thinned = chains_T_full[tempsel, :, :].copy()
      
      return chains_T, flatchain_posterior_0, lnprob_burnin, thin_steps, chains_T_full_thinned

  return chains_T, flatchain_posterior_0, lnprob_burnin, thin_steps


# ==============================================================================

def get_sigmas(best_parameters, flatchain_posterior):
  
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  sigma_perc68 = np.percentile(np.abs(flatchain_posterior - best_parameters), 68.27, axis=0, interpolation='midpoint')
  # retrieve confidence intervals of the residual distribution
  sigma_confint = np.percentile(flatchain_posterior - best_parameters, sigmas_percentiles, axis=0,
                                interpolation='midpoint')

  return sigma_perc68, sigma_confint


# ==============================================================================

def get_maxlnprob_parameters(lnprob_burnin, chains_T, flatchain_posterior):
  
  print(' -- shape(chains_T)      = ', np.shape(chains_T))
  print(' -- shape(lnprob_burnin) = ',np.shape(lnprob_burnin))
  maxlnprob_row, maxlnprob_col = get_max_indices(lnprob_burnin)
  # print(' -- maxlnprob_row, maxlnprob_col = ',maxlnprob_row, maxlnprob_col)
  maxlnprob = lnprob_burnin[maxlnprob_row, maxlnprob_col]
  # retrieve best parameters <-> best lnprobability
  # maxlnprob_parameters = chains_T[maxlnprob_col, maxlnprob_row, :]
  maxlnprob_parameters = chains_T[maxlnprob_row, maxlnprob_col, :]
  # retrieve 1sigma as 68.27th percentile of the absolute residual distribution
  # retrieve confidence intervals of the residual distribution
  maxlnprob_perc68, maxlnprob_confint = get_sigmas(maxlnprob_parameters, flatchain_posterior)

  return maxlnprob, maxlnprob_parameters, maxlnprob_perc68, maxlnprob_confint


# ==============================================================================

def get_median_parameters(flatchain_posterior):
  
  median_parameters = np.percentile(flatchain_posterior, 50.0, axis=0, interpolation='midpoint')
  median_perc68, median_confint = get_sigmas(median_parameters, flatchain_posterior)

  return median_parameters, median_perc68, median_confint


# ==============================================================================

def get_parameters_median_fitness(nwalkers, npost, nruns, lnprobability, flatchain_posterior, ln_err_const):
  
  nruns_sel = nruns - npost
  lnprob_burnin = lnprobability[:, npost:nruns]
  flat_lnprob = lnprob_burnin.T.reshape((nwalkers * nruns_sel))
  flat_fitness = -2. * (flat_lnprob - ln_err_const)
  n_med = int(nwalkers * nruns_sel * 0.5)
  id_med = np.argsort(flat_fitness)[n_med]
  # retrieve median fitness
  # median_fitness = np.percentile(flat_fitness, 50., interpolation='midpoint')
  median_fitness = flat_fitness[id_med]
  # retrieve parameters at id_med
  medfit_parameters = flatchain_posterior[id_med, :]
  medfit_perc68, medfit_confint = get_sigmas(medfit_parameters, flatchain_posterior)

  return median_fitness, medfit_parameters, medfit_perc68, medfit_confint


# ==============================================================================

# def compute_max_mean(data_vec, k):
# hist_counts, bin_edges = np.histogram(data_vec, bins=k)
# max_bin = np.argmax(np.array(hist_counts))
# max_mean = np.mean(data_vec[np.logical_and(data_vec>=bin_edges[max_bin], data_vec<=bin_edges[max_bin+1])])
# return max_mean, max_bin

# ==============================================================================

def compute_max_mean(data_vec, k):
  
  hist_counts, bin_edges = np.histogram(data_vec, bins=k)
  max_bin = np.argmax(np.array(hist_counts))
  if (max_bin == 0 or max_bin == k):
    ext_bin = 0
  elif (max_bin == 1 or max_bin == k - 1):
    ext_bin = 1
  else:
    ext_bin = 2

  data_max = data_vec[np.logical_and(data_vec >= bin_edges[max_bin - ext_bin], data_vec < bin_edges[max_bin + ext_bin])]
  if(np.shape(data_max)[0] < 1):
    print('WARNING: n(data_max) < 1!')
    return np.nan, 0
    sys.stdout.flush()
  max_mean = np.mean(data_max)

  # print 'MAX BIN: %d with %d counts' %(max_bin, np.max(hist_counts))
  # print 'Selected values in between bins: %d , %d' %(max_bin-ext_bin, max_bin+ext_bin)
  # print 'Corresponding bin_edges: %d , %d' %(bin_edges[max_bin-ext_bin], bin_edges[max_bin+ext_bin])

  return max_mean, max_bin


# ==============================================================================

def get_mode_parameters_full(flatchain_posterior, k):
  
  _, nfit = np.shape(flatchain_posterior)
  mode_parameters = np.zeros((nfit))
  mode_bin = np.zeros((nfit)).astype(int)
  
  for i_fit in range(0, nfit):
    data_vec = flatchain_posterior[:, i_fit]
    mode_parameters[i_fit], mode_bin[i_fit] = compute_max_mean(data_vec, k)

  mode_perc68, mode_confint = get_sigmas(mode_parameters, flatchain_posterior)

  return mode_bin, mode_parameters, mode_perc68, mode_confint


# ==============================================================================

def get_mode_parameters(flatchain_posterior, k):
  
  _, nfit = np.shape(flatchain_posterior)
  mode_parameters = np.zeros((nfit))
  mode_bin = np.zeros((nfit)).astype(int)
  
  print('get_mode_parameters')
  sys.stdout.flush()
  
  for i_fit in range(0, nfit):
    print('i_fit = {0:d}'.format(i_fit), end = ' ')
    data_vec = flatchain_posterior[:, i_fit]
    print('computing mode and bin ... ', end = ' ')
    mode_parameters[i_fit], mode_bin[i_fit] = compute_max_mean(data_vec, k)
    print('(%.5f , %d) done' %(mode_parameters[i_fit], mode_bin[i_fit]))
    sys.stdout.flush()
    
  return mode_bin, mode_parameters


# ==============================================================================

def get_sample_list(sample_str, parameter_names):
  
  nfit = np.shape(parameter_names)[0]
  print(' input sample to select = ', sample_str)
  str_list = sample_str.strip().split()

  name_par = parameter_names[0]
  name_excluded = None
  if (len(str_list) == 1):
    if (str_list[0].lower() == 'none'):
      name_par = parameter_names[0]
      name_excluded = None
    else:
      for ifit in range(0, nfit):
        if (str_list[0].lower() == parameter_names[ifit].lower()):
          name_par = parameter_names[ifit]
          name_excluded = None
  else:
    name_excluded = []
    for ilist in range(0, len(str_list)):
      for ifit in range(0, nfit):
        if (str_list[ilist].lower() == parameter_names[ifit].lower()):
          if (ilist == 0):
            name_par = parameter_names[ifit]
          else:
            name_excluded.append(parameter_names[ifit])

  return name_par, name_excluded


# ==============================================================================

def check_sample_parameters(parameter_names, sample_parameters, post_ci, name_excluded=None):
  
  nfit = np.shape(sample_parameters)[0]

  check_sample = True
  for ifit in range(0, nfit):
    if (name_excluded is None):
      if (sample_parameters[ifit] < post_ci[0, ifit]):
        check_sample = False
        break
      elif (sample_parameters[ifit] > post_ci[1, ifit]):
        check_sample = False
        break
    else:
      if (parameter_names[ifit] not in name_excluded):
        if (sample_parameters[ifit] < post_ci[0, ifit]):
          check_sample = False
          break
        elif (sample_parameters[ifit] > post_ci[1, ifit]):
          check_sample = False
          break

  return check_sample


# ==============================================================================

def pick_sample_parameters(posterior, parameter_names, name_par=None, name_excluded=None, post_ci=None):
  
  if (name_par is not None):
    npost, nfit = np.shape(posterior)

    if (post_ci is None):
      # post_ci = np.percentile(posterior, [15.865, 84.135], axis = 0, interpolation='midpoint')
      
      #post_ci = np.percentile(posterior, [percentile_val[2], percentile_val[3]], axis=0, interpolation='midpoint')
      post_ci = compute_hdi_full(posterior).T[0:2,:]
      
      # post_ci = np.percentile(posterior, [percentile_val[4], percentile_val[5]], axis = 0, interpolation='midpoint')
      # print np.shape(post_ci)

    sel_par = 0
    for ipar in range(0, nfit):
      if (name_par.lower() == parameter_names[ipar].lower()):
        sel_par = ipar
        break
    # print name_par, ' -> ',sel_par,': ',parameter_names[ipar]

    # use median
    # get idx sorted of the selected parameter-posterior
    idx_posterior = np.argsort(posterior[:, sel_par])
    # define the number of testing sample given the values within credible intervals
    n_test = int(0.5*np.sum(\
      np.logical_and(posterior[:, sel_par] >= post_ci[0, sel_par],
                     posterior[:, sel_par] <= post_ci[1, sel_par]).astype(int)))
    # create the list with the sorted idx, starting from the 50th and moving of +1,-1 every time
    testing = []
    testing.append(0)
    for itest in range(0, n_test):
      testing.append(itest + 1)
      testing.append(-(itest + 1))
    # print n_test
    # print testing
    n_testing = len(testing)
    # print n_testing
    sel_idx = [int(0.5 * npost) + testing[ii] for ii in range(0, n_testing)]
    # print sel_idx

    ## use mode
    # k = np.ceil(2. * posterior.shape[0]**(1./3.)).astype(int)
    # if(k>50): k=50
    # hist_counts, bin_edges = np.histogram(posterior[:,sel_par], bins=k)
    # max_bin = np.argmax(np.array(hist_counts))
    # if (max_bin == 0 or max_bin == k):
    # ext_bin = 0
    # elif (max_bin == 1 or max_bin == k-1):
    # ext_bin = 1
    # else:
    # ext_bin = 2
    # idx_posterior = np.arange(0,posterior.shape[0],1)
    # sel_post = np.logical_and(posterior[:,sel_par]>=bin_edges[max_bin-ext_bin], posterior[:,sel_par]<=bin_edges[max_bin+ext_bin])
    # sel_idx = idx_posterior[sel_post]
    # n_testing = np.shape(sel_idx)[0]

    for ii in range(0, n_testing):
      # w/ median
      idx = idx_posterior[sel_idx[ii]]
      ## w/ mode
      # idx = sel_idx[ii]

      sample_parameters = posterior[idx, :]
      check_sample = True
      check_sample = check_sample_parameters(parameter_names, sample_parameters, post_ci, name_excluded)
      if (check_sample):
        print(' found good sample parameters at index: ', idx)
        print(sample_parameters)
        return sample_parameters, idx

    return None, None


# =============================================================================

# 1) select parameters and lgllhd within ci
def select_within_all_ci(posterior, post_ci, lnprobability):
  
  npost, nfit = np.shape(posterior)
  # nrow, ncol = np.shape(post_ci)

  #if (nrow == 2 and ncol == nfit):  # check if post_ci is (nfit, 2) or (2, nfit)
    #use_ci = post_ci.T
  #else:
    #use_ci = post_ci
  use_ci = post_ci

  # that's old
  #  0       1       2       3       4       5
  # -3sigma -2sigma -1sigma +1sigma +2sigma +3sigma
  
  #print np.shape(use_ci)
  
  # use_ci should have: nfit x nci, 
  # where nci:
  # -1sigma(0) +1sigma(1) -2sigma(2) +2sigma(3) -3sigma(4) +3sigma(5)
  ok_sel = np.ones((npost)).astype(bool)
  for ifit in range(0, nfit):
    #ok_temp = np.logical_and(posterior[:, ifit] >= use_ci[ifit, 2], 
                             #posterior[:, ifit] <= use_ci[ifit, 3]
                             #)
    ok_temp = np.logical_and(posterior[:, ifit] >= use_ci[ifit, 0], 
                             posterior[:, ifit] <= use_ci[ifit, 1]
                             )
    ok_sel = np.logical_and(ok_sel, ok_temp)

  # print np.shape(posterior), np.sum(ok_sel.astype(int)), np.shape(lnprobability)
  post_sel = posterior[ok_sel, :]
  lnprob_sel = lnprobability[ok_sel]

  return post_sel, lnprob_sel


# ==============================================================================

def select_maxlglhd_with_hdi(posterior, post_ci, lnprobability):
  
  npost, _ = np.shape(posterior)
  post_sel, lnprob_sel = select_within_all_ci(posterior, 
                                              post_ci, 
                                              lnprobability.reshape((npost))
                                              )
  idmax = np.argmax(lnprob_sel)
  
  sample_par = post_sel[idmax,:]
  sample_lgllhd = lnprob_sel[idmax]
  
  return sample_par, sample_lgllhd

# ==============================================================================

# 2) determine the median lgllhd lgllhd_med of the selected pars
# 3) sort by lgllhd-lgllhd_med and take the first set of pars
def get_sample_by_sorted_lgllhd(posterior, lnprobability, post_ci=None):
  
  if (post_ci is None):
    #post_ci = np.percentile(posterior, [percentile_val[2], percentile_val[3]],
                            #axis=0, interpolation='midpoint'
                            #)
    post_ci = compute_hdi_full(posterior)
                            

  npost, _ = np.shape(posterior)

  use_lnprob = lnprobability.reshape((npost))

  post_sel, lnprob_sel = select_within_all_ci(posterior, post_ci, use_lnprob)

  lgllhd_med = np.percentile(use_lnprob, 50., interpolation='midpoint')
  abs_dlg = np.abs(lnprob_sel - lgllhd_med)

  idx_sample = np.argmin(abs_dlg)
  sample_parameters = post_sel[idx_sample, :]
  sample_lgllhd = lnprob_sel[idx_sample]

  return sample_parameters, sample_lgllhd


# ==============================================================================

def get_sample_by_par_and_lgllhd(posterior, lnprobability, parameter_names, post_ci=None, name_par=None):
  
  if (post_ci is None):
    #post_ci = np.percentile(posterior, [percentile_val[2], percentile_val[3]], axis=0, interpolation='midpoint')
    post_ci = compute_hdi_full(posterior)

  npost, nfit = np.shape(posterior)

  if (name_par is None):
    sel_par = 0
  else:
    for ifit in range(0, nfit):
      if (name_par == parameter_names[ifit].strip()):
        sel_par = ifit
        break

  use_lnprob = lnprobability.reshape((npost))

  post_sel, lnprob_sel = select_within_all_ci(posterior, post_ci, use_lnprob)

  lgllhd_med = np.percentile(use_lnprob, 50., interpolation='midpoint')
  abs_dlg = np.abs(lnprob_sel - lgllhd_med)
  lgllhd_mad = np.percentile(abs_dlg, 50., interpolation='midpoint')
  lgllhd_min = lgllhd_med - lgllhd_mad
  lgllhd_max = lgllhd_med + lgllhd_mad

  par_posterior = posterior[:, sel_par]
  par_med = np.percentile(par_posterior, 50., interpolation='midpoint')

  abs_dpar = np.abs(post_sel[:, sel_par] - par_med)
  ids_par = np.argsort(abs_dpar)
  lgllhd_check = np.logical_and(lnprob_sel >= lgllhd_min, lnprob_sel <= lgllhd_max)
  sample_parameters, sample_lgllhd = None, None
  for ii in range(0, np.shape(ids_par)[0]):
    idx = ids_par[ii]
    
    # temp_parameters = post_sel[idx, :]
    # temp_lgllhd = lnprob_sel[idx]

    # if(abs_dlg[idx] <= lgllhd_mad):
    # if(temp_lgllhd >= lgllhd_min):
    
    if (lgllhd_check[idx]):
      sample_parameters = post_sel[idx, :]
      sample_lgllhd = lnprob_sel[idx]
      break

  return sample_parameters, sample_lgllhd

# ==============================================================================

def take_n_samples(posterior, lnprob=None, post_ci=None, n_samples=100):
  
  # post_ci must have fitted parameters as cols and
  # lower and upper ci as row, i.e.:
  # post_ci(2, nfit)

  npost, nfit = np.shape(posterior)

  # idx_post = np.arange(0, npost).astype(int)

  # if(lnprob is not None and post_ci is None):
  #   idx_lnp = np.argsort(lnprob)[::-1] # sort: descending
  #   idx_sample = idx_lnp[:n_samples]

  # elif(lnprob is None and post_ci is not None):
  #   sel_within_ci = np.ones((npost)).astype(bool)
  #   for ifit in range(0,nfit):
  #     sel_par = np.logical_and(posterior[:,ifit] >= post_ci[0,ifit], posterior[:,ifit] <= post_ci[1,ifit])
  #     sel_within_ci = np.logical_and(sel_within_ci, sel_par)
    
  #   idx_within_ci = idx_post[sel_within_ci]
  #   idx_sample = np.random.choice(idx_within_ci, n_samples, replace=False)
  
  # else:
  #   idx_sample = np.random.choice(idx_post, n_samples, replace=False)
  
  # sample_parameters = posterior[idx_sample, :]

  if(lnprob is not None):
    idx_post = np.argsort(lnprob)[::-1]  # sort: descending
  else:
    idx_post = np.arange(0, npost).astype(int)

  post_sort = posterior[idx_post, :] # sort properly so then selection matches indexes
  
  if(post_ci is not None):
    # set selection withing hdi (ci) to True
    sel_within_ci = np.ones((npost)).astype(bool)
    # loop on parameters and check if within hdi (ci)
    for ifit in range(0,nfit):
      sel_par = np.logical_and(post_sort[:,ifit] >= post_ci[0,ifit], post_sort[:,ifit] <= post_ci[1,ifit])
      sel_within_ci = np.logical_and(sel_within_ci, sel_par)
  
    idx_sample = np.random.choice(idx_post[sel_within_ci], n_samples, replace=False)
  else:
    idx_sample = np.random.choice(idx_post, n_samples, replace=False)

  sample_parameters = posterior[idx_sample, :]

  return sample_parameters

# =============================================================================

def print_parameters_nolog(parameter_names, parameters, perc68, confint, par_type):
  
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  nsigmas = ['-1', '-2', '-3', '+1', '+2', '+3']

  f_w = 23
  # f_d = 16

  print()
  line_68 = '+/-1sigma(68.27|res|)'.rjust(f_w)
  line_s = '%s' % (' '.join(
    ['%2ssigma(%5.2fres)'.rjust(f_w) % (nsigmas[j], sigmas_percentiles[j]) for j in range(len(sigmas_percentiles))]))
  header = '%s %23s %s %s' % ('#'.ljust(10), par_type, line_68, line_s)
  print(header)

  for i_fit in range(0, parameters.shape[0]):
    sigma_line = ''.join(['%23.16f' % (confint[j_perc, i_fit]) for j_perc in range(len(sigmas_percentiles))])
    line = '%10s %23.16f %23.16f %s' % (parameter_names[i_fit], parameters[i_fit], perc68[i_fit], sigma_line)
    print(line)
  print()

  return


# ==============================================================================

def print_parameters_logtxt(of_run, parameter_names, parameters, perc68, confint, par_type):
  
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  nsigmas = ['-1', '-2', '-3', '+1', '+2', '+3']

  f_w = 23
  # f_d = 16

  print()
  of_run.write('\n')
  line_68 = '+/-1sigma(68.27|res|)'.rjust(f_w)
  line_s = '%s' % (' '.join(
    ['%2ssigma(%5.2fres)'.rjust(f_w) % (nsigmas[j], sigmas_percentiles[j]) for j in range(len(sigmas_percentiles))]))
  header = '%s %23s %s %s' % ('#'.ljust(10), par_type, line_68, line_s)
  print(header)
  of_run.write(header + '\n')

  for i_fit in range(0, parameters.shape[0]):
    sigma_line = ''.join(['%23.16f' % (confint[j_perc, i_fit]) for j_perc in range(len(sigmas_percentiles))])
    line = '%10s %23.16f %23.16f %s' % (parameter_names[i_fit], parameters[i_fit], perc68[i_fit], sigma_line)
    print(line)
    of_run.write(line + '\n')
  print()
  of_run.write('\n')

  return


# ==============================================================================

def print_parameters_logger(logger, parameter_names, parameters, perc68, confint, par_type):
  
  sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
  nsigmas = ['-1', '-2', '-3', '+1', '+2', '+3']

  f_w = 23
  # f_d = 16

  logger.info('')
  line_68 = '+/-1sigma(68.27|res|)'.rjust(f_w)
  line_s = '%s' % (' '.join(
    ['%2ssigma(%5.2fres)'.rjust(f_w) % (nsigmas[j], sigmas_percentiles[j]) for j in range(len(sigmas_percentiles))]))
  header = '%s %23s %s %s' % ('#'.ljust(10), par_type, line_68, line_s)
  logger.info(header)

  for i_fit in range(0, parameters.shape[0]):
    sigma_line = ''.join(['%23.16f' % (confint[j_perc, i_fit]) for j_perc in range(len(sigmas_percentiles))])
    line = '%10s %23.16f %23.16f %s' % (parameter_names[i_fit], parameters[i_fit], perc68[i_fit], sigma_line)
    logger.info(line)
  logger.info('')

  return


# =============================================================================

# OLD
def get_derived_posterior_parameters(parameter_names, chains_T, flatchain_posterior):
  
  nfit = flatchain_posterior.shape[1]
  derived_names = []
  derived_chains = []
  derived_posterior = []

  for i_fit in range(0, nfit):
    if ('ecosw' in parameter_names[i_fit]):
      body_id = parameter_names[i_fit].split('ecosw')[1]
      derived_names.append('e%s' % (body_id))
      derived_names.append('w%s' % (body_id))

      temp_e, temp_w = 0., 0.

      temp_e = np.sqrt(chains_T[:, :, i_fit] ** 2 + chains_T[:, :, i_fit + 1] ** 2)
      temp_w = np.arctan2(chains_T[:, :, i_fit + 1], chains_T[:, :, i_fit]) * 180. / np.pi
      derived_chains.append(temp_e)
      derived_chains.append(temp_w)

      temp_e, temp_w = 0., 0.

      temp_e = np.sqrt(flatchain_posterior[:, i_fit] ** 2 + flatchain_posterior[:, i_fit + 1] ** 2)
      temp_w = np.arctan2(flatchain_posterior[:, i_fit + 1], flatchain_posterior[:, i_fit]) * 180. / np.pi
      derived_posterior.append(temp_e)
      derived_posterior.append(temp_w)

    if ('icoslN' in parameter_names[i_fit]):
      body_id = parameter_names[i_fit].split('icoslN')[1]
      derived_names.append('i%s' % (body_id))
      derived_names.append('lN%s' % (body_id))

      temp_i, temp_lN = 0., 0.

      temp_i = np.sqrt(chains_T[:, :, i_fit] ** 2 + chains_T[:, :, i_fit + 1] ** 2)
      temp_lN = np.arctan2(chains_T[:, :, i_fit + 1], chains_T[:, :, i_fit]) * 180. / np.pi
      derived_chains.append(temp_i)
      derived_chains.append(temp_lN)

      temp_i, temp_lN = 0., 0.

      temp_i = np.sqrt(flatchain_posterior[:, i_fit] ** 2 + flatchain_posterior[:, i_fit + 1] ** 2)
      temp_lN = np.arctan2(flatchain_posterior[:, i_fit + 1], flatchain_posterior[:, i_fit]) * 180. / np.pi
      derived_posterior.append(temp_i)
      derived_posterior.append(temp_lN)

  nder = len(derived_names)
  derived_chains_T = np.zeros((chains_T.shape[0], chains_T.shape[1], nder))
  for i_der in range(0, nder):
    derived_chains_T[:, :, i_der] = np.array(derived_chains[i_der])

  return derived_names, derived_chains_T, np.array(derived_posterior).T


# ==============================================================================

# OLD
def save_posterior_like(out_folder, boot_id, parameter_names, flatchain_posterior):
  
  nfit = flatchain_posterior.shape[1]
  nboot = flatchain_posterior.shape[0]
  header_0000 = ' iboot %s' % (' '.join(parameter_names))
  fmt_dp = ' '.join(['%s' % ('%23.16e') for ii in range(0, nfit)])
  fmt_full = '%s %s' % ('%6d', fmt_dp)
  iboot_fake = np.arange(1, nboot + 1, 1)
  boot_fake = np.column_stack((iboot_fake, flatchain_posterior))
  boot_file = os.path.join(out_folder, '%d_posterior_sim.dat' % (boot_id))
  np.savetxt(boot_file, boot_fake, fmt=fmt_full, header=header_0000)

  return boot_file


# =============================================================================

def GelmanRubin(chains_T):
  
  n, M = np.shape(chains_T)

  theta_m = [np.mean(chains_T[:, i_m]) for i_m in range(0, M)]
  theta = np.mean(theta_m)

  d_theta2 = (theta_m - theta) ** 2
  B_n = np.sum(d_theta2) / (M - 1)

  arg_W = [np.sum((chains_T[:, i_m] - theta_m[i_m]) ** 2) / (n - 1) for i_m in range(0, M)]
  W = np.mean(arg_W)

  n_frac = (n - 1) / n
  var_plus = n_frac * W + B_n
  Var = var_plus + (B_n / M)

  Rc = np.sqrt(Var / W)

  return Rc

# ==============================================================================

# Geweke 1992 test 1
def geweke_test(chains_T, start_frac=0.05, n_sel_steps=5):
  
  n_steps, n_chains = chains_T.shape
  half = int(0.5 * n_steps)

  int_b = chains_T[half:, :]
  mean_b = np.mean(int_b, axis=0)
  var_b = np.var(int_b, axis=0, ddof=1)

  start_step = int(start_frac * n_steps)
  sel_a = np.linspace(start=start_step, stop=half, num=n_sel_steps, endpoint=False, dtype=np.int)

  zscore = np.zeros((n_sel_steps, n_chains))
  for i_a in range(0, n_sel_steps):
    int_a = chains_T[sel_a[i_a]:half, :]
    mean_a = np.mean(int_a, axis=0)
    var_a = np.var(int_a, axis=0, ddof=1)
    z_temp = (mean_a - mean_b) / np.sqrt(var_a + var_b)
    zscore[i_a, :] = z_temp

  return sel_a, zscore


# =============================================================================

# ==============================================================================

def get_units(names, mass_unit):
  
  n_names = len(names)
  units_par = []

  for i in range(0, n_names):

    if (str(names[i])[0] == 'm'):
      if ('Ms' in names[i]):
        units_par.append('[M_sun/M_star]')
      elif ('mA' in names[i]):
        units_par.append('[deg]')
      else:
        units_par.append('[%s]' % (mass_unit))

    if ('R' in str(names[i])):
      units_par.append('[R_sun]')

    if ('P' in str(names[i])):
      units_par.append('[days]')

    if (str(names[i])[0] == 'e'):  # the same for e, ecosw, esinw
      units_par.append(' ')

    if (str(names[i])[0] == 's'):  # sqrtecosw, sqrtesinw
      units_par.append(' ')

    if (str(names[i])[0] == 'w'):
      units_par.append('[deg]')

    if ('lambda' in str(names[i])):
      units_par.append('[deg]')

    if (str(names[i])[0] == 'i'):
      if ('cos' in names[i] or 'sin' in names[i]):
        units_par.append(' ')
      else:
        units_par.append('[deg]')

    if (str(names[i])[0:2] == 'lN'):
      units_par.append('[deg]')

  return units_par


# ==============================================================================

def elements(fpath, idsim, lmf=0):
  
  # kel_file = os.path.join(fpath, str(idsim) + "_" + str(lmf) + "_initialElements.dat")
  kel_file = os.path.join(fpath, '%d_%d_initialElements.dat' % (idsim, lmf))
  try:
    kep_elem = np.genfromtxt(kel_file)  # (NB-1) x (M_Msun R_Rsun P_d a_AU ecc w_deg mA_deg inc_deg lN_deg)
  except:
    print(' KEPLERIAN ELEMENTS FILE NOT FOUND %s' % (kel_file))
    sys.exit()

  return kel_file, kep_elem


# ==============================================================================

def get_case(id_body, fit_body):
  
  fit_n = np.arange(1, 11, 1)
  id_fit = [fit_n[j] for j in range(len(fit_body)) if (fit_body[j] == '1')]
  id_all = [8 * (id_body - 1) + 2 + id_fit[j] for j in range(len(id_fit))]

  if (6 not in id_fit):  # not fitting lambda
    case = [0]

  else:  # (6 in id_fit) # fitting lambda
    if (all(x in id_fit for x in [4, 5, 6, 7, 8])):  # fit ecosw, esinw, lambda, icoslN, isinlN
      case = [4]
    elif (all(x in id_fit for x in [4, 5, 6])):  # fit ecosw, esinw, lambda
      case = [2]
    elif (all(x in id_fit for x in [6, 7, 8])):  # fit lambda, icoslN, isinlN
      case = [3]
    else:
      case = [1]  # fit lambda & w || lambda & lN || lamda & w & lN

  return id_fit, id_all, case


# ==============================================================================

def get_fitted(full_path):
  
  of = open(os.path.join(full_path, 'bodies.lst'), 'r')
  lines = of.readlines()
  of.close()
  NB = len(lines)
  # n_pl = NB - 1
  bodies_file = []
  fit_string = ''
  fit_list = []
  for line in lines:
    fname = line.strip().split()[0]
    bodies_file.append(fname)
    temp_string = line.strip().split(fname)[1].split('#')[0].strip()
    if ('F' in temp_string):
      temp_string = temp_string.split('F')[0].strip()
    elif ('f' in temp_string):
      temp_string = temp_string.split('f')[0].strip()
    elif ('T' in temp_string):
      temp_string = temp_string.split('T')[0].strip()
    elif ('t' in temp_string):
      temp_string = temp_string.split('t')[0].strip()

    fit_string = '%s %s' % (fit_string, temp_string)
    fit_list.append(temp_string.split())

  nfit = np.sum(np.array(fit_string.split(), dtype=np.int))

  case = []
  id_fit = []
  id_all = []
  nfit_list = []
  cols_list = []
  nfit_cnt = 0
  for j in range(1, NB + 1):
    id_fit_j, id_all_j, case_j = get_case(j, fit_list[j - 1])
    id_fit.append(id_fit_j)
    id_all.append(id_all_j)
    case.append(case_j)
    nfit_j = len(id_fit_j)
    nfit_list.append(nfit_j)
    cols_list.append([cc for cc in range(nfit_cnt, nfit_cnt + nfit_j)])
    nfit_cnt += nfit_j
  # cols = [jj for jj in range(0, nfit)]

  return nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case


# ==============================================================================

def compute_intervals(flatchain, parameters, percentiles):
  
  sigma_par = np.percentile(np.subtract(flatchain, parameters), percentiles, axis=0,
                            interpolation='midpoint')  # (n_percentile x nfit)
  sigma_par[0] = np.percentile(np.abs(np.subtract(flatchain, parameters)), percentiles[0], axis=0,
                               interpolation='midpoint')  # 68.27th
  sigma_par[1] = np.percentile(np.abs(np.subtract(flatchain, parameters)), percentiles[1], axis=0,
                               interpolation='midpoint')  # MAD
  # sigma_par = np.percentile(flatchain, parameters, percentiles, axis=0) # (n_percentile x nfit)
  # sigma_par[0] = np.percentile(np.abs(flatchain, parameters), percentiles[0], axis=0) # 68.27th
  # sigma_par[1] = np.percentile(np.abs(flatchain, parameters), percentiles[1], axis=0) # MAD

  return sigma_par

# ==============================================================================

def compute_hdi_full(flatchains, mode_output=False):
  
  alpha=[0.3173, 0.0456, 0.0026]
  _, npar = np.shape(flatchains)
  #print '\n^^^\nIN compute_hdi_full WITH (npost , npar) = (%d , %d)' %(npost, npar)
  
  nbins = get_auto_bins(flatchains)
  #print 'nbins = ',nbins
  
  #print 'Computing hdi_full ....',
  hdi_full = [(calculate_hdi(flatchains[:,ipar],
               nbins,
               alpha,
               mode_output=mode_output
               )
               ) for ipar in range(npar)
             ]
  #print 'done with shape ', np.shape(hdi_full), ' and type ',type(hdi_full)
  #print hdi_full
  
  #print 'Reformatting hdi_full to hdi_l:'
  hdi_l = [[hdi_full[ipar][0][0], hdi_full[ipar][0][1], # -1sigma +1sigma
            hdi_full[ipar][1][0], hdi_full[ipar][1][1], # -2sigma +2sigma
            hdi_full[ipar][2][0], hdi_full[ipar][2][1]  # -3sigma +3sigma
          ] for ipar in range(npar)]
  #print 'hdi_l with shape ',np.shape(hdi_l), ' and type ',type(hdi_l)
  
  hdi_a = np.array(hdi_l, dtype=np.float64)
  #print 'hdi_a with shape ',np.shape(hdi_a), ' and type ',type(hdi_a)
  
  hdi = np.reshape(hdi_a, newshape=((npar,-1)))
  #print 'hdi   with shape ',np.shape(hdi), ' and type ',type(hdi)
  
  if(mode_output):
    mode = np.array([hdi_full[ipar][3] for ipar in range(npar)], dtype=np.float64)
    #print '^^^\n'
    return hdi, mode
  
  #print '^^^\n'
  
  return hdi

# ==============================================================================

def compute_sigma_hdi(flatchains, parameters):
  
  alpha=[0.3173, 0.0456, 0.0026]
  
  _, npar = np.shape(flatchains)
  #print '\n^^^\nIN compute_sigma_hdi WITH (npost , npar) = (%d , %d)\n^^^\n' %(npost, npar)
  nbins = get_auto_bins(flatchains)
  
  hdi_full = [(calculate_hdi(flatchains[:,ipar],
                            nbins, alpha, 
                            mode_output=False
                            )
                      ) for ipar in range(npar)]
  
  sigma_par = np.array([[hdi_full[ipar][0][0] - parameters[ipar],
                         hdi_full[ipar][0][1] - parameters[ipar],
                         hdi_full[ipar][1][0] - parameters[ipar], 
                         hdi_full[ipar][1][1] - parameters[ipar],
                         hdi_full[ipar][2][0] - parameters[ipar], 
                         hdi_full[ipar][2][1] - parameters[ipar]
                      ] for ipar in range(npar)]
                    ).reshape((npar,-1))
  
  delta_flat_T = np.abs([flatchains[:,ipar] - parameters[ipar] for ipar in range(npar)]).T
  
  sigma_r68p = np.percentile(delta_flat_T, 68.27,
                             axis=0, interpolation='midpoint')  # 68.27th
  
  sigma_mad = np.percentile(delta_flat_T, 50.,
                            axis=0, interpolation='midpoint')  # MAD
  
  sigma_par = np.column_stack((sigma_r68p, sigma_mad, sigma_par))
  
  return sigma_par

# ==============================================================================

def get_good_distribution(posterior_scale, posterior_mod):


  par_scale = np.median(posterior_scale)
  p68_scale = np.percentile(np.abs(posterior_scale - par_scale),
                            68.27, interpolation='midpoint'
                            )
  std_scale = np.std(posterior_scale, ddof=1)
  s_scale = max(p68_scale, std_scale)

  par_mod = np.median(posterior_mod)
  p68_mod = np.percentile(np.abs(posterior_mod - par_mod),
                          68.27, interpolation='midpoint'
                          )
  std_mod = np.std(posterior_mod, ddof=1)
  s_mod = max(p68_mod, std_mod)

  if(s_scale < s_mod):
    par_out = par_scale
    posterior_out = posterior_scale
  else:
    par_out = par_mod
    posterior_out = posterior_mod

  return par_out, posterior_out

# ==============================================================================

def get_proper_posterior_correlated(posterior, col=None):
  
  if (col is not None):
    posterior_scale = np.arctan2(posterior[:, col + 1], posterior[:, col]) * rad2deg
    posterior_mod = posterior_scale % 360.
  else:
    posterior_scale = np.arctan2(np.sin(posterior * deg2rad), np.cos(posterior * deg2rad)) * rad2deg
    posterior_mod = posterior_scale % 360.

  par_out, post_out = get_good_distribution(posterior_scale, posterior_mod)
  
  return par_out, post_out

# ==============================================================================
# ==============================================================================
# fix lambda only in flatchain posterior
def fix_lambda(flatchain_post, names_par):
  
  _, nfit = np.shape(flatchain_post)
  flatchain_post_out = flatchain_post.copy()
  for ifit in range(0,nfit):
    if('lambda' in names_par[ifit]):
      _, xout = \
        get_proper_posterior_correlated(flatchain_post[:,ifit])
      flatchain_post_out[:,ifit] = xout.copy()
  
  return flatchain_post_out

# ==============================================================================
# ==============================================================================

def get_trigonometric_posterior(id_NB, id_l, col, idpar,
                                posterior):  # ecosw,esinw -> (e,w) && icoslN,isinlN -> (i,lN) || No lambda

  cosboot = posterior[:, col]
  sinboot = posterior[:, col + 1]
  sum_square = (cosboot * cosboot) + (sinboot * sinboot)
  if (idpar[id_l][0:6] == 'sqrtec'):
    aboot = sum_square
  else:
    aboot = np.sqrt(sum_square)
  # aboot = np.sqrt( (cosboot*cosboot) + (sinboot*sinboot) )
  _, bboot = get_proper_posterior_correlated(posterior, col)
  if ('e' in idpar[id_l]):
    sel_ezero = aboot <= eps64bit
    bboot[sel_ezero] = 90.

  if ('ecosw' in idpar[id_l]):
    id_aa = 'e'
    id_bb = 'w'
  else:
    id_aa = 'i'
    id_bb = 'lN'

  id_a0 = '%s%s' % (id_aa, id_NB + 1)
  id_b0 = '%s%s' % (id_bb, id_NB + 1)

  return id_a0, aboot, id_b0, bboot


# ==============================================================================

def get_trigonometric_parameters(id_NB, id_l, col, idpar,
                                 parameters):  # ecosw,esinw // sqrt(e)cosw,sqrt(e)sinw -> (e,w) && icoslN,isinlN -> (i,lN) || No lambda

  cosp = parameters[col]
  sinp = parameters[col + 1]
  sum_square = (cosp * cosp) + (sinp * sinp)
  if (idpar[id_l][0:6] == 'sqrtec'):
    a_par = sum_square
  else:
    a_par = np.sqrt(sum_square)
  # a_par = np.sqrt( (cosp*cosp) + (sinp*sinp) )
  b_par = np.arctan2(sinp, cosp) * rad2deg

  if ('e' in idpar[id_l]):
    if (a_par <= eps64bit):
      b_par = 90.

  if ('ecosw' in idpar[id_l]):
    id_aa = 'e'
    id_bb = 'w'
  else:
    id_aa = 'i'
    id_bb = 'lN'

  id_a0 = '%s%s' % (id_aa, id_NB + 1)
  id_b0 = '%s%s' % (id_bb, id_NB + 1)

  return id_a0, a_par, id_b0, b_par


# ==============================================================================

def derived_posterior_case0(idpar_NB, id_fit_NB, i_NB, cols, posterior):  # not fitting lambda

  name_der = []
  der_posterior = []

  nlist = len(id_fit_NB)
  for i_l in range(nlist):
    if ('ecosw' in idpar_NB[i_l] or 'icoslN' in idpar_NB[i_l]):
      id_a0, aboot, id_b0, bboot = get_trigonometric_posterior(i_NB, i_l, cols[i_l], idpar_NB, posterior)
      name_der.append(id_a0)
      name_der.append(id_b0)
      der_posterior.append(aboot)
      der_posterior.append(bboot)

  return name_der, der_posterior


# ==============================================================================

def derived_parameters_case0(idpar_NB, id_fit_NB, i_NB, cols, parameters):  # not fitting lambda

  name_der = []
  der_par = []

  nlist = len(id_fit_NB)
  for i_l in range(nlist):
    if ('ecosw' in idpar_NB[i_l] or 'icoslN' in idpar_NB[i_l]):
      id_a0, a_par, id_b0, b_par = get_trigonometric_parameters(i_NB, i_l, cols[i_l], idpar_NB, parameters)
      name_der.append(id_a0)
      name_der.append(id_b0)
      der_par.append(a_par)
      der_par.append(b_par)

  return name_der, der_par


# ==============================================================================

def derived_posterior_case1(idpar_NB, id_fit_NB, i_NB, cols, kep_elem,
                            posterior):  # fitting lambda || lambda & w || lambda & lN || lamda & w & lN

  name_der = []
  der_posterior = []

  nlist = len(id_fit_NB)
  # for i_l in range(0,nlist):

  # if(id_fit_NB[i_l] == 6):
  # if(idpar_NB[i_l+2][0:2] != 'mA'):

  # aboot = (posterior[:,cols[i_l]] - kep_elem[5] - kep_elem[8])%360.
  # temp, aboot = get_proper_posterior_correlated(aboot)

  # name_der.append('%s%s' %('mA',i_NB+1))
  # der_posterior.append(aboot)

  for i_l in range(0, nlist):
    if (id_fit_NB[i_l] == 6):  # lambda
      lambda_post = posterior[:, cols[i_l]]
      if (i_l > 0):

        if (i_l - 1 == 5):  # id 5 == argp
          argp = posterior[:, cols[i_l - 1]]  # argp as posterior
        else:
          argp = kep_elem[5]  # argp is fixed

        if (i_l < nlist - 1 and i_l + 1 == nlist - 1):  # id 6 is not the last --> 8 is the next (longn)
          longn = posterior[:, cols[i_l + 1]]  # long node as posterior
        else:
          # if(i_l == nlist-1):
          longn = kep_elem[8]  # long of node is fixed

      else:
        argp = kep_elem[5]  # argp is fixed
        longn = kep_elem[8]  # longn is fixed

      aboot = (lambda_post - argp - longn) % 360.
      _, aboot = get_proper_posterior_correlated(aboot)
      name_der.append('%s%s' % ('mA', i_NB + 1))
      der_posterior.append(aboot)

  return name_der, der_posterior


# ==============================================================================

def derived_parameters_case1(idpar_NB, id_fit_NB, i_NB, cols, kep_elem,
                             parameters):  # fitting lambda || lambda & w || lambda & lN || lamda & w & lN
  name_der = []
  der_par = []

  nlist = len(id_fit_NB)
  # for i_l in range(0,nlist):

  # if(id_fit_NB[i_l] == 6):
  # if(idpar_NB[i_l+2][0:2] != 'mA'):

  # a_par = (parameters[cols[i_l]] - kep_elem[5] - kep_elem[8])%360.
  # name_der.append('%s%s' %('mA',i_NB+1))
  # der_par.append(a_par)

  for i_l in range(0, nlist):
    if (id_fit_NB[i_l] == 6):  # lambda
      lambda_par = parameters[cols[i_l]]
      if (i_l > 0):

        if (i_l - 1 == 5):  # id 5 == argp
          argp = parameters[cols[i_l - 1]]  # argp as parameter
        else:
          argp = kep_elem[5]  # argp is fixed

        if (i_l < nlist - 1 and i_l + 1 == nlist - 1):  # id 6 is not the last --> 8 is the next (longn)
          longn = parameters[cols[i_l + 1]]  # long node as posterior
        else:
          # if(i_l == nlist-1):
          longn = kep_elem[8]  # long of node is fixed

      else:
        argp = kep_elem[5]  # argp is fixed
        longn = kep_elem[8]  # longn is fixed

      a_par = (lambda_par - argp - longn) % 360.
      name_der.append('%s%s' % ('mA', i_NB + 1))
      der_par.append(a_par)

  return name_der, der_par


# ==============================================================================

def derived_posterior_case2(idpar_NB, id_fit_NB, i_NB, cols, kep_elem,
                            posterior):  # fit ecosw//sqrt(e)cosw, esinw//sqrt(e)sinw, lambda

  name_der = []
  der_posterior = []

  nlist = len(id_fit_NB)
  for i_l in range(0, nlist):

    if (id_fit_NB[i_l] == 4):
      id_e0, eboot, id_w0, wboot = get_trigonometric_posterior(i_NB, i_l, cols[i_l], idpar_NB, posterior)
      name_der.append(id_e0)
      name_der.append(id_w0)
      der_posterior.append(eboot)
      der_posterior.append(wboot)

      if (idpar_NB[i_l + 2][0:2] != 'mA'):
        mAboot = (posterior[:, cols[i_l + 2]] - wboot - kep_elem[8]) % 360.  # mA = lambda - w - lN
        _, mAboot = get_proper_posterior_correlated(mAboot)
        name_der.append('%s%s' % ('mA', i_NB + 1))
        der_posterior.append(mAboot)

  return name_der, der_posterior


# ==============================================================================

def derived_parameters_case2(idpar_NB, id_fit_NB, i_NB, cols, kep_elem,
                             parameters):  # fit ecosw//sqrt(e)cosw, esinw//sqrt(e)sinw, lambda

  name_der = []
  der_par = []

  nlist = len(id_fit_NB)
  for i_l in range(0, nlist):

    if (id_fit_NB[i_l] == 4):
      id_e0, e_par, id_w0, w_par = get_trigonometric_parameters(i_NB, i_l, cols[i_l], idpar_NB, parameters)
      name_der.append(id_e0)
      name_der.append(id_w0)
      der_par.append(e_par)
      der_par.append(w_par)

      if (idpar_NB[i_l + 2][0:2] != 'mA'):
        #print '*** lambda_par = ', parameters[cols[i_l + 2]]
        mA_par = (parameters[cols[i_l + 2]] - w_par - kep_elem[8]) % 360.  # mA = lambda - w - lN
        name_der.append('%s%s' % ('mA', i_NB + 1))
        der_par.append(mA_par)
        #print '*** mA_par = ', mA_par
  return name_der, der_par


# ==============================================================================

def derived_posterior_case3(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, posterior):  # fit lambda, icoslN, isinlN

  name_der = []
  der_posterior = []

  nlist = len(id_fit_NB)
  for i_l in range(0, nlist):

    if (id_fit_NB[i_l] == 7):
      id_i0, iboot, id_lN0, lNboot = get_trigonometric_posterior(i_NB, i_l, cols[i_l], idpar_NB, posterior)
      name_der.append(id_i0)
      name_der.append(id_lN0)
      der_posterior.append(iboot)
      der_posterior.append(lNboot)

      if (idpar_NB[i_l + 2][0:2] != 'mA'):
        mAboot = (posterior[:, cols[i_l - 1]] - kep_elem[5] - lNboot) % 360.  # mA = lambda - w - lN
        _, mAboot = get_proper_posterior_correlated(mAboot)
        name_der.append('%s%s' % ('mA', i_NB + 1))
        der_posterior.append(mAboot)

  return name_der, der_posterior


# ==============================================================================

def derived_parameters_case3(idpar_NB, id_fit_NB, i_NB, cols, kep_elem, parameters):  # fit lambda, icoslN, isinlN

  name_der = []
  der_par = []

  nlist = len(id_fit_NB)
  for i_l in range(0, nlist):

    if (id_fit_NB[i_l] == 7):
      id_i0, i_par, id_lN0, lN_par = get_trigonometric_parameters(i_NB, i_l, cols[i_l], idpar_NB, parameters)
      name_der.append(id_i0)
      name_der.append(id_lN0)
      der_par.append(i_par)
      der_par.append(lN_par)

      if (idpar_NB[i_l + 2][0:2] != 'mA'):
        mA_par = (parameters[cols[i_l - 1]] - kep_elem[5] - lN_par) % 360.  # mA = lambda - w - lN
        name_der.append('%s%s' % ('mA', i_NB + 1))
        der_par.append(mA_par)

  return name_der, der_par


# ==============================================================================

def derived_posterior_case4(idpar_NB, id_fit_NB, i_NB, cols, kep_elem,
                            posterior):  # fit ecosw//sqrt(e)cosw, esinw//sqrt(e)sinw, lambda, icoslN, isinlN

  name_der = []
  der_posterior = []

  nlist = len(id_fit_NB)
  for i_l in range(0, nlist):

    if (id_fit_NB[i_l] == 4):
      id_e0, eboot, id_w0, wboot = get_trigonometric_posterior(i_NB, i_l, cols[i_l], idpar_NB, posterior)
      name_der.append(id_e0)
      name_der.append(id_w0)
      der_posterior.append(eboot)
      der_posterior.append(wboot)

      id_i0, iboot, id_lN0, lNboot = get_trigonometric_posterior(i_NB, i_l + 3, cols[i_l + 3], idpar_NB, posterior)
      name_der.append(id_i0)
      name_der.append(id_lN0)
      der_posterior.append(iboot)
      der_posterior.append(lNboot)

      if (idpar_NB[i_l + 2][0:2] != 'mA'):
        mAboot = (posterior[:, cols[i_l + 2]] - wboot - lNboot) % 360.  # mA = lambda - w - lN
        _, mAboot = get_proper_posterior_correlated(mAboot)
        name_der.append('%s%s' % ('mA', i_NB + 1))
        der_posterior.append(mAboot)

  return name_der, der_posterior


# ==============================================================================

def derived_parameters_case4(idpar_NB, id_fit_NB, i_NB, cols, kep_elem,
                             parameters):  # fit ecosw//sqrt(e)cosw, esinw//sqrt(e)sinw, lambda, icoslN, isinlN

  name_der = []
  der_par = []

  nlist = len(id_fit_NB)
  for i_l in range(0, nlist):

    if (id_fit_NB[i_l] == 4):
      id_e0, e_par, id_w0, w_par = get_trigonometric_parameters(i_NB, i_l, cols[i_l], idpar_NB, parameters)
      name_der.append(id_e0)
      name_der.append(id_w0)
      der_par.append(e_par)
      der_par.append(w_par)

      id_i0, i_par, id_lN0, lN_par = get_trigonometric_parameters(i_NB, i_l + 3, cols[i_l + 3], idpar_NB, parameters)
      name_der.append(id_i0)
      name_der.append(id_lN0)
      der_par.append(i_par)
      der_par.append(lN_par)

      if (idpar_NB[i_l + 2][0:2] != 'mA'):
        mA_par = (parameters[cols[i_l + 2]] - w_par - lN_par) % 360.  # mA = lambda - w - lN
        name_der.append('%s%s' % ('mA', i_NB + 1))
        der_par.append(mA_par)

  return name_der, der_par


# ==============================================================================

def derived_posterior_check(derived_names_in, derived_posterior_in):
  # check distribution of angles...
  derived_posterior = derived_posterior_in.copy()
  derived_names = decode_list(derived_names_in)
  n_der = np.shape(derived_names)[0]
  for ider in range(0, n_der):
    if ('w' in derived_names[ider] or 'mA' in derived_names[ider] or 'lN' in derived_names[ider]):
      # print derived_names[ider]
      cosp = np.cos(derived_posterior_in[:, ider] * deg2rad)
      sinp = np.sin(derived_posterior_in[:, ider] * deg2rad)
      p_scale = np.arctan2(sinp, cosp) * rad2deg
      p_mod = (p_scale + 360.) % 360.
      _, temp_post = get_good_distribution(p_scale, p_mod)
      derived_posterior[:, ider] = temp_post

  return derived_posterior


# ==============================================================================

def derived_parameters_check(derived_names, derived_parameters_in, derived_posterior):
  derived_parameters = derived_parameters_in.copy()
  n_der = np.shape(derived_names)[0]
  for ider in range(0, n_der):
    # print derived_names[ider]
    if ('w' in derived_names[ider] or 'mA' in derived_names[ider] or 'lN' in derived_names[ider]):
      # print '[0]', derived_parameters[ider]
      if (np.min(derived_posterior[:, ider]) < 0.):
        # print np.min(derived_posterior[:,ider])
        if (np.max(derived_posterior[:, ider]) < 0.):
          # print np.max(derived_posterior[:,ider])
          derived_parameters[ider] = derived_parameters[ider] % -360.
      else:
        derived_parameters[ider] = derived_parameters[ider] % 360.
        # print '[1]', derived_parameters[ider]

  return derived_parameters


# ==============================================================================

def compute_derived_posterior(idpar, kep_elem_in, id_fit, case_list, cols_list, posterior, conv_factor=1.):
  NB = len(case_list)
  # nfit = len(id_fit)
  #print 'NB =', NB
  #print 'np.shape(kep_elem_in) = ',np.shape(kep_elem_in)
  
  if(NB == 2):
    kep_elem = np.array(kep_elem_in).reshape((1,-1))
  else:
    kep_elem = kep_elem_in
  
  #print 'np.shape(kep_elem) = ',np.shape(kep_elem)

  # posterior_out = posterior.copy()

  id_derived = []
  der_posterior = []

  cc_fit = 0

  for i_NB in range(0, NB):
    nfit_NB = len(id_fit[i_NB])

    if (nfit_NB > 0):

      idpar_NB = idpar[cols_list[i_NB][0]:cols_list[i_NB][-1] + 1]  # names of the parameter for the body i_NB
      id_fit_NB = id_fit[i_NB]  # integers that identify the proper type of the fitted parameters

      for i_fit in range(0, nfit_NB):

        if (id_fit[i_NB][i_fit] == 1):  # convert Mp and Mp/Ms into mass unit specified by the user <-> mtype

          if ('Ms' in idpar[cc_fit]):
            mboot = posterior[:, cc_fit] * conv_factor
            xid = 'm%d' % (i_NB + 1)
            id_derived.append(xid)
            der_posterior.append(mboot)

        cc_fit += 1

      id_temp, der_temp = [], []

      if (case_list[i_NB][0] == 0):
        id_temp, der_temp = derived_posterior_case0(idpar_NB, id_fit_NB, i_NB,
                                                    cols_list[i_NB],
                                                    posterior
                                                    )

      elif (case_list[i_NB][0] == 1):
        id_temp, der_temp = derived_posterior_case1(idpar_NB, id_fit_NB, i_NB, 
                                                    cols_list[i_NB],
                                                    kep_elem[i_NB - 1, :],
                                                    posterior
                                                    )

      elif (case_list[i_NB][0] == 2):
        id_temp, der_temp = derived_posterior_case2(idpar_NB, id_fit_NB, i_NB,
                                                    cols_list[i_NB],
                                                    kep_elem[i_NB - 1, :],
                                                    posterior
                                                    )

      elif (case_list[i_NB][0] == 3):
        id_temp, der_temp = derived_posterior_case3(idpar_NB, id_fit_NB, i_NB, 
                                                    cols_list[i_NB],
                                                    kep_elem[i_NB - 1, :],
                                                    posterior
                                                    )

      elif (case_list[i_NB][0] == 4):
        id_temp, der_temp = derived_posterior_case4(idpar_NB, id_fit_NB, i_NB,
                                                    cols_list[i_NB],
                                                    kep_elem[i_NB - 1, :],
                                                    posterior
                                                    )

      id_derived.append(id_temp)
      der_posterior.append(der_temp)

  n_der = 0
  n_xder = len(der_posterior)
  if (n_xder > 0):
    n_der_single = []
    derived_post = []
    names_derived = []

    for ii in range(0, n_xder):
      # temp_names = np.array(id_derived[ii],dtype=str)
      temp_names = id_derived[ii]
      ntemp = np.size(temp_names)
      nls = len(np.shape(temp_names))
      n_der_single.append(ntemp)
      n_der += ntemp
      # temp_der = np.array(der_posterior[ii],dtype=np.float64)
      temp_der = der_posterior[ii]

      if (ntemp > 0):

        if (ntemp == 1 and nls == 0):
          derived_post.append(temp_der)
          names_derived.append(temp_names)
        else:
          for jj in range(0, ntemp):
            derived_post.append(temp_der[jj])
            names_derived.append(temp_names[jj])

  return np.array(names_derived, dtype=str), np.array(derived_post, dtype=np.float64).T


# ==============================================================================

def compute_derived_parameters(idpar, kep_elem_in, id_fit, case_list, cols_list, parameters, conv_factor=1.):
  NB = len(case_list)
  # nfit = len(id_fit)
  #print 'NB =', NB
  #print 'np.shape(kep_elem_in) = ',np.shape(kep_elem_in)
  
  if(NB == 2):
    kep_elem = np.array(kep_elem_in).reshape((1,-1))
  else:
    kep_elem = kep_elem_in
  
  #print 'np.shape(kep_elem) = ',np.shape(kep_elem)

  id_derived = []
  der_par = []

  cc_fit = 0

  for i_NB in range(0, NB):
    nfit_NB = len(id_fit[i_NB])

    if (nfit_NB > 0):

      idpar_NB = idpar[cols_list[i_NB][0]:cols_list[i_NB][-1] + 1]  # names of the parameter for the body i_NB
      id_fit_NB = id_fit[i_NB]  # integers that identify the proper type of the fitted parameters

      for i_fit in range(0, nfit_NB):

        if (id_fit[i_NB][i_fit] == 1):  # convert Mp and Mp/Ms into mass unit specified by the user <-> mtype

          if ('Ms' in idpar[cc_fit]):
            m_par = parameters[cc_fit] * conv_factor
            xid = 'm%d' % (i_NB + 1)
            id_derived.append(xid)
            der_par.append(m_par)

        cc_fit += 1

      id_temp, der_temp = [], []

      if (case_list[i_NB][0] == 0):
        id_temp, der_temp = derived_parameters_case0(idpar_NB, id_fit_NB, i_NB,
                                                     cols_list[i_NB],
                                                     parameters
                                                     )

      elif (case_list[i_NB][0] == 1):
        id_temp, der_temp = derived_parameters_case1(idpar_NB, id_fit_NB, i_NB,
                                                     cols_list[i_NB],
                                                     kep_elem[i_NB - 1, :],
                                                     parameters
                                                     )

      elif (case_list[i_NB][0] == 2):
        id_temp, der_temp = derived_parameters_case2(idpar_NB, id_fit_NB, i_NB,
                                                     cols_list[i_NB],
                                                     kep_elem[i_NB - 1, :],
                                                     parameters
                                                     )

      elif (case_list[i_NB][0] == 3):
        id_temp, der_temp = derived_parameters_case3(idpar_NB, id_fit_NB, i_NB,
                                                     cols_list[i_NB],
                                                     kep_elem[i_NB - 1, :],
                                                     parameters
                                                     )

      elif (case_list[i_NB][0] == 4):
        id_temp, der_temp = derived_parameters_case4(idpar_NB, id_fit_NB, i_NB,
                                                     cols_list[i_NB],
                                                     kep_elem[i_NB - 1, :],
                                                     parameters
                                                     )

      id_derived.append(id_temp)
      der_par.append(der_temp)

  n_der = 0
  n_xder = len(der_par)
  if (n_xder > 0):
    n_der_single = []
    derived_par = []
    names_derived = []

    for ii in range(0, n_xder):
      # temp_names = np.array(id_derived[ii],dtype=str)
      temp_names = id_derived[ii]
      ntemp = np.size(temp_names)
      nls = len(np.shape(temp_names))
      n_der_single.append(ntemp)
      n_der += ntemp
      # temp_der = np.array(der_par[ii],dtype=np.float64)
      temp_der = der_par[ii]

      if (ntemp > 0):

        if (ntemp == 1 and nls == 0):
          derived_par.append(temp_der)
          names_derived.append(temp_names)
        else:
          for jj in range(0, ntemp):
            derived_par.append(temp_der[jj])
            names_derived.append(temp_names[jj])

  return np.array(names_derived, dtype=str), np.array(derived_par, dtype=np.float64).T


# ==============================================================================

def compare_par_post_angle(der_par, der_post):
  cospar = np.cos(der_par * deg2rad)
  sinpar = np.sin(der_par * deg2rad)

  par_scale = np.arctan2(sinpar, cospar) * rad2deg
  par_mod = (par_scale + 360.) % 360.

  cospos = np.cos(der_post * deg2rad)
  sinpos = np.sin(der_post * deg2rad)

  pos_scale = np.arctan2(sinpos, cospos) * rad2deg
  pos_mod = (pos_scale + 360.) % 360.

  delta_scale = np.abs(np.mean(pos_scale) - par_scale)
  delta_mod = np.abs(np.mean(pos_mod) - par_mod)

  if (np.abs(delta_scale - delta_mod) <= eps64bit):
    return par_mod, pos_mod
  else:
    if (delta_scale > delta_mod):
      return par_mod, pos_mod
    else:
      return par_scale, pos_scale
  return par_mod, pos_mod


# ==============================================================================

def adjust_derived_parameters(derived_names, derived_par, derived_post):
  n_der = np.shape(derived_par)[0]
  derived_par_adj = derived_par.copy()
  derived_post_adj = derived_post.copy()

  for i_d in range(0, n_der):
    if (derived_names[i_d][0] in ['w', 'i'] or derived_names[i_d][0:2] in ['lN', 'mA']):
      derived_par_adj[i_d], derived_post_adj[:, i_d] = compare_par_post_angle(derived_par[i_d], derived_post[:, i_d])

  return derived_par_adj, derived_post_adj


# ==============================================================================

def get_header(perc_val):

  top_header = '# %15s %20s %23s %23s %23s %23s %23s %23s %23s %23s %23s' % (' ', ' ', ' ', '+/-1sigma', 'MAD', '-1sigma', '+1sigma', '-2sigma', '+2sigma', '-3sigma', '+3sigma')
  header = '# %15s %20s %23s' % ('name', 'unit', 'parameter')
  perc_str = ' '.join(['%23s' % ('%4.2f-th' % (perc_val[i_p])) for i_p in range(len(perc_val))])
  header = '%s %s' % (header, perc_str)

  return top_header, header


# ==============================================================================

def print_parameters(top_header, header, name_parameters, unit_parameters, parameters, sigma_parameters=None, output=None):
  
  print_both(top_header, output)
  print_both(header, output)
  # n_par = parameters.shape[0]
  #n_par = len(name_parameters)
  
  if(sigma_parameters is not None):
    n_sig, n_par = np.shape(sigma_parameters)
    for i_p in range(0, n_par):
      unit = unit_parameters[i_p].strip()
      if(len(unit) == 0): unit = '-'
      sigma_line = ' '.join(['%23.16e' % (sigma_parameters[ii, i_p]) for ii in range(n_sig)])
      line = ' %15s %20s %23.16e %s' % (name_parameters[i_p],
                                        unit,
                                        parameters[i_p], sigma_line)
      print_both(line, output)
  else:
    print_both('PARAMETERS NOT AVAILABLE', output)

  return


# ==============================================================================

def print_confidence_intervals(percentiles, conf_interv=None, name_parameters=None, unit_parameters=None, output=None):
  
  if (conf_interv is not None):
    header = '# %15s %20s' % ('name', 'unit')
    perc_str = ' '.join(['%23s' %('%4.2f-th' % (percentiles[i_p])) for i_p in range(len(percentiles))])
    header = '%s %s' % (header, perc_str)
    print_both(header, output)

    n_par = len(name_parameters)
    for i_p in range(0, n_par):
      unit = unit_parameters[i_p].strip()
      if(len(unit) == 0): unit = '-'
      ci_line = ' '.join(['%23.16e' % (conf_interv[ii, i_p]) for ii in range(0, len(percentiles))])
      line = '  %15s %23s %s' % (name_parameters[i_p], unit, ci_line)
      print_both(line, output)

  else:
    print_both('EMPTY', output)

  return

# ==============================================================================

def print_hdi(conf_interv=None, name_parameters=None, unit_parameters=None, output=None):
  
  if (conf_interv is not None):
    header = '# %15s %20s %23s %23s %23s %23s %23s %23s' % ('name', 'unit',
                                                              'HDI(68.27%)lower',
                                                              'HDI(68.27%)upper',
                                                              'HDI(95.44%)lower',
                                                              'HDI(95.44%)upper',
                                                              'HDI(99.74%)lower',
                                                              'HDI(99.74%)upper'
                                                              )
    print_both(header, output)

    #n_par, n_hdi = np.shape(conf_interv)
    n_hdi, n_par = np.shape(conf_interv)
    for i_p in range(0, n_par):
      unit = unit_parameters[i_p].strip()
      if(len(unit) == 0): unit = '-'
      ci_line = ' '.join(['%23.16e' % (conf_interv[ii, i_p]) for ii in range(n_hdi)])
      line = '  %15s %23s %s' % (name_parameters[i_p], unit, ci_line)
      print_both(line, output)

  else:
    print_both('EMPTY', output)

  return

# =============================================================================

def read_fitted_file(fitted_file):
  of = open(fitted_file, 'r')
  lines = of.readlines()
  of.close()

  names, fitted_par = [], []
  for ll in lines:
    line = ll.strip()
    if (line[0] != '#'):
      lsplit = line.split()
      names.append(lsplit[0].strip())
      fitted_par.append(np.float64(lsplit[1].strip()))

  return names, np.array(fitted_par, dtype=np.float64)

# ==============================================================================
# STURGES RULE TO COMPUTE THE NUMBER OF BINS FOR HISTOGRAM
# ==============================================================================

def sturges_nbins(npost):
  '''
    Following Sturges' rule:
    nbins = log2(n) + 1
  '''
  nbins = np.ceil(np.log2(np.float64(npost))).astype(int) + 1
  
  return nbins

# ==============================================================================
# FREEDMAN-DIACONIS RULE TO COMPUTE THE NUMBER OF BINS FOR HISTOGRAM
# ==============================================================================

def freedman_diaconis_nbins(x):
  '''
    Following Freedman-Diaconis rule:
    width = 2 * IRQ / n^1/3
    nbins = [ (max-min) / width ]
  '''
  q75 = np.percentile(x, 75., interpolation='midpoint')
  q25 = np.percentile(x, 25., interpolation='midpoint')
  irq = np.abs(q75-q25)
  nx = np.float64(np.shape(x)[0])
  width = 2. * irq / np.power(nx,1./3.)
  nbins = np.ceil(np.abs(np.max(x)-np.min(x))/width).astype(int)
  
  return nbins

# ==============================================================================
# DOANE'S FORMULA TO COMPUTE THE NUMBER OF BINS FOR HISTOGRAM
# ==============================================================================

def doane_nbins(x):
  '''
   Doane's formula:
     nbins = 1 + lon2(n) + log2(1 + |g1|/s_g1)
     n: number of points/sample
     g1 = 3rd-moment-skewness
     s_g1 = sqrt( (6 x (n-2)) / ((n+1) x (n+3)) )
  '''
  nxf = np.float64(np.shape(x)[0])
  mux = np.mean(x)
  stdx = np.std(x, ddof=1)
  g1 = np.mean(((x - mux) / stdx) ** 3)  # skew
  s_g1 = np.sqrt(6. * (nxf - 2.) / ((nxf + 1.) * (nxf + 3.)))
  #print nxf
  #print mux,stdx
  #print g1, s_g1
  nbins = int(1. + np.log2(nxf) + np.log2(1. + np.abs(g1) / s_g1))

  return nbins

# ==============================================================================
# GIVEN THE POSTERIOR IT SELECT THE PROPER NBINS FOR ALL THE PARAMETER DISTRIBUTIONS
# ==============================================================================

def get_auto_bins(posterior):
  
  npost, nfit = np.shape(posterior)
  # call Freedman-Diaconis rule
  nbins_fd = [freedman_diaconis_nbins(posterior[:,ifit]) for ifit in range(nfit)]
  # doane's formula
  nbins_doa = [doane_nbins(posterior[:,ifit]) for ifit in range(nfit)]
  # and Sturges' rule
  nbins = int(np.mean([min(nbins_fd), min(nbins_doa), sturges_nbins(npost)]))
  
  #print 'fd:',nbins_fd
  #print 'do:',nbins_doa
  #print 'st:',sturges_nbins(npost)
  #print 'nb:',nbins
  
  return nbins

# ==============================================================================
# GIVEN THE POSTERIOR AND THE RULE IT COMPUTES THE NBINS
# ==============================================================================

def get_old_bins(x, rule='sq'):
  
  nx = np.shape(x)[0]
  nxf = np.float64(nx)

  # if (rule.lower() == 'sqrt'):
  if ('sq' in rule.lower()):

    nbins = int(np.sqrt(np.float64(nx)))

  # elif (rule.lower() == 'sturges'):
  elif ('stu' in rule.lower()):

    #nbins = int(1. + np.log2(np.float64(nx)))
    nbins = sturges_nbins(nxf)

  # elif (rule.lower() == 'doane'):
  elif ('doa' in rule.lower()):

    #mux = np.mean(x, axis=0)
    #stdx = np.std(x, ddof=1, axis=0)
    #g1 = np.max(np.mean(((x - mux) / stdx) ** 3, axis=0))  # skew
    ## g1 = np.min(np.mean( ((x-mux)/stdx)**3 , axis=0)) # skew
    #s_g1 = np.sqrt(6. * (nxf - 2) / ((nxf + 1.) * (nxf + 3.)))
    #nbins = int(1. + np.log2(nxf) + np.log2(1. + np.abs(g1) / s_g1))
    nbins = doane_nbins(x)

  elif('fd' in rule.lower()):
    
    nbins = freedman_diaconis_nbins(x)

  else:  # rice

    nbins = int(np.ceil(2. * np.power(nxf, 1./3.)))

  # if(nbins > 20):  nbins = 20

  return nbins

# ==============================================================================
# HIGH DENSITY INTERVALS - DOING BAYESIAN DATA ANALYSIS by Kruschke
# ==============================================================================

def calculate_hdi(x, nbins, alpha=[0.05], mode_output=False):
  
  # 68.27% (15.87th-84.13th) ==> alpha = 1. - 0.6827 = 0.3173
  # 95.44% ( 2.28th-97.72th) ==> alpha = 1. - 0.9544 = 0.0456
  # 99.74% ( 0.13th-99.87th) ==> alpha = 1. - 0.9974 = 0.0026
  
  counts, bin_edges = np.histogram(x, bins=nbins)
  #nbins = len(counts)
  
  bin_width = bin_edges[1] - bin_edges[0]
  hwidth = 0.5*bin_width
  bin_mid = [bin_edges[ibin]+hwidth for ibin in range(nbins)]
  
  # probability mass distribution = (counts / ndata)
  # probability density = probability mass distribution / bin_width
  pmd = np.asarray(counts)/np.float64(np.shape(x)[0])
  #pdd = pmd / bin_width
  spmd = np.sort(pmd)[::-1]
  hdi_ci = []
  for ialpha in alpha:
    lowestHeightIdx = np.min(np.where(np.cumsum(spmd) > (1. - ialpha)))
    lowestHeight = spmd[lowestHeightIdx]
    bin_sel = np.asarray(bin_mid)[pmd >= lowestHeight]
    hdi_l = np.min(bin_sel) - hwidth
    hdi_u = np.max(bin_sel) + hwidth
    hdi_ci.append([hdi_l, hdi_u])
  
  if(mode_output):
    max_bin = np.argmax(np.array(counts))
    #if (max_bin == 0 or max_bin == nbins):
      #ext_bin = 0
    #elif (max_bin == 1 or max_bin == nbins-1):
      #ext_bin = 1
    #else:
      #ext_bin = 2
    #sel_bin = np.logical_and(x >= bin_edges[max_bin-ext_bin],
                             #x < bin_edges[max_bin+ext_bin+1])
    sel_bin = np.logical_and(x >= bin_edges[max_bin],
                             x < bin_edges[max_bin+1])
    
    mode = np.mean(x[sel_bin])
    hdi_ci.append(mode)
  
  return hdi_ci

# ==============================================================================
# ==============================================================================

# use in emcee script sqrt(e)cos(w),sqrt(e)sin(w), while in trades use (e)cos(w),(e)sin(w)

def convert_fortran2python_strarray(array_str_in, nfit, str_len=10):
  
  # not working in python3?
  # temp = np.asarray(array_str_in, dtype=str).reshape((str_len, nfit), order='F').T
  # tested with python3...it should work...
  temp = [array_str_in[r].decode('utf-8')[c] for r in range(nfit) for c in range(str_len)]
  temp = np.asarray(temp).reshape((nfit,str_len), order='F')

  array_str_out = [''.join(temp[i, :]).strip() for i in range(0, nfit)]

  return array_str_out

# ==============================================================================

def convert_fortran_charray2python_strararray(array_str_in):
  
  ndim, str_len = np.shape(array_str_in)
  temp = np.asarray([array_str_in[r,c].decode('utf-8')
                     for r in range(ndim) for c in range(str_len)]
                   ).reshape((ndim,str_len))
  array_str_out = [''.join(temp[i, :]).strip() for i in range(0, ndim)]

  return array_str_out

# ==============================================================================

def trades_names_to_emcee(trades_names):
  nfit = np.shape(trades_names)[0]
  emcee_names = list(trades_names)
  for ifit in range(0, nfit):
    if (trades_names[ifit][:2] == 'ec'):
      emcee_names[ifit] = 'sqrt%s' % (trades_names[ifit])
      emcee_names[ifit + 1] = 'sqrt%s' % (trades_names[ifit + 1])

  return emcee_names


# ==============================================================================

def emcee_names_to_trades(emcee_names):
  nfit = np.shape(emcee_names)[0]
  trades_names = list(emcee_names)
  for ifit in range(0, nfit):
    if (emcee_names[ifit][:6] == 'sqrtec'):
      trades_names[ifit] = '%s' % (emcee_names[ifit].split('sqrt')[1])
      trades_names[ifit + 1] = '%s' % (emcee_names[ifit + 1].split('sqrt')[1])

  return trades_names


# ==============================================================================

def e_to_sqrte_boundaries(boundaries_in, names_par):
  nfit = np.shape(boundaries_in)[0]
  boundaries_out = np.array(boundaries_in, dtype=np.float64).copy()
  for ifit in range(0, nfit):
    if (names_par[ifit][:2] == 'ec'):
      ec = boundaries_in[ifit, :]
      es = boundaries_in[ifit + 1, :]
      boundaries_out[ifit, :] = np.sign(ec) * np.sqrt(np.abs(ec))
      boundaries_out[ifit + 1, :] = np.sign(es) * np.sqrt(np.abs(es))

  return boundaries_out


# ==============================================================================

def e_to_sqrte_fitting(fitting_in, names_par):
  nfit = np.shape(fitting_in)[0]
  fitting_out = np.array(fitting_in, dtype=np.float64).copy()
  for ifit in range(0, nfit):
    if (names_par[ifit][:2] == 'ec'):
      ec = fitting_in[ifit]
      es = fitting_in[ifit + 1]
      ss = ec * ec + es * es
      # sqrte = np.sqrt(np.sqrt(ss))
      sqrte = np.power(ss, 0.25)
      ww_rad = np.arctan2(es, ec)
      fitting_out[ifit] = sqrte * np.cos(ww_rad)
      fitting_out[ifit + 1] = sqrte * np.sin(ww_rad)

  return fitting_out


# ==============================================================================

def e_to_sqrte_flatchain(fitting_in, names_par):
  _, nfit = np.shape(fitting_in)
  fitting_out = np.array(fitting_in, dtype=np.float64).copy()
  for ifit in range(0, nfit):
    if (names_par[ifit][:2] == 'ec'):
      ec = fitting_in[:, ifit]
      es = fitting_in[:, ifit + 1]
      ss = ec * ec + es * es
      # sqrte = np.sqrt(np.sqrt(ss))
      sqrte = np.power(ss, 0.25)
      ww_rad = np.arctan2(es, ec)
      fitting_out[:, ifit] = sqrte * np.cos(ww_rad)
      fitting_out[:, ifit + 1] = sqrte * np.sin(ww_rad)

  return fitting_out


# ==============================================================================

def e_to_sqrte_chain(fitting_in, names_par):
  _, _, nfit = np.shape(fitting_in)
  fitting_out = np.array(fitting_in, dtype=np.float64).copy()
  for ifit in range(0, nfit):
    if (names_par[ifit][:2] == 'ec'):
      ec = fitting_in[:, :, ifit]
      es = fitting_in[:, :, ifit + 1]
      ss = ec * ec + es * es
      # sqrte = np.sqrt(np.sqrt(ss))
      sqrte = np.power(ss, np.float64(0.25))
      ww_rad = np.arctan2(es, ec)
      fitting_out[:, :, ifit] = sqrte * np.cos(ww_rad)
      fitting_out[:, :, ifit + 1] = sqrte * np.sin(ww_rad)

  return fitting_out


# ==============================================================================

def e_to_sqrte_parameters(fitting_in, names_par):
  if (len(np.shape(fitting_in)) == 1):  # fitting parameters
    fitting_out = e_to_sqrte_fitting(fitting_in, names_par)

  elif (len(np.shape(fitting_in)) == 2):  # flatchain posterior
    fitting_out = e_to_sqrte_flatchain(fitting_in, names_par)

  elif (len(np.shape(fitting_in)) == 3):  # full chain (original, nw x nr x nfit, transposed nr x nw x nfit)
    fitting_out = e_to_sqrte_chain(fitting_in, names_par)

  return fitting_out


# ==============================================================================

def sqrte_to_e_fitting(fitting_in, names_par):
  nfit = np.shape(fitting_in)[0]
  fitting_out = np.array(fitting_in, dtype=np.float64).copy()
  for ifit in range(0, nfit):
    if ('sqrtec' in names_par[ifit]):
      sec = fitting_in[ifit]
      ses = fitting_in[ifit + 1]
      ee = sec * sec + ses * ses
      ww_rad = np.arctan2(ses, sec)
      fitting_out[ifit] = ee * np.cos(ww_rad)
      fitting_out[ifit + 1] = ee * np.sin(ww_rad)

  return fitting_out


# ==============================================================================

def sqrte_to_e_flatchain(fitting_in, names_par):
  _, nfit = np.shape(fitting_in)
  fitting_out = np.array(fitting_in, dtype=np.float64).copy()
  for ifit in range(0, nfit):
    if ('sqrtec' in names_par[ifit]):
      sec = fitting_in[:, ifit]
      ses = fitting_in[:, ifit + 1]
      ee = sec * sec + ses * ses
      ww_rad = np.arctan2(ses, sec)
      fitting_out[:, ifit] = ee * np.cos(ww_rad)
      fitting_out[:, ifit + 1] = ee * np.sin(ww_rad)

  return fitting_out


# ==============================================================================

def sqrte_to_e_chain(fitting_in, names_par):
  _, _, nfit = np.shape(fitting_in)
  fitting_out = np.array(fitting_in, dtype=np.float64).copy()
  for ifit in range(0, nfit):
    if ('sqrtec' in names_par[ifit]):
      sec = fitting_in[:, :, ifit]
      ses = fitting_in[:, :, ifit + 1]
      ee = sec * sec + ses * ses
      ww_rad = np.arctan2(ses, sec)
      fitting_out[:, :, ifit] = ee * np.cos(ww_rad)
      fitting_out[:, :, ifit + 1] = ee * np.sin(ww_rad)

  return fitting_out


# ==============================================================================

def sqrte_to_e_parameters(fitting_in, names_par):
  if (len(np.shape(fitting_in)) == 1):  # fitting parameters
    fitting_out = sqrte_to_e_fitting(fitting_in, names_par)

  elif (len(np.shape(fitting_in)) == 2):  # flatchain posterior
    fitting_out = sqrte_to_e_flatchain(fitting_in, names_par)

  elif (len(np.shape(fitting_in)) == 3):  # full chain (original, nw x nr x nfit, transposed nr x nw x nfit)
    fitting_out = sqrte_to_e_chain(fitting_in, names_par)

  return fitting_out


# ==============================================================================

# ==============================================================================

# ==============================================================================

# RV semi-amplitude K in m/s
def compute_Kms(Ms_sun, Mp_jup, inc_deg, P_day, ecc):
  
  Ms_jup = Ms_sun * cst.Msjup
  sini = np.sin(inc_deg * cst.deg2rad)
  G_m_mj_s = cst.Gsi * cst.Mjup
  P_sec = P_day * cst.d2s

  P_factor = np.power((cst.dpi * G_m_mj_s / P_sec), 1.0 / 3.0)
  M_factor = (Mp_jup * sini) / np.power((Ms_jup + Mp_jup), 2.0 / 3.0)
  e_factor = 1.0 / np.sqrt(1.0 - (ecc * ecc))

  Kms = P_factor * M_factor * e_factor

  return Kms

# ==============================================================================

# ==============================================================================

# epoch or transit number for each T0 given a T0ref and a Pref
def calculate_epoch(T0, Tref, Pref):
  
  epo = np.rint((T0-Tref)/Pref).astype(int)
  
  return epo

# ==============================================================================

# 2) linear fit function with errors on y-axis.
# It returns also the errors.
# y = b_ls + m_ls * x
def lstsq_fit(x, y, yerr):
  
  A = np.vstack((np.ones_like(x), x)).T
  C = np.diag(yerr * yerr)
  cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
  b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))
  coeff_err = np.sqrt(np.diag(cov))
  
  return b_ls, m_ls, coeff_err

# ==============================================================================

def compute_lin_ephem_lstsq(T0, eT0):
  
  nT0 = np.shape(T0)[0]
  Tref0 = T0[0]
  dT = [np.abs(T0[i+1]-T0[i]) for i in range(nT0-1)]
  Pref0 = np.percentile(np.array(dT), 50., interpolation='midpoint')
  
  epo = calculate_epoch(T0, Tref0, Pref0)
  Tref, Pref, _ = lstsq_fit(epo, T0, eT0)
  epo = calculate_epoch(T0, Tref, Pref)

  return epo, Tref, Pref

# ==============================================================================

def linear_model(par, x):
  
  linmod = par[0] + x * par[1]
  
  return linmod

# ==============================================================================

def res_linear_model(par, x, y, ey=None):
  
  linmod = linear_model(par, x)
  if(ey is not None):
    wres = (y-linmod)/ey
  else:
    wres = y-linmod
  
  return wres

# ==============================================================================

def chi2r_linear_model(par, x, y, ey=None):
  
  wres = res_linear_model(par, x, y, ey)
  nfit = np.shape(par)[0]
  ndata = np.shape(x)[0]
  dof = ndata - nfit
  chi2r = np.sum(np.power(wres,2))/ np.float64(dof)
  
  return chi2r

# ==============================================================================

def compute_lin_ephem(T0, eT0=None, epoin=None, modefit='wls'):
  
  nT0 = np.shape(T0)[0]
  
  if(eT0 is None):
    errTT = np.ones((nT0))/86400.
  else:
    errTT = eT0
  
  if(epoin is None):
    # Tref0 = T0[0]
    #Tref0 = np.percentile(T0, 50., interpolation='midpoint')
    Tref0 = T0[int(0.5*nT0)]
    dT = [np.abs(T0[i+1]-T0[i]) for i in range(nT0-1)]
    Pref0 = np.percentile(np.array(dT), 50., interpolation='midpoint')
    epo = calculate_epoch(T0, Tref0, Pref0)
  else:
    epo = epoin
    Tref0, Pref0, _ = lstsq_fit(epo, T0, errTT)
  
  if(modefit in ['optimize', 'minimize']):
    # SCIPY.OPTIMIZE.MINIMIZE
    optres = sciopt.minimize(chi2r_linear_model, [Tref0, Pref0], method='nelder-mead', args=(epo, T0, errTT))
    Tref, Pref = optres.x[0], optres.x[1]
    TP_err = [0., 0.]
    epo = calculate_epoch(T0, Tref, Pref)
  
  elif(modefit in ['sklearn', 'linear_model']):
    # SKLEARN.LINEAR_MODEL
    sk_mod = sklm.LinearRegression()
    sk_mod.fit(epo.reshape((-1,1)), T0)
    Tref, Pref = np.asscalar(np.array(sk_mod.intercept_)), np.asscalar(np.array(sk_mod.coef_))
    epo = calculate_epoch(T0, Tref, Pref)

  elif(modefit == 'odr'):
    # SCIPY.ODR
    odr_lm = sodr.Model(linear_model)
    if(eT0 is not None):
      odr_data = sodr.RealData(epo, T0, sy=eT0)
    else:
      odr_data = sodr.RealData(epo, T0)
    init_coeff = [Tref0, Pref0]
    odr_mod = sodr.ODR(odr_data, odr_lm, beta0=init_coeff)
    odr_out = odr_mod.run()
    Tref, Pref = odr_out.beta[0], odr_out.beta[1]
    TP_err = odr_out.sd_beta
  
  else: # wls
    X = sm.add_constant(epo)
    if(eT0 is not None):
      wls = sm.WLS(T0, X, weights=1.0/(eT0*eT0)).fit()
    else:
      wls = sm.WLS(T0, X)
    Tref, Pref = wls.params[0], wls.params[1]
    TP_err =  wls.bse
  
  return epo, Tref, Pref, TP_err

# ==============================================================================
# ==============================================================================