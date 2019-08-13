#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
#import h5py
import sys
import constants as cst
import ancillary as anc
import glob

import matplotlib.cm as cm
from matplotlib import use as mpluse
mpluse("Agg")
#mpluse("Qt4Agg")
import matplotlib.pyplot as plt
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# ==============================================================================

# read arg.in file and retrieve t_epoch
def read_t_epoch(full_path):
  '''
    input:
    full_path [folder path]
    output:
    t_epoch [day]
  '''
  arg_file = os.path.join(full_path, 'arg.in')
  of = open(arg_file, 'r')
  lines = of.readlines()
  of.close()
  nlines = len(lines)
  epoch, start = False, False
  for i_l in range(0, nlines):
    line_split = lines[i_l].strip().split('=')
    if('tepoch' in line_split[0]):
      t_epoch = np.float64(line_split[1])
      epoch = True
      break
    elif('tstart' in line_split[0]):
      t_start = np.float64(line_split[1])
      start = True
  if(not epoch):
    if(start):
      t_epoch = t_start
    else:
      t_epoch = 1.0

  return t_epoch

# -------------------------------------

# star and initialElements file
def read_kepelem(full_path, id_sim, lm_flag, sub_folder = ''):
  '''
    input:
    full_path [folder path]
    id_sim [number of the simulation]
    lm_flag [0/1] LM flag
    sub_folder (optional): if the initialElements,
     rotorbit files are in a sub_folder
    output:
    star [[M, sigma_M],[R, sigma_R]]
    planets [n_pl x n_kepelem]
  '''
  
  # look for bodies.lst to retrieve the star name and values
  bd_lst = os.path.join(full_path, 'bodies.lst')
  of = open(bd_lst, 'r')
  line = of.readline()
  of.close()
  # star
  star_name = line.strip().split()[0]
  star_file = os.path.join(full_path, star_name)
  xx = np.genfromtxt(star_file)
  if(len(np.shape(xx)) == 2):
    star = xx
  else:
    star = np.zeros((2,2))
    star[:,0] = xx
  # initialElements of the planets
  # 
  sf_path = os.path.abspath(os.path.join(full_path, sub_folder))
  kep_file = os.path.join(sf_path, '%s_%s_initialElements.dat' %(str(id_sim), str(lm_flag)))
  # 0      1      2   3    4 5     6      7       8
  # M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg
  planets = np.genfromtxt(kep_file)
  
  return star, planets

# -------------------------------------

# read rotorbit file to retrieve RV model from la column
def read_RVmodel(full_path, id_sim, lm_flag, sub_folder = ''):
  '''
    input:
    full_path [folder path]
    id_sim [number of the simulation]
    lm_flag [0/1] LM flag
    sub_folder (optional): if the initialElements, rotorbit files are in a sub_folder
    output:
    t_model [day]
    rv_model [m/s]
  '''

  sf_path = os.path.abspath(os.path.join(full_path, sub_folder))
  rot_file = os.path.join(sf_path, '%s_%s_rotorbit.dat' %(str(id_sim), str(lm_flag)))
  # 0       1     2      3 -- -1
  # Time_JD LTE_d X_1_AU Y_1_AU Z_1_AU VX_1_AU VY_1_AU VZ_1_AU -- RV_ms^-1
  all_orbit = np.genfromtxt(rot_file, usecols=(0,1,-1))
  
  t_model = all_orbit[:,0] + all_orbit[:,1]
  rv_model = all_orbit[:,2]
  ids = np.argsort(t_model)

  return t_model[ids], rv_model[ids]

# -------------------------------------

def read_simRV(full_path, id_sim, lm_flag, sub_folder = ''):
  '''
    input:
    full_path [folder path]
    id_sim [number of the simulation]
    lm_flag [0/1] LM flag
    sub_folder (optional): if the initialElements, rotorbit files are in a sub_folder
    output:
    sim_rv [full output of simRV.dat file]
  '''
  
  sf_path = os.path.abspath(os.path.join(full_path, sub_folder))
  rv_file = os.path.join(sf_path, '%s_%s_simRV.dat' %(str(id_sim), str(lm_flag)))
  # 0  1     2      3      4      5     6       7       8
  # JD RVobs eRVobs rv_sim RV_sim gamma e_gamma RVsetID RV_stat
  sim_rv = np.genfromtxt(rv_file)
  
  n_rv = np.shape(sim_rv)[0]
  n_rvset = 1
  for irv in range(0,n_rv-1):
    if(sim_rv[irv,7] != sim_rv[irv+1,7] ):
      n_rvset += 1
  
  
  return sim_rv, n_rvset

# ==============================================================================

# given the t_ref from cli, t_epoch and numbe of planets it determines the proper t_ref
def set_t_ref(t_ref_str, t_epoch, n_pl):
  '''
    input:
    t_ref_str [str]
    t_epoch [day]
    n_pl number of planets
    output:
    t_ref [day]
  '''
  
  t_ref_list = t_ref_str.split()
  n_ref = len(t_ref_list)
  if(n_ref >= n_pl):
    try:
      t_ref = [np.float64(t_ref_list[i]) for i in range(0, n_ref)]
    except:
      t_ref = [t_epoch for i in range(0, n_pl)]
  else: #(n_ref < n_pl)
    dn = n_pl - n_ref
    try:
      t_ref = [np.float64(t_ref_list[i]) for i in range(0, n_ref)] + \
              [t_epoch for i in range(0, dn)]
    except:
      t_ref = [t_epoch for i in range(0, n_pl)]
      
  return t_ref

# ==============================================================================

# comput ther mA for each RV time from the epoch of reference and mA at epoch
def rv_mean_anomaly(t_epoch, mA_epoch, P_day, t_rv):
  '''
    input:
    t_epoch [BJD ... or JD ... day]
    mA_epoch [deg]
    P_day [day]
    t_rv [same as t_epoch]
    output:
    mA_rv [deg]
  '''
  
  mA_rv = mA_epoch - 360. * (t_epoch - t_rv) / P_day
  
  return mA_rv

# -------------------------------------

# compute the eccentri anomaly given mean anomaly and eccentricity
def mean_to_eccentric_anomaly(mA, ecc):
  '''
    input:
    mA [deg]
    ecc
    output:
    eA [deg]
  '''
  tol = 1.e-9
  eps = np.finfo(np.float64(1.0)).eps
  eA=0.
  mA_rad = (mA * cst.deg2rad)%cst.dpi
  E=mA_rad
  if (ecc > eps):
    if(ecc > 0.6):
      E = cst.pi
    imax = 1000
    for i in range(imax):
      fE = E - ecc * np.sin(E) - mA_rad
      dfE = 1. - ecc * np.cos(E)
      eA = E - (fE / dfE)
      if (np.abs(E - eA) <= tol):
        break
      E = eA
    #end do loopo
    eA=eA*cst.rad2deg
  else:
    eA = mA
  
  return eA

# -------------------------------------

# eccentric anomaly and eccentricity to true anomaly
def eA_ecc_to_true_anomaly(eA, ecc):
  '''
    input:
    eA [deg] eccentric anomaly
    ecc
    output:
    f_deg [deg] true anomaly
  '''
  
  eA_rad = eA * cst.deg2rad
  AA = np.sqrt( (1. + ecc) / (1. - ecc) )
  BB = np.tan(eA_rad * 0.5)
  tanfh = AA * BB
  f_rad = np.arctan(tanfh) * 2.
  f_deg = f_rad * cst.rad2deg
  
  return f_deg

# -------------------------------------

def mean_anomaly_to_true_anomaly(t_rv, t_epoch, mA_epoch, P_day, ecc):
  '''
    input:
    t_epoch [day]
    mA_epoch [deg]
    P_day [day]
    t_rv [day]
    ecc
    output:
    f_deg [deg]
  '''
  
  mA_rv = rv_mean_anomaly(t_epoch, mA_epoch, P_day, t_rv)
  eA_rv = np.array([mean_to_eccentric_anomaly(mA_rv[irv], ecc) for irv in range(np.shape(t_rv)[0])])
  f_deg = eA_ecc_to_true_anomaly(eA_rv, ecc)
  
  return f_deg


# -------------------------------------

# compute the phase to plot rv for each data point at t_rv and wrt t_ref
def phase_rv(t_rv, t_ref, P_day):
  '''
    input:
    t_rv [day]
    t_ref [day]
    P_day [day]
    output:
    phi_rv [-1,1] phi_rv == 0 <-> t == t_ref
  '''
  
  phi_rv = ((t_rv - t_ref)%P_day) / P_day
  
  return phi_rv

# -------------------------------------

def compute_Kms(Ms_sun,Mp_sun,inc_deg,P_day,ecc):
  '''
    input:
    Ms_sun [star mass in Msun]
    Mp_sun [planet mass in Msun]
    inc_deg [deg]
    P_day [day]
    ecc
    output:
    K [m/s]
  '''
  
  Ms_jup = Ms_sun * cst.Msjup
  Mp_jup = Mp_sun * cst.Msjup
  sini = np.sin(inc_deg*cst.deg2rad)
  G_m_mj_s = cst.Gsi * cst.Mjup
  P_sec = P_day * cst.d2s
  
  P_factor = np.power((cst.dpi*G_m_mj_s/P_sec), 1./3.)
  M_factor = (Mp_jup*sini) * np.power((Ms_jup + Mp_jup), -2./3.)
  e_factor = 1./ np.sqrt(1. - (ecc*ecc))
  
  Kms = P_factor * M_factor * e_factor
  
  return Kms

# -------------------------------------

# calculate a Keplerian for given planet
# Keplerian_1 needs: RV semi-amplitude K[m/s], argument of pericenter w[deg], true anomaly f[deg], eccentricity e
def Keplerian_1(K, w, f, e):
  '''
    input:
    K [m/s]
    w [deg]
    f [deg]
    e
    output:
    RV [m/s]
  '''

  wf = (w + np.asarray(f, dtype=np.float64)) * cst.deg2rad
  rvKep = K * ( np.cos(wf) + e * np.cos(w * cst.deg2rad) )

  return rvKep

# -------------------------------------

def keplerian_rv(Ms_sun,Mp_sun,inc_deg,P_day,ecc,w_deg,fA_deg):
  
  Kms = compute_Kms(Ms_sun, Mp_sun, inc_deg, P_day, ecc)
  rv_kep = Keplerian_1(Kms, w_deg,fA_deg, ecc)
  
  return rv_kep

# ==============================================================================

def get_labels(labels_in, n_rvset):
  
  labels_list = labels_in.split()
  n_lab = len(labels_list)
  
  if(n_lab == 1):
    if('none' in labels_in.lower()):
      labels_out = ['RVset#02d' %(iset) for iset in range(0, n_rvset)]
    else:
      labels_out = [labels_in.strip()]
  
  elif(n_lab == n_rvset):
    labels_out = labels_list
  
  elif(n_lab > 1):
    if(n_lab < n_rvset):
      labels_out = labels_list + ['RVset#02d' %(iset) for iset in range(n_lab, n_rvset)]
    else:
      labels_out = labels_list[0:n_rvset]
  
  else:
    labels_out = ['RVset#02d' %(iset) for iset in range(0, n_rvset)]
  
  
  return labels_out

# -------------------------------------

def plot_model_data(ax, t_epoch, t_model, rv_model, sim_rv, n_rvset, labels, colors, markers):
  
  t_rv = sim_rv[:,0] - t_epoch
  rv_obs = sim_rv[:,1] - sim_rv[:,5]
  rv_sim = sim_rv[:,3]
  
  xlim1, xlim2 = anc.compute_limits(t_rv, 0.05)
  miny = min(np.min(rv_obs-sim_rv[:,2]), np.min(rv_sim), np.min(rv_model))
  maxy = max(np.max(rv_obs+sim_rv[:,2]), np.min(rv_sim), np.min(rv_model))
  ylim1, ylim2 = anc.compute_limits(np.array([miny, maxy]), 0.05)
  
  ax.axhline(0., color='gray', ls='-')
  ax.plot(t_model-t_epoch, rv_model, marker='None', ls=':', lw=0.5, color='darkgray', label='model $N$-body')
  for iset in range(0, n_rvset):
    sel = sim_rv[:,7].astype(int) == iset+1
    ax.errorbar(t_rv[sel], rv_obs[sel], yerr=sim_rv[sel,2], fmt=markers[iset], mfc=colors[iset], mec='black', ecolor=colors[iset], capsize=0, label=labels[iset])
    ax.plot(t_rv[sel], rv_sim[sel], marker=markers[iset], mfc='None', mec=colors[iset], mew=1., ms=8, ls='', label='%s (TRADES)' %(labels[iset]))
  
  ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,\
           ncol=2*(n_rvset)+1, mode="expand", \
           borderaxespad=0., numpoints=1,\
           fontsize=9)
  ax.set_xlim([xlim1, xlim2])
  ax.set_xticklabels([])
  ax.set_ylim([ylim1, ylim2])
  ax.set_ylabel('RV [m/s]')
  plt.draw()
  
  return ax

# -------------------------------------

def plot_model_res(ax, t_epoch, t_model, rv_model, sim_rv, n_rvset, labels, colors, markers):
  
  t_rv = sim_rv[:,0] - t_epoch
  rv_obs = sim_rv[:,1] - sim_rv[:,5]
  rv_sim = sim_rv[:,3]
  res_rv = rv_obs - rv_sim
  
  xlim1, xlim2 = anc.compute_limits(t_rv, 0.05)
  miny = np.min(res_rv-sim_rv[:,2])
  maxy = np.max(res_rv+sim_rv[:,2])
  ylim1, ylim2 = anc.compute_limits(np.array([miny, maxy]), 0.05)
  
  ax.axhline(0., color='gray', ls='-')
  for iset in range(0, n_rvset):
    sel = sim_rv[:,7].astype(int) == iset+1
    ax.errorbar(t_rv[sel], res_rv[sel], yerr=sim_rv[sel,2], fmt=markers[iset], mfc=colors[iset], mec='black', ecolor=colors[iset], capsize=0, label=labels[iset])

  ax.set_xlabel('BJD-%.3f' %(t_epoch))
  ax.set_xlim([xlim1, xlim2])
  ax.set_ylim([ylim1, ylim2])
  ax.set_ylabel('res [m/s]')
  plt.draw()
  
  return ax

# -------------------------------------

def plot_phase_rv(ax, phi_model, rv_model, phi_rv, rv_opl, rv_spl, sim_rv, n_rvset, labels, colors, markers):
  
  #xlim1, xlim2 = anc.compute_limits(phi_rv, 0.05)
  xlim1, xlim2 = 0., 1.
  miny = np.min(rv_opl-sim_rv[:,2])
  maxy = np.max(rv_opl+sim_rv[:,2])
  ylim1, ylim2 = anc.compute_limits(np.array([miny, maxy]), 0.05)
  
  ax.axhline(0., color='gray', ls='-')
  ax.plot(phi_model, rv_model, marker='None', ls='--', lw=0.5, color='darkgray', label='model Keplerian')
  for iset in range(0, n_rvset):
    sel = sim_rv[:,7].astype(int) == iset+1
    ax.errorbar(phi_rv[sel], rv_opl[sel], yerr=sim_rv[sel,2], fmt=markers[iset], mfc=colors[iset], mec='black', ecolor=colors[iset], capsize=0)
    ax.plot(phi_rv[sel], rv_spl[sel], marker=markers[iset], mfc='None', mec=colors[iset], mew=1., ms=8, ls='') # skip simulated
  
  ax.legend(loc=0, fontsize=9, numpoints=1)
  ax.set_xlim([xlim1, xlim2])
  ax.set_xticklabels([])
  ax.set_ylim([ylim1, ylim2])
  ax.set_ylabel('RV [m/s]')
  plt.draw()
  
  return ax

# -------------------------------------

def plot_phase_res(ax, phi_rv, rv_opl, rv_spl, res_rv, sim_rv, rv_kep, n_rvset, labels, colors, markers, pl_id):
  
  res_okep = rv_opl - rv_kep
  res_skep = rv_spl - rv_kep
  
  #xlim1, xlim2 = anc.compute_limits(phi_rv, 0.05)
  xlim1, xlim2 = 0., 1.
  #miny = np.min(res_rv-sim_rv[:,2])
  #maxy = np.max(res_rv+sim_rv[:,2])
  miny = min(np.min(res_okep-sim_rv[:,2]), np.min(res_skep-sim_rv[:,2]))
  maxy = max(np.max(res_okep+sim_rv[:,2]), np.max(res_skep+sim_rv[:,2]))
  ylim1, ylim2 = anc.compute_limits(np.array([miny, maxy]), 0.05)
  
  ax.axhline(0., color='gray', ls='-')
  for iset in range(0, n_rvset):
    sel = sim_rv[:,7].astype(int) == iset+1
    ax.errorbar(phi_rv[sel], res_okep[sel], yerr=sim_rv[sel,2], fmt=markers[iset], mfc=colors[iset], mec='black', ecolor=colors[iset], capsize=0, label=labels[iset])
    ax.errorbar(phi_rv[sel], res_skep[sel], yerr=sim_rv[sel,2], fmt=markers[iset], ms= 8, mfc='None', mec=colors[iset], mew=1., ecolor=colors[iset], capsize=0, label=labels[iset])
    #ax.errorbar(phi_rv[sel], res_rv[sel], yerr=sim_rv[sel,2], fmt=markers[iset], mfc=colors[iset], mec='black', ecolor=colors[iset], capsize=0, label=labels[iset])

  ax.set_xlabel('$\phi$ [phase planet %s]' %(pl_id))
  ax.set_xlim([xlim1, xlim2])
  ax.set_ylim([ylim1, ylim2])
  ax.set_ylabel('res [m/s]')
  plt.draw()
  
  return ax

# ==============================================================================

# read command line arguments
def get_args():
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-p', '--path', '--full-path', action='store', dest='full_path', required=True, 
                      help='Folder path')
  parser.add_argument('-s', '--id', '--id-sim', '--sim', action='store', dest='id_sim', required=True, 
                      type=str, help='Simulation ID number')
  parser.add_argument('-lm', '--lm-flag', action='store', dest='lm_flag', default='0', type=str,
                      help='LM flag: 0 = not used LM (default), 1 = used LM')
  parser.add_argument('-sf', '--sub-folder', action='store', dest='sub_folder', default='', 
                      help='Sub-Folder, without path')
  parser.add_argument('-t', '--t-ref', action='store', dest='t_ref', default='t_epoch', type=str,
                      help='Value/s of the time/s needed to phase-fold each RV curve. Number of values equals the planets, otherwise the epoch time will be used for all the planets. Eg. -t "2456813. 2456816.". Default is t_epoch.')
  parser.add_argument('-la', '--label', '--labels', action='store', dest='labels', default='None', type=str,
                      help='Labels for each different RV set. Provide labels space-separated in quotes, eg. "RVset1 RVset2 ...". Number of labels should match number of RV dataset.')
  
  cli = parser.parse_args()
  
  cli.full_path = os.path.abspath(cli.full_path)
  
  return cli
# ==============================================================================

def main():
  
  # cli arguments
  cli = get_args()
  
  # get the basic parameters and data
  t_epoch = read_t_epoch(cli.full_path)
  star, planets = read_kepelem(cli.full_path, cli.id_sim, cli.lm_flag, sub_folder = cli.sub_folder)
  n_pl = np.shape(planets)[0]
  n_bd = n_pl + 1
  t_model, rv_model = read_RVmodel(cli.full_path, cli.id_sim, cli.lm_flag, sub_folder = cli.sub_folder)
  sim_rv, n_rvset = read_simRV(cli.full_path, cli.id_sim, cli.lm_flag, sub_folder = cli.sub_folder)
  n_rv = np.shape(sim_rv)[0]
  print()
  print(' t_epoch = ',t_epoch)
  print(' number of planets = ',n_pl, ' for ', n_bd, ' bodies in the system')
  print(' Star mass = ', star[0,0], ' radius = ', star[1,0])
  t_ref = set_t_ref(cli.t_ref, t_epoch, n_pl)
  print(' t_ref = ',t_ref)
  print(' labels input = ', cli.labels)
  labels = get_labels(cli.labels, n_rvset)
  print(' labels output= ', labels)
  
  #viridis_cmap = cm.viridis.colors
  #sel_col = np.random.choice(len(viridis_cmap),n_rvset)
  #colors = [viridis_cmap[isel] for isel in sel_col] # from cmap
  colors = 'blue green red orange cyan'.split()
  markers = 'o d ^ s 8 p P *'.split()
  
  rv_obs = sim_rv[:,1] - sim_rv[:,5]
  rv_sim = sim_rv[:,3]
  rv_opl, rv_spl, rv_kep = np.zeros((n_rv, n_pl)), np.zeros((n_rv, n_pl)),\
                                   np.zeros((n_rv, n_pl))

  #n_phmod = 1000
  #fA_model = np.linspace(0., 360., num=n_phmod, endpoint=True)
  #phi_model = np.linspace(0., 1., num=n_phmod, endpoint=True)
  t_smod = t_model[np.logical_and(t_model >= np.min(sim_rv[:,0]), t_model <= np.max(sim_rv[:,0]))]
  n_phmod = np.shape(t_smod)[0]
  fA_model = np.zeros((n_phmod, n_pl))
  phi_model = np.zeros((n_phmod, n_pl))
  rv_mkep = np.zeros((n_phmod, n_pl))
  
  phi_rv, fA_rv = np.zeros((n_rv, n_pl)), np.zeros((n_rv, n_pl))
  
  for ipl in range(0, n_pl):
    phi_model[:,ipl] = phase_rv(t_smod, t_ref[ipl], planets[ipl,2])
    fA_model[:,ipl] = mean_anomaly_to_true_anomaly(t_smod, t_epoch, planets[ipl,6], planets[ipl,2], planets[ipl,4])
    rv_mkep[:,ipl] = -keplerian_rv(star[0,0], planets[ipl,0], planets[ipl,7], \
                    planets[ipl,2], planets[ipl,4], planets[ipl,5]+180., \
                    fA_model[:,ipl])
    #rv_mkep[:,ipl] = -keplerian_rv(star[0,0], planets[ipl,0], planets[ipl,7], \
                    #planets[ipl,2], planets[ipl,4], planets[ipl,5], \
                    #fA_model[:,ipl])
    
    phi_rv[:,ipl] = phase_rv(sim_rv[:,0], t_ref[ipl], planets[ipl,2])
    fA_rv[:,ipl] = mean_anomaly_to_true_anomaly(sim_rv[:,0], t_epoch, planets[ipl,6], planets[ipl,2], planets[ipl,4])
    rv_kep[:,ipl] = -keplerian_rv(star[0,0], planets[ipl,0], planets[ipl,7], \
                    planets[ipl,2], planets[ipl,4], planets[ipl,5]+180., \
                    fA_rv[:,ipl])
    #rv_kep[:,ipl] = -keplerian_rv(star[0,0], planets[ipl,0], planets[ipl,7], \
                    #planets[ipl,2], planets[ipl,4], planets[ipl,5], \
                    #fA_rv[:,ipl])
    
    sel_pl = np.ones((n_pl)).astype(bool)
    sel_pl[ipl] = False
    rv_opl[:,ipl] = rv_obs - np.sum(rv_kep[:,sel_pl], axis=1)
    rv_spl[:,ipl] = rv_sim - np.sum(rv_kep[:,sel_pl], axis=1)
    #rv_opl[:,ipl] = rv_obs
    #rv_spl[:,ipl] = rv_sim
    #for ipl2 in range(0, n_pl):
      #if(ipl != ipl2):
        #rv_opl[:,ipl] -= rv_kep[:,ipl2]
        #rv_spl[:,ipl] -= rv_kep[:,ipl2]
  
  res_rv = rv_opl - rv_spl
  
  idx = np.argsort(phi_rv[:,0])
  #fig = plt.figure(figsize=(12,12))
  fig = plt.figure(figsize=(6,6))
  
  lmin, bmin = 0.07, 0.07
  rmax, umax = 1.-0.5*lmin, 1.-0.5*bmin
  print('lmin, rmax, bmin, umax:', lmin, rmax, bmin, umax)
  dhbig = 0.05 # large vertical distance, small distance between plot e res plot
  hfull = umax - bmin
  hrvres = (hfull / n_bd)- dhbig
  hrv, hres = hrvres*0.75, hrvres*0.25
  print('dhbig:', dhbig)
  print('hfull, hrvres, hrv, hres:', hfull, hrvres, hrv, hres)
  ww = rmax - lmin
  print('ww:', ww)
  
  pl_ids = 'b c d e f g h i j k l m n o p q r s t u v w x y z'.split()
  
  # plot model
  brv = umax-hrv
  ax_rvmod = fig.add_axes([lmin, brv, ww, hrv])
  ax_rvmod = plot_model_data(ax_rvmod, t_epoch, t_model, rv_model, sim_rv, n_rvset, labels, colors, markers)
  # plot res model
  bres = brv - hres
  ax_resmod = fig.add_axes([lmin, bres, ww, hres])
  ax_resmod = plot_model_res(ax_resmod, t_epoch, t_model, rv_model, sim_rv, n_rvset, labels, colors, markers)
  
  for ipl in range(0, n_pl):
    idph = np.argsort(phi_model[:,ipl])
    # plot phased-rv
    print('planet ', ipl+1)
    brv = bres - dhbig - hrv
    print('left, bottom, widht, height:', lmin, brv, ww, hrv)
    ax_rv = fig.add_axes([lmin, brv, ww, hrv])
    ax_rv = plot_phase_rv(ax_rv, phi_model[idph,ipl], rv_mkep[idph,ipl], phi_rv[:,ipl], rv_opl[:,ipl], rv_spl[:,ipl], sim_rv, n_rvset, labels, colors, markers)
    # plot phased-rv
    bres = brv - hres
    print('left, bottom, widht, height:', lmin, bres, ww, hres)
    ax_res = fig.add_axes([lmin, bres, ww, hres])
    ax_res = plot_phase_res(ax_res, phi_rv[:,ipl], rv_opl[:,ipl], rv_spl[:,ipl], res_rv[:,ipl], sim_rv, rv_kep[:,ipl], n_rvset, labels, colors, markers, pl_ids[ipl])
    print()
  
  out_folder = os.path.join(os.path.abspath(os.path.join(cli.full_path, cli.sub_folder)), 'plots')
  if(not os.path.isdir(out_folder)):
    os.makedirs(out_folder)
  out_file = os.path.join(out_folder, '%s_%s_plot_rv' %(str(cli.id_sim), str(cli.lm_flag)))
  
  fig.savefig('%s.pdf' %(out_file), bbox_inches='tight', dpi=150)
  #fig.savefig('%s.pdf' %(out_file), dpi=150)
  print(' Saved file %s.pdf' %(out_file))
  
  return

# ==============================================================================

if __name__ == "__main__":
  main()
