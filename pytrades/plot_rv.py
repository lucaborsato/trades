#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import argparse
import os #  os: operating system
import time # time: execution time
import glob # glob: globbing file...loading multiple files as *.pippa
import sys # sys: system...boh
import numpy as np # array

import h5py

import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.figsize'] = [6, 6]
plt.rcParams["figure.facecolor"] = 'white'
plt.rcParams["savefig.facecolor"] = 'white'
plt.rcParams["figure.dpi"]  = 200
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"]   = 14

# custom modules
import constants as cst
import ancillary as anc

# ==============================================================================

filled_markers = ('o', 's', 'D', '^', '*', 'v', '<', '>', \
                  '8', 'p', 'h', 'H', 'd', 'P', 'X')

# ==============================================================================

def get_sim_file(cli):
  
  main_folder = cli.full_path
  #sim_file = os.path.join(main_folder, '%s_%s_simRV.dat' %(cli.idsim, cli.lmflag))
  sim_file = os.path.join(main_folder, '%d_%d_simRV.dat' %(cli.idsim, cli.lmflag))
  
  return sim_file

# ==============================================================================
class rv_data:
  
  def __init__(self, t, RVo, eRVo, rvs, RVs, g, eg, RVset, tscale):
    self.RVset = RVset
    self.t = t
    self.time_plot = self.t - tscale
    self.RVo = RVo
    self.rvo = RVo - g
    self.eRVo = eRVo
    self.rvs = rvs
    self.RVs = RVs
    self.gamma = g
    self.egamma = eg
    self.dRVos = RVo - RVs
    self.nRV = np.shape(self.RVo)[0]

# ==============================================================================

def get_sim_data(file_in, tscale=2440000.5):
  
  # JD 0 RVobs 1 eRVobs 2 rv_sim 3 RV_sim 4 gamma 5 e_gamma 6 RVsetID 7 RV_stat 8
  sim_in = np.genfromtxt(file_in)
  rvsetid = np.unique(sim_in[:,7].astype(int))
  nset = np.shape(rvsetid)[0]
  
  print('rvsetid: ',rvsetid)
  print('nset = ', nset)
  rv_os = [] # observations and simulations data
  for i in range(0,nset):
    print('create rv set for i = ',i)
    rvsel = sim_in[:,7].astype(int) == rvsetid[i]
    print('nrvsel = ',np.sum(rvsel))
    orv = rv_data(sim_in[rvsel,0], sim_in[rvsel,1], sim_in[rvsel,2],
                  sim_in[rvsel,3], sim_in[rvsel,4],
                  sim_in[rvsel,5], sim_in[rvsel,6],
                  rvsetid[i], tscale
                 )
    rv_os.append(orv)
  
  return rv_os

# ==============================================================================

def get_sim_model(cli):
  
  orb_file = os.path.join(cli.full_path, '%s_%s_rotorbit.dat' %(cli.idsim, cli.lmflag))
  
  data = np.genfromtxt(orb_file, usecols=(0,1,-1)) # time, lte, rv
  tmod = data[:,0]+data[:,1]-cli.tscale
  rvmod = data[:,2]
  
  return tmod, rvmod

# ==============================================================================

class rv_sample:
  
  def __init__(self, tepoch, time_rv_mod, rv_mod, tscale):
    self.tepoch = tepoch
    self.time_rv_mod = time_rv_mod + tepoch
    self.tscale = tscale
    self.time_plot = time_rv_mod + tepoch - tscale
    self.rv_mod = rv_mod

def read_samples(cli, samples_file = None):
  
  if (samples_file is not None):
    
    sh5 = h5py.File(os.path.abspath(samples_file), 'r')
    # n_samples = len(list(sh5.keys()))
    samples = []
    for igr in list(sh5.keys()): # loop on groups, one for each sample
      tepoch = sh5[igr].attrs['tepoch']
      trv = sh5['%s/time_rv_mod' %(igr)][...]
      rv = sh5['%s/rv_mod' %(igr)][...]
      samples.append(rv_sample(tepoch, trv, rv, cli.tscale))
    sh5.close()
    
  else:
    
    samples = None

  return samples

# ==============================================================================

def set_axis_default(ax, ticklabel_size=10, aspect='equal', labeldown=True):
  
  ax.ticklabel_format(useOffset=False)
  ax.tick_params(direction='inout', labelsize=ticklabel_size,
                 labelbottom=labeldown
                 )
  ax.set_aspect(aspect)
  
  return

def set_symmetric_lim(ax, val):
  
  ax.set_xlim(-val, val)
  ax.set_ylim(-val, val)
  
  return

def axtitle(ax, labtitle='', fontsize=8):
  
  ax.text(0.5, 1.02, labtitle,
          horizontalalignment='center',
          verticalalignment='bottom',
          fontsize=fontsize,
          transform=ax.transAxes
          )
  
  return 

# ==============================================================================

def plot_rv(cli, samples=None):
  
  tscale = cli.tscale
  sim_file = get_sim_file(cli)
  rv_os = get_sim_data(sim_file, tscale)
  try:
    tmod, rvmod = get_sim_model(cli)
  except:
    tmod, rvmod = None, None
  
  if(tscale == 0.0):
    xlabel = 'BJD$_\\textrm{TDB}$'
  else:
    xlabel = 'BJD$_\\textrm{TDB} - %.3f$' %(tscale)
  
  fig = plt.figure()
  fig.subplots_adjust(hspace=0.05, wspace=0.25)
  
  nrows = 3
  ncols = 1
  
  lfont=10
  tfont=8
  mso=5
  mss=4
  
  # RV plot
  ax = plt.subplot2grid((nrows, ncols), (0,0), rowspan=2)
  set_axis_default(ax, ticklabel_size=tfont, aspect='auto', labeldown=False)
  ax.set_ylabel('RV (m/s)', fontsize=lfont)
  axtitle(ax, labtitle='Radial Velocities',
          fontsize=lfont
          )
  
  ax.axhline(0., color='black', ls='-',lw=0.7)
  
  lobs_a = []
  lsim_a = []
  nset = len(rv_os)
  label_data = [r'RV dataset \#%d' %(i+1) for i in range(0,nset)]
  if(cli.labels is not None):
    nlabels = len(list(cli.labels))
    label_data[0:nlabels] = [ll for ll in list(cli.labels)]
  
  x = []
  y = []
  px = 0.05
  py = 0.05

  for i in range(0, nset):
    rv = rv_os[i]
    sort_obs = np.argsort(rv.time_plot)
    x.append(np.tile(rv.time_plot[sort_obs], 3)) # for obs limits
    y.append(np.concatenate((rv.rvo[sort_obs]-rv.eRVo[sort_obs],
                             rv.rvo[sort_obs]+rv.eRVo[sort_obs],
                             rv.rvs[sort_obs]
                            ))) # for obs limits
    # data
    lobs = ax.errorbar(rv.time_plot[sort_obs], rv.rvo[sort_obs],
                       yerr=rv.eRVo[sort_obs],
                       color='black', ecolor='gray',
                       #fmt=filled_markers[i], ms=mso, mec='None',
                       fmt=filled_markers[i], ms=mso, mfc='white', mec='black',
                       ls='',
                       elinewidth=0.6,
                       capsize=0,
                       zorder=7,
                       label=label_data[i]
                       )
    lobs_a.append(lobs)
    # trades
    lsim, = ax.plot(rv.time_plot[sort_obs], rv.rvs[sort_obs],
                    #marker=filled_markers[i], ms=mss, mfc='None',
                    color='C0',
                    marker=filled_markers[i], ms=mss,
                    mec='lightgray',  mew=0.5, # mec='None',
                    ls='',
                    zorder=8,
                    label='simulations'
                    )
    lsim_a.append(lsim)

  # prepare limits
  x = np.concatenate(x)
  dx = np.max(x) - np.min(x)
  minx = np.min(x) - px*dx
  maxx = np.max(x) + px*dx
  y = np.concatenate(y)
  dy = np.max(y) - np.min(y)
  miny = np.min(y) - py*dy
  maxy = np.max(y) + py*dy
  # xlims = ax.get_xlim()

  
  if(tmod is not None):
    # plot model and samples
    sort_mod = np.argsort(tmod)
    lmod, = ax.plot(tmod[sort_mod], rvmod[sort_mod], 
                    color='black', marker='None', ls='--', lw=0.8,
                    zorder=6,
                    label='RV model'
                    )
  if(samples is not None):
    lsmp_a = []
    xx     = x
    yy     = y
    # nsmp = len(samples)
    for smp in samples:
      sort_smp = np.argsort(smp.time_plot)
      xx       = np.concatenate((xx,smp.time_plot[sort_smp]))
      yy       = np.concatenate((yy,smp.rv_mod[sort_smp]))
      lsmp,    = ax.plot(smp.time_plot[sort_smp], smp.rv_mod[sort_smp],
                      color='gray', marker='None', ls='-', lw=0.5,
                      alpha=0.33,
                      zorder=5,
                      label='RV samples'
                      )
      lsmp_a.append(lsmp)
      
    if(cli.limits == 'sam'):
      dx = np.max(xx) - np.min(xx)
      minx = np.min(xx) - px*dx
      maxx = np.max(xx) + px*dx
      dy = np.max(yy) - np.min(yy)
      miny = np.min(yy) - py*dy
      maxy = np.max(yy) + py*dy
    else: # if obs & samples not None adjust y limits
      yyy = yy[np.logical_and(xx>=minx, xx<=maxx)]
      dy  = np.max(yyy) - np.min(yyy)
      miny = np.min(yyy) - py*dy
      maxy = np.max(yyy) + py*dy

  # ax.set_xlim(xlims)
  ax.set_xlim(minx, maxx)
  ax.set_ylim(miny, maxy)
  
  if(samples is None and tmod is None):
    lhand = lobs_a + [lsim_a[0]]
  elif(samples is None and tmod is not None):
    lhand = lobs_a + [lsim_a[0], lmod]
  elif(samples is not None and tmod is None):
    lhand = lobs_a + [lsim_a[0], lsmp_a[0]]
  else:
    lhand = lobs_a + [lsim_a[0], lmod, lsmp_a[0]]
  
  # legend in a column right of the boxes
  #ax.legend(handles=lhand, bbox_to_anchor=(1., 0.5), fontsize=tfont-2)
  # legend automatic
  nlegcol = int( (len(label_data) + 3) / 2) + 1 # nRVset obs + 1 sim + 1 model + 1 samples
  ax.legend(handles=lhand, loc='best', ncol=nlegcol, fontsize=lfont)
  
  # residuals plot
  ax = plt.subplot2grid((nrows, ncols), (2,0), rowspan=1)
  set_axis_default(ax, ticklabel_size=tfont, aspect='auto')
  ax.set_ylabel('res (m/s)', fontsize=lfont)
  ax.set_xlabel(xlabel, fontsize=lfont)
  
  ax.axhline(0., color='black', ls='-',lw=0.7)
  
  for i in range(0, nset):
    rv = rv_os[i]
    sort_obs = np.argsort(rv.time_plot)
    ax.errorbar(rv.time_plot[sort_obs], rv.dRVos[sort_obs],
                yerr=rv.eRVo[sort_obs],
                #fmt=filled_markers[i], ms=4, mec='None',
                color='C0',
                fmt=filled_markers[i], ms=mss, mec='lightgray', mew=0.5,
                ls='',
                ecolor='gray',
                elinewidth=0.6,
                capsize=0,
                zorder=5,
                )
  
  ax.set_xlim(minx, maxx)

  print('full_path', cli.full_path)
  # folder_out = os.path.join(os.path.dirname(cli.full_path), 'plots')
  folder_out = os.path.join(cli.full_path, 'plots')
  print('folder_out', folder_out)
  if(not os.path.isdir(folder_out)): os.makedirs(folder_out)
  fname = os.path.basename(sim_file)
  plt_file = os.path.join(folder_out, os.path.splitext(fname)[0])
  fig.savefig('%s.png' %(plt_file), bbox_inches='tight', dpi=200)
  print('Saved plot into:')
  print('%s' %('%s.png' %(plt_file)))
  fig.savefig('%s.pdf' %(plt_file), bbox_inches='tight', dpi=72)
  print('%s' %('%s.pdf' %(plt_file)))
  plt.close(fig)
  
  return

# ==============================================================================

def get_args():
  
  parser = argparse.ArgumentParser()
  
  parser.add_argument('-p', '--p',
                      '-path', '--path',
                      '-path-folder', '--path-folder', 
                      action='store', dest='full_path',
                      required=True,
                      help='Folder path'
                      )
  parser.add_argument('-s', '--s',
                      '-sim-id', '--sim-id', 
                      action='store',
                      dest='idsim',
                      required=True,
                      help='Simulation ID number'
                      )
  parser.add_argument('-lm', '--lm', 
                      '-lm-flag', '--lm-flag',
                      action='store',
                      dest='lmflag',
                      default='0',
                      help='LM flag: 0 = not used LM (default), 1 = used LM'
                      )
  parser.add_argument('-ts', '--ts'
                      '-t-scale', '--t-scale',
                      action='store',
                      dest='tscale',
                      default=None,
                      help='Value to be subtract to the x-axis (time) in the plots. Default is None that means 2440000.5 will be subtract.'
                      )
  parser.add_argument('-labels', '--labels',
                      action='store',
                      dest='labels',
                      default='None',
                      help='Sequence of labels for RV data set, avoid spaces within 1 label, e.g.: "HARPS FIES". Put labels sorted as in the obsRV.dat file: "HARPS FIES" means 1 -> HARPS, 2 -> FIES. Default is None, meaning that it will plot "RV dataset #x" where x is an increasing number. '
                      )
  parser.add_argument('-samples-file', '--samples-file',
                      action='store',
                      dest='samples_file',
                      default='None',
                      help='HDF5 file with T0 and RV from emcee samples to overplot on O-Cs.'
                      )
  parser.add_argument('-limits', '--limits',
                      action='store',
                      dest='limits',
                      default='obs',
                      help='Axis limits based on observations (obs) or synthetic samples (sam). Default obs.'
                      )

  cli = parser.parse_args()
  
  cli.full_path = os.path.abspath(cli.full_path)
  
  cli.idsim = int(cli.idsim)
  
  cli.lmflag = int(cli.lmflag)
  
  if(str(cli.tscale).lower() != 'none'):
    try:
      cli.tscale = np.float64(cli.tscale)
    except:
      cli.tscale = 2440000.5
  else:
    cli.tscale = 2440000.5
  
  if(str(cli.labels).lower() != 'none'):
    try:
      cli.labels = cli.labels.split()
    except:
      cli.labels = None
  else:
    cli.labels = None
  
  if(str(cli.samples_file).lower() == 'none'):
    cli.samples_file = None
  else:
    if(os.path.isfile(os.path.abspath(cli.samples_file))):
      cli.samples_file = os.path.abspath(cli.samples_file)
    else:
      cli.samples_file = None
      
  if(cli.limits.lower()[:3] == 'sam'):
    cli.limits = 'sam'
  else:
    cli.limits = 'obs'

  return cli

# ==============================================================================

def main():

  cli = get_args()
  samples = read_samples(cli, samples_file = cli.samples_file)
  plot_rv(cli, samples)

  return

# ==============================================================================
# ==============================================================================

if __name__ == '__main__':
  main()
# ==============================================================================
