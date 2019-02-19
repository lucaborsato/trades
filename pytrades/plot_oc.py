#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
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
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('font',**{'family':'serif','serif':['Latin Modern Roman']})
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# custom modules
import constants as cst
import ancillary as anc

# ==============================================================================
class sim_data:
  
  def __init__(self):
    self.body_id = 0
    self.nTTs, self.ncols = 0, 0
    
    self.epo = []
    self.TTo = []
    self.eTTo = []
    self.TTs = []
    self.dTTos = []
    self.TTstat = []
    self.Teph = []
    self.Peph = []
    self.eTPeph = []
    self.TTlin = []
    self.oc_o = []
    self.oc_s = []
    
    self.T41o = []
    self.eT41o = []
    self.T41s = []
    self.dT41os = []
    self.T41stat = []
  
  def update_sim_data(self, body_id, sim_in):
    self.body_id = body_id
    self.nTTs, self.ncols = np.shape(sim_in)
    
    self.epo = sim_in[:,0].astype(int)
    self.TTo = sim_in[:,1]
    self.eTTo = sim_in[:,2]
    self.TTs = sim_in[:,3]
    self.dTTos = sim_in[:,4]
    self.TTstat = sim_in[:,5]
    
    self.epo, self.Teph, self.Peph, self.eTPeph = \
      anc.compute_lin_ephem(self.TTo, eT0=self.eTTo, epoin=self.epo,
                            modefit='odr'
                            )
    self.TTlin = self.Teph + self.epo*self.Peph
    self.oc_o = self.TTo - self.TTlin
    self.oc_s = self.TTs - self.TTlin
    
    if(self.ncols == 11):
      self.T41o = sim_in[:,6]/1440. # duration saved in min
      self.eT41o = sim_in[:,7]/1440.
      self.T41s = sim_in[:,8]/1440.
      self.dT41os = sim_in[:,9]/1440.
      self.T41stat = sim_in[:,10]
    
# ==============================================================================

class oc_sample:
  
  def __init__(self, idnum, TTs_s, T41s_s):
    self.idplanet = int(idnum)
    self.TTs = np.array(TTs_s)
    self.nTTs = np.shape(self.TTs)[0]
    self.epo = np.arange(0,self.nTTs,1)
    self.TTlin = np.zeros((self.nTTs))
    self.oc = np.zeros((self.nTTs))
    self.T41s = np.array(T41s_s)/1440. # computed in min in trades
  
  def update(self, Teph, Peph):
    self.epo = anc.calculate_epoch(self.TTs, Teph, Peph)
    self.TTlin = Teph + self.epo*Peph
    self.oc = self.TTs - self.TTlin

# ==============================================================================

def get_sim_data(file_in):
  
  # epo 0 TTo 1 eTTo 2 TTs 3 dTTos 4 TTstat 5 T41o 6 eT41o 7 T41s 8 dT41os 9 T41stat 10
  sim_in = np.genfromtxt(file_in)
  body_id = int(file_in.split('NB')[1].split('_simT0')[0])
  
  sim = sim_data()
  sim.update_sim_data(body_id, sim_in)
  
  return sim

# ==============================================================================

def get_simT0_file_list(cli):
  
  file_pattern = os.path.join(cli.full_path, '%d_%d_NB*_simT0.dat' \
                                              %(cli.idsim, cli.lmflag)
                              )
  file_list = np.sort(glob.glob(file_pattern))
  
  return file_list

# ==============================================================================

def read_samples(samples_file = None):
  
  if (samples_file is not None):
    
    sh5 = h5py.File(os.path.abspath(samples_file), 'r')
    n_samples = len(sh5.keys())
    samples = []
    for igr in sh5.keys(): # loop on groups, one for each sample
      TTfield = sh5[igr].keys()
      for kk in TTfield:
        if('TTs_' in kk):
          idpl = int(kk.split('_')[1])
          TTs = sh5['%s/%s' %(igr, kk)][...]
          try:
            T41s = sh5['%s/%s' %(igr, 'T41s_%d' %(idpl))][...]
          except:
            T41s = np.zeros(np.shape(TTs))
          samples.append(oc_sample(idpl, TTs, T41s))

    sh5.close()
    
  else:
    
    samples = None

  return samples

# ==============================================================================

def set_unit(cli):
  
  u_in = cli.ocunit
  
  try:
    if(u_in.lower()   in 's sec second seconds'.split()):
      ocu = [86400., 'sec']
    elif(u_in.lower() in 'm min minute minutes'.split()):
      ocu = [1440., 'min']
    elif(u_in.lower() in 'h hour hours'.split()        ):
      ocu = [24., 'hours']
    else:
      ocu = [1., 'days']
  except:
    ocu = [1., 'days']
  
  return ocu

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
          #verticalalignment='center',
          verticalalignment='bottom',
          fontsize=fontsize,
          transform=ax.transAxes
          )
  
  return 

# ==============================================================================

def plot_oc_T41(cli, file_in, samples):
  
  sim = get_sim_data(file_in)
  
  if(cli.tscale is None):
    tscale = 2440000.5
  else:
    tscale = np.float64(cli.tscale)
    
  ocu = set_unit(cli)
  
  letters = 'a b c d e f g h i j k l m n o p q r s t u v w x y z'.split()
  ibd = sim.body_id - 1
  planet = letters[ibd]
  
  if(samples is not None):
    samples_plt = [ss for ss in samples if ss.idplanet == sim.body_id ]
    nsmp = len(samples_plt)
  else:
    nsmp = 0
  
  lfont=12
  tfont=8
  
  xlabel = 'BJD$_\\textrm{TDB} - %.3f$' %(tscale)
  
  fig = plt.figure(figsize=(6,6))
  fig.subplots_adjust(hspace=0.05, wspace=0.25)
  
  nrows = 3
  if(sim.ncols==11):
    ncols = 2
  else:
    ncols = 1

  malpha = 0.88
  mso=5
  mss=4

  # O-C
  ax = plt.subplot2grid((nrows, ncols), (0,0), rowspan=2)
  set_axis_default(ax, ticklabel_size=tfont, aspect='auto', labeldown=False)
  ax.set_ylabel('O-C (%s)' %(ocu[1]), fontsize=lfont)
  axtitle(ax, labtitle='planet %s: transit times' %(planet),
          fontsize=lfont
          )
  
  ax.axhline(0., color='black', ls='-',lw=0.7)
  
  #x = sim.TTo - tscale
  x = sim.TTlin - tscale
  minx, maxx = np.min(x)-0.88*sim.Peph, np.max(x)+0.88*sim.Peph
  
  # data
  lobs = ax.errorbar(x, sim.oc_o*ocu[0],
                     yerr=sim.eTTo*ocu[0],
                     color='black', ecolor='gray',
                     fmt='o', ms=mso, mfc='white', mec='black',
                     ls='',
                     elinewidth=0.6,
                     capsize=0,
                     zorder=5,
                     label='observations'
                     )
  # trades
  lsim, = ax.plot(x, sim.oc_s*ocu[0],
                  marker='o', ms=mss,
                  mec='lightgray',  mew=0.5, #mec='None',
                  ls='',
                  zorder=6,
                  alpha=malpha,
                  label='simulations'
                  )
  
  # samples
  if(nsmp > 0):
    lsmp = []
    for ss in samples_plt:
      ss.update(sim.Teph, sim.Peph)
      xx = ss.TTlin - tscale
      ll, = ax.plot(xx, ss.oc*ocu[0],
                    color='gray',
                    marker='None', ls='-', lw=0.4,
                    zorder=4,
                    alpha=0.4,
                    label='samples'
                    )
      lsmp.append(ll)
  
    ax.legend(handles=[lobs, lsim, lsmp[0]], loc='best', fontsize=lfont)
  else:
    ax.legend(loc='best', fontsize=lfont)
  
  ax.set_xlim(minx, maxx)
  
  # residuals TTo-TTs
  ax = plt.subplot2grid((nrows, ncols), (2,0), rowspan=1)
  set_axis_default(ax, ticklabel_size=tfont, aspect='auto')
  ax.set_ylabel('res (%s)' %(ocu[1]), fontsize=lfont)
  ax.set_xlabel(xlabel, fontsize=lfont)
  
  ax.axhline(0., color='black', ls='-',lw=0.7)
  
  # data - trades
  ax.errorbar(x, sim.dTTos*ocu[0],
              yerr=sim.eTTo*ocu[0],
              fmt='o', ms=mss, mec='lightgray', mew=0.5,
              ls='',
              ecolor='gray',
              elinewidth=0.6,
              capsize=0,
              zorder=6,
              )
  
  ax.set_xlim(minx, maxx)
  
  if(sim.ncols == 11):
    # T41 duration

    ax = plt.subplot2grid((nrows, ncols), (0,1), rowspan=2)
    set_axis_default(ax, ticklabel_size=tfont, aspect='auto', labeldown=False)
    ax.set_ylabel('T$_{14}$ (%s)' %(ocu[1]), fontsize=lfont)
    axtitle(ax, labtitle='planet %s: durations as T$_4$-T$_1$' %(planet),
            fontsize=lfont
            )
    
    #ax.axhline(np.median(data[:,6]), color='black', ls='-',lw=0.7)
    
    # data
    ax.errorbar(x, sim.T41o*ocu[0],
                yerr=sim.eT41o*ocu[0],
                color='black', ecolor='gray',
                fmt='o', ms=mso, mfc='white', mec='black',
                ls='',
                elinewidth=0.6,
                capsize=0,
                zorder=5,
                label='observations'
                )
    # trades
    ax.plot(x, sim.T41s*ocu[0],
            marker='o', ms=mss,
            mec='lightgray',  mew=0.5, #mec='None',
            ls='',
            zorder=6,
            alpha=malpha,
            label='simulations'
            )
    
    # samples
    if(nsmp > 0):
      for ss in samples_plt:
        xx = ss.TTlin - tscale
        #print np.min(ss.T41s*ocu[0]),np.max(ss.T41s*ocu[0])
        ax.plot(xx, ss.T41s*ocu[0],
                color='gray',
                marker='None', ls='-', lw=0.4,
                zorder=4,
                alpha=0.4,
                )
        
    ax.set_xlim(minx, maxx)
    
    # residuals T41o - T41s
    ax = plt.subplot2grid((nrows, ncols), (2,1), rowspan=1)
    set_axis_default(ax, ticklabel_size=tfont, aspect='auto')
    ax.set_ylabel('res (min)', fontsize=lfont)
    ax.set_xlabel(xlabel, fontsize=lfont)
    
    ax.axhline(0., color='black', ls='-',lw=0.7)
    
    # data - trades
    ax.errorbar(x, sim.dT41os*ocu[0],
                yerr=sim.eT41o*ocu[0],
                fmt='o', ms=mss, mec='lightgray', mew=0.5,
                ls='',
                ecolor='gray',
                elinewidth=0.6,
                capsize=0,
                zorder=6,
                )
      
    ax.set_xlim(minx, maxx)
  
  folder_out = os.path.join(os.path.dirname(file_in), 'plots')
  if(not os.path.isdir(folder_out)): os.makedirs(folder_out)
  fname = os.path.basename(file_in)
  plt_file = os.path.join(folder_out, os.path.splitext(fname)[0])
  fig.savefig('%s.png' %(plt_file), bbox_inches='tight', dpi=200)
  print 'Saved plot into:'
  print '%s' %('%s.png' %(plt_file))
  fig.savefig('%s.pdf' %(plt_file), bbox_inches='tight', dpi=72)
  print '%s' %('%s.pdf' %(plt_file))
  plt.close(fig)
  
  
  return sim


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
  parser.add_argument('-u', '--u',
                      '-unit', '--unit',
                      action='store',
                      dest='ocunit',
                      default='d',
                      help='Unit of the O-C/T41 plot. Possible (lower or upper case): d, day, days or h, hour, hours or m, min, minute, minutes or s, sec, second, seconds. Default is d.'
                      )
  parser.add_argument('-samples-file', '--samples-file',
                      action='store',
                      dest='samples_file',
                      default='None',
                      help='HDF5 file with T0 and RV from emcee samples to overplot on O-Cs.'
                      )

  cli = parser.parse_args()
  
  full_path = os.path.abspath(cli.full_path)
  
  cli.idsim = int(cli.idsim)
  
  cli.lmflag = int(cli.lmflag)
  
  if(cli.samples_file.lower() == 'none'):
    cli.samples_file = None
  else:
    if(os.path.isfile(os.path.abspath(cli.samples_file))):
      cli.samples_file = os.path.abspath(cli.samples_file)
    else:
      cli.samples_file = None
      
  return cli

# ==============================================================================

def main():

  cli = get_args()
  samples = read_samples(samples_file = cli.samples_file)
  file_list = get_simT0_file_list(cli)
  sims = []
  for f in file_list:
    sims.append(plot_oc_T41(cli, f, samples))

  

  return

# ==============================================================================
# ==============================================================================

if __name__ == '__main__':
  main()
# ==============================================================================
