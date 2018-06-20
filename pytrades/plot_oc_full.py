#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# IMPORT

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import numpy as np # array
import os
import sys
import glob
import gc

import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)

# custom modules
import ancillary as anc
script_path = os.path.realpath(__file__)
subpath = os.path.abspath(os.path.join(os.path.dirname(script_path), 'cheops'))
sys.path.append(subpath)
import gls

# ==============================================================================

def get_data_from_trafile(tra_file):
  
  tra_name = os.path.basename(tra_file)
  idlmbody = tra_name.split('_')
  sim_id = idlmbody[0]
  lm_flag = idlmbody[1]
  body_id = int(idlmbody[2].split('NB')[1])
  xx = np.genfromtxt(tra_file, usecols=(0,1))
  TTs = np.sort(xx[:,0]+xx[:,1])
  
  return sim_id, lm_flag, body_id, TTs

# ==============================================================================

class Transit:
  
  def __init__(self, full_file):
    self.full_file = full_file
    self.dirname = os.path.dirname(full_file)
    self.filename = os.path.basename(full_file)
    self.sim_id, self.lm_flag, self.body_id, self.TTs = get_data_from_trafile(full_file)
    self.nTTs = np.shape(self.TTs)[0]
    #self.epo = np.arange(0,self.nTTs)
    #self.Tref = self.TTs[0]
    #self.Pref = self.TTs[1]-self.TTs[0]
    self.epo, self.Tref, self.Pref, self.TPerr = anc.compute_lin_ephem(self.TTs)
    self.linTTs = self.Tref + self.epo*self.Pref
    self.ocd = self.TTs - self.linTTs
    self.amp_ttv = 0.
    self.p_ttv = 0.
    
  def compute_TTV_gls(self):
    dur = self.TTs[-1]-self.TTs[0]
    e1sec = np.ones(np.shape(self.TTs))/86400.
    ogls = gls.Gls((self.TTs, self.ocd, e1sec),
                    Pbeg=3.*self.Pref, Pend=2.*dur,
                    ofac=10, verbose=False)
    self.amp_ttv = ogls.hpstat['amp']
    self.p_ttv = 1./ogls.hpstat['fbest']
    del ogls
    gc.collect()
    

# ==============================================================================

def get_tra_file_list(folder_path, sim_id, lm_flag):
  
  pattern = '%s_%s_NB*_tra.dat' %(str(sim_id), lm_flag)
  tra_file_list = np.sort(glob.glob(os.path.join(folder_path, pattern)))
  
  return tra_file_list

# ==============================================================================

def get_args():
  
  parser = argparse.ArgumentParser(description='PLOT OC')
  
  parser.add_argument('-f', '--f', '-folder-path','--folder-path',  
                      action='store', dest='folder_path', required=True, 
                      help='The path (absolute or relative) with simulation files for TRADES.')
  
  parser.add_argument('-s', '--s', '-sim-id','--sim-id',  
                      action='store', dest='sim_id', default='1', 
                      help='Number ID of the simulation. Default is 1.')

  parser.add_argument('-l', '--l', '-lm-flag','--lm-flag',  
                      action='store', dest='lm_flag', default='0', 
                      help='LM flag: 0 or 1. Default is 0.')

  cli = parser.parse_args()
  cli.folder_path = os.path.abspath(cli.folder_path)
  
  try:
    cli.sim_id = int(cli.sim_id)
    if(cli.sim_id <= 0): cli.sim_id = 1
  except:
    cli.sim_id = 1
  
  if(str(cli.lm_flag) == '1'):
    cli.lm_flag = '1'
  else:
    cli.lm_flag = '0'

  return cli

# ==============================================================================

def main():
  
  cli = get_args()
  
  tra_file_list = get_tra_file_list(cli.folder_path, cli.sim_id, cli.lm_flag)
  
  plot_folder = os.path.join(cli.folder_path, 'plots')
  if(not os.path.isdir(plot_folder)): os.makedirs(plot_folder)
  plot_file = os.path.join(plot_folder, '%d_%s_oc_full.png' %(cli.sim_id, cli.lm_flag))
  
  fig = plt.figure(figsize=(12,12))
  
  plt.tick_params(direction='in')

  plt.axhline(0., color='black', ls='-', lw=0.33, alpha=0.77)
  
  tra_obj = []
  cnt = 0
  for ff in tra_file_list:
    print 'READ FILE %s' %(ff)
    oo = Transit(ff)
    oo.compute_TTV_gls()
    tra_obj.append(oo)
    
    
    print 'body n. %d' %(oo.body_id)
    #print 'ephem: %20.8f + N x %20.8f' %(oo.Tref, oo.Pref)
    if(cnt == 0):
      xscale = oo.TTs[0]
      plt.xlabel('BJD - %.6f' %(xscale), fontsize=12)
      plt.ylabel('O-C (min)', fontsize=12)
    print 'plotting O-C'
    cnt += 1
    labTeph = 'TTlin $= %17.6f + N \\times %13.6f$' %(oo.Tref, oo.Pref)
    print labTeph
    labPTTV = '$P_\\textrm{TTV} = %13.6f$ days' %(oo.p_ttv)
    labATTV = '$A_\\textrm{TTV} = %10.3f$ min' %(oo.amp_ttv*1440.)
    print '%s , %s' %(labPTTV, labATTV)
    plt.plot(oo.TTs-xscale, oo.ocd*1440.,
             marker='o', ms=4, mec='None',
             ls='-', lw=0.45, alpha=0.77,
             label='body \#%d: %s ; %s , %s' %(oo.body_id, labTeph,
                                               labPTTV, labATTV)
             )
  
  print 'Saving file: %s' %(plot_file)
  plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), fontsize=12, ncol=1)
  plt.draw()
  fig.savefig(plot_file, dpi=300)
  plt.close(fig)
  print 'DONE'
  
  return

# ==============================================================================
# ==============================================================================

if __name__ == '__main__':
  main()

# ==============================================================================
