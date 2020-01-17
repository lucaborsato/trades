#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import sys
import argparse
import time
import os
import numpy as np # array
import h5py

# ==============================================================================
# custom modules
script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), '../'))
sys.path.append(module_path)
import constants as cst # local constants module
import ancillary as anc
import pytrades_lib

# ==============================================================================
def get_emcee_run(cli):

  # read emcee data
  emcee_file, _, _ = anc.get_emcee_file_and_best(cli.full_path, cli.temp_status)
  
  names_par, _, chains, _, autocor_time, lnprobability, _, completed_steps = \
      anc.get_data(emcee_file, cli.temp_status)

  nfit, nwalkers, nruns, nburnin, nruns_sel = \
      anc.get_emcee_parameters(chains, cli.temp_status, cli.nburnin, completed_steps)

  chains_T_full = np.zeros((nruns, nwalkers, nfit))
  for ii in range(0,nfit):
    chains_T_full[:,:,ii] = chains[:,:nruns,ii].T # transpose
    
  chains_T, flatchain_posterior_0, lnprob_burnin, _ = \
    anc.thin_the_chains(cli.use_thin, nburnin, nruns, nruns_sel, autocor_time, chains_T_full, lnprobability,
                        burnin_done=False
                        )
  # lambda fix
  flatchain_post = anc.fix_lambda(flatchain_posterior_0,
                                  names_par
                                  )

  return names_par, chains_T, flatchain_post, lnprob_burnin, emcee_file

class AanalysisTRADES:

  def __init__(self, cli):
    self.cli = cli
    # init trades
    pytrades_lib.pytrades.initialize_trades(os.path.join(cli.full_path, ''), '', 1)
    np.random.seed(seed=cli.seed)
    # nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case_list = anc.get_fitted(cli.full_path)
    self.nfit, self.NB, _, self.id_fit, _, _, self.cols_list, self.case_list = anc.get_fitted(cli.full_path)
    self.ndata = pytrades_lib.pytrades.ndata
    self.nfree = pytrades_lib.pytrades.nfree
    self.dof   = pytrades_lib.pytrades.dof
    self.all_parameters = pytrades_lib.pytrades.system_parameters
    self.fit_parameters = pytrades_lib.pytrades.fitting_parameters
    
    # get data from the hdf5 file
    self.names_par, self.chains_T, self.flatchain_post, self.lnprob_burnin, self.emcee_file = get_emcee_run(cli)
    
  def save_posterior(self):
    self.posterior_file = os.path.join(self.cli.full_path, 'posterior.hdf5')
    p_h5f = h5py.File(self.posterior_file, 'w')
    p_h5f.create_dataset('posterior', data=self.flatchain_post, dtype=np.float64)
    p_h5f.create_dataset('loglikelihood', data=self.lnprob_burnin.reshape((-1)), dtype=np.float64)
    p_h5f['posterior'].attrs['nfit'] = self.nfit
    p_h5f['posterior'].attrs['nposterior'] = np.shape(self.flatchain_post)[0]

    p_h5f.create_dataset('parameter_names', data=anc.encode_list(self.names_par),
          dtype="S15"
          )

    p_h5f.close()
    anc.print_both(' Saved posterior file: {}'.format(self.posterior_file))

  def fit_to_physical(self):

    self.fit_median = np.median(self.flatchain_post, axis=0)
    self.fit_mean   = np.mean(self.flatchain_post, axis=0)
    self.fit_mode   = 3*self.fit_median - 2*self.fit_mean

    print("")
    print("FIT FROM POSTERIOR")
    print("par median mean mode")
    for i,n in enumerate(self.names_par):
      print("{0:15s} {1:15.8f} {2:15.8f} {3:15.8f}".format(n, 
            self.fit_median[i], self.fit_mean[i], self.fit_mode[i])
            )


    self.phy_posterior = []
    self.names_phy     = []
    self.phy_median    = []
    self.phy_mean      = []
    self.phy_mode      = []

    # computes mass conversion factor
    MR_star = pytrades_lib.pytrades.mr_star
    m_factor_0, self.mass_unit = anc.mass_type_factor(Ms=1.0, mtype=self.cli.m_type, mscale=False)
    
    Ms_gaussian = MR_star[0,0] + np.random.normal(0.0, 1.0, size=(np.shape(self.flatchain_post)[0]))*MR_star[0,1] # if exists an error on the mass, it creates a Normal distribution for the values and use it to re-scale mp/Ms to mp.
    self.m_factor_boot = m_factor_0 * Ms_gaussian # given the factor from Msun to mass_unit it multiply it by the Normal Mstar.
    self.m_factor = m_factor_0 * MR_star[0,0]

    for ifit, name in enumerate(self.names_par):
      if('Ms' in name):
        self.flatchain_post[:,ifit] = m_factor_0 * self.flatchain_post[:,ifit]

    # kep_elem =  M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg x NB
    self.kep_elem = np.zeros((self.NB, 9))
    self.kep_elem[0,:2] = self.all_parameters[:2]
    for i in range(1,self.NB):
      self.kep_elem[i,0] = self.all_parameters[2+(i-1)*8] # M
      self.kep_elem[i,1] = self.all_parameters[3+(i-1)*8] # R
      self.kep_elem[i,2] = self.all_parameters[4+(i-1)*8] # P
      self.kep_elem[i,4] = self.all_parameters[5+(i-1)*8] # e
      self.kep_elem[i,5] = self.all_parameters[6+(i-1)*8] # w
      self.kep_elem[i,6] = self.all_parameters[7+(i-1)*8] # mA
      self.kep_elem[i,7] = self.all_parameters[8+(i-1)*8] # i
      self.kep_elem[i,8] = self.all_parameters[9+(i-1)*8] # lN
    
    self.names_phy, der_posterior = anc.compute_derived_posterior(self.names_par, self.kep_elem, 
                                                                  self.id_fit, self.case_list,
                                                                  self.cols_list, self.flatchain_post,
                                                                  conv_factor=self.m_factor_boot
                                                                  )
    self.phy_posterior = anc.derived_posterior_check(self.names_phy, der_posterior)
    self.phy_median    = np.median(self.phy_posterior, axis=0)
    self.phy_mean      = np.mean(self.phy_posterior, axis=0)
    self.phy_mode      = 3*self.phy_median-2*self.phy_mean

    print("PHY FROM POSTERIOR")
    print("par median mean mode")
    for i,n in enumerate(self.names_phy):
      print("{0:15s} {1:15.8f} {2:15.8f} {3:15.8f}".format(n, 
            self.phy_median[i], self.phy_mean[i], self.phy_mode[i])
            )

    all_median = self.all_parameters.copy()
    all_mean   = self.all_parameters.copy()
    all_mode   = self.all_parameters.copy()
    for ifit, name in enumerate(self.names_par):
      if('P' in name):
        id_pl = int(name.split('P')[1])
        pall  = 4 + (id_pl-2)*8
        all_median[pall] = self.fit_median[ifit]
        all_mean[pall]   = self.fit_mean[ifit]
        all_mode[pall]   = self.fit_mode[ifit]
      elif('i' == name[0]):
        id_pl = int(name.split('i')[1])
        iall  = 8 + (id_pl-2)*8
        all_median[iall] = self.fit_median[ifit]
        all_mean[iall]   = self.fit_mean[ifit]
        all_mode[iall]   = self.fit_mode[ifit]

    for iphy, name in enumerate(self.names_phy):
      if('m' in name and name[1] != 'A'):
        id_pl = int(name.split('m')[1])
        iall  = 2 + (id_pl-2)*8
        all_median[iall] = self.phy_median[iphy]/self.m_factor/m_factor_0
        all_mean[iall]   = self.phy_mean[iphy]/self.m_factor/m_factor_0
        all_mode[iall]   = self.phy_mode[iphy]/self.m_factor/m_factor_0
      elif('e' in name):
        id_pl = int(name.split('e')[1])
        eall  = 5 + (id_pl-2)*8
        wall  = 6 + (id_pl-2)*8
        all_median[eall] = self.phy_median[iphy]
        all_mean[eall]   = self.phy_mean[iphy]
        all_mode[eall]   = self.phy_mode[iphy]
        all_median[wall] = self.phy_median[iphy+1]
        all_mean[wall]   = self.phy_mean[iphy+1]
        all_mode[wall]   = self.phy_mode[iphy+1]
      elif('mA' in name):
        id_pl = int(name.split('mA')[1])
        iall  = 7 + (id_pl-2)*8
        all_median[iall] = self.phy_median[iphy]
        all_mean[iall]   = self.phy_mean[iphy]
        all_mode[iall]   = self.phy_mode[iphy]
      elif('i' in name):
        id_pl = int(name.split('m')[1])
        iall  = 8 + (id_pl-2)*8
        all_median[iall] = self.phy_median[iphy]
        all_mean[iall]   = self.phy_mean[iphy]
        all_mode[iall]   = self.phy_mode[iphy]
      elif('lN' in name):
        id_pl = int(name.split('m')[1])
        iall  = 9 + (id_pl-2)*8
        all_median[iall] = self.phy_median[iphy]
        all_mean[iall]   = self.phy_mean[iphy]
        all_mode[iall]   = self.phy_mode[iphy]

    # POI USO QUESTO PER PASSARE DA all_PARAMETRO A fit_PARAMETRO
     #WARNING: QUAN NON COVERTE IN VALORI GIUSTI (E,W) -> (SQRTECOSW, SQRTESINW)
    self.fit_median = pytrades_lib.pytrades.init_fit_parameters(all_median, self.nfit)
    self.fit_mean   = pytrades_lib.pytrades.init_fit_parameters(all_mean, self.nfit)
    self.fit_mode   = pytrades_lib.pytrades.init_fit_parameters(all_mode, self.nfit)

    print("FIT FROM PHY")
    print("par median mean mode")
    for i,n in enumerate(self.names_par):
      print("{0:15s} {1:15.8f} {2:15.8f} {3:15.8f}".format(n, 
            self.fit_median[i], self.fit_mean[i], self.fit_mode[i])
            )


# ==============================================================================

# ==============================================================================

# main

def main():
  print()
  print(' TRADES: EMCEE confidence intervals')
  print()

  cli = anc.get_args()

  analysis = AanalysisTRADES(cli)
  analysis.fit_to_physical()

  pytrades_lib.pytrades.deallocate_variables()

  return

# ==============================================================================
# ==============================================================================

if __name__ == "__main__":
  main()
