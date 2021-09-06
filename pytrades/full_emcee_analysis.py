#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import sys
import argparse
import os
import multiprocessing as mp
import subprocess as sp
import time
# ==============================================================================
import ancillary as anc
# ==============================================================================
# ==============================================================================

script_dir = os.path.dirname(os.path.realpath(__file__))

def geweke(cli):

  log_folder = os.path.join(cli.full_path, 'logs')
  log_file = os.path.join(log_folder, 'log_geweke.txt')
  olog = open(log_file, 'w')
  
  print('RUN geweke.py')
  # pyscript = os.path.join(os.path.abspath(cli.pyscript), 'geweke.py')
  pyscript = os.path.join(script_dir, 'geweke.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', '0', '-m', str(cli.m_type), '-t', str(cli.temp_status), '-s', str(cli.sel_steps)], 
                 #stdout = sp.PIPE, stderr = sp.PIPE,
                 stdout = olog, stderr = olog,
                 bufsize=1
                 )
  # sout, serr = run.communicate()
  _, _ = run.communicate()
  print('COMPLETED geweke.py')
  olog.close()
  
  return

# ==============================================================================

def gelman_rubin(cli):
  
  log_folder = os.path.join(cli.full_path, 'logs')
  log_file = os.path.join(log_folder, 'log_gelman_rubin.txt')
  olog = open(log_file, 'w')
  
  print('RUN gelman_rubin.py')
  # pyscript = os.path.join(os.path.abspath(cli.pyscript), 'gelman_rubin.py')
  pyscript = os.path.join(script_dir, 'gelman_rubin.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', '0', '-m', str(cli.m_type), '-t', str(cli.temp_status), '-s', str(cli.sel_steps)],
                 #stdout = sp.PIPE, stderr = sp.PIPE,
                 stdout = olog, stderr = olog,
                 bufsize=1
                 )
  
  _, _ = run.communicate()  
  print('COMPLETED gelman_rubin.py')
  olog.close()
  
  return

# ==============================================================================

def chains(cli):
  
  log_folder = os.path.join(cli.full_path, 'logs')
  log_file = os.path.join(log_folder, 'log_chains.txt')
  olog = open(log_file, 'w')
  
  print('RUN chains_summary_plot.py')
  # pyscript = os.path.join(os.path.abspath(cli.pyscript), 'chains_summary_plot.py')
  pyscript = os.path.join(script_dir, 'chains_summary_plot.py')
  run = sp.Popen(['python', pyscript, '-p',cli.full_path, '-nb', str(cli.nburnin), '-m', str(cli.m_type), '-t', str(cli.temp_status), '-u', str(cli.use_thin), '--overplot', str(cli.overplot)],
                 #stdout = sp.PIPE, stderr = sp.PIPE,
                 stdout = olog, stderr = olog,
                 bufsize=1
                 )
  _, _ = run.communicate()  
  print('COMPLETED chains_summary_plot.py')
  olog.close()
  
  return

# ==============================================================================

def correlation(cli):
  
  log_folder = os.path.join(cli.full_path, 'logs')
  log_file = os.path.join(log_folder, 'log_corr_plot_fitted.txt')
  olog = open(log_file, 'w')
  
  print('RUN correlation_triangle_plot.py')
  # pyscript = os.path.join(os.path.abspath(cli.pyscript), 'correlation_triangle_plot.py')
  pyscript = os.path.join(script_dir, 'correlation_triangle_plot.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', str(cli.nburnin), '-m', str(cli.m_type), '-t', str(cli.temp_status), '-u', str(cli.use_thin), '-c', str(cli.cumulative), '--overplot', str(cli.overplot) ],
                 #stdout = sp.PIPE, stderr = sp.PIPE,
                 stdout = olog, stderr = olog,
                 bufsize=1
                 )
  _, _ = run.communicate()  
  print('COMPLETED correlation_triangle_plot.py')
  olog.close()
  
  return

# ==============================================================================

def derived_correlation(cli):
  
  log_folder = os.path.join(cli.full_path, 'logs')
  log_file = os.path.join(log_folder, 'log_corr_plot_derived.txt')
  olog = open(log_file, 'w')
  
  print('RUN derived_correlation_plot.py')
  # pyscript = os.path.join(os.path.abspath(cli.pyscript), 'derived_correlation_plot.py')
  pyscript = os.path.join(script_dir, 'derived_correlation_plot.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', str(cli.nburnin), '-m', str(cli.m_type), '-t', str(cli.temp_status), '-u', str(cli.use_thin), '--overplot', str(cli.overplot)],
                 #stdout = sp.PIPE, stderr = sp.PIPE,
                 stdout = olog, stderr = olog,
                 bufsize=1
                 )
  _, _ = run.communicate()  
  print('COMPLETED derived_correlation_plot.py')
  olog.close()
  
  return

# ==============================================================================

def ci_chains_correlation(cli):
  
  log_folder = os.path.join(cli.full_path, 'logs')
  log_file = os.path.join(log_folder, 'log_confidence_intervals.txt')
  olog = open(log_file, 'w')
  
  print('RUN confidence_intervals.py')
  # pyscript = os.path.join(os.path.abspath(cli.pyscript), 'confidence_intervals.py')
  pyscript = os.path.join(script_dir, 'confidence_intervals.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', str(cli.nburnin), '-m', str(cli.m_type), '-t', str(cli.temp_status), '-u', str(cli.use_thin), '--sample', str(cli.sample_str), '--seed', str(cli.seed), '--overplot', str(cli.overplot), '--ad-hoc', str(cli.adhoc), '--n-samples', str(cli.n_samples) ],
                 #stdout = sp.PIPE, stderr = sp.PIPE,
                 stdout = olog, stderr = olog,
                 bufsize=1
                 )
  _, _ = run.communicate()  
  print('COMPLETED confidence_intervals.py')
  olog.close()
  
  sub_th = []
  
  # run chains plot
  sub_pp = mp.Process(target=chains, args=(cli,))
  sub_pp.start()
  sub_th.append(sub_pp)
  
  # run correlation plot
  sub_pp = mp.Process(target=correlation, args=(cli,))
  sub_pp.start()
  sub_th.append(sub_pp)
  
  # run derived correlation plot
  sub_pp = mp.Process(target=derived_correlation, args=(cli,))
  sub_pp.start()
  sub_th.append(sub_pp)
  
  for sub_pp in sub_th:
    sub_pp.join()
    
  # Check processes status
  for sub_pp in sub_th:
    print("Process ",sub_pp, " @ ", sub_pp.pid, " is ", status(sub_pp))
  # Wait processes to finish
  while len(sub_th) != 0:
    for sub_pp in sub_th:
      if status(sub_pp) == "dead":
        sub_th.pop(sub_th.index(sub_pp))
    # loose 10 seconds before checking again
    time.sleep(5)
  
  return

# ==============================================================================

# Define the function for checking the process status
def status(proc):
  #print proc.is_alive(), type(proc.is_alive())
  if (proc.is_alive()==True):
    return 'alive'
  elif (proc.is_alive()==False):
    return 'dead'
  else:
    return proc.is_alive()

# ==============================================================================

def main():
  cli = anc.get_args()
  
  log_folder = os.path.join(cli.full_path, 'logs')
  if(not os.path.isdir(log_folder)):
    os.makedirs(log_folder)

  threads = []
  
  # run geweke
  pp = mp.Process(target=geweke, args=(cli,))
  pp.start()
  threads.append(pp)
  
  # run gelman-rubin
  pp = mp.Process(target=gelman_rubin, args=(cli,))
  pp.start()
  threads.append(pp)
  
  ## run chains plot
  #pp = mp.Process(target=chains, args=(cli,))
  #pp.start()
  #threads.append(pp)
  
  ## run correlation plot
  #pp = mp.Process(target=correlation, args=(cli,))
  #pp.start()
  #threads.append(pp)

  # run credible intervals, chains, and derived correlation plot
  pp = mp.Process(target=ci_chains_correlation, args=(cli,))
  pp.start()
  threads.append(pp)


  for pp in threads:
    pp.join()
  
  # Check processes status
  for pp in threads:
    print("Process ",pp, " @ ", pp.pid, " is ", status(pp))
  # Wait processes to finish
  while len(threads) != 0:
    for pp in threads:
      if status(pp) == "dead":
        threads.pop(threads.index(pp))
    # loose 10 seconds before checking again
    time.sleep(2.5)
  
  print('COMPLETED ALL SCRIPTS')
  
  return

# ==============================================================================
# ==============================================================================
if __name__ == "__main__":
  main()
# ==============================================================================
# ==============================================================================
  
  
