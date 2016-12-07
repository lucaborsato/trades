#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import sys
import argparse
import os
import ancillary as anc
import multiprocessing as mp
import subprocess as sp
import time


def geweke(cli):

  print 'RUN geweke.py'
  pyscript = os.path.join(os.path.abspath('./local_trades/trades-dev/pytrades'), 'geweke.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', '0', '-m', str(cli.m_type), '-t', str(cli.temp_status), '-s', str(cli.sel_steps)], stdout = sp.PIPE, stderr = sp.PIPE, bufsize=1)
  sout, serr = run.communicate()  
  print 'COMPLETED geweke.py'
  
  return

def gelman_rubin(cli):
  
  print 'RUN gelman_rubin.py'
  pyscript = os.path.join(os.path.abspath('./local_trades/trades-dev/pytrades'), 'gelman_rubin.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', '0', '-m', str(cli.m_type), '-t', str(cli.temp_status), '-s', str(cli.sel_steps)], stdout = sp.PIPE, stderr = sp.PIPE, bufsize=1)
  sout, serr = run.communicate()  
  print 'COMPLETED gelman_rubin.py'
  
  return

def chains(cli):
  
  print 'RUN chains_summary_plot.py'
  pyscript = os.path.join(os.path.abspath('./local_trades/trades-dev/pytrades'), 'chains_summary_plot.py')
  run = sp.Popen(['python', pyscript, '-p',cli.full_path, '-nb', str(cli.npost), '-m', str(cli.m_type), '-t', str(cli.temp_status), '-u', str(cli.use_thin), '--sample', str(cli.sample_str) ], stdout = sp.PIPE, stderr = sp.PIPE, bufsize=1)
  sout, serr = run.communicate()  
  print 'COMPLETED chains_summary_plot.py'
  
  return

def correlation(cli):
  
  print 'RUN correlation_triangle_plot.py'
  pyscript = os.path.join(os.path.abspath('./local_trades/trades-dev/pytrades'), 'correlation_triangle_plot.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', str(cli.npost), '-m', str(cli.m_type), '-t', str(cli.temp_status), '-u', str(cli.use_thin), '-c', str(cli.cumulative), '--sample', str(cli.sample_str) ], stdout = sp.PIPE, stderr = sp.PIPE, bufsize=1)
  sout, serr = run.communicate()  
  print 'COMPLETED correlation_triangle_plot.py'
  
  return


def ci_and_derived(cli):
  
  print 'RUN confidence_intervals.py'
  pyscript = os.path.join(os.path.abspath('./local_trades/trades-dev/pytrades'), 'confidence_intervals.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', str(cli.npost), '-m', str(cli.m_type), '-t', str(cli.temp_status), '-u', str(cli.use_thin), '--sample', str(cli.sample_str), '--seed', str(cli.seed) ], stdout = sp.PIPE, stderr = sp.PIPE, bufsize=1)
  sout, serr = run.communicate()  
  print 'COMPLETED confidence_intervals.py'
  
  print 'RUN derived_correlation_plot.py'
  pyscript = os.path.join(os.path.abspath('./local_trades/trades-dev/pytrades'), 'derived_correlation_plot.py')
  run = sp.Popen(['python', pyscript, '-p', cli.full_path, '-nb', str(cli.npost), '-m', str(cli.m_type), '-t', str(cli.temp_status), '-u', str(cli.use_thin)], stdout = sp.PIPE, stderr = sp.PIPE, bufsize=1)
  sout, serr = run.communicate()  
  print 'COMPLETED derived_correlation_plot.py'
  
  return

# Define the function for checking the process status
def status(proc):
  #print proc.is_alive(), type(proc.is_alive())
  if (proc.is_alive()==True):
    return 'alive'
  elif (proc.is_alive()==False):
    return 'dead'
  else:
    return proc.is_alive()


def main():
  cli = anc.get_args()

  threads = []
  
  # run geweke
  pp = mp.Process(target=geweke, args=(cli,))
  pp.start()
  threads.append(pp)
  
  # run gelman-rubin
  pp = mp.Process(target=gelman_rubin, args=(cli,))
  pp.start()
  threads.append(pp)
  
  # run chains plot
  pp = mp.Process(target=chains, args=(cli,))
  pp.start()
  threads.append(pp)
  
  # run correlation plot
  pp = mp.Process(target=correlation, args=(cli,))
  pp.start()
  threads.append(pp)

  # run credible intervals and derived correlation plot
  pp = mp.Process(target=ci_and_derived, args=(cli,))
  pp.start()
  threads.append(pp)


  for pp in threads:
    pp.join()
  
  # Check processes status
  for pp in threads:
    print "Process ",pp, " @ ", pp.pid, " is ", status(pp)
  # Wait processes to finish
  while len(threads) != 0:
    for pp in threads:
      if status(pp) == "dead":
        threads.pop(threads.index(pp))
    # loose 10 seconds before checking again
    time.sleep(5)
  
  print 'COMPLETED ALL SCRIPTS'
  
  return


if __name__ == "__main__":
  main()
  
  
  
