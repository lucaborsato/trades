# TRADES
  

**`TRADES`: TRAnsits and Dynamics of Exoplanetary Systems**  
I have developed a computer program for determining 
the possible physical and dynamical configurations of extra-solar planetary 
systems from observational data.  
`TRADES` models the dynamics of multiple planet systems and
reproduces the observed transit times (T0s, or mid-transit times),
full transit durations (T14 or T41),
 and radial velocities (RVs).
These T0s and RVs are computed during the integration of the planetary orbits.  

Download and compile `fortran90` and `python` libraries:  

```bash
git clone https://github.com/lucaborsato/trades.git
cd trades/src
make cleanall
make full_parallel_release
```

To use with `python` as it would be in [PyORBIT](https://github.com/LucaMalavolta/PyORBIT) 
see notebook [import_trades_for_PyORBIT](trades_example/python_examples/import_trades_for_PyORBIT.ipynb) 
in [trades_example/python_examples/](trades_example/python_examples/).  
In this way you will integrates the orbits of the planets and it will return the Radial Velocities (RVs) if present, and the Transit Times (T0s), if present, with the Transit Durations (T14s) and the orbital elements at each T0s.  

To use as it is you can:  

1. copy a folder in [trades_example](trades_example) based on the number of planets (i.e. 2p for 2, 3 for 3p, and so on)  
2. copy in it the [configuration.yml](trades_example/configuration.yml) and adapt the `run` section, change `PyDE` and `emcee` sections etc  
3. from terminal type:  
  > `python /path/to/pytrades/trades_emcee.py --input configuration.yml`  
4. after it finished you have to modify the `analysis`, `OC` and `RV` sections of the `configuration.yml` file and you will have the analysis results with  
  > `python /path/to/pytrades/trades_emcee_analysis.py --input configuration.yml`  

You can look at the `trades_emcee_analysis.py` to learn how to use the different part of `pytrades` and there are a few other notebooks to plot `PyDE` run, `emcee` chains, RV and OCs.  

---
**`TRADES` v2.19.0 by Luca Borsato - 2016-2023**  

Long README and description in [README_long](README_long.md).
