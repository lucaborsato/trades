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

I suggest to create an `anaconda environment`, install all the dependencies and then do `make full_parallel_release`.  
**WARNING**: not working with `python>=3.11s` and `numpy>=1.26` due to deprecated `distutils` for `f2py`, now using `meson`.  
**SUGGESTION**: `conda create --name trades_env python=3.10 numpy=1.23.5` and other packages should be installed properly.  

Possible way to create an environment and install all dependencies:  

create env with python and numpy and matplotlib
`conda create --name trades_env python=3.10 numpy=1.23.5 matplotlib`  
activate env  
`conda activate trades_env`  
optional `cython`:  
`conda install cython`  
install `h5py`, `pyyaml`, `tqdm`, `emcee`, `scipy`, `pandas`, `pygtc`:  

`conda install h5py`  
`pip install pyyaml`  
`conda install -c conda-forge tqdm` or `pip install tqdm`  
`pip install corner`
`conda install -c conda-forge emcee`  
`conda install pandas` or `conda install -c conda-forge pandas`  
`conda install scipy` or `conda install -c conda-forge scipy`  
`conda install -c conda-forge astropy`  
`conda install -c conda-forge pygtc`  

some packages needs `celerite`:  
`conda install -c conda-forge pybind11`  
`conda install -c conda-forge celerite`  

for photo-dynamical approach it is needed to install `pytransit`:  
`conda install numba cudatoolkit`  
`pip install semantic-version`  
`conda install -c conda-forge arviz`  
`git clone https://github.com/hpparvi/PyTransit.git`  
`cd PyTransit`  
`python setup.py install`  

Still not fully tested and in development, `dynesty` and `ultranest`:  
`pip install dynesty` 
`conda install -c conda-forge ultranest`  
(could use `mpi4py`, I was not able to install it)  


To use with `python` as it would be in [PyORBIT](https://github.com/LucaMalavolta/PyORBIT) 
see notebook [import_trades_for_python](trades_example/python_examples/import_trades_for_python.ipynb) 
in [trades_example/python_examples/](trades_example/python_examples/).  
In this way you will integrates the orbits of the planets and it will return the Radial Velocities (RVs) if present, and the Transit Times (T0s), if present, with the Transit Durations (T14s) and the orbital elements at each T0s.  
Currently, `TRADES` is ready to be used as a photo-dynamical (`photoTRADES`) code,
the example, as python notebook, is still in development.  

To use as it was initially intended you can:  

1. copy a folder in [trades_example](trades_example) based on the number of planets (i.e. 2p for 2, 3 for 3p, and so on). The `base_2p_grid` is a version for the old `Fortran` grid version.  
2. copy inside the folder, if not present, the [configuration.yml](trades_example/configuration.yml) and adapt the `run` section, change `PyDE` and `emcee` sections etc (`ultranest` and `dynesty` in development, not fully tested).  
3. from terminal type:  
  > `python /path/to/pytrades/trades_emcee.py --input configuration.yml`  
4. after it finished you have to modify the `analysis`, `OC` and `RV` sections of the `configuration.yml` file and you will have the analysis results with  
  > `python /path/to/pytrades/trades_emcee_analysis.py --input configuration.yml`  

You can look at the `trades_emcee_analysis.py` to learn how to use the different part of `pytrades` and there are a few other notebooks to plot `PyDE` run, `emcee` chains, RV and OCs.  

---
**`TRADES` v2.20.0 by Luca Borsato - 2016-2023**  

`TRADES` v2.19.0 by Luca Borsato - 2016-2023  

Long README and description in [README_long](README_long.md).
