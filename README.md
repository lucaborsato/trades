# TRADES
  

**`TRADES`: TRAnsits and Dynamics of Exoplanetary Systems**  

Keywords: _exoplanets, multi-planet systems, transit times, transit timing variations, TTV, radial velocities, planetary orbits, N-body, numerical integration_  

I have developed a computer program for determining 
the possible physical and dynamical configurations of extra-solar planetary 
systems from observational data.  
`TRADES` models the dynamics of multiple planet systems and
reproduces the observed transit times (T0s, or mid-transit times),
full transit durations (T14 or T41),
 and radial velocities (RVs).
These T0s and RVs are computed during the integration of the planetary orbits.  

### References

[Borsato et al., 2014](https://ui.adsabs.harvard.edu/abs/2014A%26A...571A..38B/abstract): main reference and description of the initial development of `TRADES`.  

[Borsato et al., 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.3233B/abstract): updated version with python interface.  

[Borsato et al., 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...689A..52B/abstract): updated version with photo-dynamical approach.


## Download `trades` source

```bash
git clone https://github.com/lucaborsato/trades.git
```

## Create conda environment


I suggest to create an `anaconda environment`, install all the dependencies and then do `make full_parallel_release`.  
**WARNING**: not working with `python>3.10`, `numpy>1.23.5`, and `setuptools>65.6.3` due to deprecated `distutils` for `f2py`, now using `meson`.  

**Strongly suggested** way to create an environment and install all dependencies:  

create env with python and numpy and matplotlib  

`conda create --name trades_env python=3.10 numpy=1.23.5 matplotlib`  

activate env  
`conda activate trades_env`  

install `cython` (optional, could be used from other packages):  
`conda install cython`  

install `h5py`, `pyyaml`, `tqdm`, `emcee`, `scipy`, `pandas`, `pygtc`:  

`conda install h5py`  
`pip install pyyaml`  
`conda install -c conda-forge tqdm`  
`pip install corner`  
`conda install -c conda-forge emcee`  
`conda install -c conda-forge pandas`  
`conda install -c conda-forge scipy`  
`conda install -c conda-forge astropy`  
`conda install -c conda-forge pygtc`  

some packages needs `celerite`:  
`conda install -c conda-forge pybind11`  
`conda install -c conda-forge celerite`  

for photo-dynamical approach it is needed to install `pytransit`:  
`conda install -c conda-forge numba cudatoolkit`  
`pip install semantic-version`  
`conda install -c conda-forge arviz`  
`git clone https://github.com/hpparvi/PyTransit.git`  
`cd PyTransit`  
`pip install .`  
Test with `import pytransit` to check if some requirements is still missing for `pytransit`.

### Optional - `dynesty` and `ultranest`
Still not fully tested and in development, `dynesty` and `ultranest`:  
`pip install dynesty`  
`conda install -c conda-forge ultranest`  
(could use `mpi4py`, install could lead to some issues)  

### IMPORTANT

Due to upgrades of the different packages, `setuptools` will be upgraded breaking the compilation through `f2py`.
Force `setuptools` to the version `65.6.3` and it should compile fine.  
`conda install setuptools=65.6.3`  

## Compile the `fortran` sources and build the `python` library

Download and compile `fortran90` and `python` libraries:  

If not done:  
```bash
git clone https://github.com/lucaborsato/trades.git
conda activate trades_env
```
Enter `TRADES` `fortran90` source directory, clean from previous compilations, and compile everything:  
```bash
cd trades/src
make cleanall
make full_parallel_release
```

## Use `pyTRADES`

To use with `python` see notebook [import_trades_for_python](trades_example/python_examples/import_trades_for_python.ipynb) 
in [trades_example/python_examples/](trades_example/python_examples/).  
In this way you will integrates the orbits of the planets and it will return the Radial Velocities (RVs) if present, and the Transit Times (T0s), if present, with the Transit Durations (T14s) and the orbital elements at each T0s.  
It will be implemented soon in [PyORBIT](https://github.com/LucaMalavolta/PyORBIT)  
Currently, `TRADES` is ready to be used as a photo-dynamical (`photoTRADES`) code,
the example, as python notebook, is still in development (but it is working).  

An example of how to call `pytrades` and other modules:  
```python
import numpy as np
import os
from pytrades import constants as cst
from pytrades import ancillary as anc
from pytrades import pytrades
...
```

or you can access nested functions as:

```python
from pytrades.ancillary import set_rcParams
from pytrades.pytrades import linear_fit
```

Set number of bodies in the system and times:  

```python
t_start = 0.0  # start time of the integration
t_epoch = t_start  # reference time of the orbital parameters, could be equal to, earlier or later than `t_start`
t_int = 1000.0  # total integration time from the start time in days

n_body = 3  # number of bodies (NOT PLANETS) in the system, that is star + planets
duration_check = 1  # to check or not the duration in this case it is a return flag, 1 means return it

pytrades.args_init(
    n_body,
    duration_check,
    t_epoch=t_epoch,
    t_start=t_start,
    t_int=t_int,
    do_hill_check=False,
    amd_hill_check=False,
    rv_res_gls=False,
)
# n_body (int): Number of bodies.
# duration_check (int): Duration check value.
# t_epoch (int, optional): Optional epoch time value.
# t_start (int, optional): Optional start time value.
# t_int (int, optional): Optional interval time value.
# do_hill_check (bool, optional): Optional flag for Hill check.
# amd_hill_check (bool, optional): Optional flag for AMD Hill check.
```


Set masses, radii, and orbital parameters of the bodies for forward modelling:
```python
# parameter = array([star, planet b, planet c]) = array([body_1, body_2, body_3])
mass   = np.array([1.022, 43.40 * cst.Mears,  29.90 * cst.Mears])  # Masses in Solar unit
radius = np.array([0.958,  8.29 * cst.Rears,   8.08 * cst.Rears])  # Radii in Solar unit
period = np.array([0.0,   19.23891,           38.9853])  # Periods in days
ecc    = np.array([0.0,    0.0609,             0.06691])  # eccentricities
argp   = np.array([0.0,  357.0,              167.5])  # argument of pericenters in degrees
meana  = np.array([0.0,    2.6,              307.4])  # mean anonalies in degrees
inc    = np.array([0.0,   88.9,               89.188])  # inclinations in degrees
longn  = np.array([0.0,  180.0,              179.0])  # longitude of ascending nodes in degrees
```

### Load datasets and get the RVs and T0s to be used in a $\log \mathcal{L}$

Data folder:  
```python
data_folder = os.path.abspath("./trades_example/Kepler-9_example")
```


Read and load into shared memory the Radial Velocity (RV) data:  
```python
# read the RVs -- as you want, here an example
(
    time_rv, # times of observed RVs, same unit and frame system of t_epoch
    rv_mps, # RV in m/s
    erv_mps, # error on RV in m/s
    rv_setid, # RV set IDs, 1, 2, etc depending of the RV source, for example all RVs from HARPS set to 1, RVs from ESPRESSO set to 2
) = np.genfromtxt(
    os.path.join(
      data_folder,
      "obsRV.dat"
    ),
    unpack=True
)

set_rv = np.unique(rv_setid)
n_rvset = len(set_rv)

# deallocation, this is a test,this function deallocates RV data before re-allocate them
pytrades.deallocate_rv_dataset()
# add the RV dataset to the common variable of TRADES
pytrades.set_rv_dataset(time_rv, rv_mps, erv_mps, rv_setid=rv_setid, n_rvset=n_rvset)
```

Read and load into shared memory the Transit times (T0s) data and (optional) durations (T14s in minutes) data:
```python
# read T0s -- as you want, here an example
body_id = 2 # planet b has number 2 in this system with fortran indexing
(
    epoch_b,
    T0s_b,
    eT0s_b,
    t14s_b,
    et14s_b,
    b_sources_id,
) = np.genfromtxt(
    os.path.join(
      data_folder, 
      "NB{}_observations.dat".format(body_id),
    )
    usecols=(0,1,2,3,4,5), 
    unpack=True
)

# for example set only transit times, no durations
pytrades.set_t0_dataset(body_id, epoch_b, T0s_b, eT0s_b, sources_id=b_sources_id)

# deallocate current T0s dataset of b
pytrades.deallocate_t0_dataset(body_id)
# ... and load T0s and T14s
pytrades.set_t0_dataset(body_id, epoch_b, T0s_b, eT0s_b, t14=t14s_b, et14=et14s_b, sources_id=b_sources_id)

body_id = 3 # planet c has number 3 in this system with fortran indexing
(
    epoch_c,
    T0s_c,
    eT0s_c,
    t14s_c,
    et14s_c,
    c_sources_id,
) = np.genfromtxt(
    os.path.join(
      data_folder, 
      "NB{}_observations.dat".format(body_id),
    )
    usecols=(0,1,2,3,4,5), 
    unpack=True
)

pytrades.set_t0_dataset(body_id, epoch_c, T0s_c, eT0s_c, t14=t14s_c, et14=et14s_c, sources_id=b_sources_id)
```

Required to define a flag vector to identify which bodies has to transit or not.
The `star` does not, so first element has to be flagged to `False`.
If you don't know if they transit or not set `star` to `False`, all `other bodies` to `True`.  
```python
transit_flag = np.array([False, True, True])
```

Define the reference time of the orbital parameters, and integration time, based on the data:  
```python
t_epoch = 2455088.212 # reference time of the orbital parameters
t_start = 2454964.00 # start time, before the first data point
t_int   = 1977.00 # integration time (days) spanning all data
```
**BE CAREFUL** always double check these three times with your data. Most of the time some issues like `cored dumped` or very long running are due to mismatch between these times and the data.  


Let's run the orbital integration and get the simulated RVs and T0s (with orbital elements):  
```python
# **input:**
#
# - `t_start`: start time of the integration
# - `t_epoch`: reference time of the orbital parameters, it could be equal to, earlier, or later than `t_start`
# - `t_int`: integration time in days from t_start, remember that `t_start + t_int` has to cover all your observations!
# - `mass`: masses of the bodies (star, planet 1, planet 2, etc) in Msun
# - `radius`: radii of the bodies (star, planet 1, planet 2, etc) in Rsun
# - `period`: periods of the bodies (0, planet 1, planet 2, etc) in days
# - `ecc`: eccentricities of the bodies (0, planet 1, planet 2, etc)
# - `argp`: argument of pericenter/periastron of the bodies (0, planet 1, planet 2, etc) in deg, if ecc=0 set it to 90°
# - `meana`: mean anamolies of the bodies (0, planet 1, planet 2, etc) in deg
# - `inc`: inclination of the bodies (0, planet 1, planet 2, etc) in deg, remember inc=90° means the planet pass exact at the center of the star
# - `longn`: longitude of the ascending nodes of the bodies (0, planet 1, planet 2, etc) in deg, set for a reference body to 180°, if unknown set 180° for all the planets
# - `transit_flag`: flag of the transiting bodies, suggested default: `[False, True, True, ...]`
#
# **output:**
#
# - `rv_sim(n_rv)` in m/s
# - `body_tra_flag_sim(n_T0s)` id of the body, from 2 to n_body, of each transit time
# - `epo_sim(n_T0s)` epoch number of the transit w.r.t. the linear ephemeris of the corresponding body
# - `transits_sim(n_T0s)` simulated transit time in days corresponding to epoch and body at the same row
# - `durations_sim(n_T0s)` simulated Total Duration as T4 - T1 in minutes
# - `lambda_rm_sim` simulated spin-orbit misaligned at each transit time in degrees
# - `kep_elem_sim(n_T0s, 8)` keplerian elements for each transit time, each column an orbital elements:
#     - `kep_elem_sim(, 0)` period in days
#     - `kep_elem_sim(, 1)` semi-major axis in au
#     - `kep_elem_sim(, 2)` eccentricity
#     - `kep_elem_sim(, 3)` inclination in degrees
#     - `kep_elem_sim(, 4)` mean anomaly in degrees
#     - `kep_elem_sim(, 5)` argument of pericenter in degrees
#     - `kep_elem_sim(, 6)` true anomaly in degrees
#     - `kep_elem_sim(, 7)` longitude of ascending node in degrees
# - `stable` simulation stable is True/1 else is False/0

(
    rv_sim,
    body_tra_flag_sim,
    epo_sim,
    transits_sim,
    durations_sim,
    lambda_rm_sim,
    kep_elem_sim,
    stable,
) = pytrades.kelements_to_rv_and_t0s(
    t_start,
    t_epoch,
    t_int,
    mass,
    radius,
    period,
    ecc,
    argp,
    meana,
    inc,
    longn,
    transit_flag,
)
```

`body_tra_flag_sim` is a vector of integers with values corresponding to the bodies in `transit_flag` (`Fortran` indexing, star is 1, b is 2, and so on).  
To get only the simulate transits of planet b:
```python
body_id = 2 # planet b
epo_b_sim = epo_sim[body_tra_flag_sim==body_id] # ordered epochs as in the observed data
t0s_b_sim = transits_sim[body_tra_flag_sim==body_id] # each transit corresponds to each epoch and observe transit
# you can do the same for the durations etc.
```


Testing only RVs output:

```python
rv_sim = pytrades.kelements_to_rv(
    t_start,
    t_epoch,
    t_int,
    mass,
    radius,
    period,
    ecc,
    argp,
    meana,
    inc,
    longn,
)
```

Or only T0s output:

```python
(
    body_tra_flag_sim,
    epo_sim,
    transits_sim,
    durations_sim,
    lambda_rm_sim,
    kep_elem_sim,
    stable,
) = pytrades.kelements_to_t0s(
    t_start,
    t_epoch,
    t_int,
    mass,
    radius,
    period,
    ecc,
    argp,
    meana,
    inc,
    longn,
    transit_flag,
)

```

If you want to set-up a `log-Likelihood` function, I suggest to define:  
1. a function that converts fitted parameters `fit_pars` into physical parameters;
2. a function that computes the log-likelihood with `fit_pars` as input.

```python
def fit_to_physical(fit_pars):

    # here convert the fit_pars into:
    # `mass` of all bodies (star, b, c, ...) in Msun,
    # `radius` in Rsun,
    # `period` in day,
    # `ecc`, # if you have circular orbits set ecc = np.zeros((n_body))
    # `argp` in deg,
    # `meana` in deg,
    # `inc` in deg, # if fixed inclinations: inc = inc_deg
    # `longn` in deg,
    # use GLOBAL!

    return mass, radius, period, ecc, argp, meana, inc, longn


def loglike_function(fit_pars):

    mass, radius, period, ecc, argp, meana, inc, longn = fit_to_physical(fit_pars)

    (
        rv_sim,
        body_tra_flag_sim,
        epo_sim,
        transits_sim,
        durations_sim,
        lambda_rm_sim,
        kep_elem_sim,
        stable,
    ) = pytrades.kelements_to_rv_and_t0s(
        t_start,
        t_epoch,
        t_int,
        mass,
        radius,
        period,
        ecc,
        argp,
        meana,
        inc,
        longn,
        transit_flag,  # GLOBAL
    )

    ## RV
    lnL_rv = 0
    # add RVgamma, trends, gp to rv_sim
    rv = rv_sim + rv_gamma + trend_rv + gp_rv
    res_rv = rv_obs - rv  # rv_obs GLOBAL
    # lnL_rv = f(res_rv, jitter, etc)

    ## T0s
    lnL_T0 = 0
    res_T0b = t0_b - transits_sim[body_tra_flag_sim == 2]  # t0_b GLOBAL
    res_T0c = t0_c - transits_sim[body_tra_flag_sim == 3]  # t0_c GLOBAL
    # lnL_T0 = g(res_T0b, res_T0c)

    lnL = lnL_rv + lnL_T0

    return lnL
```

### Forward modelling for full RVs and T0s

Computes orbits:  
```python
time_steps, orbits, stable = pytrades.kelements_to_orbits_full(
    t_epoch,
    t_start,
    t_int,
    mass,
    radius,
    period,
    ecc,
    argp,
    meana,
    inc,
    longn,
    specific_times=None,  # add additional times to compute orbits, e.g. time of RV observations
    step_size=None,  # define the output stepsize
    n_steps_smaller_orbits=10,  # number of output steps of the inner planet, Default: 10
)
# sort them
n_steps = len(time_steps)
sort_time = np.argsort(time_steps)
initial_idx = np.arange(0, n_steps, 1).astype(int)[sort_time] # if integration forward and backward, the initial index will be present twice, in this way you could remove/skip it
time_steps = time_steps[sort_time]
orbits = orbits[sort_time, :] # orbits in astrocentric cartesian coordinates: cols (x, y, z, vx, vy, vz) x n_body, row for each step
```

Let's plot the orbits:  
```python
figsize = (4, 4)
body_names = ["star", "b", "c"]

fig = pytrades.base_plot_orbits(
    time_steps,
    orbits,
    radius,
    n_body,
    body_names,
    figsize=figsize,
    sky_scale="star",
    side_scale="star",
    title=None,
    show_plot=True,
)
plt.close(fig)
```

You can have the barycentric orbits and the barycenter:  
```python
bary_orbits, barycentre = pytrades.astrocentric_to_barycentric_orbits(mass, orbits)
```

Get the RV in m/s for each time step:  
```python
rvs = pytrades.orbits_to_rvs(mass, orbits)
```

Get all possible Transit Times with associated durations, mis-alignments, and Keplerian elements:
```python
transiting_body = 1  # all planets, use fortran indexing, so the first planet has 2, the second 3, and so on
n_transits = (t_int / period[1:]).astype(int)
n_all_transits = np.sum(n_transits) + (n_body - 1)  # star has no transits by definition
# you could set `n_all_transits = len(time_steps)*(n_body-1)` but it would use more memory.

transits, durations, lambda_rm, kep_elem, body_flag = pytrades.orbits_to_transits(
    n_all_transits, time_steps, mass, radius, orbits, transiting_body
)

# let's put info in a dictionary
planets_transits = {}

for i, pl_name in enumerate(body_names[1:]): # start from the first planet
    pl_num = i+1 # first planet will have i = 1, pl_num = 2
    sel = body_flag == pl_num
    planets_transits[pl_name] = {
      "planet_num": pl_num,
      "transits": transits[sel],
      "durations": durations[sel],
      "lambda_rm": lambda_rm[sel],
      "kep_elem": kep_elem[sel],
    }
```


### From script and folder with `yml` configurations

To use as it was _initially_ intended (in `python`) you can:  

1. copy a folder in [trades_example](trades_example) based on the number of planets (i.e. 2p for 2, 3 for 3p, and so on). The `base_2p_grid` is a version for the old `Fortran` grid version.  
2. copy inside the folder, if not present, the [configuration.yml](trades_example/configuration.yml) and adapt the `run` section, change `PyDE` and `emcee` sections etc (`ultranest` and `dynesty` in development, not fully tested).  
3. from terminal type:  
  > `python /path/to/pytrades/trades_emcee.py --input configuration.yml`  
4. after it finished you have to modify the `analysis`, `OC` and `RV` sections of the `configuration.yml` file and you will have the analysis results with  
  > `python /path/to/pytrades/trades_emcee_analysis.py --input configuration.yml`  

You can look at the `trades_emcee.py` and `trades_emcee_analysis.py` to learn how to use the different part of `pytrades` and there are a few other notebooks to plot `PyDE` run, `emcee` chains, RV and OCs.  

---
**`TRADES` v2.21.0 by Luca Borsato - 2016-2024**  
It is now possible to add flag to distinguish different telescopes for each Transit times.  

`TRADES` v2.20.0 by Luca Borsato - 2016-2023  

`TRADES` v2.19.0 by Luca Borsato - 2016-2023  

Long README and description in [README_long](README_long.md).

