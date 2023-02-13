#!/usr/bin/env python
# coding: utf-8

# # TRADES FROM PYTHON

# In[ ]:


import numpy as np  # array
import os
import sys


# Add `pytrades` folder

# In[ ]:


trades_path = os.path.abspath("/path/to/trades/")
pytrades_path = os.path.join(trades_path, "pytrades")
sys.path.append(pytrades_path)
import constants as cst
from pytrades_lib import pytrades as f90trades


# ## integration parameters

# In[ ]:


t_start = 2454965.0 # start time of the integration
t_epoch = 2455088.212 # reference time of the orbital parameters
t_int   = 500.0 # total integration time from the start time

n_body = 3 # number of bodies (NOT PLANETS) in the system, that is star + planets
duration_check = 1 # to check or not the duration in this case it is a return flag, 1 means return it


# set the integration arguments:  
# - `n_body` **mandatory**  
# - `duration_check` **mandatory**  
# - `t_start` **optional**  
# - `t_epoch` **optional**  
# - `t_int` **optional**  

# In[ ]:


# f90trades.args_init(n_body, duration_check, t_start, t_epoch, t_int)
f90trades.args_init(n_body, duration_check)


# ## RV  
# provides radial velocities info: time, rv (m/s), and error_rv (m/s)

# In[ ]:


t_rv = np.array([
    2455342.947836,
    2455344.061419,
    2455351.037415,
    2455352.011289,
    2455367.083543,
    2455376.930532,
])
rv_obs = np.array([
    17.15,
    2.23,
    -7.89,
    -4.61,
    -25.45,
    17.39,
])
erv_obs = np.array([
    2.51,
    2.74,
    2.36,
    2.21,
    2.49,
    2.61,
])
# rv_data = np.column_stack((t_rv, rv_obs, erv_obs))
n_rv = len(t_rv)
print("Defined RVs with n_rv = ", n_rv)


# deallocation, this is a test, but the next fuction should deallocate RV data before re-allocate them

# In[ ]:


# f90trades.deallocate_rv_dataset()


# add the RV dataset to the common variable of TRADES

# In[ ]:


f90trades.set_rv_dataset(t_rv, rv_obs, erv_obs)


# ## TOs of planet b
# define transit times (T0s) of planet b, that will be the `2nd` body, where the `1st` body is the star!

# In[ ]:


epo_b = np.array([
    -5,
    -4,
    -2,
    -1,
    0,
    2,
    3,
    4,
    5,
])

t0_b = np.array([
    2454977.24875,
    2454996.48240,
    2455034.95437,
    2455054.19058,
    2455073.43412,
    2455111.92589,
    2455131.17167,
    2455150.42951,
    2455169.68103,
])

et0_b = np.array([
    0.00059,
    0.00062,
    0.00052,
    0.00052,
    0.00072,
    0.00050,
    0.00048,
    0.00048,
    0.00046,
])

# t0b_data = np.column_stack((epo_b, t0_b, et0_b))
n_t0b = len(t0_b)

body_id = 2 # fortran numbering starting with 1 == star


# add transits of planet b to the common variable of TRADES

# In[ ]:


f90trades.set_t0_dataset(body_id, epo_b, t0_b, et0_b)


# ## T0s planet c
# define T0s for planet c, that is `3rd` body

# In[ ]:


epo_c = np.array([
    -3,
    -2,
    -1,
    0,
    1,
    2,
])

t0_c = np.array([
    2454969.30577,
    2455008.33086,
    2455047.33560,
    2455086.31251,
    2455125.26284,
    2455164.18168,
])
et0_c = np.array([
    0.00070,
    0.00061,
    0.00058,
    0.00059,
    0.00053,
    0.00055,
])
# t0c_data = np.column_stack((epo_c, t0_c, et0_c))
n_t0c = len(t0_c)

body_id = 3 # fortran numbering starting with 1 == star


# add transits of planet c to the common variable of TRADES

# In[ ]:


f90trades.set_t0_dataset(body_id, epo_c, t0_c, et0_c)


# if you try to repeat the `set_t0_dataset` it will deallocate and reallocate the data ... hopefully :)

# In[ ]:


f90trades.set_t0_dataset(body_id, epo_c, t0_c, et0_c) # yes it works...


# ## define the orbital parameter of system at the `t_epoch`

# In[ ]:


# parameter = array([star, planet b, planet c]) = array([body_1, body_2, body_3])
M_msun   = np.array([1.0, 0.136985 * cst.Mjups, 0.094348 * cst.Mjups]) # Masses in Solar unit
R_rsun   = np.array([1.1, 0.842 * cst.Rjups, 0.823 * cst.Rjups]) # Radii in Solar unit
P_day    = np.array([0.0, 19.238756, 38.986097]) # Periods in days
ecc      = np.array([0.0, 0.058, 0.068]) # eccentricities
argp_deg = np.array([0.0, 356.1, 167.6]) # argument of pericenters in degrees
mA_deg   = np.array([0.0, 3.8, 307.4]) # mean anonalies in degrees
inc_deg  = np.array([0.0, 88.55, 89.12]) # inclinations in degrees
lN_deg   = np.array([0.0, 0.0, 0.0])# longitude of ascending nodes in degrees


# it is needed to define a flag vector to identify which bodies has to transit or not, the `star` does not, so first element has to be flagged to `False`.  
# If you don't know if they transit or not set `star` to `False`, all `other bodies` to `True`.

# In[ ]:


transit_flag = np.array([False, True, True])


# It necessary to prepare output arrays

# In[ ]:


# RV sim dataset
rv_sim = np.zeros((n_rv))
# T0s
n_T0s = n_t0b + n_t0c
n_kep = 8 # number of keplerian elements in output for each T0s
body_T0_sim, epo_sim = np.zeros((n_T0s)).astype(int), np.zeros((n_T0s)).astype(int)
t0_sim, t14_sim = np.zeros((n_T0s)), np.zeros((n_T0s))
kel_sim = np.zeros((n_T0s, n_kep))


# ### OUTPUT:
# - `rv_sim(n_rv)` in m/s  
# - `body_T0_sim(n_T0s)` id of the body, from 2 to n_body, of each transit time    
# - `epo_sim(n_T0s)` epoch number of the transit w.r.t. the linear ephemeris of the corresponding body  
# - `t0_sim(n_T0s)` simulated transit time in days corresponding to epoch and body at the same row  
# - `t14_sim(n_T0s)` simulated Total Duration as T4 - T1 in minutes  
# - `kel_sim(n_T0s, 8)` keplerian elements for each transit time, each column an orbital elements:  
#     - `kel_sim(, 0)` period in days  
#     - `kel_sim(, 1)` semi-major axis in au  
#     - `kel_sim(, 2)` eccentricity  
#     - `kel_sim(, 3)` inclination in degrees  
#     - `kel_sim(, 4)` mean anomaly in degrees  
#     - `kel_sim(, 5)` argument of pericenter in degrees  
#     - `kel_sim(, 6)` true anomaly in degrees  
#     - `kel_sim(, 7)` longitude of ascending node in degrees  

# let's run the orbital integration and get the simulated RVs and T0s (with orbital elements)

# In[ ]:


rv_sim, body_T0_sim, epo_sim, t0_sim, t14_sim, kel_sim = f90trades.kelements_to_rv_and_t0s(
    t_start,
    t_epoch,
    t_int,
    M_msun,
    R_rsun,
    P_day,
    ecc,
    argp_deg,
    mA_deg,
    inc_deg,
    lN_deg,
    transit_flag,
    n_rv,
    n_T0s,
    n_kep
)


# In[ ]:


print("time RVobs RVsim")
for to, rvo, rvs in zip(t_rv, rv_obs, rv_sim):
    print(to, rvo, rvs)


# In[ ]:


print("body epoch T0 T14")
for ibd, epo, t0s, t14 in zip(body_T0_sim, epo_sim, t0_sim, t14_sim):
    print(ibd, epo, t0s, t14)


# ### try to deallocate a T0s dataset: planet c

# In[ ]:


f90trades.deallocate_t0_dataset(3)


# so it is needed to re-define the output variables

# In[ ]:


n_T0s = n_t0b 
n_kep = 8
body_T0_sim, epo_sim = np.zeros((n_T0s)).astype(int), np.zeros((n_T0s)).astype(int)
t0_sim, t14_sim = np.zeros((n_T0s)), np.zeros((n_T0s))
kel_sim = np.zeros((n_T0s, n_kep))


# In[ ]:


rv_sim, body_T0_sim, epo_sim, t0_sim, t14_sim, kel_sim = f90trades.kelements_to_rv_and_t0s(
    t_start,
    t_epoch,
    t_int,
    M_msun,
    R_rsun,
    P_day,
    ecc,
    argp_deg,
    mA_deg,
    inc_deg,
    lN_deg,
    transit_flag,
    n_rv,
    n_T0s,
    n_kep
)


# In[ ]:


print("time RVobs RVsim")
for to, rvo, rvs in zip(t_rv, rv_obs, rv_sim):
    print(to, rvo, rvs)


# In[ ]:


print("body epoch T0 T14")
for ibd, epo, t0s, t14 in zip(body_T0_sim, epo_sim, t0_sim, t14_sim):
    print(ibd, epo, t0s, t14)


# ### testing only RVs output

# In[ ]:


rv_sim = np.zeros((n_rv))
rv_sim = f90trades.kelements_to_rv(
    t_start,
    t_epoch,
    t_int,
    M_msun,
    R_rsun,
    P_day,
    ecc,
    argp_deg,
    mA_deg,
    inc_deg,
    lN_deg,
    n_rv,
)


# In[ ]:


print("time RVobs RVsim")
for to, rvo, rvs in zip(t_rv, rv_obs, rv_sim):
    print(to, rvo, rvs)


# ### testing only T0s output  
# re-adding T0s of planet c

# In[ ]:


f90trades.set_t0_dataset(body_id, epo_c, t0_c, et0_c)


# resetting size and output variables

# In[ ]:


n_T0s = n_t0b + n_t0c
n_kep = 8
body_T0_sim, epo_sim = np.zeros((n_T0s)).astype(int), np.zeros((n_T0s)).astype(int)
t0_sim, t14_sim = np.zeros((n_T0s)), np.zeros((n_T0s))
kel_sim = np.zeros((n_T0s, n_kep))


# In[ ]:


body_T0_sim, epo_sim, t0_sim, t14_sim, kel_sim = f90trades.kelements_to_t0s(
    t_start,
    t_epoch,
    t_int,
    M_msun,
    R_rsun,
    P_day,
    ecc,
    argp_deg,
    mA_deg,
    inc_deg,
    lN_deg,
    transit_flag,
    n_T0s,
    n_kep
)


# In[ ]:


print("body epoch T0 T14")
for ibd, epo, t0s, t14 in zip(body_T0_sim, epo_sim, t0_sim, t14_sim):
    print(ibd, epo, t0s, t14)


# In[ ]:




