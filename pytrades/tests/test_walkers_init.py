import numpy as np
import os
import sys

# import `pytrades` and `ancillary` modules

pytrades_dir = os.path.abspath("../")
sys.path.append(pytrades_dir)

import pytrades
import ancillary as anc

# init `pytrades` with 2p-system

exosystem_folder = os.path.abspath("../../trades_example/base_2p/")
sim = pytrades.TRADES(
    exosystem_folder,
    sub_folder="",
    nthreads=1,
    seed=42,
    m_type="e"
)

for k, v in vars(sim).items():
    print(k, v)

sim.init_trades()


