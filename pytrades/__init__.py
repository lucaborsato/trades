#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .__version__ import __version__

# Import specific modules and assign them to the package namespace
from . import constants
from . import ancillary
from . import pytrades
from . import plot_oc
from . import plot_rv
from . import gelman_rubin
from . import geweke
from . import convergence
from . import chains_summary_plot
from . import de
from . import pyde
from . import evolution_pso_plot
from . import fitted_correlation_plot
from . import physical_correlation_plot
from . import gls
from . import plot_oc_full


# Optionally, you can import specific functions or classes from these modules
# For example, if ancillary.py has a function called `some_function`:
# from .ancillary import some_function

# If constants.py has a constant called `SOME_CONSTANT`:
# from .constants import SOME_CONSTANT

# If pytrades.py has a class called `Trade`:
# from .pytrades import Trade

# Expose them in the package namespace if needed
# __all__ = ['anc', 'cst', 'trades', 'some_function', 'SOME_CONSTANT', 'Trade']
		