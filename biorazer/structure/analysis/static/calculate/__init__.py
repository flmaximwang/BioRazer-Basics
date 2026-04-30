"""
This package contains functions that return a scalar value for a given structure, such as
= the number of unsatisfied hydrogen bonds or
- the solvent accessible surface area.

These functions are used in the ``calculate`` module to calculate features for a structure.
"""

from .array import sasa_array, buried_unsat_hbond, hbond
from .scaler import sasa_value
