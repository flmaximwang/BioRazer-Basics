# The batch analyzer subpackage helps you to analyze a batch of structure files.
# Lots of metrics are provided to help you understand your structures.
# - Multiple structures are defined as a project, and the summary of the project is saved as a .csv table.
# - Every structure in the project is recorded as an entry, with files used to analyze generated and stored as separate files to avoid recalculation.

from .abstracts import *
from . import utils