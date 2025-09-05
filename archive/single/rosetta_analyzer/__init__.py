# This package requires to have Rosetta installed and 
# the environment variable ROSETTA_BIN set to the path of the Rosetta bin directory.

import os

ROSETTA_BIN = os.environ["ROSETTA_BIN"]
ROSETTA_PLATFORM = os.environ["ROSETTA_PLATFORM"]

def set_rosetta_platform(platform):
    ROSETTA_PLATFORM = platform