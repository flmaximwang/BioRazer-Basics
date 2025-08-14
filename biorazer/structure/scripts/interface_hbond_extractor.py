import pyrosetta

pyrosetta_option_list = [
    "-mute core basic",
    "-ignore_unrecognized_res",
    "ignore_zero_occupancy",
    "-holes:dalphaball ~/Applications/bin/DAlphaBall.gcc"
]
pyrosetta.init(
    " ".join(pyrosetta_option_list)
)

def 