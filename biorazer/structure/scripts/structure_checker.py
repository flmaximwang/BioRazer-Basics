import Bio.PDB as biop_PDB
from Bio.PDB.Structure import Structure as biop_PDB_Structure
import numpy as np

# Atomic radii for various atom types. 
# You can comment out the ones you don't care about or add new ones
atom_radii = {
#    "H": 1.20,  # Who cares about hydrogen??
    "C": 1.70, 
    "N": 1.55, 
    "O": 1.52,
    "S": 1.80,
    "F": 1.47, 
    "P": 1.80, 
    "CL": 1.75, 
    "MG": 1.73,
}

def count_clashes(structure: biop_PDB_Structure, clash_cutoff=0.63):
    # Set what we count as a clash for each pair of atoms
    clash_cutoffs = {i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j])) for i in atom_radii for j in atom_radii}

    # Extract atoms for which we have a radii
    atoms = [x for x in structure.get_atoms() if x.element in atom_radii]
    coords = np.array([a.coord for a in atoms], dtype="d")

    # Build a KDTree (speedy!!!)
    kdt = biop_PDB.kdtrees.KDTree(coords)

    # Initialize a list to hold clashes
    clashes = []

    # Iterate through all atoms
    for atom_1 in atoms:
        # Find atoms that could be clashing
        kdt_search = kdt.search(np.array(atom_1.coord, dtype="d"), max(clash_cutoffs.values()))

        # Get index and distance of potential clashes
        potential_clash = [(a.index, a.radius) for a in kdt_search]

        for ix, atom_distance in potential_clash:
            atom_2 = atoms[ix]

            # Exclude clashes from atoms in the same residue
            if atom_1.parent.id == atom_2.parent.id:
                continue

            # Exclude clashes from peptide bonds
            elif (atom_2.name == "C" and atom_1.name == "N") or (atom_2.name == "N" and atom_1.name == "C"):
                continue

            # Exclude clashes from disulphide bridges
            elif (atom_2.name == "SG" and atom_1.name == "SG") and atom_distance > 1.88:
                continue

            if atom_distance < clash_cutoffs[atom_2.element + "_" + atom_1.element]:
                clashes.append((atom_1, atom_2))

    return len(clashes) // 2