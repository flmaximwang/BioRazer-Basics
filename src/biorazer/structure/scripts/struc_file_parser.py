from biotite.structure.io import pdb, pdbx
import biotite.structure as bio_struc
import biotite.sequence as bio_seq

def read_struc_seq(struc_file, protein_chains=[]):
    '''
    This function reads a structure file and returns a dictionary with the chain IDs as keys and the sequences as values.
    '''
    if struc_file.endswith(".pdb"):
        my_atom_array = pdb.get_structure(pdb.PDBFile.read(struc_file))[0]
    elif struc_file.endswith(".cif"):
        my_atom_array = pdbx.get_structure(pdbx.CIFFile.read(struc_file))[0]
    res = {}
    if not protein_chains:
        protein_chains = list(bio_struc.get_chains(my_atom_array))
    for i in protein_chains:
        seq_i = "".join(
            list(
                map(
                    lambda x: bio_seq.ProteinSequence.convert_letter_3to1(x),
                    bio_struc.get_residues(my_atom_array[my_atom_array.chain_id == i])[1]
                )
            )
        )
        res[i] = seq_i
    return res