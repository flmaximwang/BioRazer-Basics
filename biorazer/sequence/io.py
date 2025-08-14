def sequence_to_fasta(sequence_dict, fasta_file):
    """
    This function writes a dictionary of sequences to a FASTA file.

    Parameters:
    sequence_dict (dict): A dictionary with sequence IDs as keys and sequences as values.
    fasta_file (str): The path to the output FASTA file.
    """
    with open(fasta_file, "w") as f:
        for seq_id, sequence in sequence_dict.items():
            f.write(f">{seq_id}\n")
            f.write(f"{sequence}\n")
