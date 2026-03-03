from ..util.dictionary.codon import *
from ..util.dictionary import codon as codon_dicts


def full_reverse_translate(
    amino_acid_sequence: str,
    codon_dict: dict,
    nucleotide_sequences: list = [],
    prob_cutoff: float = 0.000001,
) -> str:
    """
    Perform a full reverse translation of an amino acid sequence into a nucleotide sequence
    using the provided codon dictionary.

    Parameters
    ----------
    amino_acid_sequence : str
        The amino acid sequence to be reverse translated.
    codon_dict : dict
        A dictionary mapping amino acids to their corresponding codons and frequencies.

    Returns
    -------
    dict
        All possible nucleotide sequences for the given amino acid sequence and their probabilities.
    """

    # Base case: If the amino acid sequence is empty, return the current nucleotide sequences
    if len(amino_acid_sequence) == 0:
        return nucleotide_sequences

    # Recursive case: Process the first amino acid and extend the nucleotide sequences
    current_aa = amino_acid_sequence[0]
    if current_aa not in codon_dict:
        raise ValueError(f"Amino acid '{current_aa}' not found in codon dictionary.")

    possible_codons = codon_dict[current_aa]
    new_nucleotide_sequences = []
    if nucleotide_sequences == []:
        for codon, freq in possible_codons:
            new_nucleotide_sequences.append((codon, freq))
    else:
        for nucleotide_seq, prob in nucleotide_sequences:
            for codon, freq in possible_codons:
                new_seq = nucleotide_seq + codon
                new_prob = prob * freq
                if new_prob >= prob_cutoff:
                    new_nucleotide_sequences.append((new_seq, new_prob))
    return full_reverse_translate(
        amino_acid_sequence[1:],
        codon_dict,
        new_nucleotide_sequences,
        prob_cutoff=prob_cutoff,
    )
