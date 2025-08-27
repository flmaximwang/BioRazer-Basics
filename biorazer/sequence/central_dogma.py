from .util.dictionary.codon import E_COLI_RT


def transcribe(dna: str):

    rna = ""
    return rna


def reverse_transcribe(rna: str):

    dna = ""
    return dna


def translate(nucleotide: str, codon_table: dict) -> str:
    """Translate a nucleotide sequence into an amino acid sequence using the provided codon table.

    Args:
        nucleotide (str): The nucleotide sequence to be translated.
        codon_table (dict): A dictionary mapping codons (3-letter nucleotide sequences) to amino acids.

    Returns:
        str: The resulting amino acid sequence.
    """
    amino_acid_sequence = []
    # Process the nucleotide sequence in chunks of 3 (codons)
    for i in range(0, len(nucleotide) - 2, 3):
        codon = nucleotide[i : i + 3]
        amino_acid = codon_table.get(codon, "X")  # Use 'X' for unknown codons
        amino_acid_sequence.append(amino_acid)
    return "".join(amino_acid_sequence)


def reverse_translate(protein: str, rt_table: dict) -> str:
    res = ""
    for aa in protein:
        res += rt_table[aa][0][0]
    return res
