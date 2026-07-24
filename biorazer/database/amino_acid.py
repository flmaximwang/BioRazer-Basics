import numpy as np
from biotite.sequence import LetterAlphabet

# --- from aa_list.py ---
AMINO_ACIDS_1LETTER = tuple("ACDEFGHIKLMNPQRSTVWY")
AMINO_ACIDS_3LETTER = (
    "ALA",
    "CYS",
    "ASP",
    "GLU",
    "PHE",
    "GLY",
    "HIS",
    "ILE",
    "LYS",
    "LEU",
    "MET",
    "ASN",
    "PRO",
    "GLN",
    "ARG",
    "SER",
    "THR",
    "VAL",
    "TRP",
    "TYR",
)

# --- from aa_types.py ---
TYPES2AA = {
    "hydrophobic": "AFILMPVWY",
    "polar": "CDEGHKNRST",
    "charged": "DEHKR",
    "positive": "HKR",
    "negative": "DE",
    "aromatic": "FHWY",
    "aliphatic": "AILMPV",
    "small": "ACDGNPSTV",
    "large": "EFHIKLRWY",
    "all": "ACDEFGHIKLMNPQRSTVWY",
}

# --- from protein.py ---
COMMON_AA_NAME1 = "ACDEFGHIKLMNPQRSTVWY"
COMMON_AA_ALPHABET = LetterAlphabet(list(COMMON_AA_NAME1))
HYDROPHOBIC_AA_NAME1 = "AILMFWYV"
HYDROPHOBIC_AA_ALPHABET = LetterAlphabet(list(HYDROPHOBIC_AA_NAME1))
POLAR_AA_NAME1 = "NQSTC"
POLAR_AA_ALPHABET = LetterAlphabet(list(POLAR_AA_NAME1))
POSITIVE_AA_NAME1 = "KRH"
POSITIVE_AA_ALPHABET = LetterAlphabet(list(POSITIVE_AA_NAME1))
NEGATIVE_AA_NAME1 = "DE"
NEGATIVE_AA_ALPHABET = LetterAlphabet(list(NEGATIVE_AA_NAME1))
APOLAR_AA_NAME1 = "AILMFWYVGP"
APOLAR_AA_ALPHABET = LetterAlphabet(list(APOLAR_AA_NAME1))


def sequences_to_symbols(sequences: list[str], AA_seq: str) -> str:
    symbols = np.zeros((len(sequences[0]), len(AA_seq)), dtype=int)
    for seq in sequences:
        for i, aa in enumerate(seq):
            symbols[i, AA_seq.index(aa)] += 1
    return symbols
