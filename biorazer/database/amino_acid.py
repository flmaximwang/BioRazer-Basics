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

# The variants differ in the VALUE case only -- keys are always the
# standard upper-case form. 1TO3 has 3 value forms (MET/met/Met),
# 3TO1 has 2 (M/m). Use AMINO_ACIDS_1TO3_INITIAL_CAPITAL by default.
AMINO_ACIDS_1TO3_UPPER = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}

AMINO_ACIDS_1TO3_LOWER = {
    "A": "ala",
    "C": "cys",
    "D": "asp",
    "E": "glu",
    "F": "phe",
    "G": "gly",
    "H": "his",
    "I": "ile",
    "K": "lys",
    "L": "leu",
    "M": "met",
    "N": "asn",
    "P": "pro",
    "Q": "gln",
    "R": "arg",
    "S": "ser",
    "T": "thr",
    "V": "val",
    "W": "trp",
    "Y": "tyr",
}

AMINO_ACIDS_1TO3_INITIAL_CAPITAL = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
}

AMINO_ACIDS_3TO1_UPPER = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}

AMINO_ACIDS_3TO1_LOWER = {
    "ALA": "a",
    "CYS": "c",
    "ASP": "d",
    "GLU": "e",
    "PHE": "f",
    "GLY": "g",
    "HIS": "h",
    "ILE": "i",
    "LYS": "k",
    "LEU": "l",
    "MET": "m",
    "ASN": "n",
    "PRO": "p",
    "GLN": "q",
    "ARG": "r",
    "SER": "s",
    "THR": "t",
    "VAL": "v",
    "TRP": "w",
    "TYR": "y",
}

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
