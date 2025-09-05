import numpy as np
from biotite.sequence import Alphabet, LetterAlphabet

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
