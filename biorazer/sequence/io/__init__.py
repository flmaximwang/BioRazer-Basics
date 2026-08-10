from .string import Fasta_Profile, Fasta_StrDict, StrDict_Fasta
from .protein import (
    Fasta_ProteinSequence,
    Fasta_ProteinSequenceDict,
    ProteinSequenceDict_Fasta,
)
from . import nucleotide

# 用户侧别名: 与既有脚本 (seqlogo.py / mutation_pair.py) 及文档一致
FASTA_PROFILE = Fasta_Profile
FASTA_DICT = Fasta_StrDict