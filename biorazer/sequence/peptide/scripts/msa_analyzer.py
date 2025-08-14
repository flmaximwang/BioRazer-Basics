from Bio import SeqIO

LEGAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

def summarize_msa_aa_count_per_position(msa_list):
    '''
    Return: a list of dictionaries, where each dictionary contains the count of each amino acid at that position
    '''
    seq_len = len(msa_list[0])
    for i in range(1, len(msa_list)):
        if len(msa_list[i]) != seq_len:
            raise ValueError(f"Sequence {i} length {len(msa_list[i])} does not match first sequence length {seq_len}.")
    aa_counts_per_position = []
    for i in range(seq_len):
        aa_counts_per_position.append({})
        for aa in LEGAL_AA:
            aa_counts_per_position[i][aa] = 0
        for j in range(len(msa_list)):
            aa = msa_list[j][i].upper()
            if aa in LEGAL_AA:
                aa_counts_per_position[i][aa] += 1
            else:
                continue
    
    return aa_counts_per_position

def aa_count_to_aa_freq(aa_counts_per_position):
    aa_freq_per_position = []
    for i in range(len(aa_counts_per_position)):
        total_count = sum(aa_counts_per_position[i].values())
        aa_freq_per_position.append({aa: count / total_count for aa, count in aa_counts_per_position[i].items()})
    return aa_freq_per_position

def find_most_frequent_aas(aa_counts_per_position):
    '''
    Return: a list of tuples (amino_acid, frequency) for the most frequent amino acid at each position
    '''
    aa_freq_per_position = aa_count_to_aa_freq(aa_counts_per_position)
    most_frequent_aa = []
    
    for aa_freq in aa_freq_per_position:
        max_freq = 0
        max_aa = None
        for aa, freq in aa_freq.items():
            if freq > max_freq:
                max_freq = freq
                max_aa = aa
        most_frequent_aa.append((max_aa, max_freq))
        
    return most_frequent_aa
