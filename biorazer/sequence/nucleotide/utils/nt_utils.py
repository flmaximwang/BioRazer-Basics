def get_mRNA_from_transcript(transcript_seq: str, exons: list[tuple[int, int]]):
    '''
    Joining exons from a transcript sequence to get the mRNA sequence.
    
    :param transcript_seq: The transcript sequence from which the mRNA sequence 
    will be extracted.
    
    :param exons: The exons of the transcript sequence. Every exon marked with a tuple, 
    where the first element is the start position and the second element is the end position.
    
    :return: The mRNA sequence.
    '''
    print(transcript_seq)
    print(exons)
    mRNA_seq = ""
    for start, end in exons:
        mRNA_seq += transcript_seq[start:end+1]
    return mRNA_seq