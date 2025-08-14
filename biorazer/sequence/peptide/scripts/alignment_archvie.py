
# def align_two_fasta_sequences(
#     seq1_path: str,
#     seq2_path: str,
#     label1: str,
#     label2: str,
#     title: str,
#     figsize=(10, 10),
#     sequence_slice=None
# ):
#     seq1 = fasta.get_sequence(fasta.FastaFile.read(seq1_path))
#     seq2 = fasta.get_sequence(fasta.FastaFile.read(seq2_path))
#     matrix = bio_seq_align.SubstitutionMatrix.std_protein_matrix()
#     alignment = bio_seq_align.align_optimal(seq1, seq2, matrix)[0]
#     fig, ax = plt.subplots(figsize=figsize)
#     if not sequence_slice:
#         bio_seq_graphics.plot_alignment_similarity_based(
#             ax, alignment, labels=[label1, label2], show_numbers=True
#         )
#     else:
#         bio_seq_graphics.plot_alignment_similarity_based(
#             ax, alignment[sequence_slice], labels=[label1, label2], show_numbers=True
#         )
#     ax.set_title(title)
#     return alignment



# def align_fasta_sequences(
#     fasta_paths: list[str] = [],
#     names: list[str] = ["all"],
#     labels: list[str] = [],
#     title: str = "Default Title",
#     figsize=(10, 10),
#     gap_penalty=-5,
#     terminal_penalty=True,
#     sequence_slice=None
# ):
#     '''
#     This helper function aligns multiple fasta sequences and plots the alignment.
    
#     Args:
#     - fasta_paths: list of paths to fasta files.
#         - When only one fasta file is provided, the function will align all sequences in the file.
#         - When multiple fasta files are provided, the function will align the sequences in the same order as the fasta_paths list.
#     '''
#     res = read_multiple_fasta_files(fasta_paths, names, labels)
#     sequences = res["sequences"]
#     labels = res["labels"]
                    
#     matrix = bio_seq_align.SubstitutionMatrix.std_protein_matrix()
#     alignment, order, tree, distance_matrix = bio_seq_align.align_multiple(sequences, matrix, gap_penalty=gap_penalty, terminal_penalty=terminal_penalty)
#     fig, ax = plt.subplots(figsize=figsize)
#     if not sequence_slice:
#         bio_seq_graphics.plot_alignment_similarity_based(
#             ax, alignment, labels=labels, show_numbers=True
#         )
#     else:
#         bio_seq_graphics.plot_alignment_similarity_based(
#             ax, alignment[sequence_slice], labels=labels, show_numbers=True
#         )
#     ax.set_title(title)
#     return alignment

# def clustalo_fasta_sequences(
#     fasta_paths: list[str] = [],
#     names: list[str] = ["all"],
#     labels: list[str] = [],
#     title: str = "Default Title",
#     figsize=(10, 10),
#     sequence_slice=None
# ):
#     res = read_multiple_fasta_files(fasta_paths, names, labels)
#     sequences = res["sequences"]
#     labels = res["labels"]

#     app = clustalo(sequences)
#     app.start()
#     app.join()
#     alignment = app.get_alignment()
#     fig, ax = plt.subplots(figsize=figsize)
#     if not sequence_slice:
#         bio_seq_graphics.plot_alignment_similarity_based(
#             ax, alignment, labels=labels, show_numbers=True
#         )
#     else:
#         bio_seq_graphics.plot_alignment_similarity_based(
#             ax, alignment[sequence_slice], labels=labels, show_numbers=True
#         )
#     ax.set_title(title)
#     return alignment