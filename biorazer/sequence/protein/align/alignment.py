from biotite.sequence.align import Alignment


class AlignmentHelper:

    @staticmethod
    def get(alignment: Alignment, key):
        return Alignment(alignment.sequences[key], alignment.trace[:, key])

    @staticmethod
    def concat_alignments(alignments: list[Alignment]):
        """
        Merge multiple alignments into a single alignment.

        The alignments are merged by concatenating the sequences
        in the order they are provided.

        Parameters
        ----------
        alignments : list of Alignment
            A list of `Alignment` objects to be merged.

        Returns
        -------
        Alignment
            A new `Alignment` object containing the merged sequences.
        """
        if not alignments:
            raise ValueError("The list of alignments is empty.")

        # Ensure all alignments have traces of the same length
        trace_length = alignments[0].trace.shape[0]
        for aln in alignments[1:]:
            if aln.trace.shape[0] != trace_length:
                raise ValueError("All alignments must have the same trace length.")

        trace_concat = np.concatenate([aln.trace for aln in alignments], axis=1)
        sequences_concat = []
        for aln in alignments:
            sequences_concat.extend(aln.sequences)
        alignment_concat = Alignment(sequences_concat, trace_concat)

        return alignment_concat
