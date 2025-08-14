import biotite.sequence as bio_seq
import biotite.sequence.io.fasta as fasta
import biotite.sequence.align as bio_seq_align
from biotite.application.clustalo import ClustalOmegaApp as clustalo
import biotite.sequence.graphics as bio_seq_graphics
import matplotlib.pyplot as plt

def read_multiple_fasta_files(
    fasta_paths: list[str] = [],
    names: list[str] = ["all"],
    labels: list[str] = []
):
    '''
    This helper function reads multiple fasta files and returns the sequences as a dictionary.
    
    Args:
    - fasta_paths: list of paths to fasta files.
        - When only one fasta file is provided, the function will read all sequences in the file.
        - When multiple fasta files are provided, the function will read the sequences in the same order as the fasta_paths list.
    '''
    assert isinstance(fasta_paths, list), "fasta_paths must be a list"
    if len(fasta_paths) ==0:
        raise ValueError("fasta_paths must contain at least one fasta file")
    else:
        sequence_dict = {}
        for name in fasta_paths:
            sequence_dict.update(fasta.get_sequences(fasta.FastaFile.read(name)))
        if "all" in names:
            labels_temp = list(sequence_dict.keys())
            sequences = list(sequence_dict.values())
            if len(labels) == 0:
                labels = labels_temp
            else:
                if len(labels) != len(labels_temp):
                    raise ValueError("labels must have the same length as the number of sequences")
        else:
            sequences = []
            for name in names:
                sequences.append(sequence_dict[name])
            if len(labels) == 0:
                for name in names:
                    labels.append(name)
    return {
        "labels": labels,
        "sequences": sequences
    }

class Alignment:
    
    @staticmethod
    def from_fasta_files(
        fasta_paths: list[str] = [],
        names: list[str] = ["all"],
        labels: list[str] = [],
        title: str = "Default Title",
        plot_options={},
        align_options={},
    ):
        '''
        This helper function reads multiple fasta files and returns the sequences as a dictionary.
        
        Args:
        - fasta_paths: list of paths to fasta files.
            - When only one fasta file is provided, the function will read all sequences in the file.
            - When multiple fasta files are provided, the function will read the sequences in the same order as the fasta_paths list.
        '''
        assert isinstance(fasta_paths, list), "fasta_paths must be a list"
        if len(fasta_paths) ==0:
            raise ValueError("fasta_paths must contain at least one fasta file")
        else:
            sequence_dict = {}
            for name in fasta_paths:
                sequence_dict.update(fasta.get_sequences(fasta.FastaFile.read(name)))
            if "all" in names:
                labels_temp = list(sequence_dict.keys())
                sequences = list(sequence_dict.values())
                if len(labels) == 0:
                    labels = labels_temp
                else:
                    if len(labels) != len(labels_temp):
                        raise ValueError("labels must have the same length as the number of sequences")
            else:
                sequences = []
                for name in names:
                    sequences.append(sequence_dict[name])
                if len(labels) == 0:
                    for name in names:
                        labels.append(name)
        
        return Alignment(sequences, labels, title=title, plot_options=plot_options, align_options=align_options)
    
    @staticmethod
    def from_strings(
        sequences: list[str],
        labels: list[str],
        title: str = "Default Title",
        method="clustalo",
        plot_options={},
        align_options={}
    ):
        '''
        Options used for plot_alignment_similarity_based is passed by plot_options.
        '''
        protein_sequences = []
        for sequence in sequences:
            protein_sequences.append(bio_seq.NucleotideSequence(sequence))
        return Alignment(protein_sequences, labels, title=title, method=method, plot_options=plot_options, align_options=align_options)
    
    def __init__(self,
            sequences: list[bio_seq.NucleotideSequence],
            labels: list[str],
            title: str = "Default Title",
            method="clustalo",
            align_options={},
            plot_options={},
        ):
        '''
        Options used for plot_alignment_similarity_based is passed by plot_options.
        
        (function) def plot_alignment_similarity_based(
            axes: Any,
            alignment: Any,
            symbols_per_line: int = 50,
            show_numbers: bool = False,
            number_size: Any | None = None,
            number_functions: Any | None = None,
            labels: Any | None = None,
            label_size: Any | None = None,
            show_line_position: bool = False,
            spacing: int = 1,
            color: Any | None = None,
            cmap: Any | None = None,
            matrix: Any | None = None,
            color_symbols: bool = False,
            symbol_spacing: Any | None = None,
            symbol_size: Any | None = None,
            symbol_param: Any | None = None
        ) -> None:
        '''
        self.sequences = sequences
        self.labels = labels
        self.title = title
        self.multiple_method = method
        self.sequence_slice = None
        self.figsize = (10, 10)
        
        self.align_options = {
            "gap_penalty": -10,
            "terminal_penalty": True,
            "local": False,
            "max_number": 1000,
        }
        for key in align_options:
            self.align_options[key] = align_options[key]
        
        self.plot_options = plot_options
        
        if len(self.sequences) < 2:
            raise ValueError("Alignment requires at least two sequences")
        if len(self.sequences) != len(self.labels):
            raise ValueError("Number of sequences must match the number of labels")
        
        if len(self.sequences) == 2:
            matrix = bio_seq_align.SubstitutionMatrix.std_nucleotide_matrix()
            my_align_options = {}
            for i in ["gap_penalty", "terminal_penalty", "local", "max_number"]:
                my_align_options[i] = self.align_options[i]
            self.alignments = bio_seq_align.align_optimal(
                self.sequences[0], self.sequences[1], matrix, 
                **my_align_options
            )
            self.alignment = self.alignments[-1]
        else:
            
            if self.multiple_method == "clustalo":
                app = clustalo(self.sequences)
                app.start()
                app.join()
                self.alignment = app.get_alignment()
            else:
                matrix = bio_seq_align.SubstitutionMatrix.std_nucleotide_matrix()
                my_align_options = {}
                for i in ["gap_penalty", "terminal_penalty"]:
                    my_align_options[i] = self.align_options[i]
                self.alignment = bio_seq_align.align_multiple(
                    self.sequences, matrix,
                    **my_align_options
                )
    
    def set_slice(self, sequence_slice):
        self.sequence_slice = sequence_slice
        
    def set_method(self, method):
        '''
        Set method to clustalo to use clustalo for alignment.
        Otherwise, the sequences will be aligned using bio_seq_align.align_multiple.
        '''
        self.multiple_method = method
    
    def plot(self):
        fig, ax = plt.subplots(figsize=self.figsize)
        if not self.sequence_slice:
            bio_seq_graphics.plot_alignment_similarity_based(
                ax, self.alignment, labels=self.labels, show_numbers=True, **self.plot_options
            )
        else:
            bio_seq_graphics.plot_alignment_similarity_based(
                ax, self.alignment[self.sequence_slice], labels=self.labels, show_numbers=True, **self.plot_options
            )
        ax.set_title(self.title)
        return fig, ax

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