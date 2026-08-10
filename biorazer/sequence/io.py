import re
from collections.abc import Iterable

import numpy as np
from biotite.sequence import ProteinSequence
from biotite.sequence.io import fasta
from biotite.sequence.profile import SequenceProfile

from biorazer.io import Converter


class StrDict_Fasta(Converter):
    def write(self, tmp: dict):
        """
        Write dictionary of sequences to a FASTA file.

        Parameters
        ----------
        tmp: dict[str, str]
            A dictionary with sequence IDs as keys and sequences as values.
        """

        fasta_file = fasta.FastaFile()
        for key, value in tmp.items():
            fasta_file[key] = value
        fasta_file.write(self.output_io)


class Fasta_StrDict(Converter):

    def read(self):
        """
        Read sequences from a FASTA file into a dictionary.

        Returns
        -------
        dict[str, str]
            A dictionary with sequence IDs as keys and sequences as values.
        """

        fasta_file = fasta.FastaFile.read(self.input_io)
        seq_dict = {}
        for key in fasta_file.keys():
            seq_dict[key] = str(fasta_file[key])
        return seq_dict


class Fasta_Profile(Converter):
    """
    Read a FASTA file of equal-length sequences into per-chain
    SequenceProfile objects.

    Sequences may contain multiple chains separated by ``':'`` (e.g.
    LigandMPNN output). FASTA records are counted 0-based for
    ``ignore_seqs``; chains are counted 0-based for ``chains``.
    """

    @staticmethod
    def _parse_index_spec(indices, n, what, ordered=False):
        """
        Coerce a raw index argument into validated 0-based indices in 0..n-1.

        The raw argument may be ``None`` (select nothing), a single int
        (select one), or an iterable of ints (select several). It is
        returned as a set, or as a list preserving first-occurrence order
        with duplicates removed when ``ordered=True``.

        Parameters
        ----------
        indices : int, iterable of int, or None
            The requested 0-based indices; ``None`` means "none".
        n : int
            The number of items (FASTA records or chains); valid indices
            are 0..n-1.
        what : str
            Label of the indexed thing, used in error messages
            (e.g. "sequence", "chain").
        ordered : bool, optional
            If True, return a list that preserves first-occurrence order
            and removes duplicates; otherwise return a set.

        Returns
        -------
        set[int] or list[int]
            Empty when ``indices`` is None: ``[]`` for ``ordered=True``,
            otherwise ``set()``.

        Raises
        ------
        ValueError
            If any index is not an int or is outside 0..n-1.

        Examples
        --------
        >>> FASTA_PROFILE._parse_index_spec(2, 5, "sequence")
        {2}
        >>> FASTA_PROFILE._parse_index_spec([2, 0, 2], 5, "sequence")
        {0, 2}
        >>> FASTA_PROFILE._parse_index_spec([2, 0, 2], 5, "chain", ordered=True)
        [2, 0]
        >>> FASTA_PROFILE._parse_index_spec(None, 5, "chain", ordered=True)
        []
        """
        if indices is None:
            return [] if ordered else set()
        if isinstance(indices, int):
            indices = [indices]
        indices = list(indices)
        bad_type = [i for i in indices if not isinstance(i, int)]
        if bad_type:
            raise ValueError(f"{what} index must be an int, got {bad_type}")
        bad = [i for i in indices if not 0 <= i < n]
        if bad:
            raise ValueError(f"{what} index out of range 0..{n - 1}: {bad}")
        if ordered:
            seen, out = set(), []
            for i in indices:
                if i not in seen:
                    seen.add(i)
                    out.append(i)
            return out
        return set(indices)

    @staticmethod
    def _make_profile(chain_seqs: list[str]) -> SequenceProfile:
        """Per-column symbol counts -> biotite SequenceProfile."""
        alphabet = ProteinSequence.alphabet
        length = len(chain_seqs[0])
        codes = np.array(
            [[alphabet.encode(ch) for ch in s] for s in chain_seqs], dtype=int
        )
        symbols = np.zeros((length, len(alphabet)), dtype=int)
        for i in range(length):
            symbols[i] = np.bincount(codes[:, i], minlength=len(alphabet))
        return SequenceProfile(symbols, np.zeros(length, dtype=bool), alphabet)

    def read(
        self,
        ignore_seqs: int | Iterable[int] | None = None,
        chains: int | Iterable[int] | None = None,
    ) -> dict[int, SequenceProfile]:
        """
        Parameters
        ----------
        ignore_seqs : int or iterable of int, optional
            0-based indices of FASTA records to skip.
        chains : int or iterable of int, optional
            0-based indices of chains to keep, in the given order
            (default: all chains, in file order).

        Returns
        -------
        dict[int, SequenceProfile]
            Selected chain index (0-based) -> column-count profile
            (symbols shape (L, 24), gap mask all False).
        """
        raw = Fasta_StrDict(self.input_io, self.output_io).read()
        items = list(raw.items())
        n_seqs = len(items)
        if n_seqs == 0:
            raise ValueError(f"FASTA file contains no sequences: {self.input_io}")

        ignore = self._parse_index_spec(ignore_seqs, n_seqs, "sequence")
        seqs = {}
        for idx, (name, s) in enumerate(items):
            if idx in ignore:
                continue
            seqs[name] = re.sub(r"\s+", "", s).upper()
        if not seqs:
            raise ValueError("ignore_seqs skipped all sequences")

        # Split into chains, requiring a consistent chain count
        n_chain_sizes = {len(s.split(":")) for s in seqs.values()}
        if len(n_chain_sizes) > 1:
            raise ValueError(
                "inconsistent chain counts across sequences: "
                f"{sorted(n_chain_sizes)}"
            )
        n_chains = n_chain_sizes.pop()

        chain_seqs = [[] for _ in range(n_chains)]
        for name, s in seqs.items():
            parts = s.split(":")
            for c, part in enumerate(parts):
                if not part:
                    raise ValueError(f"sequence {name!r} has an empty chain")
                chain_seqs[c].append(part)

        # Validate per-chain length and alphabet membership
        symbols = set(ProteinSequence.alphabet.get_symbols())
        for c, parts in enumerate(chain_seqs):
            lens = {len(x) for x in parts}
            if len(lens) > 1:
                raise ValueError(
                    f"chain {c} has unequal sequence lengths: {sorted(lens)}"
                )
            bad = sorted({ch for x in parts for ch in x} - symbols)
            if bad:
                raise ValueError(
                    f"chain {c} has non-protein symbols: {bad}"
                )

        if chains is None:
            selected = list(range(n_chains))
        else:
            selected = self._parse_index_spec(
                chains, n_chains, "chain", ordered=True
            )

        return {c: self._make_profile(chain_seqs[c]) for c in selected}
