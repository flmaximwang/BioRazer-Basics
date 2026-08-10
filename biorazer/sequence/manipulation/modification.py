import re

from biotite.sequence import AlphabetError, ProteinSequence

_MUTATION_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")


def mutate(
    sequence: ProteinSequence | str,
    mutation_spec: str,
    first_res_id: int = 1,
) -> ProteinSequence:
    """
    Apply mutations to a protein sequence and return the mutated copy.

    Each mutation in ``mutation_spec`` has the form
    ``<wild_type><position><new_type>`` (e.g. ``A16T``), where
    ``wild_type`` is the 1-letter residue expected at the position and
    ``new_type`` is the 1-letter residue to substitute. The mutation is
    rejected if the residue at the position does not match
    ``wild_type``.

    Parameters
    ----------
    sequence : ProteinSequence or str
        The sequence to mutate. A str is converted to a
        ProteinSequence (alphabet membership is validated). The input
        is not modified.
    mutation_spec : str
        Comma-separated mutations, e.g. ``"A16T,A18L"``.
    first_res_id : int, optional
        The residue number of the first residue of ``sequence``, e.g.
        the PDB author numbering of the first residue of a fragment.
        Positions in ``mutation_spec`` are interpreted as absolute
        residue numbers: the residue at position ``p`` is at index
        ``p - first_res_id``. The default ``1`` gives standard 1-based
        numbering relative to the start of ``sequence``.

    Returns
    -------
    ProteinSequence
        A copy of ``sequence`` with all mutations applied. Always a
        ProteinSequence, even when ``sequence`` was given as a str.

    Raises
    ------
    TypeError
        If ``sequence`` is neither a ProteinSequence nor a str, or
        ``first_res_id`` is not an int.
    ValueError
        If ``sequence`` is a str containing symbols outside the protein
        alphabet, or a mutation is malformed, refers to a position
        outside the sequence, repeats a position, or the residue at the
        position does not match the wild-type letter in the spec.
    """
    if isinstance(sequence, str):
        try:
            sequence = ProteinSequence(sequence)
        except AlphabetError as exc:
            raise ValueError(
                f"sequence is not a valid protein sequence: {exc}"
            ) from exc
    if not isinstance(sequence, ProteinSequence):
        raise TypeError(
            f"sequence must be a ProteinSequence or str, got {type(sequence).__name__}"
        )
    if not isinstance(first_res_id, int):
        raise TypeError(
            f"first_res_id must be an int, got {type(first_res_id).__name__}"
        )
    if not mutation_spec.strip():
        raise ValueError("mutation_spec must not be empty")

    mutated = sequence.copy()
    applied = set()
    for token in mutation_spec.split(","):
        token = token.strip()
        match = _MUTATION_RE.match(token)
        if match is None:
            raise ValueError(
                f"malformed mutation {token!r}: expected "
                "<wild_type><position><new_type>, e.g. A16T"
            )
        wild_type, position_str, new_type = match.groups()
        wild_type = wild_type.upper()
        new_type = new_type.upper()
        position = int(position_str)

        index = position - first_res_id
        if not 0 <= index < len(sequence):
            raise ValueError(
                f"mutation {token!r} at position {position} is outside "
                f"the sequence of length {len(sequence)}"
            )
        if index in applied:
            raise ValueError(
                f"duplicate mutation position {position} in {mutation_spec!r}"
            )
        actual = sequence[index]
        if actual != wild_type:
            raise ValueError(
                f"mutation {token!r}: expected {wild_type} at position "
                f"{position}, found {actual}"
            )
        mutated[index] = new_type
        applied.add(index)
    return mutated
