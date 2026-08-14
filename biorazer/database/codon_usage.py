"""Codon-usage frequency tables, backed by the python-codon-tables package.

Exposes the species codon-frequency tables used e.g. by dnachisel's
CodonOptimize. python-codon-tables stores them per amino acid
({aa: {codon: freq}}); ``get_codon_usage_table`` flattens that to
{codon: freq} for direct use, ``get_codon_usage_table_by_aa`` returns the
native view. Table names follow "<species>_<taxid>" (e.g. e_coli_316407);
short species names like "e_coli" are resolved by prefix match.
"""

from python_codon_tables import available_codon_tables_names, get_codons_table

AVAILABLE_CODON_TABLES = tuple(available_codon_tables_names)


def _resolve_table_name(species: str) -> str:
    """Accept a full table name ('e_coli_316407') or a short species name
    ('e_coli'); the short form must match exactly one table prefix."""
    if species in AVAILABLE_CODON_TABLES:
        return species
    matches = [name for name in AVAILABLE_CODON_TABLES
               if name.startswith(species + "_")]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(
        f"unknown codon-usage species {species!r}; available: "
        f"{', '.join(sorted(AVAILABLE_CODON_TABLES))}"
    )


def get_codon_usage_table(species: str) -> dict:
    """Return the flattened {codon: frequency} usage table for a species."""
    by_aa = get_codons_table(_resolve_table_name(species))
    flat = {}
    for codons in by_aa.values():
        flat.update(codons)
    return flat


def get_codon_usage_table_by_aa(species: str) -> dict:
    """Return the native python-codon-tables {aa: {codon: freq}} view."""
    return get_codons_table(_resolve_table_name(species))