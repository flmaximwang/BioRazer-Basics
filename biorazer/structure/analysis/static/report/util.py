"""Shared helpers for structure report modules.

``_normalize_fmt`` lives here (instead of in a report package) so that
report modules in any package -- including the lower-level
``selection`` package -- can reuse it without inverting the package
dependency direction.
"""

from biorazer.structure.util.report import _normalize_fmt
