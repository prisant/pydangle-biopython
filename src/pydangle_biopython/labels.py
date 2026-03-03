"""Residue classification labels for pydangle.

Provides per-residue categorical labels such as cis/trans peptide bond
classification, pre-proline detection, and Ramachandran category
assignment.  These are implemented as a new ``label`` function type
alongside the existing ``distance``, ``angle``, and ``dihedral`` types.

Label functions receive the residue list, the current residue index,
and return a string classification or the unknown marker.

Ramachandran categories follow the wwPDB Validation Task Force
conventions (Read et al., Structure 19:1395-1412, 2011):

    General    -- not Gly, Pro, Ile/Val, and not followed by Pro
    IleVal     -- Ile or Val (beta-branched)
    Gly        -- Glycine
    TransPro   -- Proline with |omega| >= 30 degrees
    CisPro     -- Proline with |omega| < 30 degrees
    PrePro     -- any residue (except Gly/Pro) followed by Proline
"""

from __future__ import annotations

import math
from typing import Any

from Bio.PDB.vectors import calc_dihedral

# ---------------------------------------------------------------------------
# Omega computation helper
# ---------------------------------------------------------------------------

def _compute_omega(residue_list: list[Any], index: int) -> float | None:
    """Compute the omega dihedral angle for residue *index* in degrees.

    Omega is defined as: i-1 CA, i-1 C, i N, i CA (IUPAC convention).
    Returns None if the required atoms are missing or index is 0.
    """
    if index < 1:
        return None
    prev = residue_list[index - 1]
    curr = residue_list[index]
    try:
        v_ca0 = prev['CA'].get_vector()
        v_c0 = prev['C'].get_vector()
        v_n1 = curr['N'].get_vector()
        v_ca1 = curr['CA'].get_vector()
    except (KeyError, TypeError):
        return None
    radians = calc_dihedral(  # type: ignore[no-untyped-call]
        v_ca0, v_c0, v_n1, v_ca1,
    )
    return float(math.degrees(radians))


# ---------------------------------------------------------------------------
# Primitive label functions
# ---------------------------------------------------------------------------

def label_is_cis(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if the peptide bond to this residue is cis.

    A peptide bond is cis when |omega| < 30 degrees.  Returns *unknown*
    for the first residue in a fragment or if omega cannot be computed.
    """
    omega = _compute_omega(residue_list, index)
    if omega is None:
        return unknown
    return str(abs(omega) < 30.0)


def label_is_trans(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if the peptide bond to this residue is trans.

    A peptide bond is trans when |omega| >= 30 degrees.  Returns *unknown*
    for the first residue in a fragment or if omega cannot be computed.
    """
    omega = _compute_omega(residue_list, index)
    if omega is None:
        return unknown
    return str(abs(omega) >= 30.0)


def label_is_gly(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if the residue is glycine."""
    return str(residue_list[index].get_resname() == 'GLY')


def label_is_pro(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if the residue is proline."""
    return str(residue_list[index].get_resname() == 'PRO')


def label_is_ileval(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if the residue is Ile or Val (beta-branched)."""
    return str(residue_list[index].get_resname() in ('ILE', 'VAL'))


def label_is_prepro(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if the next residue is proline.

    Returns *unknown* for the last residue in a fragment.
    """
    if index + 1 >= len(residue_list):
        return unknown
    return str(residue_list[index + 1].get_resname() == 'PRO')


# ---------------------------------------------------------------------------
# Composite label: Ramachandran category
# ---------------------------------------------------------------------------

def label_rama_category(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Assign the wwPDB Ramachandran category for the residue.

    Decision tree (Read et al. 2011, Richardson & Richardson):

    1. Gly -> ``'Gly'``
    2. Ile or Val -> ``'IleVal'``
    3. Pro -> measure omega -> ``'CisPro'`` (|omega| < 30) or ``'TransPro'``
    4. Next residue is Pro -> ``'PrePro'``
    5. Otherwise -> ``'General'``

    Returns *unknown* if the residue identity cannot be determined.
    """
    resname = residue_list[index].get_resname()

    # 1. Glycine
    if resname == 'GLY':
        return 'Gly'

    # 2. Ile / Val
    if resname in ('ILE', 'VAL'):
        return 'IleVal'

    # 3. Proline -- subdivide by cis/trans
    if resname == 'PRO':
        omega = _compute_omega(residue_list, index)
        if omega is None:
            # Can't determine cis/trans (e.g. first residue).
            # Default to TransPro since trans is overwhelmingly common.
            return 'TransPro'
        if abs(omega) < 30.0:
            return 'CisPro'
        return 'TransPro'

    # 4. Pre-proline
    if index + 1 < len(residue_list):
        if residue_list[index + 1].get_resname() == 'PRO':
            return 'PrePro'

    # 5. General
    return 'General'


# ---------------------------------------------------------------------------
# Label dispatch registry
# ---------------------------------------------------------------------------

#: Callable type for label functions.
LabelFunc = Any

#: Maps label names to their computation functions.
#: Each function has signature:
#:     (residue_list: list, index: int, unknown_str: str) -> str
LABEL_REGISTRY: dict[str, LabelFunc] = {
    'is_cis':        label_is_cis,
    'is_trans':       label_is_trans,
    'is_gly':         label_is_gly,
    'is_pro':         label_is_pro,
    'is_ileval':      label_is_ileval,
    'is_prepro':      label_is_prepro,
    'rama_category':  label_rama_category,
}
