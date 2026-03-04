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
# Atom completeness labels
# ---------------------------------------------------------------------------

#: Mainchain atom names required for all amino acid residues.
_MAINCHAIN_ATOMS: frozenset[str] = frozenset({'N', 'CA', 'C', 'O'})


def label_has_all_mc(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if all mainchain heavy atoms are present.

    Checks for N, CA, C, and O.
    """
    residue = residue_list[index]
    atom_names = {atom.get_fullname().strip() for atom in residue}
    return str(_MAINCHAIN_ATOMS.issubset(atom_names))


#: Expected sidechain heavy atoms for each standard amino acid.
#: Atom names are stripped (no space padding).
_SIDECHAIN_ATOMS: dict[str, frozenset[str]] = {
    'GLY': frozenset(),
    'ALA': frozenset({'CB'}),
    'VAL': frozenset({'CB', 'CG1', 'CG2'}),
    'LEU': frozenset({'CB', 'CG', 'CD1', 'CD2'}),
    'ILE': frozenset({'CB', 'CG1', 'CG2', 'CD1'}),
    'PRO': frozenset({'CB', 'CG', 'CD'}),
    'PHE': frozenset({'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'}),
    'TYR': frozenset({
        'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH',
    }),
    'TRP': frozenset({
        'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3',
        'CZ2', 'CZ3', 'CH2',
    }),
    'SER': frozenset({'CB', 'OG'}),
    'THR': frozenset({'CB', 'OG1', 'CG2'}),
    'CYS': frozenset({'CB', 'SG'}),
    'MET': frozenset({'CB', 'CG', 'SD', 'CE'}),
    'ASP': frozenset({'CB', 'CG', 'OD1', 'OD2'}),
    'GLU': frozenset({'CB', 'CG', 'CD', 'OE1', 'OE2'}),
    'ASN': frozenset({'CB', 'CG', 'OD1', 'ND2'}),
    'GLN': frozenset({'CB', 'CG', 'CD', 'OE1', 'NE2'}),
    'LYS': frozenset({'CB', 'CG', 'CD', 'CE', 'NZ'}),
    'ARG': frozenset({
        'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2',
    }),
    'HIS': frozenset({'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'}),
}


def label_has_all_sc(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if all expected sidechain heavy atoms are present.

    Uses a lookup table of expected atoms for each standard amino acid.
    Returns *unknown* for non-standard residues not in the table.
    GLY always returns ``'True'`` (no sidechain atoms expected).
    """
    residue = residue_list[index]
    resname = residue.get_resname()
    expected = _SIDECHAIN_ATOMS.get(resname)
    if expected is None:
        return unknown
    if not expected:
        return 'True'  # GLY
    atom_names = {atom.get_fullname().strip() for atom in residue}
    return str(expected.issubset(atom_names))



# ---------------------------------------------------------------------------
# Chirality labels
# ---------------------------------------------------------------------------

def _compute_ca_chirality(residue: Any) -> float | None:
    """Compute the improper dihedral N-CA-C-CB in degrees.

    Negative for L-amino acids (~-33), positive for D-amino acids.
    Returns None if any required atom is missing (e.g. GLY).
    """
    try:
        v_cb = residue['CB'].get_vector()
        v_n = residue['N'].get_vector()
        v_c = residue['C'].get_vector()
        v_ca = residue['CA'].get_vector()
    except (KeyError, TypeError):
        return None
    radians = calc_dihedral(  # type: ignore[no-untyped-call]
        v_cb, v_n, v_c, v_ca,
    )
    return float(math.degrees(radians))


def label_is_left(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if the residue has L-amino acid chirality.

    Based on the improper dihedral CB-N-C-CA being negative.
    Returns *unknown* for GLY or if CB is missing.
    """
    chi = _compute_ca_chirality(residue_list[index])
    if chi is None:
        return unknown
    return str(chi < 0.0)


def label_is_right(
    residue_list: list[Any], index: int, unknown: str,
) -> str:
    """Return ``'True'`` if the residue has D-amino acid chirality.

    Based on the improper dihedral CB-N-C-CA being positive.
    Returns *unknown* for GLY or if CB is missing.
    """
    chi = _compute_ca_chirality(residue_list[index])
    if chi is None:
        return unknown
    return str(chi > 0.0)

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
    'has_all_mc':     label_has_all_mc,
    'has_all_sc':     label_has_all_sc,
    'is_left':        label_is_left,
    'is_right':       label_is_right,
    'rama_category':  label_rama_category,
}
