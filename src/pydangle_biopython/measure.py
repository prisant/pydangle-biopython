"""Measurement computation for pydangle.

Computes distances, angles, and dihedral angles from parsed command
specifications and BioPython structure objects.
"""

from __future__ import annotations

import math
import random
from typing import Any

from Bio.PDB.Polypeptide import CaPPBuilder
from Bio.PDB.vectors import Vector, calc_angle, calc_dihedral

from pydangle_biopython.parser import ParsedCommand, command_string_parser

# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def angle_to_string(radians: float) -> str:
    """Convert angle from radians to a formatted degree string.

    Parameters
    ----------
    radians : float
        Angle in radians.

    Returns
    -------
    str
        Angle in degrees, formatted to 3 decimal places with trailing
        zeros stripped.
    """
    return f"{math.degrees(radians):.3f}".rstrip('0').rstrip('.')


def number_to_string(number: float) -> str:
    """Convert a number to a formatted string.

    Parameters
    ----------
    number : float

    Returns
    -------
    str
        Number formatted to 3 decimal places with trailing zeros stripped.
    """
    return f"{number:.3f}".rstrip('0').rstrip('.')


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def calc_dist(v1: Any, v2: Any) -> float:
    """Compute the Euclidean distance between two BioPython Vector objects.

    Parameters
    ----------
    v1, v2 : Bio.PDB.vectors.Vector

    Returns
    -------
    float
    """
    return float((v1 - v2).norm())


def _is_origin(cartesian_vector: Any) -> bool:
    """Check whether a vector is at the origin (all components zero)."""
    return all(abs(c) < 1e-12 for c in cartesian_vector)


def _add_jitter(cartesian_vector: Any, jitter_range: float = 0.0001) -> Any:
    """Add small random jitter to a vector to avoid degenerate geometry.

    This prevents division-by-zero errors when atoms happen to sit at
    the coordinate origin.

    Parameters
    ----------
    cartesian_vector : Bio.PDB.vectors.Vector
    jitter_range : float

    Returns
    -------
    Bio.PDB.vectors.Vector
    """
    return Vector(  # type: ignore[no-untyped-call]
        cartesian_vector[0] + random.uniform(-jitter_range, jitter_range),
        cartesian_vector[1] + random.uniform(-jitter_range, jitter_range),
        cartesian_vector[2] + random.uniform(-jitter_range, jitter_range),
    )


# ---------------------------------------------------------------------------
# Calculation dispatch
# ---------------------------------------------------------------------------

def calc_wrapper(function_key: str, args: list[Any], unknown_str: str) -> str:
    """Dispatch a geometry calculation based on the function key.

    Parameters
    ----------
    function_key : str
        One of ``'distance'``, ``'angle'``, ``'dihedral'``.
    args : list[Bio.PDB.vectors.Vector]
        Atomic coordinate vectors.
    unknown_str : str
        String to return if the calculation cannot be performed.

    Returns
    -------
    str
        Formatted measurement result, or *unknown_str* on failure.
    """
    # Guard against degenerate origin vectors
    for i, point in enumerate(args):
        if _is_origin(point):
            args[i] = _add_jitter(point)

    n = len(args)
    try:
        if function_key == 'dihedral' and n == 4:
            return angle_to_string(
                calc_dihedral(  # type: ignore[no-untyped-call]
                    args[0], args[1], args[2], args[3],
                )
            )
        elif function_key == 'angle' and n == 3:
            return angle_to_string(
                calc_angle(  # type: ignore[no-untyped-call]
                    args[0], args[1], args[2],
                )
            )
        elif function_key == 'distance' and n == 2:
            return number_to_string(calc_dist(args[0], args[1]))
    except Exception:
        # Degenerate geometry, missing atoms, etc.
        return unknown_str

    return unknown_str


# ---------------------------------------------------------------------------
# Per-residue measurement
# ---------------------------------------------------------------------------

def compute_measurement(
    command: ParsedCommand,
    residue_list: list[Any],
    residue_index: int,
    unknown_str: str,
) -> str:
    """Compute a single measurement for one residue.

    Parameters
    ----------
    command : tuple
        Parsed command tuple ``(function_key, label, arg_lists)``.
    residue_list : list[Bio.PDB.Residue.Residue]
        Contiguous list of residues (e.g. a polypeptide fragment).
    residue_index : int
        Index of the current residue in *residue_list*.
    unknown_str : str
        String to return if the measurement cannot be computed.

    Returns
    -------
    str
        Formatted measurement result, or *unknown_str*.
    """
    output = unknown_str
    list_length = len(residue_list)
    function_key = command[0]

    # Loop over alternative argument lists (separated by | in the command)
    for arg_list in command[2]:
        vector_list: list[Any] = []
        for offset, name_regex in arg_list:
            pos = residue_index + offset
            if pos < 0 or pos >= list_length:
                break  # Out of bounds → can't compute this alternative
            residue = residue_list[pos]
            matched = False
            for atom in residue:
                if name_regex.match(atom.get_fullname()):
                    vector_list.append(atom.get_vector())
                    matched = True
                    break  # Use the first matching atom
            if not matched:
                break  # Required atom not found → try next alternative

        result = calc_wrapper(function_key, vector_list, unknown_str)
        # If any alternative produces a valid result, use it
        if result != unknown_str:
            output = result
            break  # No need to try further alternatives

    return output


# ---------------------------------------------------------------------------
# Polymer residue filter
# ---------------------------------------------------------------------------

# Protein residue names accepted by Java Dangle's isProtOrNucAcid().
# Includes the 20 standard amino acids plus common modified residues.
_PROTEIN_NAMES: frozenset[str] = frozenset({
    # 20 standard amino acids
    'GLY', 'ALA', 'VAL', 'PHE', 'PRO', 'MET', 'ILE', 'LEU',
    'ASP', 'GLU', 'LYS', 'ARG', 'SER', 'THR', 'TYR', 'HIS',
    'CYS', 'ASN', 'GLN', 'TRP',
    # Ambiguous / modified
    'ASX', 'GLX',           # ambiguous Asp/Asn, Glu/Gln
    'ACE', 'FOR',           # acetyl, formyl caps
    'NH2', 'NME',           # amide caps
    'MSE',                  # selenomethionine
    'AIB', 'ABU',           # alpha-aminoisobutyric, alpha-aminobutyric acid
    'PCA',                  # pyroglutamic acid
    'MLY', 'M3L',           # methylated lysines
    'CYO', 'CSD',           # oxidized cysteine, cysteine sulfonic acid
    'DGN',                  # diaminoglutaric acid
})

# Nucleic acid residue names accepted by Java Dangle's isProtOrNucAcid().
_NUCLEIC_ACID_NAMES: frozenset[str] = frozenset({
    '  C', '  G', '  A', '  T', '  U',
    'CYT', 'GUA', 'ADE', 'THY', 'URA', 'URI',
    'GSP', 'H2U', 'PSU', '4SU', '1MG', '2MG', 'M2G',
    '5MC', '5MU', 'T6A', '1MA', 'RIA', 'OMC', 'OMG',
    ' YG', '  I', '7MG', 'YYG', 'YG ', 'A2M', '5FU',
    'G7M', 'OMU', 'PR5', 'FHU', 'XUG', 'A23', 'UMS',
    'FMU', 'UR3', 'CFL', 'UD5', 'CSL', 'UFT', '5IC',
    '5BU', 'BGM', 'CBR', 'U34', 'CCC', 'AVC', 'TM2',
    'AET', ' IU',
    'C  ', 'G  ', 'A  ', 'T  ', 'U  ', 'I  ',
    'C', 'G', 'A', 'T', 'U', 'I',
    ' rC', ' rG', ' rA', ' rT', ' rU',
    ' dC', ' dG', ' dA', ' dT', ' dU',
    ' DC', ' DG', ' DA', ' DT', ' DU',
})


def _is_polymer_residue(residue: Any) -> bool:
    """Return True if *residue* is a protein or nucleic acid residue.

    Uses the same residue-name whitelist as Java Dangle's
    ``isProtOrNucAcid()`` method in ``Measurement.java``.
    """
    resname = residue.get_resname()
    return resname in _PROTEIN_NAMES or resname in _NUCLEIC_ACID_NAMES


# ---------------------------------------------------------------------------
# Formatting helpers for output lines
# ---------------------------------------------------------------------------

def _format_residue_id(residue: Any) -> str:
    """Format a residue identifier string.

    Returns
    -------
    str
        e.g. ``"  42: :ALA"``
    """
    res_id = residue.get_id()
    seq_num = str(res_id[1]).rjust(4)
    icode = res_id[2]
    resname = residue.get_resname()
    return f"{seq_num}:{icode}:{resname}"


# ---------------------------------------------------------------------------
# Polypeptide fragment builder
# ---------------------------------------------------------------------------

def _build_fragments(model: Any) -> list[tuple[str, list[Any]]]:
    """Split a model into contiguous polypeptide fragments.

    Uses BioPython's ``CaPPBuilder`` with ``aa_only=0`` to detect
    chain breaks based on CA-CA distance (< 4.6 Angstrom), then
    filters each fragment through the Java Dangle residue-name
    whitelist.  For chains with no peptide fragments (e.g. nucleic
    acids), falls back to the raw residue list.

    Parameters
    ----------
    model : Bio.PDB.Model.Model

    Returns
    -------
    list of (chain_id, residue_list) tuples
        Each tuple contains the chain identifier and a list of
        contiguous residues forming a polypeptide fragment.
    """
    ppb = CaPPBuilder()  # type: ignore[no-untyped-call]
    fragments: list[tuple[str, list[Any]]] = []
    for chain in model:
        pp_list = ppb.build_peptides(  # type: ignore[no-untyped-call]
            chain, aa_only=0,
        )
        if pp_list:
            for pp in pp_list:
                # Filter fragment to only keep whitelisted residues
                residue_list = [
                    r for r in pp if _is_polymer_residue(r)
                ]
                if residue_list:
                    fragments.append((chain.id, residue_list))
        else:
            # No peptide fragments (e.g. nucleic acid chain).
            # Fall back to raw residue list.
            residue_list = list(chain.get_residues())
            if residue_list:
                fragments.append((chain.id, residue_list))
    return fragments


# ---------------------------------------------------------------------------
# Top-level processing
# ---------------------------------------------------------------------------

def process_measurement_for_residue(
    label: str,
    residue_list: list[Any],
    residue_index: int,
    command_list: list[ParsedCommand],
    unknown_str: str = "__?__",
) -> str | None:
    """Compute all measurements for a single residue.

    Parameters
    ----------
    label : str
        Prefix for the output line (typically filename:model:chain).
    residue_list : list[Bio.PDB.Residue.Residue]
        Contiguous polypeptide fragment.
    residue_index : int
    command_list : list[tuple]
        Parsed command list from :func:`command_string_parser`.
    unknown_str : str

    Returns
    -------
    str or None
        Formatted output line, or ``None`` if no measurements were computed.
    """
    residue = residue_list[residue_index]
    if not _is_polymer_residue(residue):
        return None

    res_str = _format_residue_id(residue)
    out_str = label + res_str
    hit_count = 0

    for command in command_list:
        result = compute_measurement(
            command, residue_list, residue_index, unknown_str,
        )
        if result != unknown_str:
            hit_count += 1
        out_str += ":" + result

    if hit_count > 0:
        return out_str
    return None


def process_measurement_commands(
    label: str,
    structure: Any,
    commands: str,
    unknown_str: str = "__?__",
) -> list[str]:
    """Process measurement commands on an entire structure.

    Parameters
    ----------
    label : str
        Identifier for the structure (typically the filename).
    structure : Bio.PDB.Structure.Structure
        Parsed BioPython structure object.
    commands : str
        User-supplied measurement command string.
    unknown_str : str
        String used when a measurement cannot be computed.

    Returns
    -------
    list[str]
        Formatted output lines, one per residue that has at least one
        valid measurement.
    """
    command_list = command_string_parser(commands)
    output_lines: list[str] = []
    output_lines.append(f"# pydangle: Query string: {commands}")

    filename_prefix = label + ":"
    for model in structure:
        model_str = filename_prefix + str(model.id + 1) + ":"
        for chain_id, residue_list in _build_fragments(model):
            chain_str = model_str + str(chain_id) + ":"
            for i in range(len(residue_list)):
                line = process_measurement_for_residue(
                    chain_str, residue_list, i,
                    command_list, unknown_str,
                )
                if line is not None:
                    output_lines.append(line)

    return output_lines
