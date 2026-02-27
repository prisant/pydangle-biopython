"""Measurement computation for pydangle.

Computes distances, angles, and dihedral angles from parsed command
specifications and BioPython structure objects.
"""

import math
import random
import re

from Bio.PDB.vectors import Vector, calc_angle, calc_dihedral

from pydangle_biopython.parser import command_string_parser

# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def angle_to_string(radians):
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


def number_to_string(number):
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

def calc_dist(v1, v2):
    """Compute the Euclidean distance between two BioPython Vector objects.

    Parameters
    ----------
    v1, v2 : Bio.PDB.vectors.Vector

    Returns
    -------
    float
    """
    return (v1 - v2).norm()


def _is_origin(cartesian_vector):
    """Check whether a vector is at the origin (all components zero)."""
    return all(abs(c) < 1e-12 for c in cartesian_vector)


def _add_jitter(cartesian_vector, jitter_range=0.0001):
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
    return Vector(
        cartesian_vector[0] + random.uniform(-jitter_range, jitter_range),
        cartesian_vector[1] + random.uniform(-jitter_range, jitter_range),
        cartesian_vector[2] + random.uniform(-jitter_range, jitter_range),
    )


# ---------------------------------------------------------------------------
# Calculation dispatch
# ---------------------------------------------------------------------------

def calc_wrapper(function_key, args, unknown_str):
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
                calc_dihedral(args[0], args[1], args[2], args[3])
            )
        elif function_key == 'angle' and n == 3:
            return angle_to_string(
                calc_angle(args[0], args[1], args[2])
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

def compute_measurement(command, chain, residue_index, unknown_str):
    """Compute a single measurement for one residue.

    Parameters
    ----------
    command : tuple
        Parsed command tuple ``(function_key, label, arg_lists)``.
    chain : Bio.PDB.Chain.Chain
        The chain containing the residue.
    residue_index : int
        Index of the current residue in ``chain.child_list``.
    unknown_str : str
        String to return if the measurement cannot be computed.

    Returns
    -------
    str
        Formatted measurement result, or *unknown_str*.
    """
    output = unknown_str
    chain_length = len(chain.child_list)
    function_key = command[0]

    # Loop over alternative argument lists (separated by | in the command)
    for arg_list in command[2]:
        vector_list = []
        for offset, name_regex in arg_list:
            pos = residue_index + offset
            if pos < 0 or pos >= chain_length:
                break  # Out of bounds → can't compute this alternative
            residue = chain.child_list[pos]
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
# Heteroatom filter
# ---------------------------------------------------------------------------

_HET_RE = re.compile(r'^H_')


def _is_het_residue(residue):
    """Return True if *residue* is a heteroatom group (water, ligand, etc.)."""
    het_flag = residue.get_id()[0]
    return het_flag != ' '


# ---------------------------------------------------------------------------
# Formatting helpers for output lines
# ---------------------------------------------------------------------------

def _format_residue_id(residue):
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
# Top-level processing
# ---------------------------------------------------------------------------

def process_measurement_for_residue(
    label, chain, residue_index, command_list, unknown_str="__?__"
):
    """Compute all measurements for a single residue.

    Parameters
    ----------
    label : str
        Prefix for the output line (typically filename:model:chain).
    chain : Bio.PDB.Chain.Chain
    residue_index : int
    command_list : list[tuple]
        Parsed command list from :func:`command_string_parser`.
    unknown_str : str

    Returns
    -------
    str or None
        Formatted output line, or ``None`` if no measurements were computed.
    """
    residue = chain.child_list[residue_index]
    if _is_het_residue(residue):
        return None

    res_str = _format_residue_id(residue)
    out_str = label + res_str
    hit_count = 0

    for command in command_list:
        result = compute_measurement(command, chain, residue_index, unknown_str)
        if result != unknown_str:
            hit_count += 1
        out_str += ":" + result

    if hit_count > 0:
        return out_str
    return None


def process_measurement_commands(label, structure, commands, unknown_str="__?__"):
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
    output_lines = []
    output_lines.append(f"# pydangle: Query string: {commands}")

    filename_prefix = label + ":"
    for model in structure:
        model_str = filename_prefix + str(model.id + 1) + ":"
        for chain in model:
            chain_str = model_str + str(chain.id) + ":"
            for i in range(len(chain.child_list)):
                line = process_measurement_for_residue(
                    chain_str, chain, i, command_list, unknown_str
                )
                if line is not None:
                    output_lines.append(line)

    return output_lines
