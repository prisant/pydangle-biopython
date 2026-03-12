"""DSSP secondary structure assignment for pydangle.

Calls the external ``mkdssp`` binary to compute per-residue secondary
structure assignments using the Kabsch & Sander hydrogen-bond energy
method (Biopolymers 22:2577-2637, 1983).

Requires ``mkdssp`` to be installed and available on ``$PATH``.
Install via: ``apt install dssp`` (Debian/Ubuntu),
``brew install brewsci/bio/dssp`` (macOS), or
``conda install -c conda-forge dssp``.

Provides two label functions:

    dssp   -- full 8-state DSSP code (H, B, E, G, I, T, S, or C for coil)
    dssp3  -- reduced 3-state (H = helix, E = strand, C = coil)
"""

from __future__ import annotations

import shutil
import subprocess
import warnings
from typing import Any

# ---------------------------------------------------------------------------
# mkdssp binary discovery
# ---------------------------------------------------------------------------

def find_mkdssp() -> str | None:
    """Return the path to the ``mkdssp`` executable, or None."""
    return shutil.which("mkdssp")


# ---------------------------------------------------------------------------
# Run mkdssp and parse output
# ---------------------------------------------------------------------------

def run_dssp(filepath: str) -> str | None:
    """Run mkdssp on *filepath* and return the classic-format output.

    Returns None if mkdssp is not installed or the command fails.
    """
    exe = find_mkdssp()
    if exe is None:
        return None
    try:
        result = subprocess.run(
            [exe, "--output-format=dssp", filepath],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            return None
        return result.stdout
    except (subprocess.TimeoutExpired, OSError):
        return None


def parse_dssp_output(
    text: str,
) -> dict[tuple[str, int, str], str]:
    """Parse classic DSSP output into a per-residue assignment dict.

    Returns a dict mapping ``(chain_id, residue_number, insertion_code)``
    to the single-character DSSP assignment code.  A space (coil) is
    returned as ``'C'`` for clarity.

    The classic DSSP format has fixed-width columns::

        columns  6-10: residue number (right-justified)
        column     12: chain ID
        column     14: amino acid one-letter code
        column     17: DSSP secondary structure code (space = coil)
    """
    assignments: dict[tuple[str, int, str], str] = {}
    in_data = False

    for line in text.splitlines():
        # Data section starts after the header line containing
        # "#  RESIDUE AA STRUCTURE"
        if line.strip().startswith("#  RESIDUE AA"):
            in_data = True
            continue
        if not in_data:
            continue
        if len(line) < 17:
            continue

        # Skip chain-break markers (indicated by '!' in column 14)
        aa = line[13:14]
        if aa == "!":
            continue

        # Parse fields from fixed-width columns
        try:
            resnum = int(line[5:10].strip())
        except ValueError:
            continue

        chain_id = line[11:12]
        icode = line[10:11].strip()  # insertion code (usually space)
        dssp_code = line[16:17]

        # Normalize: space means coil -> 'C'
        if dssp_code == " ":
            dssp_code = "C"

        assignments[(chain_id, resnum, icode)] = dssp_code

    return assignments


# ---------------------------------------------------------------------------
# 8-state to 3-state reduction
# ---------------------------------------------------------------------------

#: Map 8-state DSSP codes to 3-state (H=helix, E=strand, C=coil).
#: H, G, I (alpha, 3-10, pi helix) -> H
#: E, B (strand, bridge) -> E
#: T, S, C (turn, bend, coil) -> C
_DSSP8_TO_DSSP3: dict[str, str] = {
    "H": "H",
    "G": "H",
    "I": "H",
    "E": "E",
    "B": "E",
    "T": "C",
    "S": "C",
    "C": "C",
    "P": "C",
}


# ---------------------------------------------------------------------------
# Module-level state: current DSSP assignments
#
# This is set by the measurement pipeline before processing residues
# and cleared afterwards.  It avoids changing the label function
# signature while providing DSSP data to the label functions.
# ---------------------------------------------------------------------------

_current_assignments: dict[tuple[str, int, str], str] | None = None


def set_dssp_assignments(
    assignments: dict[tuple[str, int, str], str] | None,
) -> None:
    """Set the current DSSP assignments for label dispatch."""
    global _current_assignments  # noqa: PLW0603
    _current_assignments = assignments


def get_dssp_assignments_for_file(filepath: str) -> (
    dict[tuple[str, int, str], str] | None
):
    """Run mkdssp on *filepath* and return assignments.

    Issues a warning if mkdssp is not available.
    """
    text = run_dssp(filepath)
    if text is None:
        warnings.warn(
            "mkdssp not found or failed. "
            "Install with: apt install dssp",
            stacklevel=2,
        )
        return None
    return parse_dssp_output(text)


# ---------------------------------------------------------------------------
# Label functions
# ---------------------------------------------------------------------------

def label_dssp(
    residue_list: list[Any],
    index: int,
    unknown: str,
) -> str:
    """Return the 8-state DSSP secondary structure code.

    Codes: H (alpha helix), B (beta bridge), E (extended strand),
    G (3-10 helix), I (pi helix), T (turn), S (bend), C (coil).

    Returns *unknown* if mkdssp is not available or the residue
    is not found in the DSSP output.
    """
    if _current_assignments is None:
        return unknown

    residue = residue_list[index]
    chain_id = residue.get_parent().get_id()
    res_id = residue.get_id()
    resnum = res_id[1]
    icode = res_id[2].strip()

    key = (chain_id, resnum, icode)
    code = _current_assignments.get(key)
    if code is None:
        return unknown
    return code


def label_dssp3(
    residue_list: list[Any],
    index: int,
    unknown: str,
) -> str:
    """Return the 3-state DSSP secondary structure code.

    Codes: H (helix: H/G/I), E (strand: E/B), C (coil: T/S/C).

    Returns *unknown* if mkdssp is not available or the residue
    is not found in the DSSP output.
    """
    code8 = label_dssp(residue_list, index, unknown)
    if code8 == unknown:
        return unknown
    return _DSSP8_TO_DSSP3.get(code8, unknown)
