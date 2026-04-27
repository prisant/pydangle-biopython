"""DSSP secondary structure assignment for pydangle.

Calls the external ``mkdssp`` binary to compute per-residue secondary
structure assignments using the Kabsch & Sander hydrogen-bond energy
method (Biopolymers 22:2577-2637, 1983).

Requires ``mkdssp`` to be installed and available on ``$PATH``.
Install via: ``apt install dssp`` (Debian/Ubuntu),
``brew install brewsci/bio/dssp`` (macOS), or
``conda install -c conda-forge dssp``.

Provides two label functions:

    dssp   -- DSSP code (H, B, E, G, I, T, S, P) or null for loop/coil
    dssp3  -- reduced 3-state (H = helix, E = strand, C = coil/loop)
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
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

# Record prefixes that cause mkdssp 4.x failures
_STRIP_PREFIXES = ("ANISOU", "SIGATM", "SIGUIJ")


def _needs_pdb_cleanup(filepath: str) -> bool:
    """Check if a PDB file needs cleanup before mkdssp."""
    with open(filepath, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(_STRIP_PREFIXES):
                return True
            if line.startswith("HETATM"):
                return True
            if (
                line.startswith(("ATOM  ", "TER   "))
                and len(line) > 21
                and line[21] == " "
            ):
                return True
            if (
                line.startswith("SEQRES")
                and len(line) > 11
                and line[11] == " "
            ):
                return True
    return False


def _clean_pdb_for_dssp(filepath: str) -> str:
    """Write a temporary copy of *filepath* cleaned for mkdssp.

    Strips ANISOU/SIGATM/SIGUIJ records and selectively strips HETATM
    records: only HETATM whose residue name is not part of the
    protein chain (i.e., not in SEQRES and not in ATOM resnames) are
    removed.  This drops ligand / water / cofactor heteroatoms (HOH,
    SO4, GOL, NAG, ZN, ...) that would collide with the SEQRES
    protein alignment and silently null mkdssp's output, while
    preserving modified amino-acid residues (MSE, PCA, MEN, PTR, LOV,
    ...) that some old PDBs deposit as HETATM rather than ATOM.

    Also fills blank chain IDs in ATOM/HETATM/TER records (column 22)
    and SEQRES records (column 12) with the first non-blank chain ID
    found in ATOM records, falling back to ``'A'``.

    Returns the path to the temporary file.  Caller is responsible
    for cleanup.
    """
    # Single first pass: collect default chain + protein residue names
    # from SEQRES and ATOM records.
    default_chain = "A"
    chain_seen = False
    protein_residues: set[str] = set()
    with open(filepath, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM  ") and len(line) > 21:
                if not chain_seen:
                    ch = line[21]
                    if ch != " ":
                        default_chain = ch
                        chain_seen = True
                if len(line) >= 20:
                    protein_residues.add(line[17:20].strip())
            elif line.startswith("SEQRES") and len(line) > 19:
                # SEQRES residue names are 3-char fields starting at
                # column 20 (0-indexed 19), spaced every 4 columns.
                idx = 19
                while idx + 3 <= len(line):
                    name = line[idx:idx + 3].strip()
                    if name:
                        protein_residues.add(name)
                    idx += 4

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".pdb", delete=False,
    ) as tmp:
        with open(filepath, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(_STRIP_PREFIXES):
                    continue
                if line.startswith("HETATM"):
                    resname = line[17:20].strip() if len(line) >= 20 else ""
                    if resname not in protein_residues:
                        continue
                out = line
                if (
                    line.startswith(("ATOM  ", "HETATM", "TER   "))
                    and len(line) > 21
                    and line[21] == " "
                ):
                    out = line[:21] + default_chain + line[22:]
                elif (
                    line.startswith("SEQRES")
                    and len(line) > 11
                    and line[11] == " "
                ):
                    out = line[:11] + default_chain + line[12:]
                tmp.write(out)
        tmp_name = tmp.name
    return tmp_name


def run_dssp(filepath: str) -> str | None:
    """Run mkdssp on *filepath* and return the classic-format output.

    Returns None if mkdssp is not installed or the command fails.

    Before invoking mkdssp, the input is cleaned to work around
    mkdssp 4.x issues: ANISOU/SIGATM/SIGUIJ records are stripped
    (they cause incomplete or zero-residue output), HETATM records
    are stripped (they cause silent zero-residue output when ligand /
    water / cofactor heteroatoms collide with SEQRES protein
    alignment), and blank chain IDs are filled (they cause parse
    failures on old PDB depositions).

    mkdssp 4.x exits 1 with the message "This file contains data that
    won't fit in the original DSSP format" for structures whose HEADER
    or COMPND fields exceed legacy fixed-width column widths (typically
    multi-chain structures with >9 chains and long molecule descriptors,
    e.g. 1h64, 1ryp, 2zzs, 3mt6 in top8000).  The residue-by-residue
    secondary-structure data in stdout is nonetheless complete and
    correct in this case (verified residue-identical against mkdssp's
    file-output mode which exits 0 for the same input).  This wrapper
    treats stdout as authoritative when it contains the residue-table
    marker, regardless of returncode.
    """
    exe = find_mkdssp()
    if exe is None:
        return None

    # Clean PDB if needed (mkdssp 4.x workarounds)
    tmp_path = None
    try:
        needs_clean = _needs_pdb_cleanup(filepath)
    except FileNotFoundError:
        return None
    if needs_clean:
        tmp_path = _clean_pdb_for_dssp(filepath)
        run_path = tmp_path
    else:
        run_path = filepath

    try:
        result = subprocess.run(
            [exe, "--output-format=dssp", run_path],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        # mkdssp 4.x exits 1 with a "won't fit in original DSSP format"
        # warning for structures whose header fields overflow legacy
        # fixed-width columns; the residue table in stdout is still
        # complete and correct.  Trust stdout when it contains the
        # residue-table marker, regardless of returncode.
        if "#  RESIDUE" not in result.stdout:
            return None
        return result.stdout
    except (subprocess.TimeoutExpired, OSError):
        return None
    finally:
        if tmp_path is not None:
            os.unlink(tmp_path)


def parse_dssp_output(
    text: str,
) -> dict[tuple[str, int, str], str | None]:
    """Parse classic DSSP output into a per-residue assignment dict.

    Returns a dict mapping ``(chain_id, residue_number, insertion_code)``
    to the single-character DSSP assignment code, or ``None`` for
    residues that DSSP classified as loop (space in column 17).

    DSSP never outputs a literal ``'C'`` character.  The valid codes
    are: H, B, E, G, I, T, S, P, or space (loop).

    The classic DSSP format has fixed-width columns::

        columns  6-10: residue number (right-justified)
        column     12: chain ID
        column     14: amino acid one-letter code
        column     17: DSSP secondary structure code (space = loop)
    """
    assignments: dict[tuple[str, int, str], str | None] = {}
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
        dssp_code: str | None = line[16:17]

        # Space means loop — preserve as None (DSSP never outputs 'C')
        if dssp_code == " ":
            dssp_code = None

        assignments[(chain_id, resnum, icode)] = dssp_code

    return assignments


# ---------------------------------------------------------------------------
# 8-state to 3-state reduction
# ---------------------------------------------------------------------------

#: Map DSSP codes to 3-state (H=helix, E=strand, C=coil).
#: H, G, I (alpha, 3-10, pi helix) -> H
#: E, B (strand, bridge) -> E
#: T, S, P (turn, bend, PPII) -> C
#: None (loop/space) -> C
_DSSP8_TO_DSSP3: dict[str | None, str] = {
    "H": "H",
    "G": "H",
    "I": "H",
    "E": "E",
    "B": "E",
    "T": "C",
    "S": "C",
    "P": "C",
    None: "C",
}


# ---------------------------------------------------------------------------
# Module-level state: current DSSP assignments
#
# This is set by the measurement pipeline before processing residues
# and cleared afterwards.  It avoids changing the label function
# signature while providing DSSP data to the label functions.
# ---------------------------------------------------------------------------

_current_assignments: dict[tuple[str, int, str], str | None] | None = None


def set_dssp_assignments(
    assignments: dict[tuple[str, int, str], str | None] | None,
) -> None:
    """Set the current DSSP assignments for label dispatch."""
    global _current_assignments  # noqa: PLW0603
    _current_assignments = assignments


def get_dssp_assignments_for_file(filepath: str) -> (
    dict[tuple[str, int, str], str | None] | None
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

def _lookup_dssp(
    residue_list: list[Any],
    index: int,
) -> tuple[bool, str | None]:
    """Look up the DSSP code for a residue.

    Returns ``(found, code)`` where *found* is False if mkdssp
    data is unavailable or the residue is absent from the output,
    and *code* is the DSSP character or None for loop.
    """
    if _current_assignments is None:
        return False, None

    residue = residue_list[index]
    chain_id = residue.get_parent().get_id()
    res_id = residue.get_id()
    resnum = res_id[1]
    icode = res_id[2].strip()

    key = (chain_id, resnum, icode)
    if key not in _current_assignments:
        return False, None
    return True, _current_assignments[key]


def label_dssp(
    residue_list: list[Any],
    index: int,
    unknown: str,
) -> str:
    """Return the DSSP secondary structure code.

    Codes: H (alpha helix), B (beta bridge), E (extended strand),
    G (3-10 helix), I (pi helix), T (turn), S (bend), P (PPII helix).

    Returns *unknown* if mkdssp is not available, the residue is not
    found in the DSSP output, or DSSP classified the residue as loop
    (space in the native output).
    """
    found, code = _lookup_dssp(residue_list, index)
    if not found or code is None:
        return unknown
    return code


def label_dssp3(
    residue_list: list[Any],
    index: int,
    unknown: str,
) -> str:
    """Return the 3-state DSSP secondary structure code.

    Codes: H (helix: H/G/I), E (strand: E/B), C (coil: T/S/P/loop).

    Returns *unknown* if mkdssp is not available or the residue
    is not found in the DSSP output.
    """
    found, code = _lookup_dssp(residue_list, index)
    if not found:
        return unknown
    return _DSSP8_TO_DSSP3.get(code, unknown)
