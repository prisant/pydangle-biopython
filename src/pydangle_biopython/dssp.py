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
from typing import IO, Any

# ---------------------------------------------------------------------------
# mkdssp binary discovery
# ---------------------------------------------------------------------------

def find_mkdssp() -> str | None:
    """Return the path to the ``mkdssp`` executable, or None."""
    return shutil.which("mkdssp")


def mkdssp_version(exe: str | None = None) -> str | None:
    """Return the version reported by ``mkdssp --version`` (e.g. ``"4.2.2"``).

    Uses *exe* or, if omitted, the ``mkdssp`` found on ``$PATH``.
    Returns None if mkdssp is not found or reports no version.  Recorded
    in the JSONL metadata record because DSSP assignments depend on the
    mkdssp version (e.g. 4.2.2 vs 4.6.1 differ in beta-bridge calls).
    """
    exe = exe or find_mkdssp()
    if exe is None:
        return None
    try:
        result = subprocess.run(
            [exe, "--version"], capture_output=True, text=True,
            timeout=30, check=False,
        )
    except (subprocess.TimeoutExpired, OSError):
        return None
    for line in (result.stdout + result.stderr).splitlines():
        words = line.split()
        if "version" in words and words.index("version") + 1 < len(words):
            return words[words.index("version") + 1]
    return None


# ---------------------------------------------------------------------------
# Run mkdssp and parse output
# ---------------------------------------------------------------------------

# Record prefixes that cause mkdssp 4.x failures
_STRIP_PREFIXES = ("ANISOU", "SIGATM", "SIGUIJ")

#: mkdssp 4.x classic-format output holds at most this many distinct
#: chain IDs.  Beyond it mkdssp exits 1 with "This file contains data
#: that won't fit in the original DSSP format" and silently omits the
#: surplus chains from the residue table.  Structures with more chains
#: are run in batches of this size (see ``_run_dssp_chain_batched``),
#: following richardson-dataset-curation ``pipeline/dssp_chain_batch_fix.py``.
CHAIN_BATCH_SIZE = 26

#: Timeouts (seconds) for a whole-structure mkdssp run and for one batch
#: (the latter matches ``dssp_chain_batch_fix.py``).
_RUN_TIMEOUT = 120
_BATCH_TIMEOUT = 300


class DsspWarning(UserWarning):
    """A structure (or part of one) got no DSSP assignments.

    Issued with the file name and reason whenever mkdssp is missing,
    times out, fails, or a chain batch fails.  The command-line
    interface counts these per file and exits nonzero if any occurred,
    so affected residues (``dssp`` = null) are never silently produced.
    """


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


def _chains_in_file(filepath: str) -> list[str]:
    """Distinct chain IDs from ATOM records, in file order.

    Same rule as ``chains_in_pdb`` in ``dssp_chain_batch_fix.py``.
    """
    seen: set[str] = set()
    ordered: list[str] = []
    with open(filepath, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM  ") and len(line) > 21:
                ch = line[21]
                if ch not in seen:
                    seen.add(ch)
                    ordered.append(ch)
    return ordered


def _write_chain_batch(src: str, keep_chains: set[str], out: IO[str]) -> None:
    """Write the records of *src* belonging to *keep_chains* to *out*.

    Same record filter as ``write_clean_batch_pdb`` in
    ``dssp_chain_batch_fix.py``: HEADER and CRYST1 always; SEQRES when
    its column-12 chain is kept; ATOM/TER/HETATM when the column-22
    chain is kept.  *src* is expected to be already cleaned.
    """
    with open(src, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("HEADER", "CRYST1")):
                out.write(line)
            elif line.startswith("SEQRES"):
                if len(line) > 11 and line[11] in keep_chains:
                    out.write(line)
            elif line.startswith(("ATOM  ", "TER   ", "HETATM")):
                if len(line) > 21 and line[21] in keep_chains:
                    out.write(line)
    out.write("END\n")


def _run_mkdssp(exe: str, path: str, timeout: int) -> tuple[str | None, str]:
    """Run mkdssp once.

    Returns ``(stdout, "")`` if stdout holds a residue table, else
    ``(None, reason)``.  stdout is accepted even when mkdssp exits
    nonzero, as long as it contains the residue-table marker (see
    ``run_dssp``).
    """
    try:
        result = subprocess.run(
            [exe, "--output-format=dssp", path],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return None, f"timed out after {timeout} s"
    except OSError as exc:
        return None, f"could not run mkdssp: {exc}"
    if "#  RESIDUE" not in result.stdout:
        last = result.stderr.strip().splitlines()[-1] if result.stderr.strip() else ""
        return None, f"no residue table (exit {result.returncode}) {last}".rstrip()
    return result.stdout, ""


def _run_dssp_chain_batched(
    exe: str, path: str, chains: list[str], display_name: str,
) -> str | None:
    """Run mkdssp in batches of ``CHAIN_BATCH_SIZE`` chains and merge.

    Mirrors ``run_dssp_chain_batched`` in ``dssp_chain_batch_fix.py``:
    chains in file order are split into consecutive batches, each batch
    is written with ``_write_chain_batch`` and run separately.  The
    returned text is the first successful batch's complete output
    followed by the residue-table lines of the remaining batches.
    Unlike the pipeline module, sequential indices and BP1/BP2 columns
    are not renumbered: pydangle only reads each residue's chain,
    number, insertion code and structure code.

    Hydrogen bonds between residues in different batches are not seen
    by mkdssp, so a strand or helix that spans a batch boundary can
    lose its assignment there (documented in ``dssp_chain_batch_fix.py``).

    A failed batch issues a ``DsspWarning`` naming *display_name* (the
    original input file) and the batch's chains; the other batches are
    still used.  Returns None only if every batch fails.
    """
    batches = [
        chains[i:i + CHAIN_BATCH_SIZE]
        for i in range(0, len(chains), CHAIN_BATCH_SIZE)
    ]
    merged: list[str] = []
    for i, batch in enumerate(batches):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False,
        ) as tmp:
            _write_chain_batch(path, set(batch), tmp)
            tmp_name = tmp.name
        try:
            text, reason = _run_mkdssp(exe, tmp_name, _BATCH_TIMEOUT)
        finally:
            os.unlink(tmp_name)
        if text is None:
            warnings.warn(
                f"mkdssp failed on {display_name}, chain batch "
                f"{i + 1}/{len(batches)} (chains {''.join(batch)}): {reason}; "
                f"those residues get no DSSP code",
                DsspWarning,
                stacklevel=2,
            )
            continue
        if not merged:
            merged.append(text if text.endswith("\n") else text + "\n")
        else:
            lines = text.splitlines(keepends=True)
            start = next(
                (j + 1 for j, line in enumerate(lines)
                 if line.lstrip().startswith("#  RESIDUE")),
                len(lines),
            )
            merged.extend(lines[start:])
    return "".join(merged) if merged else None


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

    Structures with more than ``CHAIN_BATCH_SIZE`` (26) distinct chain
    IDs exceed mkdssp 4.x's classic-format limit: mkdssp exits 1 with
    "This file contains data that won't fit in the original DSSP
    format" and its residue table omits every chain after the 26th.
    (An earlier version of this docstring attributed that message to
    header-field overflow and called the table complete; for top8000
    1h64, 1ryp, 2zzs and 3mt6 -- 28, 28, 32 and 56 chains -- it was
    truncated.)  Such structures are therefore run in chain batches by
    ``_run_dssp_chain_batched``; the decision depends only on the
    chain count of the input, not on timing.

    For a whole-structure run, stdout is accepted whenever it contains
    the residue-table marker, regardless of returncode.

    Every failure issues a ``DsspWarning`` naming *filepath* and the
    reason, and the function returns None (or, for a partly failed
    chain-batched run, the assignments of the batches that succeeded).
    """
    def fail(reason: str) -> None:
        warnings.warn(
            f"no DSSP assignments for {filepath}: {reason}",
            DsspWarning,
            stacklevel=3,
        )

    exe = find_mkdssp()
    if exe is None:
        fail("mkdssp not found on PATH (install with: apt install dssp)")
        return None

    # Clean PDB if needed (mkdssp 4.x workarounds)
    tmp_path = None
    try:
        needs_clean = _needs_pdb_cleanup(filepath)
    except FileNotFoundError:
        fail("file not found")
        return None
    if needs_clean:
        tmp_path = _clean_pdb_for_dssp(filepath)
        run_path = tmp_path
    else:
        run_path = filepath

    try:
        chains = _chains_in_file(run_path)
        if len(chains) > CHAIN_BATCH_SIZE:
            batched = _run_dssp_chain_batched(exe, run_path, chains, filepath)
            if batched is None:
                fail(f"all chain batches failed ({len(chains)} chains)")
            return batched
        text, reason = _run_mkdssp(exe, run_path, _RUN_TIMEOUT)
        if text is None:
            fail(f"mkdssp {reason}")
        return text
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

    Returns None if no assignments could be computed; ``run_dssp`` has
    then already issued a ``DsspWarning`` naming the file and reason.
    """
    text = run_dssp(filepath)
    if text is None:
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
