#!/usr/bin/env python3
"""Build ersatz full-structure PDB files for the top2018 dataset.

For each entry in the dataset, creates an "ersatz" PDB by starting with
the Reduce-processed FH file (which has correct NQH flips on all chains)
and stripping hydrogens and problematic headers.  Also generates a
residue mask file from the pruned file.

The ersatz files preserve multi-chain context needed for correct DSSP
secondary structure assignments, with NQH corrections applied to all
chains.

Prerequisites:
    Run scripts/run_reduce.py first to generate {pdbid}FH.pdb files.

Usage:
    python scripts/build_ersatz_pdbs.py SRC_DIR [-j JOBS]

Input per entry (e.g., 1a/1a2z/):
    1a2z.pdb              - original full multi-chain PDB (no H, no flips)
    1a2zFH.pdb            - Reduce output: full PDB with H + NQH flips
    1a2z_C_pruned_all.pdb - pruned chain C (for mask generation only)

Output per entry (in-place):
    1a2z_ersatz.pdb       - full PDB with NQH corrections, no H, clean headers
    1a2z_C.mask           - pruned residue list (chain resnum icode)
"""

from __future__ import annotations

import argparse
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# PDB record types to keep in ersatz files (same as clean_pdb_headers.py)
_KEEP_PREFIXES = (
    "HEADER", "CRYST1",
    "SCALE1", "SCALE2", "SCALE3",
    "ORIGX1", "ORIGX2", "ORIGX3",
    "MTRIX1", "MTRIX2", "MTRIX3",
    "MODEL ", "ENDMDL", "END   ",
    "ATOM  ", "HETATM", "ANISOU", "TER   ",
)


def _build_mask(pruned_path: Path) -> list[tuple[str, int, str]]:
    """Extract residue identifiers from pruned file for mask output.

    Returns list of (chain_id, resnum, icode) tuples, one per residue.
    """
    seen: set[tuple[str, int, str]] = set()
    residues: list[tuple[str, int, str]] = []
    with open(pruned_path) as fh:
        for line in fh:
            if not line.startswith("ATOM  "):
                continue
            chain_id = line[21]
            resnum = int(line[22:26])
            icode = line[26]
            key = (chain_id, resnum, icode)
            if key not in seen:
                seen.add(key)
                residues.append(key)
    return residues


def _is_hydrogen(line: str) -> bool:
    """Check if an ATOM/HETATM line is a hydrogen atom."""
    # Element symbol is in columns 77-78 (0-indexed 76-77)
    if len(line) >= 78:
        element = line[76:78].strip()
        return element in ("H", "D")
    # Fallback: check atom name (col 13-16) for hydrogen naming conventions
    atom_name = line[12:16]
    first_nonspace = atom_name.lstrip()[0] if atom_name.strip() else ""
    return first_nonspace in ("H", "D")


def process_entry(entry_dir: Path) -> tuple[str, bool, str, int]:
    """Process a single entry directory.

    Creates {pdbid}_ersatz.pdb and {pdbid}_{chain}.mask in the entry dir.

    Returns (pdb_id, success, message, n_mask_residues).
    """
    pdb_id = entry_dir.name

    # Find input files
    fh_pdb = entry_dir / f"{pdb_id}FH.pdb"
    pruned_files = list(entry_dir.glob(f"{pdb_id}_*_pruned_all.pdb"))

    if not fh_pdb.exists():
        return (pdb_id, False, "missing FH file", 0)
    if not pruned_files:
        return (pdb_id, False, "missing pruned PDB", 0)

    pruned_pdb = pruned_files[0]
    # Extract target chain from filename: {pdb}_{chain}_pruned_all.pdb
    chain_match = re.match(
        rf"{re.escape(pdb_id)}_(\w)_pruned_all\.pdb", pruned_pdb.name
    )
    if not chain_match:
        return (pdb_id, False, f"cannot parse chain from {pruned_pdb.name}", 0)
    target_chain = chain_match.group(1)

    try:
        # Build ersatz PDB: FH file with H stripped and headers cleaned
        ersatz_path = entry_dir / f"{pdb_id}_ersatz.pdb"
        atom_count = 0

        with open(fh_pdb) as fin, open(ersatz_path, "w") as fout:
            for line in fin:
                # Header filtering
                if not line.startswith(_KEEP_PREFIXES):
                    continue

                # Skip hydrogen atoms
                if line.startswith(("ATOM  ", "HETATM")) and _is_hydrogen(line):
                    continue

                fout.write(line)
                if line.startswith("ATOM  "):
                    atom_count += 1

        # Build mask file from pruned PDB
        mask_residues = _build_mask(pruned_pdb)
        mask_path = entry_dir / f"{pdb_id}_{target_chain}.mask"
        with open(mask_path, "w") as fout:
            for chain_id, resnum, icode in mask_residues:
                fout.write(f"{chain_id} {resnum} {icode}\n")

        msg = (
            f"ok: chain={target_chain}, "
            f"{atom_count} heavy atoms, "
            f"{len(mask_residues)} mask residues"
        )
        return (pdb_id, True, msg, len(mask_residues))

    except Exception as exc:
        return (pdb_id, False, str(exc), 0)


def collect_entries(src_dir: Path) -> list[Path]:
    """Find all entry directories (two-level: prefix/pdbid/)."""
    entries: list[Path] = []
    for prefix_dir in sorted(src_dir.iterdir()):
        if not prefix_dir.is_dir():
            continue
        for entry_dir in sorted(prefix_dir.iterdir()):
            if not entry_dir.is_dir():
                continue
            entries.append(entry_dir)
    return entries


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build ersatz PDB files from Reduce FH output. "
            "Strips hydrogens and problematic headers, generates mask files."
        )
    )
    parser.add_argument("src_dir", type=Path, help="Top2018 dataset directory")
    parser.add_argument(
        "-j", "--jobs", type=int, default=8,
        help="Number of parallel workers (default: 8)",
    )
    args = parser.parse_args()

    if not args.src_dir.is_dir():
        print(f"Error: {args.src_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning {args.src_dir} ...")
    entries = collect_entries(args.src_dir)
    print(f"Found {len(entries)} entries")

    if not entries:
        print("Nothing to do.")
        return

    done = 0
    errors = 0
    total_mask = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(process_entry, entry): entry
            for entry in entries
        }
        for future in as_completed(futures):
            pdb_id, success, msg, n_mask = future.result()
            done += 1
            if not success:
                errors += 1
                print(f"  ERROR {pdb_id}: {msg}", file=sys.stderr)
            else:
                total_mask += n_mask
            if done % 500 == 0 or done == len(entries):
                print(
                    f"  {done}/{len(entries)} entries processed",
                    flush=True,
                )

    print(f"\nDone: {done - errors} succeeded, {errors} errors")
    print(f"Total mask residues: {total_mask}")
    print(f"Output: in-place alongside originals in {args.src_dir}")


if __name__ == "__main__":
    main()
