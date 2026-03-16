#!/usr/bin/env python3
"""Build ersatz PDB files from Reduce FH output (no pruned files needed).

For datasets without Richardson Lab pruned files (e.g., top8000), this script
builds ersatz PDB files directly from the FH files by stripping hydrogens
and problematic headers. Mask generation is handled separately by
build_bfactor_masks.py.

Usage:
    python scripts/build_ersatz_only.py SRC_DIR [-j JOBS]

Output per entry:
    {pdbid}_ersatz.pdb — full PDB with NQH corrections, no H, clean headers
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_KEEP_PREFIXES = (
    "HEADER", "CRYST1",
    "SCALE1", "SCALE2", "SCALE3",
    "ORIGX1", "ORIGX2", "ORIGX3",
    "MTRIX1", "MTRIX2", "MTRIX3",
    "MODEL ", "ENDMDL", "END   ",
    "ATOM  ", "HETATM", "ANISOU", "TER   ",
)


def _is_hydrogen(line: str) -> bool:
    """Check if an ATOM/HETATM line is a hydrogen atom."""
    if len(line) >= 78:
        element = line[76:78].strip()
        return element in ("H", "D")
    atom_name = line[12:16]
    first_nonspace = atom_name.lstrip()[0] if atom_name.strip() else ""
    return first_nonspace in ("H", "D")


def process_entry(entry_dir: Path) -> tuple[str, bool, str]:
    """Build ersatz PDB from FH file.

    Returns (pdb_id, success, message).
    """
    pdb_id = entry_dir.name
    fh_pdb = entry_dir / f"{pdb_id}FH.pdb"

    if not fh_pdb.exists():
        return (pdb_id, False, "missing FH file")

    try:
        ersatz_path = entry_dir / f"{pdb_id}_ersatz.pdb"
        atom_count = 0

        with open(fh_pdb) as fin, open(ersatz_path, "w") as fout:
            for line in fin:
                if not line.startswith(_KEEP_PREFIXES):
                    continue
                if line.startswith(("ATOM  ", "HETATM")) and _is_hydrogen(line):
                    continue
                fout.write(line)
                if line.startswith("ATOM  "):
                    atom_count += 1

        return (pdb_id, True, f"ok: {atom_count} heavy atoms")

    except Exception as exc:
        return (pdb_id, False, str(exc))


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
            "Strips hydrogens and problematic headers. "
            "Does not require pruned files."
        )
    )
    parser.add_argument(
        "src_dir", type=Path,
        help="Dataset directory (two-level: prefix/pdbid/)",
    )
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

    done = 0
    errors = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(process_entry, entry): entry
            for entry in entries
        }
        for future in as_completed(futures):
            pdb_id, success, msg = future.result()
            done += 1
            if not success:
                errors += 1
                print(f"  ERROR {pdb_id}: {msg}", file=sys.stderr)
            if done % 500 == 0 or done == len(entries):
                print(
                    f"  {done}/{len(entries)} entries processed",
                    flush=True,
                )

    print(f"\nDone: {done - errors} succeeded, {errors} errors")


if __name__ == "__main__":
    main()
