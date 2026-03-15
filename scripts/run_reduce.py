#!/usr/bin/env python3
"""Run Reduce on full PDB files to add hydrogens and flip NQH sidechains.

Processes each {pdbid}.pdb in the top2018 dataset directory and writes
{pdbid}FH.pdb alongside it (Richardson Lab naming convention: Flipped +
Hydrogens).

Usage:
    python scripts/run_reduce.py SRC_DIR [-j JOBS] [--het-db PATH]

Example:
    python scripts/run_reduce.py ~/Desktop/top2018_pdbs_full_filtered_hom70 -j 8
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_DEFAULT_HET_DB = os.path.expanduser(
    "~/bin/reduce_wwPDB_het_dict.txt"
)


def run_reduce_on_file(
    pdb_path: Path, het_db: str
) -> tuple[str, bool, str]:
    """Run Reduce on a single PDB file.

    Returns (pdb_id, success, message).
    """
    pdb_id = pdb_path.stem  # e.g., "1a2z"
    output_path = pdb_path.parent / f"{pdb_id}FH.pdb"

    try:
        result = subprocess.run(
            [
                "reduce",
                "-db", het_db,
                "-allalt",
                "-build",
                "-norotmet",
                str(pdb_path),
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0 and not result.stdout:
            return (pdb_id, False, f"reduce failed: {result.stderr[:200]}")

        with open(output_path, "w") as fh:
            fh.write(result.stdout)

        return (pdb_id, True, "ok")

    except subprocess.TimeoutExpired:
        return (pdb_id, False, "timeout (120s)")
    except Exception as exc:
        return (pdb_id, False, str(exc))


def collect_pdb_files(src_dir: Path) -> list[Path]:
    """Find original (non-pruned, non-FH) PDB files in two-level structure."""
    files: list[Path] = []
    for prefix_dir in sorted(src_dir.iterdir()):
        if not prefix_dir.is_dir():
            continue
        for entry_dir in sorted(prefix_dir.iterdir()):
            if not entry_dir.is_dir():
                continue
            pdb_id = entry_dir.name
            pdb_file = entry_dir / f"{pdb_id}.pdb"
            if pdb_file.exists():
                files.append(pdb_file)
    return files


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run Reduce on full PDB files (add H + flip NQH)"
    )
    parser.add_argument("src_dir", type=Path, help="Top2018 dataset directory")
    parser.add_argument(
        "-j", "--jobs", type=int, default=8,
        help="Number of parallel workers (default: 8)",
    )
    parser.add_argument(
        "--het-db", type=str, default=_DEFAULT_HET_DB,
        help=f"Path to Reduce het dictionary (default: {_DEFAULT_HET_DB})",
    )
    args = parser.parse_args()

    if not args.src_dir.is_dir():
        print(f"Error: {args.src_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.het_db):
        print(f"Error: het dictionary not found: {args.het_db}", file=sys.stderr)
        sys.exit(1)

    # Verify reduce is available
    try:
        result = subprocess.run(
            ["reduce", "-version"], capture_output=True, text=True
        )
        lines = result.stderr.strip().splitlines()
        version_line = lines[-1] if lines else "unknown"
        print(f"Reduce version: {version_line}")
    except FileNotFoundError:
        print("Error: reduce not found in PATH", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning {args.src_dir} ...")
    pdb_files = collect_pdb_files(args.src_dir)
    print(f"Found {len(pdb_files)} PDB files")

    if not pdb_files:
        print("Nothing to do.")
        return

    done = 0
    errors = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(run_reduce_on_file, pdb, args.het_db): pdb
            for pdb in pdb_files
        }
        for future in as_completed(futures):
            pdb_id, success, msg = future.result()
            done += 1
            if not success:
                errors += 1
                print(f"  ERROR {pdb_id}: {msg}", file=sys.stderr)
            if done % 500 == 0 or done == len(pdb_files):
                print(
                    f"  {done}/{len(pdb_files)} files processed",
                    flush=True,
                )

    print(f"\nDone: {done - errors} succeeded, {errors} errors")


if __name__ == "__main__":
    main()
