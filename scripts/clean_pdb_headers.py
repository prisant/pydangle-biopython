#!/usr/bin/env python3
"""Preprocess PDB files for mkdssp compatibility.

Strips non-essential header records (USER, TITLE, SOURCE, COMPND, etc.)
that can cause mkdssp 4.x to reject files due to strict format parsing.
Keeps only records needed for structure analysis:

    HEADER, CRYST1, SCALEn, ORIGXn, MTRIXn,
    ATOM, HETATM, ANISOU, TER, MODEL, ENDMDL, END

Preserves the source directory structure in the output directory.

Usage:
    python scripts/clean_pdb_headers.py SRC_DIR DST_DIR [-j JOBS]
"""

from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# PDB record types that mkdssp needs or that contain coordinate data
_KEEP_PREFIXES = (
    "HEADER", "CRYST1",
    "SCALE1", "SCALE2", "SCALE3",
    "ORIGX1", "ORIGX2", "ORIGX3",
    "MTRIX1", "MTRIX2", "MTRIX3",
    "MODEL ", "ENDMDL", "END   ",
    "ATOM  ", "HETATM", "ANISOU", "TER   ",
)


def clean_pdb(src: Path, dst: Path) -> tuple[str, bool, str]:
    """Clean a single PDB file.  Returns (relpath, success, message)."""
    relpath = str(src)
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        with open(src) as fin, open(dst, "w") as fout:
            for line in fin:
                if line.startswith(_KEEP_PREFIXES):
                    fout.write(line)
        return (relpath, True, "ok")
    except Exception as exc:
        return (relpath, False, str(exc))


def collect_pdb_files(src_dir: Path) -> list[tuple[Path, Path]]:
    """Walk src_dir and return (src, dst) pairs preserving structure."""
    pairs: list[tuple[Path, Path]] = []
    for root, _dirs, files in os.walk(src_dir):
        for fname in files:
            if fname.endswith(".pdb"):
                src = Path(root) / fname
                pairs.append((src, src.relative_to(src_dir)))
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Clean PDB headers for mkdssp compatibility"
    )
    parser.add_argument("src_dir", type=Path, help="Source PDB directory")
    parser.add_argument("dst_dir", type=Path, help="Output directory")
    parser.add_argument(
        "-j", "--jobs", type=int, default=8,
        help="Number of parallel workers (default: 8)",
    )
    args = parser.parse_args()

    if not args.src_dir.is_dir():
        print(f"Error: {args.src_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    # Collect files
    print(f"Scanning {args.src_dir} ...")
    pairs = collect_pdb_files(args.src_dir)
    print(f"Found {len(pairs)} PDB files")

    if not pairs:
        print("Nothing to do.")
        return

    args.dst_dir.mkdir(parents=True, exist_ok=True)

    # Process
    done = 0
    errors = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(clean_pdb, src, args.dst_dir / rel): rel
            for src, rel in pairs
        }
        for future in as_completed(futures):
            relpath, success, msg = future.result()
            done += 1
            if not success:
                errors += 1
                print(f"  ERROR: {relpath}: {msg}", file=sys.stderr)
            if done % 1000 == 0 or done == len(pairs):
                print(f"  {done}/{len(pairs)} files processed", flush=True)

    print(f"\nDone: {done - errors} cleaned, {errors} errors")
    print(f"Output: {args.dst_dir}")

    # Quick size comparison
    src_size = sum(
        os.path.getsize(args.src_dir / rel)
        for _, rel in pairs
    )
    dst_size = sum(
        os.path.getsize(args.dst_dir / rel)
        for _, rel in pairs
        if (args.dst_dir / rel).exists()
    )
    print(f"Source: {src_size / 1e9:.2f} GB -> Cleaned: {dst_size / 1e9:.2f} GB "
          f"({dst_size / src_size * 100:.0f}%)")


if __name__ == "__main__":
    main()
