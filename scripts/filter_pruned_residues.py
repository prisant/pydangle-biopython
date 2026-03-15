#!/usr/bin/env python3
"""Filter pydangle JSONL output to keep only pruned residues.

Reads mask files and USER INC fragment records from pruned PDB files to
filter a JSONL file to only quality-filtered residues, and nulls out
measurements that depend on neighbor residues not in the quality set.

The USER INC records from the Richardson Lab pruned files define the
uninterrupted fragments that survived quality filtering.  A residue at
the start of a fragment has no valid i-1 neighbor; at the end, no valid
i+1 neighbor.  Measurements depending on those missing neighbors are
set to null, since the neighbor coordinates are unreliable (they failed
quality filtering).

Neighbor-dependent measurements:
    phi, omega, is_cis, is_trans  — require residue i-1
    psi                           — requires residue i+1
    tau, rama_category            — require both i-1 and i+1
    chi1, dssp                    — use only residue i (unaffected)

Usage:
    python scripts/filter_pruned_residues.py INPUT.jsonl SRC_DIR > OUTPUT.jsonl

SRC_DIR is the top2018 dataset directory containing both .mask files and
the *_pruned_all.pdb files (for USER INC fragment records).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


def _parse_inc_fragments(
    pruned_path: Path,
) -> list[tuple[str, int, str, int, str]]:
    """Parse USER INC fragment records from a pruned PDB file.

    USER INC format:
        USER  INC: C:   2: : C:   6: :3
        (chain1:resseq1:icode1:chain2:resseq2:icode2:length)

    Returns list of (chain, start_resnum, start_icode, end_resnum, end_icode).
    """
    fragments: list[tuple[str, int, str, int, str]] = []
    with open(pruned_path) as fh:
        for line in fh:
            if not line.startswith("USER  INC:"):
                continue
            # Parse: "USER  INC: C:   2: : C:   6: :3"
            # Split on ":" gives 7 parts:
            #   [0] chain1  [1] resnum1  [2] icode1
            #   [3] chain2  [4] resnum2  [5] icode2  [6] length
            parts = line[10:].split(":")
            if len(parts) >= 7:
                chain = parts[0].strip()
                start_resnum = int(parts[1].strip())
                start_icode = parts[2].strip() or " "
                end_resnum = int(parts[4].strip())
                end_icode = parts[5].strip() or " "
                fragments.append(
                    (chain, start_resnum, start_icode,
                     end_resnum, end_icode)
                )
    return fragments


def load_masks_and_fragments(
    src_dir: Path,
) -> tuple[
    dict[str, set[tuple[str, int, str]]],
    dict[str, set[tuple[str, int, str]]],
    dict[str, set[tuple[str, int, str]]],
]:
    """Load masks and parse fragment boundaries from pruned files.

    Returns:
        masks: {pdb_id: set of (chain, resnum, icode)}
        frag_starts: {pdb_id: set of (chain, resnum, icode)} — first
            residue of each fragment (no valid i-1)
        frag_ends: {pdb_id: set of (chain, resnum, icode)} — last
            residue of each fragment (no valid i+1)
    """
    masks: dict[str, set[tuple[str, int, str]]] = {}
    frag_starts: dict[str, set[tuple[str, int, str]]] = {}
    frag_ends: dict[str, set[tuple[str, int, str]]] = {}

    for mask_file in src_dir.rglob("*.mask"):
        stem = mask_file.stem  # e.g., "1a2z_C"
        pdb_id = stem.rsplit("_", 1)[0]
        residues: set[tuple[str, int, str]] = set()
        with open(mask_file) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) >= 3:
                    residues.add((parts[0], int(parts[1]), parts[2]))
                elif len(parts) == 2:
                    residues.add((parts[0], int(parts[1]), " "))
        masks[pdb_id] = residues

        # Find matching pruned file in same directory
        entry_dir = mask_file.parent
        pruned_files = list(entry_dir.glob(f"{pdb_id}_*_pruned_all.pdb"))
        starts: set[tuple[str, int, str]] = set()
        ends: set[tuple[str, int, str]] = set()
        if pruned_files:
            fragments = _parse_inc_fragments(pruned_files[0])
            for chain, s_res, s_ic, e_res, e_ic in fragments:
                starts.add((chain, s_res, s_ic))
                ends.add((chain, e_res, e_ic))
        frag_starts[pdb_id] = starts
        frag_ends[pdb_id] = ends

    return masks, frag_starts, frag_ends


# Measurements that require the previous residue (i-1)
_NEEDS_PREV = {"phi", "omega", "is_cis", "is_trans"}

# Measurements that require the next residue (i+1)
_NEEDS_NEXT = {"psi"}

# Measurements that require both neighbors (phi needs i-1, psi needs i+1)
_NEEDS_BOTH = {"rama_category"}


def null_neighbor_dependent(
    record: dict,
    is_frag_start: bool,
    is_frag_end: bool,
) -> dict:
    """Null out measurements that depend on missing neighbor residues.

    A residue at a fragment start has no valid i-1 neighbor (it was
    quality-filtered).  A residue at a fragment end has no valid i+1.
    """
    if is_frag_start:
        for col in _NEEDS_PREV:
            if col in record:
                record[col] = None
    if is_frag_end:
        for col in _NEEDS_NEXT:
            if col in record:
                record[col] = None
    if is_frag_start or is_frag_end:
        for col in _NEEDS_BOTH:
            if col in record:
                record[col] = None
    return record


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter JSONL to pruned residues using mask files"
    )
    parser.add_argument(
        "input_jsonl", type=Path, help="Input JSONL file"
    )
    parser.add_argument(
        "src_dir", type=Path,
        help="Dataset directory with .mask and pruned PDB files",
    )
    args = parser.parse_args()

    if not args.input_jsonl.exists():
        print(f"Error: {args.input_jsonl} not found", file=sys.stderr)
        sys.exit(1)
    if not args.src_dir.is_dir():
        print(
            f"Error: {args.src_dir} is not a directory", file=sys.stderr
        )
        sys.exit(1)

    masks, frag_starts, frag_ends = load_masks_and_fragments(args.src_dir)
    print(
        f"Loaded masks for {len(masks)} PDB entries", file=sys.stderr
    )

    kept = 0
    skipped = 0
    no_mask = 0
    nulled_start = 0
    nulled_end = 0

    # Compile filename normalization pattern
    ersatz_re = re.compile(r"_ersatz\.pdb$|\.pdb$")

    with open(args.input_jsonl) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            filename = record.get("file", "")
            pdb_id = ersatz_re.sub("", filename)

            mask = masks.get(pdb_id)
            if mask is None:
                no_mask += 1
                continue

            chain_id = record.get("chain", "")
            resnum = record.get("resnum", 0)
            icode = record.get("ins", " ")
            res_key = (chain_id, resnum, icode)

            if res_key not in mask:
                skipped += 1
                continue

            # Check fragment boundaries from USER INC records
            is_start = res_key in frag_starts.get(pdb_id, set())
            is_end = res_key in frag_ends.get(pdb_id, set())
            if is_start:
                nulled_start += 1
            if is_end:
                nulled_end += 1
            record = null_neighbor_dependent(record, is_start, is_end)

            print(json.dumps(record, separators=(",", ":")))
            kept += 1

    print(
        f"Kept {kept}, skipped {skipped}, no mask {no_mask}",
        file=sys.stderr,
    )
    print(
        f"Fragment starts (no i-1): {nulled_start}, "
        f"fragment ends (no i+1): {nulled_end}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
