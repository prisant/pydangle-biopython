#!/usr/bin/env python3
"""Compare two pydangle JSONL output files for regression testing.

Reports per-column match/difference counts and samples of differing
records for manual inspection.  Designed to compare the ersatz-based
top2018 output against the previous pruned-file output.

Usage:
    python scripts/regression_compare.py OLD.jsonl.gz NEW.jsonl \
        [--sample N] [--columns COL1,COL2,...]

Expected results for the top2018 ersatz comparison:
  - Backbone geometry (phi, psi, omega, tau) should match exactly
  - rama_category, is_cis, is_trans, is_left, is_right should match exactly
  - dssp may differ (multi-chain context changes H-bond assignments)
  - chi1 may differ for NQH-flipped Asn/Gln/His residues
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path


def open_jsonl(path: Path):
    """Open a JSONL file, handling .gz compression."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def _normalize_filename(fname: str) -> str:
    """Normalize a filename to just the PDB ID for cross-dataset matching.

    Handles patterns like:
        1a2z_ersatz.pdb      -> 1a2z
        1a2z_C_pruned_all.pdb -> 1a2z
        1a2z.pdb             -> 1a2z
    """
    import re

    basename = fname.replace(".pdb", "")
    basename = basename.replace("_ersatz", "")
    basename = re.sub(r"_[A-Za-z]_pruned_(all|mc|full)", "", basename)
    return basename


def record_key(rec: dict) -> tuple:
    """Generate a unique key for a JSONL record."""
    return (
        _normalize_filename(rec.get("file", "")),
        rec.get("model", 0),
        rec.get("chain", ""),
        rec.get("resnum", 0),
        rec.get("ins", " "),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare two pydangle JSONL files for regression testing"
    )
    parser.add_argument("old_file", type=Path, help="Old/reference JSONL file")
    parser.add_argument("new_file", type=Path, help="New JSONL file to compare")
    parser.add_argument(
        "--sample", type=int, default=5,
        help="Number of sample differences to show per column (default: 5)",
    )
    parser.add_argument(
        "--columns", type=str, default="",
        help="Comma-separated list of columns to compare (default: all)",
    )
    args = parser.parse_args()

    # Load old records into dict by key
    print(f"Loading {args.old_file} ...", file=sys.stderr)
    old_records: dict[tuple, dict] = {}
    with open_jsonl(args.old_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("_meta"):
                continue
            # Normalize filename: strip _ersatz suffix for matching
            fname = rec.get("file", "")
            rec["file"] = _normalize_filename(fname)
            old_records[record_key(rec)] = rec

    print(f"Loading {args.new_file} ...", file=sys.stderr)
    new_records: dict[tuple, dict] = {}
    with open_jsonl(args.new_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("_meta"):
                continue
            fname = rec.get("file", "")
            rec["file"] = _normalize_filename(fname)
            new_records[record_key(rec)] = rec

    print(f"Old: {len(old_records)} records", file=sys.stderr)
    print(f"New: {len(new_records)} records", file=sys.stderr)

    # Find common keys
    old_keys = set(old_records.keys())
    new_keys = set(new_records.keys())
    common = old_keys & new_keys
    only_old = old_keys - new_keys
    only_new = new_keys - old_keys

    print("\n=== Record counts ===")
    print(f"Common:   {len(common)}")
    print(f"Only old: {len(only_old)}")
    print(f"Only new: {len(only_new)}")

    if only_old:
        print("\nSample records only in old:")
        for key in sorted(only_old)[:args.sample]:
            print(f"  {key}")

    if only_new:
        print("\nSample records only in new:")
        for key in sorted(only_new)[:args.sample]:
            print(f"  {key}")

    # Compare common records column by column
    # Identify all measurement columns (exclude metadata fields)
    meta_fields = {"file", "model", "chain", "resnum", "ins", "resname"}
    if args.columns:
        compare_cols = args.columns.split(",")
    else:
        all_cols: set[str] = set()
        for rec in list(old_records.values())[:100]:
            all_cols.update(rec.keys())
        for rec in list(new_records.values())[:100]:
            all_cols.update(rec.keys())
        compare_cols = sorted(all_cols - meta_fields)

    print(f"\n=== Per-column comparison ({len(common)} common records) ===")
    col_diffs: dict[str, list[tuple]] = {}
    match_counts: Counter[str] = Counter()
    diff_counts: Counter[str] = Counter()
    null_counts: Counter[str] = Counter()

    for key in common:
        old_rec = old_records[key]
        new_rec = new_records[key]
        for col in compare_cols:
            old_val = old_rec.get(col)
            new_val = new_rec.get(col)
            if old_val == new_val:
                match_counts[col] += 1
            elif old_val is None and new_val is None:
                null_counts[col] += 1
            else:
                diff_counts[col] += 1
                if col not in col_diffs:
                    col_diffs[col] = []
                if len(col_diffs[col]) < args.sample:
                    col_diffs[col].append((key, old_val, new_val))

    for col in compare_cols:
        m = match_counts[col]
        d = diff_counts[col]
        n = null_counts[col]
        pct = d / len(common) * 100 if common else 0
        status = "MATCH" if d == 0 else "DIFF"
        print(f"  {col:20s}: {status:5s}  match={m}  diff={d} ({pct:.2f}%)  null={n}")

    # Show sample differences
    for col in compare_cols:
        if col in col_diffs:
            print(f"\n--- Sample differences for '{col}' ---")
            for key, old_val, new_val in col_diffs[col]:
                file_label, model, chain, resnum, ins = key
                resname_old = old_records.get(key, {}).get("resname", "?")
                print(
                    f"  {file_label} {chain}{resnum}{ins.strip()} "
                    f"{resname_old}: {old_val!r} -> {new_val!r}"
                )

    # Summary
    total_diffs = sum(diff_counts.values())
    total_comparisons = len(common) * len(compare_cols)
    print("\n=== Summary ===")
    print(f"Total comparisons: {total_comparisons}")
    print(f"Total differences: {total_diffs}")
    if total_comparisons > 0:
        print(f"Overall match rate: {(1 - total_diffs/total_comparisons)*100:.4f}%")

    # Check if all diffs are explainable
    expected_diff_cols = {"dssp", "dssp3", "chi1"}
    unexpected_diffs = {
        col: diff_counts[col]
        for col in compare_cols
        if diff_counts[col] > 0 and col not in expected_diff_cols
    }
    if unexpected_diffs:
        print(f"\nWARNING: Unexpected differences in: {unexpected_diffs}")
    elif total_diffs > 0:
        print(f"\nAll differences are in expected columns: "
              f"{', '.join(c for c in compare_cols if diff_counts[c] > 0)}")
    else:
        print("\nPerfect match across all columns.")


if __name__ == "__main__":
    main()
