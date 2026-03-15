#!/usr/bin/env python3
"""Compare NQH flip decisions between our Reduce run and the Zenodo pruned files.

For each entry in the top2018 dataset, parses USER MOD records from both
the locally-generated FH file and the downloaded pruned file, extracts
all NQH flip/keep decisions, and reports matches and discrepancies.

Usage:
    python scripts/compare_flips.py SRC_DIR [-j JOBS] [--verbose]

Example:
    python scripts/compare_flips.py ~/Desktop/top2018_pdbs_full_filtered_hom70

Output:
    Summary statistics plus per-entry details for any mismatches.
"""

from __future__ import annotations

import argparse
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Regex to parse NQH evaluation records from Reduce USER MOD lines.
# Captures: chain, resnum, restype, decision (F or K)
#
# Examples matched:
#   USER  MOD Single : C  34 ASN     :FLIP  amide:sc=  0.597  F(o=-0.26,f=0.6)
#   USER  MOD Set 2.2: D 174 GLN     :      amide:sc=  0.898  K(o=3.3,f=-2.4!)
#   USER  MOD Set 7.2: D 152 HIS     :FLIP no HD1:sc=  0.898  F(o=-4.5!,f=2.1)
#   USER  MOD Set12.2: C 152 HIS     :     no HD1:sc=  0.662  K(o=3.1,f=-6.2!)
#
# The pattern:
#   - "USER  MOD" prefix with Set/Single identifier ending in ":"
#   - chain (single letter), resnum (digits), restype (ASN/GLN/HIS)
#   - eventually "F(" or "K(" near end of line
_NQH_RE = re.compile(
    r"USER  MOD .+?:\s+(\w)\s+(\d+)\s+(ASN|GLN|HIS)\s+:.+\b([FK])\("
)


def parse_nqh_decisions(
    pdb_path: Path,
) -> dict[tuple[str, int, str], str]:
    """Parse NQH flip/keep decisions from a PDB file's USER MOD records.

    Returns {(chain, resnum, restype): decision} where decision is 'F' or 'K'.
    """
    decisions: dict[tuple[str, int, str], str] = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("USER  MOD"):
                continue
            m = _NQH_RE.match(line)
            if m:
                chain, resnum_str, restype, decision = m.groups()
                key = (chain, int(resnum_str), restype)
                decisions[key] = decision
    return decisions


def compare_entry(
    entry_dir: Path,
) -> tuple[str, bool, str, dict[str, int]]:
    """Compare flip decisions for one entry.

    Returns (pdb_id, success, message, counts_dict).
    counts_dict keys: match_F, match_K, mismatch_FK, mismatch_KF,
                      only_ours, only_zenodo
    """
    pdb_id = entry_dir.name
    counts: dict[str, int] = {
        "match_F": 0,
        "match_K": 0,
        "mismatch_FK": 0,
        "mismatch_KF": 0,
        "only_ours": 0,
        "only_zenodo": 0,
    }

    fh_pdb = entry_dir / f"{pdb_id}FH.pdb"
    pruned_files = list(entry_dir.glob(f"{pdb_id}_*_pruned_all.pdb"))

    if not fh_pdb.exists():
        return (pdb_id, False, "missing FH file", counts)
    if not pruned_files:
        return (pdb_id, False, "missing pruned file", counts)

    pruned_pdb = pruned_files[0]

    try:
        ours = parse_nqh_decisions(fh_pdb)
        zenodo = parse_nqh_decisions(pruned_pdb)

        all_keys = set(ours.keys()) | set(zenodo.keys())
        mismatches: list[str] = []

        for key in sorted(all_keys):
            our_dec = ours.get(key)
            zen_dec = zenodo.get(key)

            if our_dec is not None and zen_dec is not None:
                if our_dec == zen_dec:
                    if our_dec == "F":
                        counts["match_F"] += 1
                    else:
                        counts["match_K"] += 1
                else:
                    if our_dec == "F" and zen_dec == "K":
                        counts["mismatch_FK"] += 1
                    else:
                        counts["mismatch_KF"] += 1
                    chain, resnum, restype = key
                    mismatches.append(
                        f"  {chain}{resnum} {restype}: "
                        f"ours={our_dec} zenodo={zen_dec}"
                    )
            elif our_dec is not None:
                counts["only_ours"] += 1
            else:
                counts["only_zenodo"] += 1

        if mismatches:
            detail = "\n".join(mismatches)
            return (pdb_id, True, f"mismatches:\n{detail}", counts)
        return (pdb_id, True, "ok", counts)

    except Exception as exc:
        return (pdb_id, False, str(exc), counts)


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
        description="Compare NQH flip decisions: our Reduce vs Zenodo pruned"
    )
    parser.add_argument(
        "src_dir", type=Path, help="Top2018 dataset directory"
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=8,
        help="Number of parallel workers (default: 8)",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Show details for every mismatched entry",
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

    # Aggregate counts
    totals: dict[str, int] = {
        "match_F": 0,
        "match_K": 0,
        "mismatch_FK": 0,
        "mismatch_KF": 0,
        "only_ours": 0,
        "only_zenodo": 0,
    }
    done = 0
    errors = 0
    mismatch_entries = 0
    mismatch_details: list[tuple[str, str]] = []

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(compare_entry, entry): entry
            for entry in entries
        }
        for future in as_completed(futures):
            pdb_id, success, msg, counts = future.result()
            done += 1

            if not success:
                errors += 1
                print(f"  ERROR {pdb_id}: {msg}", file=sys.stderr)
            else:
                for k in totals:
                    totals[k] += counts[k]
                if "mismatches" in msg:
                    mismatch_entries += 1
                    mismatch_details.append((pdb_id, msg))

            if done % 2000 == 0 or done == len(entries):
                print(
                    f"  {done}/{len(entries)} entries compared",
                    flush=True,
                )

    # Report
    total_decisions = sum(totals.values())
    total_matches = totals["match_F"] + totals["match_K"]

    print(f"\n{'=' * 60}")
    print("NQH Flip Decision Comparison")
    print(f"{'=' * 60}")
    print(f"Entries compared:    {done - errors}")
    print(f"Entries with errors: {errors}")
    print()
    print(f"Total NQH decisions: {total_decisions}")
    print(f"  Both agree FLIP:   {totals['match_F']}")
    print(f"  Both agree KEEP:   {totals['match_K']}")
    print(f"  Ours=F, Zenodo=K:  {totals['mismatch_FK']}")
    print(f"  Ours=K, Zenodo=F:  {totals['mismatch_KF']}")
    print(f"  Only in ours:      {totals['only_ours']}")
    print(f"  Only in Zenodo:    {totals['only_zenodo']}")
    print()
    if total_decisions > 0:
        pct = total_matches / total_decisions * 100
        print(f"Agreement rate: {pct:.2f}% ({total_matches}/{total_decisions})")
    print(f"Entries with mismatches: {mismatch_entries}")

    if mismatch_details:
        # Sort and show
        mismatch_details.sort()
        show = mismatch_details if args.verbose else mismatch_details[:20]
        print(f"\n--- Mismatch details "
              f"(showing {len(show)}/{len(mismatch_details)}) ---")
        for pdb_id, detail in show:
            print(f"\n{pdb_id}: {detail}")
        if not args.verbose and len(mismatch_details) > 20:
            print(f"\n... {len(mismatch_details) - 20} more entries "
                  f"with mismatches (use --verbose to see all)")


if __name__ == "__main__":
    main()
