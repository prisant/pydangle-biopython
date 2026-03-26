#!/usr/bin/env python3
"""Build residue masks and fragment records using mainchain B-factor filtering.

For datasets without Richardson Lab pruned files (e.g., top8000), this script
implements the simple Ramachandran quality filter described in Williams et al.
(2018, doi:10.1002/pro.3330): exclude any residue where a mainchain atom has
B-factor > threshold (default 30).

The output is compatible with filter_pruned_residues.py: per-chain .mask files
and synthetic _pruned_bfilt.pdb files containing USER INC fragment records
derived from contiguous runs of passing residues.

Usage:
    python scripts/build_bfactor_masks.py SRC_DIR CHAIN_LIST [-j JOBS] [-b BMAX]

Arguments:
    SRC_DIR      Dataset directory with two-level structure (prefix/pdbid/)
                 containing {pdbid}_ersatz.pdb files
    CHAIN_LIST   CSV file with columns pdb_id,chain
    -j, --jobs   Number of parallel workers (default: 8)
    -b, --bmax   Maximum mainchain B-factor (default: 30.0)

Output per chain (in entry directory):
    {pdbid}_{chain}.mask           — residue list (chain resnum icode)
    {pdbid}_{chain}_pruned_bfilt.pdb — synthetic file with USER INC records

Example:
    python scripts/build_bfactor_masks.py \\
        ~/Desktop/Data/top8000_pdbs_hom70 \\
        ~/Desktop/Data/top8000_pdbs_hom70/Top8000-best_hom70_pdb_chain_list.csv \\
        -j 8 -b 30
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Mainchain atom names (PDB columns 13-16, space-padded)
_MC_ATOMS = {" N  ", " CA ", " C  ", " O  "}


def _parse_residues(
    ersatz_path: Path, chain_id: str, bmax: float,
) -> tuple[
    list[tuple[str, int, str]],   # all residues in chain
    list[tuple[str, int, str]],   # passing residues
]:
    """Parse ersatz PDB and classify residues by mainchain B-factor.

    Returns (all_residues, passing_residues) where each is a list of
    (chain, resnum, icode) tuples in file order.
    """
    # Collect per-residue mainchain B-factors
    # Key: (chain, resnum, icode) -> list of B-factors
    mc_bfactors: dict[tuple[str, int, str], list[float]] = defaultdict(list)
    seen_order: list[tuple[str, int, str]] = []
    seen_set: set[tuple[str, int, str]] = set()

    with open(ersatz_path) as fh:
        for line in fh:
            if not line.startswith("ATOM  "):
                continue
            ch = line[21]
            if ch != chain_id:
                continue
            resnum = int(line[22:26])
            icode = line[26] if line[26] != " " else " "
            atom_name = line[12:16]
            key = (ch, resnum, icode)

            if key not in seen_set:
                seen_set.add(key)
                seen_order.append(key)

            if atom_name in _MC_ATOMS:
                try:
                    bfactor = float(line[60:66])
                    mc_bfactors[key].append(bfactor)
                except (ValueError, IndexError):
                    pass

    # Classify: pass if ALL mainchain atoms have B <= bmax
    # Also require at least N, CA, C (3 mainchain atoms) to be present
    passing = []
    for key in seen_order:
        bvals = mc_bfactors.get(key, [])
        if len(bvals) >= 3 and all(b <= bmax for b in bvals):
            passing.append(key)

    return seen_order, passing


def _find_fragments(
    all_residues: list[tuple[str, int, str]],
    passing: set[tuple[str, int, str]],
) -> list[tuple[str, int, str, int, str]]:
    """Find contiguous fragments of passing residues.

    A fragment is a maximal run of consecutive residues (by position in
    the PDB file) that all pass the B-factor filter.

    Returns list of (chain, start_resnum, start_icode, end_resnum, end_icode).
    """
    fragments: list[tuple[str, int, str, int, str]] = []

    # Walk through all residues in file order, tracking runs of passing
    current_start = None
    current_end = None

    for key in all_residues:
        if key in passing:
            if current_start is None:
                current_start = key
            current_end = key
        else:
            # Gap: emit current fragment if any
            if current_start is not None:
                fragments.append((
                    current_start[0], current_start[1], current_start[2],
                    current_end[1], current_end[2],
                ))
                current_start = None
                current_end = None

    # Emit final fragment
    if current_start is not None:
        fragments.append((
            current_start[0], current_start[1], current_start[2],
            current_end[1], current_end[2],
        ))

    return fragments


def _write_synthetic_pruned(
    output_path: Path,
    fragments: list[tuple[str, int, str, int, str]],
    passing: list[tuple[str, int, str]],
) -> None:
    """Write synthetic pruned file with USER INC fragment records.

    Format matches Richardson Lab convention:
        USER  INC: C:   2: : C:   6: :5
    """
    with open(output_path, "w") as fh:
        fh.write(
            "USER  MOD synthetic B-factor filter "
            "(no Reduce, no quality pruning)\n"
        )
        for chain, s_res, s_ic, e_res, e_ic in fragments:
            # Count residues in this fragment
            in_frag = False
            count = 0
            for key in passing:
                if not in_frag:
                    if (key[0] == chain and key[1] == s_res
                            and key[2] == s_ic):
                        in_frag = True
                        count = 1
                elif key[0] == chain and key[1] == e_res and key[2] == e_ic:
                    count += 1
                    break
                else:
                    count += 1

            s_ic_str = s_ic if s_ic.strip() else " "
            e_ic_str = e_ic if e_ic.strip() else " "
            fh.write(
                f"USER  INC: {chain}:{s_res:>4}:{s_ic_str}:"
                f" {chain}:{e_res:>4}:{e_ic_str}:{count}\n"
            )


def process_entry(
    entry_dir: Path,
    pdb_id: str,
    chains: list[str],
    bmax: float,
) -> tuple[str, bool, str, int]:
    """Process one entry: build masks and synthetic pruned files.

    Returns (pdb_id, success, message, n_mask_residues).
    """
    ersatz = entry_dir / f"{pdb_id}_ersatz.pdb"
    if not ersatz.exists():
        return (pdb_id, False, "no ersatz file", 0)

    total_mask = 0
    chain_info: list[str] = []

    try:
        for chain_id in chains:
            all_res, passing = _parse_residues(ersatz, chain_id, bmax)

            if not passing:
                chain_info.append(f"{chain_id}:0/{len(all_res)}")
                continue

            # Write mask file
            mask_path = entry_dir / f"{pdb_id}_{chain_id}.mask"
            with open(mask_path, "w") as fh:
                for ch, resnum, icode in passing:
                    fh.write(f"{ch} {resnum} {icode}\n")

            # Find fragments and write synthetic pruned file
            passing_set = set(passing)
            fragments = _find_fragments(all_res, passing_set)

            pruned_path = (
                entry_dir / f"{pdb_id}_{chain_id}_pruned_bfilt.pdb"
            )
            _write_synthetic_pruned(pruned_path, fragments, passing)

            total_mask += len(passing)
            chain_info.append(
                f"{chain_id}:{len(passing)}/{len(all_res)}"
                f"({len(fragments)} frags)"
            )

        msg = f"ok: {', '.join(chain_info)}"
        return (pdb_id, True, msg, total_mask)

    except Exception as exc:
        return (pdb_id, False, str(exc), 0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build residue masks using mainchain B-factor filter. "
            "Creates .mask files and synthetic pruned files with "
            "USER INC fragment records."
        )
    )
    parser.add_argument(
        "src_dir", type=Path,
        help="Dataset directory (two-level: prefix/pdbid/)",
    )
    parser.add_argument(
        "chain_list", type=Path,
        help="CSV file with columns pdb_id,chain",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=8,
        help="Number of parallel workers (default: 8)",
    )
    parser.add_argument(
        "-b", "--bmax", type=float, default=30.0,
        help="Maximum mainchain B-factor (default: 30.0)",
    )
    args = parser.parse_args()

    if not args.src_dir.is_dir():
        print(f"Error: {args.src_dir} is not a directory", file=sys.stderr)
        sys.exit(1)
    if not args.chain_list.exists():
        print(f"Error: {args.chain_list} not found", file=sys.stderr)
        sys.exit(1)

    # Load chain list: group chains by pdb_id
    entries: dict[str, list[str]] = defaultdict(list)
    with open(args.chain_list) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pdb_id = row["pdb_id"].lower().strip()
            chain = row["chain"].strip()
            entries[pdb_id].append(chain)

    total_chains = sum(len(v) for v in entries.values())
    print(f"Chain list: {total_chains} chains from {len(entries)} structures")
    print(f"B-factor threshold: {args.bmax}")
    print(f"Source: {args.src_dir}")

    done = 0
    errors = 0
    total_mask = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {}
        for pdb_id, chains in entries.items():
            prefix = pdb_id[:2]
            entry_dir = args.src_dir / prefix / pdb_id
            futures[pool.submit(
                process_entry, entry_dir, pdb_id, chains, args.bmax
            )] = pdb_id

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
    print(f"Total passing residues: {total_mask}")


if __name__ == "__main__":
    main()
