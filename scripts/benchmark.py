#!/usr/bin/env python3
"""Performance and correctness benchmark: pydangle-biopython vs Java Dangle.

Runs both tools on the top100 PDB benchmark set and produces a report
with timing comparison, per-file breakdown, and correctness diff.
"""

import argparse
import glob
import os
import statistics
import subprocess
import sys
import time

DEFAULT_PDB_DIR = os.path.expanduser(
    "~/Desktop/MyProgs/PolarWork/top100pdb"
)
PYDANGLE_MEASUREMENTS = "phi; psi; omega; tau"
DANGLE_MEASUREMENTS = "phi, psi, omega, tau"
DEFAULT_TOLERANCE = 0.01
BATCH_REPEATS = 3


# ---------------------------------------------------------------------------
# Running tools
# ---------------------------------------------------------------------------


def run_pydangle_batch(pdb_dir: str) -> tuple[float, str]:
    """Run pydangle on entire directory, return (elapsed_seconds, stdout)."""
    t0 = time.perf_counter()
    result = subprocess.run(
        ["pydangle-biopython", "-c", PYDANGLE_MEASUREMENTS, "-d", pdb_dir],
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - t0
    if result.returncode != 0:
        print(f"pydangle batch failed: {result.stderr}", file=sys.stderr)
    return elapsed, result.stdout


def run_dangle_batch(pdb_files: list[str]) -> tuple[float, str]:
    """Run dangle on all files, return (elapsed_seconds, stdout)."""
    t0 = time.perf_counter()
    result = subprocess.run(
        ["dangle", DANGLE_MEASUREMENTS] + pdb_files,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - t0
    if result.returncode != 0:
        print(f"dangle batch failed: {result.stderr}", file=sys.stderr)
    return elapsed, result.stdout


def run_single(tool: str, filepath: str) -> tuple[float, str]:
    """Run a single file through one tool, return (elapsed, stdout)."""
    if tool == "pydangle":
        cmd = ["pydangle-biopython", "-c", PYDANGLE_MEASUREMENTS, filepath]
    else:
        cmd = ["dangle", DANGLE_MEASUREMENTS, filepath]
    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    return elapsed, result.stdout


# ---------------------------------------------------------------------------
# Parsing output
# ---------------------------------------------------------------------------


def parse_output(raw: str) -> dict[tuple[str, ...], list[str]]:
    """Parse colon-separated output into {row_key: [values]}.

    Row key is (filename, model, chain, resnum, ins, resname) with
    chain ID whitespace-stripped for normalization.
    """
    rows: dict[tuple[str, ...], list[str]] = {}
    for line in raw.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split(":")
        if len(fields) < 7:
            continue
        # fields: file, model, chain, num, ins, resname, val1, val2, ...
        key = (
            fields[0].strip(),           # filename
            fields[1].strip(),           # model
            fields[2].strip(),           # chain (normalize whitespace)
            fields[3].strip(),           # resnum
            fields[4],                   # ins (preserve single space)
            fields[5].strip(),           # resname
        )
        values = [f.strip() for f in fields[6:]]
        rows[key] = values
    return rows


# ---------------------------------------------------------------------------
# Correctness comparison
# ---------------------------------------------------------------------------


def compare_outputs(
    pydangle_rows: dict[tuple[str, ...], list[str]],
    dangle_rows: dict[tuple[str, ...], list[str]],
    tolerance: float,
) -> dict[str, object]:
    """Compare parsed outputs and return comparison statistics."""
    pd_keys = set(pydangle_rows.keys())
    dg_keys = set(dangle_rows.keys())

    common = pd_keys & dg_keys
    only_pydangle = pd_keys - dg_keys
    only_dangle = dg_keys - pd_keys

    matches = 0
    mismatches: list[tuple[tuple[str, ...], list[str], list[str]]] = []

    for key in sorted(common):
        pd_vals = pydangle_rows[key]
        dg_vals = dangle_rows[key]

        row_ok = True
        n = min(len(pd_vals), len(dg_vals))
        for i in range(n):
            pv, dv = pd_vals[i], dg_vals[i]
            if pv == "__?__" and dv == "__?__":
                continue
            if pv == "__?__" or dv == "__?__":
                row_ok = False
                break
            try:
                if abs(float(pv) - float(dv)) > tolerance:
                    row_ok = False
                    break
            except ValueError:
                if pv != dv:
                    row_ok = False
                    break

        if row_ok:
            matches += 1
        else:
            mismatches.append((key, pd_vals, dg_vals))

    return {
        "common": len(common),
        "matches": matches,
        "mismatches": mismatches,
        "only_pydangle": sorted(only_pydangle),
        "only_dangle": sorted(only_dangle),
    }


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def print_report(
    batch_times: dict[str, list[float]],
    per_file_times: dict[str, dict[str, float]],
    comparison: dict[str, object],
    tolerance: float,
) -> None:
    """Print the full benchmark report."""
    sep = "=" * 70

    # -- Batch timing --
    print(sep)
    print("BATCH TIMING (full top100 set)")
    print(sep)
    for tool in ["pydangle", "dangle"]:
        times = batch_times[tool]
        print(
            f"  {tool:12s}  "
            f"mean={statistics.mean(times):7.3f}s  "
            f"min={min(times):7.3f}s  "
            f"max={max(times):7.3f}s  "
            f"(n={len(times)})"
        )
    pd_mean = statistics.mean(batch_times["pydangle"])
    dg_mean = statistics.mean(batch_times["dangle"])
    if dg_mean > 0:
        ratio = pd_mean / dg_mean
        print(f"\n  pydangle / dangle = {ratio:.2f}x")
    print()

    # -- Per-file timing --
    print(sep)
    print("PER-FILE TIMING")
    print(sep)
    for tool in ["pydangle", "dangle"]:
        times_list = list(per_file_times[tool].values())
        if not times_list:
            continue
        print(f"\n  {tool}:")
        print(
            f"    mean={statistics.mean(times_list):.4f}s  "
            f"median={statistics.median(times_list):.4f}s  "
            f"stdev={statistics.stdev(times_list):.4f}s"
            if len(times_list) > 1
            else f"    mean={statistics.mean(times_list):.4f}s"
        )
        ranked = sorted(
            per_file_times[tool].items(), key=lambda x: x[1], reverse=True
        )
        print("    5 slowest:")
        for fname, t in ranked[:5]:
            print(f"      {os.path.basename(fname):20s} {t:.4f}s")
    print()

    # -- Correctness --
    print(sep)
    print(f"CORRECTNESS (tolerance = {tolerance}°)")
    print(sep)
    common = comparison["common"]
    matches = comparison["matches"]
    mismatches = comparison["mismatches"]
    only_pd = comparison["only_pydangle"]
    only_dg = comparison["only_dangle"]

    assert isinstance(common, int)
    assert isinstance(matches, int)
    assert isinstance(mismatches, list)
    assert isinstance(only_pd, list)
    assert isinstance(only_dg, list)

    print(f"  Rows in common:       {common}")
    print(f"  Rows agreeing:        {matches} / {common}")
    print(f"  Rows only in pydangle: {len(only_pd)}")
    print(f"  Rows only in dangle:   {len(only_dg)}")

    if mismatches:
        print(f"\n  VALUE MISMATCHES ({len(mismatches)}):")
        for key, pd_vals, dg_vals in mismatches[:20]:
            label = ":".join(key)
            print(f"    {label}")
            print(f"      pydangle: {':'.join(pd_vals)}")
            print(f"      dangle:   {':'.join(dg_vals)}")
        if len(mismatches) > 20:
            print(f"    ... and {len(mismatches) - 20} more")

    if only_pd:
        print(f"\n  ROWS ONLY IN PYDANGLE ({len(only_pd)}):")
        for key in only_pd[:10]:
            print(f"    {':'.join(key)}")
        if len(only_pd) > 10:
            print(f"    ... and {len(only_pd) - 10} more")

    if only_dg:
        print(f"\n  ROWS ONLY IN DANGLE ({len(only_dg)}):")
        for key in only_dg[:10]:
            print(f"    {':'.join(key)}")
        if len(only_dg) > 10:
            print(f"    ... and {len(only_dg) - 10} more")

    print()
    print(sep)
    if matches == common and not mismatches:
        print("RESULT: All common rows agree within tolerance.")
    else:
        print(f"RESULT: {len(mismatches)} mismatches found.")
    print(sep)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Benchmark pydangle-biopython vs Java Dangle.",
    )
    ap.add_argument(
        "--pdb-dir",
        default=DEFAULT_PDB_DIR,
        help=f"Directory with PDB files (default: {DEFAULT_PDB_DIR})",
    )
    ap.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help=f"Numeric comparison tolerance in degrees (default: {DEFAULT_TOLERANCE})",
    )
    ap.add_argument(
        "--skip-per-file",
        action="store_true",
        help="Skip per-file timing (much faster).",
    )
    args = ap.parse_args()

    pdb_dir = args.pdb_dir
    if not os.path.isdir(pdb_dir):
        print(f"Error: directory not found: {pdb_dir}", file=sys.stderr)
        sys.exit(1)

    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
    if not pdb_files:
        print(f"Error: no .pdb files in {pdb_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Benchmark: {len(pdb_files)} PDB files from {pdb_dir}")
    print(f"Measurements: {PYDANGLE_MEASUREMENTS}")
    print()

    # -- Batch timing --
    batch_times: dict[str, list[float]] = {"pydangle": [], "dangle": []}
    last_output: dict[str, str] = {}

    for i in range(BATCH_REPEATS):
        print(f"  Batch run {i + 1}/{BATCH_REPEATS} ...", end="", flush=True)
        t, out = run_pydangle_batch(pdb_dir)
        batch_times["pydangle"].append(t)
        last_output["pydangle"] = out
        print(f" pydangle={t:.3f}s", end="", flush=True)

        t, out = run_dangle_batch(pdb_files)
        batch_times["dangle"].append(t)
        last_output["dangle"] = out
        print(f"  dangle={t:.3f}s")

    print()

    # -- Per-file timing --
    per_file_times: dict[str, dict[str, float]] = {
        "pydangle": {},
        "dangle": {},
    }
    if not args.skip_per_file:
        total = len(pdb_files)
        for idx, f in enumerate(pdb_files):
            if sys.stderr.isatty():
                print(
                    f"\r  Per-file timing: {idx + 1}/{total}",
                    end="",
                    file=sys.stderr,
                )
            t, _ = run_single("pydangle", f)
            per_file_times["pydangle"][f] = t
            t, _ = run_single("dangle", f)
            per_file_times["dangle"][f] = t
        if sys.stderr.isatty():
            print("", file=sys.stderr)
        print()

    # -- Correctness --
    pd_rows = parse_output(last_output["pydangle"])
    dg_rows = parse_output(last_output["dangle"])
    comparison = compare_outputs(pd_rows, dg_rows, args.tolerance)

    # -- Report --
    print_report(batch_times, per_file_times, comparison, args.tolerance)


if __name__ == "__main__":
    main()
