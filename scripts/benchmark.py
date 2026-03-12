#!/usr/bin/env python3
"""Performance and correctness benchmark: pydangle-biopython vs Java Dangle.

Runs both tools on the top100 PDB benchmark set and produces a report
with timing comparison, per-file breakdown, and correctness diff.
"""

import argparse
import glob
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

DEFAULT_PDB_DIR = os.path.expanduser(
    "~/Desktop/MyProgs/PolarWork/top100pdb"
)
PYDANGLE_MEASUREMENTS = "phi; psi; omega; tau"
DANGLE_MEASUREMENTS = "phi, psi, omega, tau"
DEFAULT_TOLERANCE = 0.01
BATCH_REPEATS = 3


# ---------------------------------------------------------------------------
# PDB preprocessing for dangle compatibility
# ---------------------------------------------------------------------------


def _needs_dangle_fix(filepath: str) -> bool:
    """Return True if PDB file needs preprocessing for dangle.

    Checks for blank chain IDs or hydrogen atoms, either of which
    prevents dangle from computing backbone dihedrals.
    """
    with open(filepath) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")):
                if len(line) > 21 and line[21] == " ":
                    return True
                atom_name = line[12:16].strip()
                if atom_name.startswith(("H", "1H", "2H", "3H")):
                    return True
    return False


def _is_hydrogen(line: str) -> bool:
    """Return True if an ATOM/HETATM line is a hydrogen."""
    atom_name = line[12:16].strip()
    return atom_name.startswith(("H", "1H", "2H", "3H"))


def _fix_for_dangle(filepath: str, dest: str) -> None:
    """Copy PDB file to dest, fixing chain IDs and stripping hydrogens."""
    with open(filepath) as fin, open(dest, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM  ", "HETATM")):
                if _is_hydrogen(line):
                    continue
                if len(line) > 21 and line[21] == " ":
                    line = line[:21] + "A" + line[22:]
            fout.write(line)


def prepare_dangle_files(
    pdb_files: list[str],
) -> tuple[list[str], str | None]:
    """Prepare PDB files for dangle, fixing chain IDs and stripping hydrogens.

    Dangle cannot compute backbone dihedrals when chain IDs are blank
    (common in Reduce-processed files) or when hydrogen atoms are present
    (they confuse dangle's residue/atom matching).

    Returns (file_list_for_dangle, temp_dir_or_None).
    The caller should clean up temp_dir when done.
    """
    needs_fix = [f for f in pdb_files if _needs_dangle_fix(f)]
    if not needs_fix:
        return pdb_files, None

    tmp_dir = tempfile.mkdtemp(prefix="pydangle_bench_")
    result: list[str] = []
    fixed_set = set(needs_fix)

    for f in pdb_files:
        if f in fixed_set:
            dest = os.path.join(tmp_dir, os.path.basename(f))
            _fix_for_dangle(f, dest)
            result.append(dest)
        else:
            result.append(f)

    print(
        f"  (Preprocessed {len(needs_fix)} files for dangle compatibility:"
        " fixed blank chain IDs and stripped hydrogens)"
    )
    return result, tmp_dir


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


def _chainless_key(key: tuple[str, ...]) -> tuple[str, ...]:
    """Return a row key with the chain field removed (index 2)."""
    return (key[0], key[1], key[3], key[4], key[5])


def _match_rows(
    pydangle_rows: dict[tuple[str, ...], list[str]],
    dangle_rows: dict[tuple[str, ...], list[str]],
) -> tuple[
    list[tuple[tuple[str, ...], list[str], tuple[str, ...], list[str]]],
    list[tuple[str, ...]],
    list[tuple[str, ...]],
]:
    """Match rows between pydangle and dangle, tolerating chain ID mismatches.

    First tries exact key match.  For unmatched rows, falls back to
    matching by (filename, model, resnum, ins, resname) — ignoring
    chain ID — to handle files where dangle uses the PDB segment ID
    as chain while pydangle reports an empty chain.

    Returns (matched_pairs, only_pydangle, only_dangle) where each
    matched pair is (pd_key, pd_vals, dg_key, dg_vals).
    """
    pd_keys = set(pydangle_rows)
    dg_keys = set(dangle_rows)

    # Phase 1: exact key match
    exact_common = pd_keys & dg_keys
    pd_remaining = pd_keys - exact_common
    dg_remaining = dg_keys - exact_common

    matched: list[
        tuple[tuple[str, ...], list[str], tuple[str, ...], list[str]]
    ] = []
    for key in sorted(exact_common):
        matched.append((key, pydangle_rows[key], key, dangle_rows[key]))

    # Phase 2: chainless fallback for remaining rows
    if pd_remaining and dg_remaining:
        dg_by_chainless: dict[tuple[str, ...], tuple[str, ...]] = {}
        for key in dg_remaining:
            cl = _chainless_key(key)
            dg_by_chainless[cl] = key  # last one wins if duplicates

        still_pd: set[tuple[str, ...]] = set()
        consumed_dg: set[tuple[str, ...]] = set()
        for pd_key in sorted(pd_remaining):
            cl = _chainless_key(pd_key)
            if cl in dg_by_chainless:
                dg_key = dg_by_chainless[cl]
                matched.append((
                    pd_key,
                    pydangle_rows[pd_key],
                    dg_key,
                    dangle_rows[dg_key],
                ))
                consumed_dg.add(dg_key)
            else:
                still_pd.add(pd_key)
        dg_remaining = dg_remaining - consumed_dg
        pd_remaining = still_pd

    return matched, sorted(pd_remaining), sorted(dg_remaining)


# ---------------------------------------------------------------------------
# Correctness comparison
# ---------------------------------------------------------------------------


def _compare_values(
    pd_vals: list[str],
    dg_vals: list[str],
    tolerance: float,
) -> tuple[bool, bool]:
    """Compare value lists, return (agree, is_coverage_diff).

    agree: True if all comparable values match within tolerance.
    is_coverage_diff: True if the only differences are one tool
    reporting __?__ where the other has a value (coverage gap,
    not a computational disagreement).
    """
    agree = True
    coverage_only = True
    n = min(len(pd_vals), len(dg_vals))
    for i in range(n):
        pv, dv = pd_vals[i], dg_vals[i]
        if pv == "__?__" and dv == "__?__":
            continue
        if pv == "__?__" or dv == "__?__":
            agree = False
            continue  # coverage diff — keep checking other fields
        try:
            if abs(float(pv) - float(dv)) > tolerance:
                agree = False
                coverage_only = False
        except ValueError:
            if pv != dv:
                agree = False
                coverage_only = False
    return agree, coverage_only


def compare_outputs(
    pydangle_rows: dict[tuple[str, ...], list[str]],
    dangle_rows: dict[tuple[str, ...], list[str]],
    tolerance: float,
) -> dict[str, object]:
    """Compare parsed outputs and return comparison statistics."""
    matched, only_pydangle, only_dangle = _match_rows(
        pydangle_rows, dangle_rows
    )

    agrees = 0
    coverage_diffs: list[
        tuple[tuple[str, ...], list[str], list[str]]
    ] = []
    value_mismatches: list[
        tuple[tuple[str, ...], list[str], list[str]]
    ] = []

    for pd_key, pd_vals, _dg_key, dg_vals in matched:
        ok, coverage_only = _compare_values(pd_vals, dg_vals, tolerance)
        if ok:
            agrees += 1
        elif coverage_only:
            coverage_diffs.append((pd_key, pd_vals, dg_vals))
        else:
            value_mismatches.append((pd_key, pd_vals, dg_vals))

    return {
        "matched": len(matched),
        "agrees": agrees,
        "coverage_diffs": coverage_diffs,
        "value_mismatches": value_mismatches,
        "only_pydangle": only_pydangle,
        "only_dangle": only_dangle,
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
    matched = comparison["matched"]
    agrees = comparison["agrees"]
    coverage_diffs = comparison["coverage_diffs"]
    value_mismatches = comparison["value_mismatches"]
    only_pd = comparison["only_pydangle"]
    only_dg = comparison["only_dangle"]

    assert isinstance(matched, int)
    assert isinstance(agrees, int)
    assert isinstance(coverage_diffs, list)
    assert isinstance(value_mismatches, list)
    assert isinstance(only_pd, list)
    assert isinstance(only_dg, list)

    print(f"  Rows matched:          {matched}")
    print(f"  Rows agreeing:         {agrees} / {matched}")
    print(f"  Coverage diffs:        {len(coverage_diffs)}"
          "  (one tool has __?__ where the other computed a value)")
    print(f"  Value mismatches:      {len(value_mismatches)}"
          "  (both computed values but they disagree)")
    print(f"  Rows only in pydangle: {len(only_pd)}")
    print(f"  Rows only in dangle:   {len(only_dg)}")

    if value_mismatches:
        print(f"\n  VALUE MISMATCHES ({len(value_mismatches)}):")
        for key, pd_vals, dg_vals in value_mismatches[:20]:
            label = ":".join(key)
            print(f"    {label}")
            print(f"      pydangle: {':'.join(pd_vals)}")
            print(f"      dangle:   {':'.join(dg_vals)}")
        if len(value_mismatches) > 20:
            print(f"    ... and {len(value_mismatches) - 20} more")

    if coverage_diffs:
        print(f"\n  COVERAGE DIFFS ({len(coverage_diffs)}):")
        for key, pd_vals, dg_vals in coverage_diffs[:20]:
            label = ":".join(key)
            print(f"    {label}")
            print(f"      pydangle: {':'.join(pd_vals)}")
            print(f"      dangle:   {':'.join(dg_vals)}")
        if len(coverage_diffs) > 20:
            print(f"    ... and {len(coverage_diffs) - 20} more")

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
    if agrees == matched and not value_mismatches:
        print("RESULT: All matched rows agree within tolerance.")
    else:
        print(f"RESULT: {len(value_mismatches)} value mismatches found.")
    if coverage_diffs:
        print(f"        {len(coverage_diffs)} coverage diffs"
              " (pydangle computed values where dangle reported __?__).")
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

    # Preprocess files with blank chain IDs for dangle compatibility
    dangle_files, tmp_dir = prepare_dangle_files(pdb_files)

    try:
        _run_benchmark(pdb_dir, pdb_files, dangle_files, args)
    finally:
        if tmp_dir:
            shutil.rmtree(tmp_dir, ignore_errors=True)


def _run_benchmark(
    pdb_dir: str,
    pdb_files: list[str],
    dangle_files: list[str],
    args: argparse.Namespace,
) -> None:
    """Run the benchmark after preprocessing."""
    # -- Batch timing --
    batch_times: dict[str, list[float]] = {"pydangle": [], "dangle": []}
    last_output: dict[str, str] = {}

    for i in range(BATCH_REPEATS):
        print(f"  Batch run {i + 1}/{BATCH_REPEATS} ...", end="", flush=True)
        t, out = run_pydangle_batch(pdb_dir)
        batch_times["pydangle"].append(t)
        last_output["pydangle"] = out
        print(f" pydangle={t:.3f}s", end="", flush=True)

        t, out = run_dangle_batch(dangle_files)
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
        for idx in range(total):
            if sys.stderr.isatty():
                print(
                    f"\r  Per-file timing: {idx + 1}/{total}",
                    end="",
                    file=sys.stderr,
                )
            f = pdb_files[idx]
            t, _ = run_single("pydangle", f)
            per_file_times["pydangle"][f] = t
            t, _ = run_single("dangle", dangle_files[idx])
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
