#!/usr/bin/env python3
"""Analyze pydangle JSONL output files.

Produces summary statistics including residue composition, Ramachandran
category distribution, cis peptide bond rates, chirality, DSSP secondary
structure, geometric distributions, and data completeness.

Supports both human-readable text reports and JSON output.

Note on non-standard residues: selenomethionine (MSE) is grouped with
methionine (MET) for analysis purposes, as MSE is a standard
crystallographic isomorph used for phasing.  All other non-standard
residues are reported separately.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from typing import Any

# -------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------

#: The 20 standard amino acid 3-letter codes.
STANDARD_AA: frozenset[str] = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
})

#: Non-standard residues treated as standard equivalents.
#: MSE (selenomethionine) is a crystallographic isomorph of MET.
NONSTANDARD_MAPPING: dict[str, str] = {"MSE": "MET"}


# -------------------------------------------------------------------
# Data loading
# -------------------------------------------------------------------

def load_records(path: str) -> list[dict[str, Any]]:
    """Load all records from a JSONL file."""
    records: list[dict[str, Any]] = []
    with open(path) as f:
        for line in f:
            record = json.loads(line)
            if record.get("_meta"):
                continue
            records.append(record)
    return records


def effective_resname(resname: str) -> str:
    """Map non-standard residues to their standard equivalents."""
    return NONSTANDARD_MAPPING.get(resname, resname)


# -------------------------------------------------------------------
# Statistics helpers
# -------------------------------------------------------------------

class RunningStats:
    """Welford's online algorithm for mean and std deviation."""

    def __init__(self) -> None:
        self.n = 0
        self.mean = 0.0
        self._m2 = 0.0

    def add(self, x: float) -> None:
        self.n += 1
        delta = x - self.mean
        self.mean += delta / self.n
        delta2 = x - self.mean
        self._m2 += delta * delta2

    def std(self) -> float:
        if self.n < 2:
            return 0.0
        return math.sqrt(self._m2 / (self.n - 1))

    def as_dict(self) -> dict[str, Any]:
        return {
            "n": self.n,
            "mean": round(self.mean, 3),
            "std": round(self.std(), 3),
        }


class CircularRunningStats:
    """Online circular statistics for angular data (degrees).

    Uses sum-of-sines / sum-of-cosines accumulation.
    Circular mean = atan2(mean_sin, mean_cos).
    Circular std  = sqrt(-2 * ln(R)) where R is the mean resultant length.
    """

    def __init__(self) -> None:
        self.n = 0
        self._sum_sin = 0.0
        self._sum_cos = 0.0

    def add(self, x: float) -> None:
        self.n += 1
        rad = math.radians(x)
        self._sum_sin += math.sin(rad)
        self._sum_cos += math.cos(rad)

    def mean(self) -> float:
        if self.n == 0:
            return 0.0
        return math.degrees(
            math.atan2(self._sum_sin / self.n, self._sum_cos / self.n)
        )

    def std(self) -> float:
        if self.n < 2:
            return 0.0
        r = math.sqrt(
            (self._sum_sin / self.n) ** 2
            + (self._sum_cos / self.n) ** 2
        )
        # Clamp R to avoid log(0) when data is uniformly distributed
        if r < 1e-15:
            return 180.0  # maximum dispersion
        return math.degrees(math.sqrt(-2.0 * math.log(r)))

    def as_dict(self) -> dict[str, Any]:
        return {
            "n": self.n,
            "mean": round(self.mean(), 3),
            "std": round(self.std(), 3),
        }


def pct(count: int, total: int) -> float:
    """Percentage with safe division."""
    if total > 0:
        return round(100.0 * count / total, 4)
    return 0.0


def frac(count: int, total: int) -> float:
    """Fraction with safe division."""
    if total > 0:
        return round(count / total, 6)
    return 0.0


def sorted_counter(
    counter: Counter[str],
) -> list[tuple[str, int]]:
    """Return counter items sorted by count descending."""
    return sorted(counter.items(), key=lambda x: (-x[1], x[0]))


# -------------------------------------------------------------------
# Analysis
# -------------------------------------------------------------------

def analyze(  # noqa: C901 — intentionally monolithic
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    """Compute all summary statistics from JSONL records."""
    results: dict[str, Any] = {}
    total = len(records)

    # --- Detect available columns ---
    if not records:
        return {"error": "No records found"}
    available = set(records[0].keys())

    # --- Dataset overview ---
    files: set[str] = set()
    chains: set[tuple[str, str]] = set()
    eff_resname_counts: Counter[str] = Counter()
    standard_count = 0
    nonstandard_counts: Counter[str] = Counter()
    mapped_counts: Counter[str] = Counter()

    for rec in records:
        files.add(rec["file"])
        chains.add((rec["file"], rec.get("chain", "")))
        rn = rec["resname"]
        eff_resname_counts[effective_resname(rn)] += 1
        if rn in STANDARD_AA:
            standard_count += 1
        elif rn in NONSTANDARD_MAPPING:
            mapped_counts[rn] += 1
            standard_count += 1
        else:
            nonstandard_counts[rn] += 1

    # Count structures with multiple chains
    chains_per_file: Counter[str] = Counter()
    for file_id, chain_id in chains:
        chains_per_file[file_id] += 1
    multi_chain_structures = sum(
        1 for n in chains_per_file.values() if n > 1
    )

    results["overview"] = {
        "total_residues": total,
        "unique_structures": len(files),
        "unique_chains": len(chains),
        "multi_chain_structures": multi_chain_structures,
        "standard_residues": standard_count,
        "nonstandard_residues": sum(nonstandard_counts.values()),
        "mapped_residues": {
            k: {
                "count": v,
                "mapped_to": NONSTANDARD_MAPPING[k],
            }
            for k, v in sorted_counter(mapped_counts)
        },
        "nonstandard_types": dict(
            sorted_counter(nonstandard_counts)
        ),
    }

    # --- Residue composition (using effective names) ---
    composition: list[dict[str, Any]] = []
    for rn, count in sorted_counter(eff_resname_counts):
        composition.append({
            "resname": rn,
            "count": count,
            "fraction": frac(count, total),
            "percent": pct(count, total),
            "is_standard": rn in STANDARD_AA,
        })
    results["residue_composition"] = composition

    # --- Ramachandran category distribution ---
    if "rama_category" in available:
        rama_overall: Counter[str] = Counter()
        rama_by_rt: dict[str, Counter[str]] = defaultdict(Counter)
        for rec in records:
            cat = rec["rama_category"]
            if cat is not None:
                eff = effective_resname(rec["resname"])
                rama_overall[cat] += 1
                rama_by_rt[eff][cat] += 1

        rama_total = sum(rama_overall.values())
        results["rama_category"] = {
            "overall": {
                cat: {
                    "count": n,
                    "percent": pct(n, rama_total),
                }
                for cat, n in sorted_counter(rama_overall)
            },
            "by_residue_type": {
                rn: {
                    cat: {
                        "count": n,
                        "percent": pct(n, sum(cats.values())),
                    }
                    for cat, n in sorted_counter(cats)
                }
                for rn, cats in sorted(rama_by_rt.items())
            },
        }

    # --- Cis peptide bond analysis ---
    if "is_cis" in available:
        _analyze_cis(records, results)

    # --- Chirality analysis ---
    if "is_left" in available:
        _analyze_chirality(records, results)

    # --- DSSP distribution ---
    if "dssp" in available:
        _analyze_dssp(records, results)

    # --- Geometric distributions ---
    geom_fields = ["phi", "psi", "omega", "tau", "chi1", "chi2", "chi3", "chi4"]
    geom_avail = [f for f in geom_fields if f in available]
    if geom_avail:
        _analyze_geometry(records, geom_avail, available, results)

    # --- Data completeness ---
    check_fields = [
        "phi", "psi", "omega", "tau", "chi1",
        "is_cis", "is_left", "dssp", "rama_category",
    ]
    completeness: dict[str, dict[str, Any]] = {}
    for f in check_fields:
        if f in available:
            present = sum(
                1 for rec in records if rec.get(f) is not None
            )
            completeness[f] = {
                "present": present,
                "missing": total - present,
                "percent_present": pct(present, total),
            }
    results["data_completeness"] = completeness

    return results


def _analyze_cis(
    records: list[dict[str, Any]],
    results: dict[str, Any],
) -> None:
    """Cis peptide bond analysis."""
    cis_total = 0
    cis_count = 0
    cis_by_rt: Counter[str] = Counter()
    total_by_rt: Counter[str] = Counter()
    cis_by_dssp: Counter[str] = Counter()
    total_by_dssp: Counter[str] = Counter()

    for rec in records:
        ic = rec.get("is_cis")
        if ic is None:
            continue
        eff = effective_resname(rec["resname"])
        is_cis = (ic == "True") if isinstance(ic, str) else ic
        cis_total += 1
        if is_cis:
            cis_count += 1
            cis_by_rt[eff] += 1
        total_by_rt[eff] += 1

        dssp_val = rec.get("dssp")
        if dssp_val is not None:
            if is_cis:
                cis_by_dssp[dssp_val] += 1
            total_by_dssp[dssp_val] += 1

    cis_pro = cis_by_rt.get("PRO", 0)
    total_pro = total_by_rt.get("PRO", 0)
    cis_nonpro = cis_count - cis_pro
    total_nonpro = cis_total - total_pro

    results["cis_peptide"] = {
        "overall": {
            "cis": cis_count,
            "total": cis_total,
            "percent": pct(cis_count, cis_total),
        },
        "proline": {
            "cis": cis_pro,
            "total": total_pro,
            "percent": pct(cis_pro, total_pro),
        },
        "non_proline": {
            "cis": cis_nonpro,
            "total": total_nonpro,
            "percent": pct(cis_nonpro, total_nonpro),
        },
        "by_residue_type": {
            rn: {
                "cis": cis_by_rt.get(rn, 0),
                "total": n,
                "percent": pct(cis_by_rt.get(rn, 0), n),
            }
            for rn, n in sorted_counter(total_by_rt)
            if cis_by_rt.get(rn, 0) > 0
        },
        "by_dssp": {
            ds: {
                "cis": cis_by_dssp.get(ds, 0),
                "total": n,
                "percent": pct(cis_by_dssp.get(ds, 0), n),
            }
            for ds, n in sorted_counter(total_by_dssp)
        },
    }


def _analyze_chirality(
    records: list[dict[str, Any]],
    results: dict[str, Any],
) -> None:
    """Chirality analysis."""
    chiral_total = 0
    d_count = 0
    d_by_rt: Counter[str] = Counter()
    total_by_rt: Counter[str] = Counter()

    for rec in records:
        il = rec.get("is_left")
        if il is None:
            continue
        eff = effective_resname(rec["resname"])
        is_d = (il == "False") if isinstance(il, str) else (not il)
        chiral_total += 1
        total_by_rt[eff] += 1
        if is_d:
            d_count += 1
            d_by_rt[eff] += 1

    results["chirality"] = {
        "overall": {
            "L": chiral_total - d_count,
            "D": d_count,
            "total": chiral_total,
            "D_percent": pct(d_count, chiral_total),
        },
        "by_residue_type": {
            rn: {
                "D": d_by_rt.get(rn, 0),
                "total": n,
                "D_percent": pct(d_by_rt.get(rn, 0), n),
            }
            for rn, n in sorted_counter(total_by_rt)
            if d_by_rt.get(rn, 0) > 0
        },
    }


def _analyze_dssp(
    records: list[dict[str, Any]],
    results: dict[str, Any],
) -> None:
    """DSSP secondary structure analysis."""
    dssp_overall: Counter[str] = Counter()
    dssp_by_rama: dict[str, Counter[str]] = defaultdict(Counter)

    for rec in records:
        ds = rec.get("dssp")
        if ds is None:
            continue
        dssp_overall[ds] += 1
        rama = rec.get("rama_category")
        if rama is not None:
            dssp_by_rama[rama][ds] += 1

    dssp_total = sum(dssp_overall.values())
    results["dssp"] = {
        "overall": {
            ds: {
                "count": n,
                "percent": pct(n, dssp_total),
            }
            for ds, n in sorted_counter(dssp_overall)
        },
        "by_rama_category": {
            rama: {
                ds: {
                    "count": n,
                    "percent": pct(
                        n, sum(cats.values())
                    ),
                }
                for ds, n in sorted_counter(cats)
            }
            for rama, cats in sorted(dssp_by_rama.items())
        },
    }


def _analyze_geometry(
    records: list[dict[str, Any]],
    geom_avail: list[str],
    available: set[str],
    results: dict[str, Any],
) -> None:
    """Geometric distribution analysis."""
    # Torsion angles use circular statistics; tau (bond angle) uses linear.
    circular_fields = {"phi", "psi", "omega", "chi1", "chi2", "chi3", "chi4"}

    def _make_stats(fields: list[str]) -> dict[str, RunningStats | CircularRunningStats]:
        return {
            f: CircularRunningStats() if f in circular_fields else RunningStats()
            for f in fields
        }

    geom_overall = _make_stats(geom_avail)
    geom_by_rama: dict[str, dict[str, RunningStats | CircularRunningStats]] = (
        defaultdict(lambda: _make_stats(geom_avail))
    )
    omega_dev = RunningStats()

    for rec in records:
        rama = rec.get("rama_category")
        for f in geom_avail:
            val = rec.get(f)
            if val is not None:
                geom_overall[f].add(val)
                if rama is not None:
                    geom_by_rama[rama][f].add(val)

        # Omega deviation for trans peptides
        if "omega" in available and "is_cis" in available:
            omega_val = rec.get("omega")
            cis_val = rec.get("is_cis")
            if omega_val is not None and cis_val is not None:
                is_trans = (
                    (cis_val == "False")
                    if isinstance(cis_val, str)
                    else (not cis_val)
                )
                if is_trans:
                    dev = abs(abs(omega_val) - 180.0)
                    omega_dev.add(dev)

    results["geometry"] = {
        "overall": {
            f: geom_overall[f].as_dict()
            for f in geom_avail
        },
        "by_rama_category": {
            rama: {
                f: stats[f].as_dict()
                for f in geom_avail
                if stats[f].n > 0
            }
            for rama, stats in sorted(geom_by_rama.items())
        },
    }
    if omega_dev.n > 0:
        results["geometry"]["omega_deviation_trans"] = (
            omega_dev.as_dict()
        )


# -------------------------------------------------------------------
# Text report formatting
# -------------------------------------------------------------------

def _fmt_ratio(
    num: int, denom: int, pct_val: float,
) -> str:
    """Format 'num / denom  (pct%)' string."""
    return f"{num:>6,} / {denom:>10,}  ({pct_val:.4f}%)"


def format_text_report(stats: dict[str, Any]) -> str:
    """Format analysis results as a human-readable text report."""
    lines: list[str] = []
    sep = "=" * 70

    def section(title: str) -> None:
        lines.append("")
        lines.append(sep)
        lines.append(f"  {title}")
        lines.append(sep)
        lines.append("")

    def subsection(title: str) -> None:
        lines.append("")
        lines.append(f"  --- {title} ---")
        lines.append("")

    # --- Overview ---
    section("DATASET OVERVIEW")
    ov = stats["overview"]
    lines.append(f"  Total residues:       {ov['total_residues']:>12,}")
    n_struct = ov.get("unique_structures", ov.get("unique_files", 0))
    n_chain_ids = ov.get("unique_chains", n_struct)
    n_multi = ov.get("multi_chain_structures", 0)
    lines.append(f"  Unique structures:    {n_struct:>12,}")
    lines.append(f"  Output chain IDs:     {n_chain_ids:>12,}  [1]")
    if n_multi > 0:
        lines.append(
            f"  Multi-chain structs:  {n_multi:>12,}"
            f"  ({n_chain_ids - n_struct} extra chains)"
        )
    lines.append(f"  Standard residues:    {ov['standard_residues']:>12,}")
    lines.append(f"  Non-standard:         {ov['nonstandard_residues']:>12,}")
    if ov["mapped_residues"]:
        lines.append("")
        lines.append(
            "  Mapped non-standard residues"
            " (treated as standard):"
        )
        for rn, info in ov["mapped_residues"].items():
            mapped = info["mapped_to"]
            cnt = info["count"]
            lines.append(f"    {rn:5s} -> {mapped:5s}  ({cnt:,})")
    if ov["nonstandard_types"]:
        lines.append("")
        lines.append("  Remaining non-standard residues:")
        for rn, count in ov["nonstandard_types"].items():
            lines.append(f"    {rn:5s}  {count:>8,}")

    # --- Residue composition ---
    section("RESIDUE COMPOSITION")
    hdr = (
        f"  {'Res':>7s}  {'Count':>10s}"
        f"  {'Pct':>8s}  {'Std':>4s}"
    )
    lines.append(hdr)
    lines.append(
        f"  {'---':>7s}  {'-----':>10s}"
        f"  {'---':>8s}  {'---':>4s}"
    )
    for entry in stats["residue_composition"]:
        flag = "yes" if entry["is_standard"] else "no"
        lines.append(
            f"  {entry['resname']:>7s}"
            f"  {entry['count']:>10,}"
            f"  {entry['percent']:>7.2f}%"
            f"  {flag:>4s}"
        )

    # --- Rama category ---
    if "rama_category" in stats:
        section("RAMACHANDRAN CATEGORY DISTRIBUTION")
        rc = stats["rama_category"]
        subsection("Overall")
        lines.append(
            f"  {'Category':>10s}"
            f"  {'Count':>10s}  {'Percent':>8s}"
        )
        lines.append(
            f"  {'--------':>10s}"
            f"  {'-----':>10s}  {'-------':>8s}"
        )
        for cat, info in rc["overall"].items():
            lines.append(
                f"  {cat:>10s}"
                f"  {info['count']:>10,}"
                f"  {info['percent']:>7.2f}%"
            )

        subsection("By residue type (top 20)")
        shown = 0
        for rn in sorted(rc["by_residue_type"].keys()):
            if rn not in STANDARD_AA or shown >= 20:
                continue
            shown += 1
            cats = rc["by_residue_type"][rn]
            parts = [
                f"{c}:{info['percent']:.1f}%"
                for c, info in cats.items()
            ]
            lines.append(f"  {rn:5s}  {'  '.join(parts)}")

    # --- Cis peptide ---
    if "cis_peptide" in stats:
        section("CIS PEPTIDE BONDS")
        cp = stats["cis_peptide"]
        ov = cp["overall"]
        pr = cp["proline"]
        np_ = cp["non_proline"]
        lines.append(
            f"  Overall:      {_fmt_ratio(ov['cis'], ov['total'], ov['percent'])}"
        )
        lines.append(
            f"  Proline:      {_fmt_ratio(pr['cis'], pr['total'], pr['percent'])}"
        )
        lines.append(
            f"  Non-proline:  {_fmt_ratio(np_['cis'], np_['total'], np_['percent'])}"
        )

        subsection("Cis by residue type (non-zero only)")
        lines.append(
            f"  {'Res':>7s}  {'Cis':>6s}"
            f"  {'Total':>10s}  {'Pct':>8s}"
        )
        lines.append(
            f"  {'---':>7s}  {'---':>6s}"
            f"  {'-----':>10s}  {'---':>8s}"
        )
        for rn, info in cp["by_residue_type"].items():
            lines.append(
                f"  {rn:>7s}  {info['cis']:>6,}"
                f"  {info['total']:>10,}"
                f"  {info['percent']:>7.4f}%"
            )

        if "by_dssp" in cp:
            subsection("Cis by DSSP state")
            lines.append(
                f"  {'DSSP':>6s}  {'Cis':>6s}"
                f"  {'Total':>10s}  {'Pct':>8s}"
            )
            lines.append(
                f"  {'----':>6s}  {'---':>6s}"
                f"  {'-----':>10s}  {'---':>8s}"
            )
            for ds, info in cp["by_dssp"].items():
                lines.append(
                    f"  {ds:>6s}  {info['cis']:>6,}"
                    f"  {info['total']:>10,}"
                    f"  {info['percent']:>7.4f}%"
                )

    # --- Chirality ---
    if "chirality" in stats:
        section("CHIRALITY")
        ch = stats["chirality"]
        ov = ch["overall"]
        lines.append(f"  L-amino acids:  {ov['L']:>10,}")
        lines.append(f"  D-amino acids:  {ov['D']:>10,}")
        lines.append(f"  Total:          {ov['total']:>10,}")
        lines.append(f"  D-amino acid %: {ov['D_percent']:>10.4f}%")
        if ch["by_residue_type"]:
            subsection("D-amino acids by residue type")
            for rn, info in ch["by_residue_type"].items():
                lines.append(
                    f"  {rn:>7s}  "
                    f"{_fmt_ratio(info['D'], info['total'], info['D_percent'])}"
                )

    # --- DSSP ---
    if "dssp" in stats:
        section("DSSP SECONDARY STRUCTURE")
        ds = stats["dssp"]
        subsection("Overall")
        lines.append(
            f"  {'State':>6s}"
            f"  {'Count':>10s}  {'Percent':>8s}"
        )
        lines.append(
            f"  {'-----':>6s}"
            f"  {'-----':>10s}  {'-------':>8s}"
        )
        for state, info in ds["overall"].items():
            lines.append(
                f"  {state:>6s}"
                f"  {info['count']:>10,}"
                f"  {info['percent']:>7.2f}%"
            )

        subsection("By Ramachandran category")
        for rama, states in ds["by_rama_category"].items():
            parts = [
                f"{s}:{info['percent']:.1f}%"
                for s, info in states.items()
            ]
            lines.append(
                f"  {rama:10s}  {'  '.join(parts)}"
            )

    # --- Geometry ---
    if "geometry" in stats:
        section("GEOMETRIC DISTRIBUTIONS")
        lines.append(
            "  Note: Torsion angles (phi, psi, omega, chi) use circular"
        )
        lines.append(
            "  statistics (atan2 mean, sqrt(-2 ln R) std). Tau is linear."
        )
        lines.append("")
        geom = stats["geometry"]
        subsection("Overall")
        lines.append(
            f"  {'Measure':>8s}  {'N':>10s}"
            f"  {'Mean':>10s}  {'Std':>10s}"
        )
        lines.append(
            f"  {'-------':>8s}  {'-':>10s}"
            f"  {'----':>10s}  {'---':>10s}"
        )
        for f, info in geom["overall"].items():
            lines.append(
                f"  {f:>8s}  {info['n']:>10,}"
                f"  {info['mean']:>10.3f}"
                f"  {info['std']:>10.3f}"
            )

        if "omega_deviation_trans" in geom:
            od = geom["omega_deviation_trans"]
            lines.append("")
            lines.append(
                "  Omega deviation from 180 (trans only):"
            )
            lines.append(
                f"    N={od['n']:,}"
                f"  Mean={od['mean']:.3f} deg"
                f"  Std={od['std']:.3f} deg"
            )

        subsection("By Ramachandran category")
        for rama, measures in geom["by_rama_category"].items():
            lines.append(f"  {rama}:")
            for f, info in measures.items():
                lines.append(
                    f"    {f:>8s}"
                    f"  N={info['n']:>10,}"
                    f"  Mean={info['mean']:>10.3f}"
                    f"  Std={info['std']:>10.3f}"
                )

    # --- Data completeness ---
    if "data_completeness" in stats:
        section("DATA COMPLETENESS")
        lines.append(
            f"  {'Field':>15s}  {'Present':>10s}"
            f"  {'Missing':>10s}  {'% Present':>10s}"
        )
        lines.append(
            f"  {'-----':>15s}  {'-------':>10s}"
            f"  {'-------':>10s}  {'---------':>10s}"
        )
        for f, info in stats["data_completeness"].items():
            lines.append(
                f"  {f:>15s}"
                f"  {info['present']:>10,}"
                f"  {info['missing']:>10,}"
                f"  {info['percent_present']:>9.2f}%"
            )

    # --- Notes ---
    n_chain_ids = stats.get("overview", {}).get("unique_chains", 0)
    n_struct = stats.get("overview", {}).get(
        "unique_structures",
        stats.get("overview", {}).get("unique_files", 0),
    )
    if n_chain_ids != n_struct:
        lines.append("")
        lines.append(
            f"  [1] The output contains {n_chain_ids:,} unique"
            f" (file, chain) combinations."
        )
        lines.append(
            "      This may differ from the source chain list count"
        )
        lines.append(
            "      if RCSB PDB remediation renamed or split chains"
        )
        lines.append(
            "      relative to the original dataset. See issues"
        )
        lines.append(
            "      for per-entry details."
        )

    lines.append("")
    return "\n".join(lines)


# -------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze pydangle JSONL output files.",
    )
    parser.add_argument("input", help="Path to JSONL file")
    parser.add_argument(
        "-o", "--output-format",
        choices=["text", "json"],
        default="text",
        help="Output format (default: text)",
    )
    parser.add_argument(
        "--output", "-O",
        help="Output file (default: stdout)",
    )
    args = parser.parse_args()

    print(f"Loading {args.input}...", file=sys.stderr)
    records = load_records(args.input)
    print(f"Loaded {len(records):,} records.", file=sys.stderr)

    print("Analyzing...", file=sys.stderr)
    stats = analyze(records)

    if args.output_format == "json":
        output = json.dumps(stats, indent=2)
    else:
        output = format_text_report(stats)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
            f.write("\n")
        print(f"Written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()
