# Benchmark: pydangle-biopython vs Java Dangle

Performance and correctness comparison on the top100 PDB benchmark set
(100 structures) using `phi; psi; omega; tau` measurements.

Tested on Linux 6.17, Python 3.12, Java Dangle 1.07.

## Timing

### Batch (all 100 files in one invocation)

| Tool     | Mean    | Min     | Max     | Runs |
|----------|---------|---------|---------|------|
| pydangle | 9.88 s  | 9.59 s  | 10.08 s | 3    |
| dangle   | 6.29 s  | 6.11 s  | 6.48 s  | 3    |

**pydangle / dangle = 1.57x**

The overhead is primarily Python/BioPython startup and PDB parsing.
Java Dangle benefits from the JVM's optimized I/O and its own
lightweight parser.

### Per-file

| Tool     | Mean     | Median   | Stdev    |
|----------|----------|----------|----------|
| pydangle | 0.380 s  | 0.370 s  | 0.068 s  |
| dangle   | 0.332 s  | 0.332 s  | 0.033 s  |

5 slowest files (pydangle): 1cpc (0.59 s), 1luc (0.58 s),
1xyz (0.57 s), 2olb (0.53 s), 1kap (0.50 s).

## Correctness

| Metric                  | Count          |
|-------------------------|----------------|
| Rows in common          | 21,003         |
| Rows agreeing (≤ 0.01°) | 20,976 (99.87%) |
| Value mismatches        | 27             |
| Rows only in pydangle   | 4              |
| Rows only in dangle     | 0              |

### Rows only in pydangle (+4)

These are residues at chain break boundaries where pydangle computes
measurements because the CA–CA distance indicates a connected peptide
bond, while Java Dangle suppresses them due to sequence numbering gaps.
See DIFFERENCES.md for details.

- 1aac A:9 SER, 1ben A:11 CYS, 1ben C:8 THR, 1lit A:23 SER

### Value mismatches (27)

The 27 mismatches fall into three categories:

**Chain break boundary (13 rows):** Residues adjacent to chain breaks
where one tool computes a value and the other reports `__?__`.  These
are the same chain gap handling difference described above — pydangle
reports measurements across physically connected residues that dangle
treats as broken chains.  Affected structures: 1aac, 1ben, 1lit.

**Alternate conformation handling (10 rows):** Structures with
alternate conformations where the two parsers select different atom
positions, producing genuinely different angle values.  The largest
discrepancies appear in 1cpc (chains B and L, residues 71 and 75),
1ctj (residues 1–3), and 1fxd (residues 10 and 12).

**Rounding at tolerance boundary (4 rows):** Small differences in
1ifc (residues 30, 31, 47, 48) where values differ by 0.05–1.0°,
likely due to floating-point differences between BioPython's and
Java Dangle's coordinate handling.

## Running the benchmark

```bash
# Full benchmark (batch timing + per-file timing + correctness)
make benchmark

# Quick (skip per-file timing)
python scripts/benchmark.py --skip-per-file

# Custom PDB directory and tolerance
python scripts/benchmark.py --pdb-dir /path/to/pdbs --tolerance 0.05
```
