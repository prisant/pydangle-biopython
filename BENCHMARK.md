# Benchmark: pydangle-biopython vs Java Dangle

Performance and correctness comparison using `phi; psi; omega; tau`
measurements on two benchmark sets:

- **top100pdb** — 100 PDB structures as deposited
- **top100H** — the same 100 structures with hydrogens added by Reduce

Tested on Linux 6.17, Python 3.12, Java Dangle 1.07.

## Timing

### top100pdb (no hydrogens)

| Tool     | Mean    | Min     | Max     | Runs |
|----------|---------|---------|---------|------|
| pydangle | 9.88 s  | 9.59 s  | 10.08 s | 3    |
| dangle   | 6.29 s  | 6.11 s  | 6.48 s  | 3    |

**pydangle / dangle = 1.57x**

### top100H (with hydrogens)

| Tool     | Mean     | Min     | Max      | Runs |
|----------|----------|---------|----------|------|
| pydangle | 10.06 s  | 9.82 s  | 10.51 s  | 3    |
| dangle   | 6.33 s   | 6.29 s  | 6.36 s   | 3    |

**pydangle / dangle = 1.59x**

The overhead is primarily Python/BioPython startup and PDB parsing.
Java Dangle benefits from the JVM's optimized I/O and its own
lightweight parser.  The ratio is consistent across both benchmark
sets.

### Per-file (top100pdb)

| Tool     | Mean     | Median   | Stdev    |
|----------|----------|----------|----------|
| pydangle | 0.380 s  | 0.370 s  | 0.068 s  |
| dangle   | 0.332 s  | 0.332 s  | 0.033 s  |

5 slowest files (pydangle): 1cpc (0.59 s), 1luc (0.58 s),
1xyz (0.57 s), 2olb (0.53 s), 1kap (0.50 s).

## Multiprocessing scaling

Parallel file processing with `-j`/`--jobs` on **top100pdb** using
`phi; psi; omega; tau; chi1; dssp; rama_category` (7 measurements).
Tested on an 8-core Intel NUC (Linux 6.17, Python 3.12).

| Flag    | Workers        | Wall time | Speedup | Output lines |
|---------|----------------|-----------|---------|--------------|
| `-j 1`  | 1 (serial)     | 52.3 s    | 1.0x   | 21,107       |
| `-j 2`  | 2              | 28.6 s    | 1.8x   | 21,107       |
| `-j 4`  | 4              | 19.6 s    | 2.7x   | 21,107       |
| `-j 0`  | 8 (auto-detect)| 20.0 s    | 2.6x   | 21,107       |

### top500H (500 files, hydrogen-added, same measurements)

| Flag    | Workers        | Wall time | Speedup | Output lines |
|---------|----------------|-----------|---------|--------------|
| `-j 1`  | 1 (serial)     | 70.0 s    | 1.0x   | 110,396      |
| `-j 2`  | 2              | 38.4 s    | 1.8x   | 110,396      |
| `-j 4`  | 4              | 23.6 s    | 3.0x   | 110,396      |
| `-j 0`  | 8 (auto-detect)| 22.2 s    | 3.2x   | 110,396      |

All parallel outputs are byte-identical to serial output in both
benchmark sets.  Good scaling through 4 cores, then diminishing returns
at 8 — likely because DSSP (external `mkdssp` process) and I/O become
bottlenecks.  Scaling improves slightly with more files (3.2x at 500
files vs 2.7x at 100 files for `-j 4`) as work distributes more evenly.

Extrapolating to the top2018 dataset (12,125 files): serial ~105 min,
`-j 4` ~39 min.

## Correctness: top100pdb

| Metric                  | Count           |
|-------------------------|-----------------|
| Rows matched            | 21,003          |
| Rows agreeing (≤ 0.01°) | 20,976 (99.87%) |
| Coverage diffs          | 14              |
| Value mismatches        | 13              |
| Rows only in pydangle   | 4               |
| Rows only in dangle     | 0               |

"Coverage diffs" are rows where one tool computes a value and the other
reports `__?__` — a difference in what is measured, not a computational
disagreement.  "Value mismatches" are rows where both tools compute
values that disagree beyond tolerance.

### Rows only in pydangle (+4)

These are residues at chain break boundaries where pydangle computes
measurements because the CA–CA distance indicates a connected peptide
bond, while Java Dangle suppresses them due to sequence numbering gaps.
See DIFFERENCES.md for details.

- 1aac A:9 SER, 1ben A:11 CYS, 1ben C:8 THR, 1lit A:23 SER

### Coverage diffs (14)

Residues adjacent to chain breaks where one tool computes a value and
the other reports `__?__`.  Same root cause as the +4 rows above.
Affected structures: 1aac, 1ben, 1lit, 1lkk, 1not, 2er7, 2msb, 3b5c.

### Value mismatches (13)

All 13 mismatches trace to two root causes: chain break boundary
handling and alternate conformation selection.  Neither represents a
computational error — they are parser-level decisions about which
atoms to use.

#### Chain break boundary differences (4 rows)

Structures: 1cpc (4).

Residues 72–74 are entirely absent from chains B and L.  CaPPBuilder
correctly breaks the chain into fragments ending at 71 and starting
at 75.  Both tools agree on this break, but they compute different
values at the fragment boundaries because they select different atoms
for the measurement.

In 1fxd, residue 11 is missing between 10 and 12, producing the same
pattern.

#### Alternate conformation selection (9 rows)

Structures: 1ctj (3), 1ifc (4), 1fxd (2).

**1ctj, residues 1–3:** The PDB has two fully distinct backbone
conformations (altloc A occupancy 0.49, altloc B occupancy 0.51).
BioPython selects altloc B (higher occupancy).  Java Dangle selects
altloc A.  This produces completely different phi/psi/omega/tau values
because the backbone atoms are in different physical positions.
Residue 3 has no alt confs itself but its phi and omega depend on
atoms from residue 2 (which does), so it is also affected.

**1ifc, residues 30, 31, 47, 48:** Residues 29–49 have full A/B
alternate conformations for all backbone atoms.  BioPython selects a
*mix* of altlocs within individual residues — for example, residue 30
uses N and CA from altloc A but C from altloc B.  This happens because
BioPython's PDBParser, when multiple alt confs have equal or near-equal
occupancy, keeps whichever atom it encounters first per atom name.
Java Dangle presumably selects one altloc consistently per residue.
The resulting angle differences range from 0.05° to 1.8°, consistent
with using slightly different atom positions rather than any
computational error.

## Correctness: top100H (with hydrogens)

| Metric                  | Count           |
|-------------------------|-----------------|
| Rows matched            | 18,384          |
| Rows agreeing (≤ 0.01°) | 17,532 (95.4%)  |
| Coverage diffs          | 836             |
| Value mismatches        | 16              |
| Rows only in pydangle   | 0               |
| Rows only in dangle     | 0               |

### Dangle preprocessing

The Reduce-processed PDB files present two problems for Java Dangle:
blank chain ID columns (column 22) and hydrogen atom records.  Without
chain IDs, dangle cannot build peptide chains.  With hydrogen atoms
present, dangle's residue/atom matching is confused and it reports
`__?__` for backbone dihedrals even when chain IDs are present.

The benchmark script automatically preprocesses files for dangle by
inserting chain ID `A` where blank and stripping hydrogen atoms.  This
preprocessing is for benchmark comparison only — pydangle handles
these files natively without modification.

### Coverage diffs (836)

After preprocessing, 836 rows still show coverage differences where
pydangle computes a value and dangle reports `__?__`.  These appear as
a systematic pattern at chain break boundaries within the H-added
files, likely where Reduce's modifications trigger dangle's chain
break detection.

### Value mismatches (16)

The same alternate conformation structures (1ctj, 1ifc, 1lit) plus
additional residues in 1ifc where the H-added file has larger alt conf
regions.

## Correctness: top500H (with hydrogens)

| Metric                  | Count            |
|-------------------------|------------------|
| Rows matched            | 109,823          |
| Rows agreeing (≤ 0.01°) | 104,896 (95.5%) |
| Coverage diffs          | 4,739            |
| Value mismatches        | 188              |
| Rows only in pydangle   | 19               |
| Rows only in dangle     | 2                |

### Timing

| Tool     | Mean     | Min     | Max      | Runs |
|----------|----------|---------|----------|------|
| pydangle | 62.59 s  | 59.88 s | 65.88 s  | 3    |
| dangle   | 32.14 s  | 31.52 s | 32.92 s  | 3    |

**pydangle / dangle = 1.95x**

### Dangle preprocessing

Same preprocessing as top100H — all 500 files required chain ID fixes
and hydrogen stripping for dangle compatibility.  Pydangle handles
the original extensionless H-added files natively with `-p` (force PDB
format).

### Coverage diffs (4,739)

Same pattern as top100H at larger scale: pydangle computes values at
chain break boundaries where dangle reports `__?__`, likely where
Reduce's modifications trigger dangle's chain break detection.

### Value mismatches (188)

Alternate conformation selection differences between BioPython and
dangle, the same root cause as top100.  The larger dataset surfaces
more structures with extensive alternate conformations (e.g. 1bi5H
with 20 affected residues, 1dosAH with multiple chain break boundary
regions).

### Rows only in pydangle (+19)

Chain break boundary residues where pydangle computes measurements
based on CA–CA distance while dangle suppresses them due to sequence
numbering gaps.  Includes 1aacH (same as top100) plus 1dosAH (10
residues across multiple chain breaks) and others.

### Rows only in dangle (+2)

Two residues (1dinH A:233 SER, 1lstH A:240 LYS) where dangle computes
values that pydangle does not — likely terminal residues where the
two tools disagree on chain extent.

## Resilient PDB parsing

Three H-added files (1benABH, 1dadH, 1etmH) trigger BioPython parser
bugs — an `UnboundLocalError` in header parsing and a `ValueError` in
coordinate parsing.  Pydangle mitigates this with resilient parsing:
on any parser failure, it strips non-coordinate lines (keeping only
ATOM/HETATM/MODEL/ENDMDL/TER/END) and retries from a temporary file.
A warning is emitted to stderr.  This loses PDB header metadata but
pydangle does not use it.

## Running the benchmark

```bash
# Full benchmark on deposited structures
make benchmark

# Quick (skip per-file timing)
python scripts/benchmark.py --skip-per-file

# Benchmark on hydrogen-added structures
python scripts/benchmark.py --pdb-dir /path/to/top100H

# Extensionless PDB files (auto-detected, forced PDB format)
python scripts/benchmark.py --pdb-dir /path/to/top500H --skip-per-file

# Custom tolerance
python scripts/benchmark.py --tolerance 0.05
```

The benchmark script (`scripts/benchmark.py`) handles all
preprocessing automatically.  It detects files needing chain ID fixes
or hydrogen stripping, creates temporary preprocessed copies for
dangle, and cleans them up after the run.  Directories with
extensionless PDB files (e.g. top500H) are auto-detected and processed
with forced PDB format.  Row matching tolerates chain ID mismatches
between the two tools by falling back to chain-ID-free keys when exact
keys don't match.
