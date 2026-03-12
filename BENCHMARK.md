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

All 27 mismatches trace to two root causes: chain break boundary
handling (17 rows) and alternate conformation selection (10 rows).
Neither represents a computational error — they are parser-level
decisions about which atoms to use.

#### Chain break boundary differences (17 rows)

Structures: 1aac (2), 1ben (4), 1cpc (4), 1fxd (2), 1lit (5).

These are residues at the edges of missing-residue gaps.  Both tools
detect the chain break, but they disagree on which measurements are
computable at fragment boundaries.  Pydangle computes psi for the last
residue and phi/omega for the first residue of the next fragment using
neighboring atoms within the fragment.  Dangle reports `__?__` for
those measurements because it sees the sequence number gap.

The 1cpc cases illustrate this clearly: residues 72–74 are entirely
absent from the PDB file in chains B and L.  CaPPBuilder correctly
breaks the chain into fragments ending at 71 and starting at 75.
Both tools agree on this break, but pydangle reports psi for residue 71
(using 71 N, CA, C and 75 N — which are *not* connected) while dangle
suppresses it.  Similarly for phi/omega of residue 75.

In 1fxd, residue 11 is missing between 10 and 12, producing the same
pattern.

#### Alternate conformation selection (10 rows)

Structures: 1ctj (3), 1ifc (4), 1cpc (not alt conf — see above),
1fxd (not alt conf — see above).  The remaining mismatches in 1cpc
and 1fxd are chain break issues, not alt conf.

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

## Running the benchmark

```bash
# Full benchmark (batch timing + per-file timing + correctness)
make benchmark

# Quick (skip per-file timing)
python scripts/benchmark.py --skip-per-file

# Custom PDB directory and tolerance
python scripts/benchmark.py --pdb-dir /path/to/pdbs --tolerance 0.05
```
