# Differences from Java Dangle

Pydangle aims for scientific compatibility with Java Dangle while improving
on some design choices and adding new capabilities.

## Measurement compatibility

On the top100 PDB benchmark set with `phi; psi; omega`, pydangle matches
Java Dangle output to within rounding precision on all residues, with the
intentional differences noted below.  See ``BENCHMARK.md`` for full
performance and correctness statistics.

### Chain gap handling

Java Dangle treats discontinuities in residue sequence numbering as chain
breaks, suppressing measurements across the gap even when coordinates are
present and the peptide bond is intact.  Pydangle uses BioPython's
`CaPPBuilder` to detect chain breaks based on CA–CA distance (< 4.6 Å),
so measurements are computed whenever residues are physically connected.
This adds a small number of correct measurements that Java Dangle omits
(+4 lines across the top100 set).

### Modified amino acids

Both programs use the same residue-name whitelist (derived from Java Dangle's
`isProtOrNucAcid()`) to identify polymer residues.  Standard amino acids and
common modifications such as PCA (pyroglutamic acid), MSE (selenomethionine),
and CSD (cysteine sulfonic acid) are included.  Non-polymer het groups
(chromophores, cofactors, etc.) are excluded.

### Chain ID formatting

Java Dangle space-pads chain identifiers (e.g. ` A:`), following legacy PDB
conventions.  Pydangle uses compact chain IDs (e.g. `A:`) for easier
programmatic parsing.

### Omega definition

Both programs define omega as the dihedral `i-1 CA, i-1 C, i N, i CA`,
following the IUPAC-IUB 1970 convention.  Omega is reported on residue *i*
(the residue containing the N atom of the peptide bond).

## Features in pydangle but not in Java Dangle

- **Residue classification labels** — categorical per-residue annotations
  including Ramachandran category (rama3–6), cis/trans peptide bond
  (is_cis, is_trans, peptide_bond), residue identity (is_gly, is_pro,
  is_ileval, is_prepro), atom completeness (has_all_mc, has_all_sc),
  and Cα chirality (is_left, is_right, chirality).
- **DSSP secondary structure** — dssp (8-state) and dssp3 (3-state)
  labels via external mkdssp, providing hydrogen-bond-based secondary
  structure assignment independent of phi/psi geometry.
- **Multiple input files** — process many structures in a single command.
- **Bulk input methods** — file lists (`-f`), stdin piping (`-f -`),
  Python glob patterns (`-g`), and recursive directory search (`-d`)
  for processing large numbers of structures without shell ARG_MAX limits.
- **Parallel processing** — multiprocessing support (`-j`/`--jobs`) for
  processing large datasets across multiple CPU cores.
- **JSONL output format** — `-o jsonl` produces one JSON object per line with
  typed values (floats for angles, strings for labels, null for unknown),
  suitable for `jq` pipelines and programmatic consumption.
- **mmCIF input support** — reads both PDB and mmCIF formats via BioPython.
- **Distance-based chain break detection** — uses CA–CA distance (< 4.6 Å)
  rather than residue sequence numbering to identify chain breaks.
- **Resilient PDB parsing** — recovers from BioPython header parsing
  failures by retrying with coordinate-only lines, enabling processing
  of files (e.g. Reduce hydrogen-added PDBs) that would otherwise crash.
- **Hydrogen-added PDB support** — handles PDB files with hydrogen atoms
  and blank chain IDs (common in Reduce-processed files) that Java Dangle
  cannot compute backbone dihedrals on.

## Features in Java Dangle not yet in pydangle

### Protein validation measurements

- **cbdev** — Cbeta deviation: distance from observed CB to ideal CB
  position constructed from backbone geometry (N, CA, C) using
  amino-acid-specific parameters (Ala, Pro, Val/Thr/Ile, Gly, and
  a general class).
- **hadev** — HA deviation: distance from observed HA to ideal position.
- **nhdev** — N-H deviation: distance from observed H to ideal position
  using Reduce-style N-H placement parameters.
- **codev** — C=O deviation: distance from observed O to ideal position
  using Engh & Huber parameters.

### Disulfide bond measurements

- **ss** (disulfides) — composite measurement of phi, psi, chi1, chi2,
  chi3, and S-S distance for disulfide-bonded cysteines.

### Nucleic acid measurements

- **pperp** — base-phosphate perpendicular distance for sugar pucker
  diagnostics.
- **PuckerAmp** — sugar pucker amplitude.
- **PuckerAng** — sugar pucker phase angle.
- **suitefit** — nine bond distances used for RNA suite fitting.
- **rnabb** / **dnabb** — composite RNA/DNA backbone measurement shortcuts.

### Other measurement types

- **VectorAngle** — angle between two vectors.
- **MaxB** — maximum B-factor for a set of atoms.
- **MinQ** — minimum occupancy for a set of atoms.
- **Planarity** — deviation from planarity for a set of atoms.

### Input

(No remaining input gaps — stdin piping, file lists, glob patterns, and
directory recursion are now supported in pydangle.)
