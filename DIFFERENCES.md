# Differences from Java Dangle

Pydangle aims for scientific compatibility with Java Dangle while improving
on some design choices and adding new capabilities.

## Measurement compatibility

On the top100 PDB benchmark set with `phi; psi; omega`, pydangle matches
Java Dangle output to within rounding precision on all residues, with the
intentional differences noted below.  See ``BENCHMARK.md`` for full
performance and correctness statistics.

### Chain gap handling

Both programs decide that residue *i* and *i+1* are joined by testing for a
peptide bond: C(*i*)–N(*i+1*) ≤ 2.0 Å.  Java Dangle steps to the
neighbouring residue in the model and applies this test
(`AtomSpec.checkConnection`, `sqDistance <= 4.00`, "from Molikin"; O3′–P
≤ 2.2 Å for nucleic acids).  Pydangle builds fragments with BioPython's
`PPBuilder` using the same 2.0 Å test (`_PEPTIDE_BOND_MAX` in
`measure.py`; since 2026-11-08 — earlier versions used `CaPPBuilder`, a
CA–CA distance test, and then removed non-whitelisted residues from the
fragment, which joined the residues on either side of a modified residue;
see `CHANGELOG.md`, 5a22669).

Remaining differences:

- **Altloc labels.**  Pydangle joins two residues only through a
  compatible pair of alternates (same altloc, or one of them blank).
  Java Dangle ignores altloc labels in the bond test.  Example: 3rdo
  A−6 exists only as altloc A and A−5 only as altloc B (C–N 1.33 Å);
  Java links them, pydangle does not.  cctbx (MolProbity) does not link
  them either.
- **Isolated residues are not reported.**  A residue with no peptide bond
  to either neighbour — a free amino acid modelled as a ligand, or a
  residue stranded between two gaps — has no φ, ψ or ω.  Pydangle leaves
  it out of the output entirely; Java Dangle reports it with τ (and χ)
  only.  Decision 2026-11-08 (M. Prisant): drop them; they contribute
  nothing to backbone analysis and are not polypeptide residues.  About
  6–10 such records per dataset (top500: 3std C2; top8000/top2018: e.g.
  3gfv B304 ASN, 1i9c B1305 GLU, 1qgj B235–236).  Exception: a chain with
  no peptide bond at all (e.g. nucleic acid) falls back to its raw
  residue list.
- **C-terminal stubs.**  A residue bonded to its predecessor but modelled
  with only an N atom is reported (all angles null); Java Dangle omits it.

### Alternate conformations

By default (`--altloc first`, since v0.6) pydangle measures the first
conformer, as Java Dangle and cctbx (the source of the Richardson Lab's
reference angles from Top8000 onward) do.  "First" is per residue: among
the alternates of the atoms concerned, the label written first in the file
(lowest atom serial; normally A).  Ranking labels over the whole file
instead fails where the file's first alternate residue has no A (3zeu A1
has only B and C; the lab measured 3zeu B168 in A).  Each residue is
measured in the first conformer shared by every atom of its φ/ψ window,
C(i−1), N, CA, C, N(i+1); a measurement whose atoms lack that conformer
uses the first conformer its own atoms share, and failing that each atom's
first conformer.  So a residue present only as altloc B, and the residue
after it, are measured in B, as in the lab's Top2018 database (6b6w
B301–B302, 5j1s B431–B432).  Over all 9,095 Top2018 entries with
alternates, pydangle and that database agree within 0.01° on 451,148 of
451,186 φ/ψ values near an alternate.  Of the rest, 26 are in two entries
the database numbers differently from the PDB files (4oe9, 4z54), and 12,
in 7 entries, are values where the database uses conformer B and pydangle
A.  The code is in `altloc.py`.

At microheterogeneity (two residue types at one position) pydangle keeps
one residue type, so each position gives one record: the first conformer
whose type is a standard residue (the 20 amino acids or MSE), or the first
conformer if none is.  2zs0 D31 Ile/Val gives the Ile, as Java Dangle
measures it; 5k0a A157 CSO/Cys (CSO is conformer A) gives the Cys.  The
lab's Top2018 database measures the standard conformer at all 13 such
sites (26 φ/ψ values agree) but names the residue as the modified type.
BioPython represents such sites either as one disordered residue or, when
one type is written as HETATM (1din A123 CSD/Cys), as two residues with
the same number; both are handled.  Where the first type has an atom
without an altloc label, BioPython's parser keeps only that type (1ejg
A22 Pro/Ser: the Pro; Java Dangle measured the Ser).

`--altloc occupancy` restores the behaviour before v0.6, BioPython's
default: for each atom separately the alternate with the highest occupancy
(ties: the first listed), so one dihedral can mix atoms of two conformers,
and split microheterogeneity sites give two records with the same residue
number.  It reproduces earlier pydangle output byte for byte on the five
Richardson datasets.  The JSONL `_meta` record states the mode
(`"altloc"`).  Decisions 2026-11-09 (first conformer) and 2026-11-10
(per-residue order, standard residue type).

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
  is_ileval, is_prepro, is_standard), atom completeness (has_all_mc, has_all_sc),
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
