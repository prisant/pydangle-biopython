# pydangle-biopython

Michael G. Prisant, Ph.D.
Prisant Scientific
info@prisantscientific.org

Compute distances, angles, dihedrals, per-residue classifications,
and DSSP secondary structure assignments in protein and nucleic acid
structures.

This package uses **BioPython** as the structure-parsing backend.
Sibling packages `pydangle-cctbx` and `pydangle-gemmi` will provide the same
functionality with alternative backends.

## Development story

[Pydangle](https://github.com/prisant/pydangle-biopython) is a Python reimagining
of the Java [Dangle](https://github.com/rlabduke/javadev) program
(`chiropraxis.dangle.Dangle`).

The original Java Dangle was written in the Richardson Laboratory at Duke
University by Ian Davis and Vincent B. Chen as part of the
[MolProbity](https://github.com/rlabduke/MolProbity) suite.

This new Python reimagining began with foundation code hand-coded by the author while still working as a senior research scientist in the Richardson Laboratory at the Duke University Medical School's Department of Biochemistry. 
The foundation code was back-engineered from the original program's
documentation and output rather than a literal translation of the original Java.

Development continued after the author's retirement from Duke in July 2025 with assistance from Anthropic's Claude.
The new project used the hello project (also developed jointly with Claude) as a template.
The original code was first backported and refactored to this template.

Subsequent early development added nucleic acid builtins and was tested against archived versions of older dangle output and documentation.
Later extensions included the ability to run on multiple files.
Extensive testing and debugging followed using the top100 pdb files as well as a newer version of Java Dangle -- output, documentation, and source code -- for comparison.

This project is intended as a base tool for the author's collaborative work with Ezra Miller's group in the Duke Math Department on application of Topological Data Analysis (TDA) to problems in statistical analysis of macromolecular structural data starting with protein Ramachandran maps.

## Installation

```bash
git clone https://github.com/prisant/pydangle-biopython.git
cd pydangle-biopython

python -m venv .venv
source .venv/bin/activate

make install          # editable install + dev deps + git hooks
```

## Usage

```bash
# Default measurements (phi, psi, chi1–chi4)
pydangle-biopython structure.pdb

# Protein backbone
pydangle-biopython -c 'tau; phi; psi; omega' structure.pdb

# Nucleic acid backbone
pydangle-biopython -c 'alpha; beta; gamma; delta; epsilon; zeta' rna.cif

# Custom distance
pydangle-biopython -c 'distance: Ca_Ca: i-1 _CA_, i _CA_' structure.pdb

# Multiple input files
pydangle-biopython -c 'phi; psi; omega' *.pdb

# Bulk input: read paths from a file list
pydangle-biopython -c 'phi; psi' -f file_list.txt

# Bulk input: pipe paths from find/ls
find /data -name '*.pdb' | pydangle-biopython -c 'phi; psi' -f -

# Bulk input: Python glob (bypasses shell ARG_MAX)
pydangle-biopython -c 'phi; psi' -g '**/*.pdb'

# Bulk input: recursive directory search
pydangle-biopython -c 'phi; psi' -d /data/structures/

# Parallel processing (auto-detect CPU cores)
pydangle-biopython -j -c 'phi; psi' -d /data/structures/

# Parallel processing (4 workers)
pydangle-biopython -j 4 -c 'phi; psi; omega; tau' -d /data/structures/

# Residue classification labels
pydangle-biopython -c 'phi; psi; rama_category; is_cis' structure.pdb

# Mixed measurements, labels, and completeness checks
pydangle-biopython -c 'phi; psi; rama_category; has_all_mc; has_all_sc' structure.pdb

# DSSP secondary structure (requires mkdssp on PATH)
pydangle-biopython -c 'phi; psi; dssp; dssp3; rama_category' structure.pdb
```

## Output format

Output is colon-separated with one line per residue:

```
filename:model:chain:number:ins:resname:measurement1:measurement2:...
```

For example:

```
1ubq.pdb:1:A:2: :GLN:-93.066:132.565:General:False
```

Lines beginning with `#` are header comments.  Unknown or unmeasurable values
are reported as `__?__`.

## Measurement types

Pydangle supports four types of per-residue output:

**Geometric measurements** compute distances (Å), angles (°), and dihedral
angles (°) from atom coordinates.  These can be specified as builtins
(e.g. `phi`, `psi`, `chi1`) or with custom atom specifications
(e.g. `distance: CaCa: i _CA_, i+1 _CA_`).

**Residue classification labels** return categorical strings for each residue.
Labels use a two-field syntax with no atom specifications
(e.g. `rama_category`, `is_cis`).

**DSSP secondary structure** labels call the external ``mkdssp`` binary
to assign per-residue secondary structure using the Kabsch & Sander
hydrogen-bond energy method.  This requires ``mkdssp`` to be installed
on ``$PATH`` (``apt install dssp`` on Debian/Ubuntu, or
``conda install -c conda-forge dssp``).  If ``mkdssp`` is not found,
DSSP labels return ``__?__``.

### Builtin labels

`rama_category` (also `rama6`) assigns wwPDB Validation Task Force Ramachandran
categories (Read et al., Structure 19:1395–1412, 2011): General, Gly, IleVal,
TransPro, CisPro, or PrePro.  Reduced schemes are available as `rama5`
(merges cis/trans Pro), `rama4` (also merges IleVal into General), and
`rama3` (General, Gly, Pro only, following Lovell et al. 2003).

`is_cis` / `is_trans` classify the peptide bond based on the omega dihedral
angle (cis when |ω| < 30°).

`is_gly`, `is_pro`, `is_ileval` identify residue type.  `is_prepro` flags
residues immediately preceding proline.

`has_all_mc` / `has_all_sc` check whether all expected mainchain (N, CA, C, O)
or sidechain heavy atoms are present for each residue.

`is_left` / `is_right` determine Cα chirality from the improper dihedral
CB–N–C–CA (L-amino acids are negative, D-amino acids positive).  GLY returns
unknown since it has no CB.

See `examples/README.md` for the complete builtin reference table covering
protein backbone, sidechain, nucleic acid, and label builtins.

## Differences from Java Dangle

Pydangle matches Java Dangle output to within rounding precision on the
top100 PDB benchmark set.  Key differences include distance-based chain
break detection (rather than sequence numbering gaps) and compact chain ID
formatting.

Pydangle extends Java Dangle with residue classification labels (rama3–6,
cis/trans, chirality, atom completeness), DSSP secondary structure
integration, mmCIF input support, multi-file processing, flexible
bulk input (file lists, stdin piping, Python glob patterns, recursive
directory search), multiprocessing for parallel file processing
(``-j``/``--jobs``), and resilient PDB parsing that handles files
(e.g. Reduce hydrogen-added PDBs) that crash or confuse other parsers.
Java Dangle includes validation measurements (cbdev, hadev, nhdev, codev),
disulfide analysis, and nucleic acid diagnostics (pperp, pucker, suitefit)
not yet in pydangle.

See ``DIFFERENCES.md`` for a comprehensive comparison and ``BENCHMARK.md``
for performance and correctness statistics.

## Development

```bash
make test        # run tests
make test-cov    # tests with coverage
make lint        # ruff check
make typecheck   # mypy strict
make all         # lint + typecheck + test
make benchmark   # performance/correctness comparison with Java Dangle
make bump        # bump version, commit, and tag
make clean       # remove build artifacts and editor backups
make docs        # serve docs locally
```

See [DEVELOPMENT.md](DEVELOPMENT.md) for the full developer guide including
git conventions, branching, releasing, and applying patches.

## License

MIT
