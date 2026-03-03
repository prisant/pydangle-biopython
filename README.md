# pydangle-biopython

Michael G. Prisant, Ph.D.
Prisant Scientific
info@prisantscientific.com

Compute distances, angles, and dihedral angles in protein and nucleic acid
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
```

## Output format

Output is colon-separated with one line per residue:

```
filename:model:chain:number:ins:resname:phi:psi:omega
```

For example:

```
1ubq.pdb:1:A:2: :GLN:-93.066:132.565:179.137
```

Lines beginning with `#` are header comments.  Unknown or unmeasurable values
are reported as `__?__`.

## Differences from Java Dangle

Pydangle aims for scientific compatibility with Java Dangle while improving
on some design choices.  On the top100 PDB benchmark set with `phi; psi; omega`,
pydangle matches Java Dangle output to within rounding precision on all
residues, with the following intentional differences.

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

## Development

```bash
make test        # run tests
make test-cov    # tests with coverage
make lint        # ruff check
make typecheck   # mypy strict
make all         # lint + typecheck + test
make docs        # serve docs locally
```

## License

MIT
