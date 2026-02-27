# pydangle-biopython

Compute distances, angles, and dihedral angles in protein and nucleic acidstructures.  

This package uses **BioPython** as the structure-parsing backend.
Sibling packages `pydangle-cctbx` and `pydangle-gemmi` will provide the same functionality with alternative backends.

[Pydangle](https://github.com/prisant/pydangle-byopython) is a Python rewrite of the Java [Dangle](https://github.com/rlabduke/javadev) program (`chiropraxis.dangle.Dangle.`) 
The original Java [Dangle] was written in the Richardson Laboratory at Duke University by Ian Davis and Vincent B. Chen as part of the [MolProbity](https://github.com/rlabduke/MolProbity) suite.
This new Python version [Pydangle] was written based on the original program's documentation and output as opposed to a translation of code.
The first version was handcoded by Michael G. Prisant but was then subsequently revised and extended with the assistence of Anthropic Claude Code. 

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
```

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
