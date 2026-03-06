# Getting Started

## Prerequisites

- Python 3.10 or later
- pip
- git (for changelog automation)

## Installation

Clone the repository and install in a virtual environment:

```bash
git clone https://github.com/prisant/pydangle-biopython.git
cd pydangle-biopython

python -m venv .venv
source .venv/bin/activate    # Linux / macOS
# .venv\Scripts\activate     # Windows

make install
```

`make install` does two things: installs the package with dev dependencies
(`pip install -e ".[dev]"`) and configures the git changelog hook.

## Usage

### Command line

After installation, the `pydangle-biopython` command is available:

```bash
# Show version
pydangle-biopython --version

# Default measurements (phi, psi, chi1–chi4)
pydangle-biopython structure.pdb

# Specify measurements
pydangle-biopython -c 'phi; psi; omega' structure.pdb

# Nucleic acid backbone
pydangle-biopython -c 'alpha; beta; gamma; delta; epsilon; zeta' rna.cif

# Force mmCIF format
pydangle-biopython --mmcif structure.cif

# Custom distance measurement
pydangle-biopython -c 'distance: Ca_Ca: i-1 _CA_, i _CA_' structure.pdb
```

### As a library

```python
from Bio.PDB import PDBParser
from pydangle_biopython import process_measurement_commands

parser = PDBParser(QUIET=True)
structure = parser.get_structure('X', 'structure.pdb')

lines = process_measurement_commands('structure.pdb', structure, 'phi; psi; omega')
for line in lines:
    print(line)
```

### Running without installing

```bash
PYTHONPATH=src python -m pydangle_biopython.cli structure.pdb
```

## Running tests

```bash
# Basic test run
pytest

# With coverage report
make test-cov

# Lint, type check, and test in one shot
make all
```

## Changelog automation

Every `git commit` automatically regenerates `CHANGELOG.md` from your commit
history.  The hook groups entries by tag.

To create a release:

```bash
git tag -a v0.1.0 -m "First release"
git commit --allow-empty -m "chore: trigger changelog for v0.1.0"
```

Commits prefixed with `chore:` or `wip:` are excluded from the changelog.

To skip the hook for a single commit:

```bash
SKIP_CHANGELOG=1 git commit -m "quick fix"
```
