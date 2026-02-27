# pydangle-biopython

Compute distances, angles, and dihedral angles in protein and nucleic acid
structures using BioPython.

## Overview

`pydangle-biopython` is a Python translation of the Java
[Dangle](https://github.com/rlabduke/javadev) program from the Richardson
Laboratory at Duke University.  It uses BioPython as the structure-parsing
backend.  Sibling packages `pydangle-cctbx` and `pydangle-gemmi` will
provide the same functionality with alternative backends.

## Features

- **Protein measurements** — phi, psi, omega, tau, chi1–chi4, backbone
  bond lengths and angles, virtual Cα geometry
- **Nucleic acid measurements** — alpha through zeta backbone dihedrals,
  eta/theta pseudo-torsions, glycosidic chi, sugar pucker (nu0–nu4)
- **Custom measurements** — specify arbitrary distances, angles, and
  dihedrals using atom name patterns and relative residue offsets
- **PDB and mmCIF input** — auto-detected from file extension or forced
  via command-line flags
- **Fully tested** — pytest suite with sample protein and RNA structures

## Quick start

```bash
pip install -e "."
pydangle-biopython structure.pdb
```

Output:

```
# pydangle: Query string: phi; psi; chi1; chi2; chi3; chi4
structure.pdb:1:A:   2: :GLY:-166.362:174.478:__?__:__?__:__?__:__?__
structure.pdb:1:A:   3: :PRO:-59.004:-27.599:...
```

See [Getting Started](getting-started.md) for full setup instructions.

## Project layout

```
pydangle-biopython/
├── pyproject.toml
├── Makefile
├── mkdocs.yml
├── CHANGELOG.md
├── docs/
│   ├── index.md
│   ├── getting-started.md
│   ├── api.md
│   └── contributing.md
├── scripts/
│   ├── generate_changelog.sh
│   └── hooks/
│       └── post-commit
├── src/
│   └── pydangle_biopython/
│       ├── __init__.py
│       ├── py.typed
│       ├── builtins.py
│       ├── cli.py
│       ├── exceptions.py
│       ├── measure.py
│       └── parser.py
├── tests/
│   ├── __init__.py
│   ├── conftest.py
│   ├── test_builtins.py
│   ├── test_cli.py
│   ├── test_measure.py
│   └── test_parser.py
└── examples/
    └── AGPQVS.pdb
```
