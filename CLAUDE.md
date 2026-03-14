# CLAUDE.md — Project context for Claude Code

## Project

pydangle-biopython: Python reimagining of Java Dangle for computing distances,
angles, dihedrals, per-residue classifications, and DSSP secondary structure
in protein and nucleic acid structures. Uses BioPython as the structure-parsing
backend.

**Author:** Michael G. Prisant, Ph.D. (info@prisantscientific.org)
**Repository:** https://github.com/prisant/pydangle-biopython
**License:** MIT

## Development workflow

- **Branch:** Work directly on `main` for now.
- **Quality gate:** Always run `make all` (lint + typecheck + test) before committing.
- **Lint fix:** Run `ruff check --fix src/ tests/` before `make all` if needed.
- **Commit style:** Use conventional commits (feat:, fix:, docs:, refactor:, chore:).
- **Push:** Push to origin after each commit.

## Build commands

```bash
make all          # ruff check + mypy + pytest (the full quality gate)
make test         # pytest only
make lint         # ruff check only
make typecheck    # mypy only
make clean        # remove build artifacts and editor backups
make install      # editable install + dev deps + git hooks
```

## Architecture

### Source layout

```
src/pydangle_biopython/
├── cli.py          # Command-line interface (argparse)
├── fileinput.py    # File input collection (file lists, glob, directory)
├── parser.py       # Measurement command string parser
├── builtins.py     # Built-in measurement definitions
├── measure.py      # Measurement computation and output pipeline
├── labels.py       # Residue classification labels + LABEL_REGISTRY
├── dssp.py         # DSSP secondary structure via external mkdssp
└── exceptions.py   # Custom exceptions
```

### Measurement types

- **distance / angle / dihedral** — geometric measurements from atom coordinates
- **label** — categorical per-residue classifications (no atom specs, 2-field parser syntax)
- **property** (planned) — numeric per-residue values (cbdev, B-factors)

### Adding a new label

1. Write the function in `labels.py` with signature `(residue_list, index, unknown) -> str`
2. Add to `LABEL_REGISTRY` in `labels.py`
3. Add builtin entry in `builtins.py` `RESIDUE_LABELS` dict
4. Add tests in `tests/test_labels.py`
5. Update `examples/README.md` builtin reference table
6. Update `README.md` and `DIFFERENCES.md` as needed

### Test fixtures

- **hexapeptide** — real coordinates from ubiquitin residues 34-39 (GLU GLY ILE PRO PRO ASP), defined in `tests/conftest.py`
- **1ubq.pdb** — full ubiquitin structure in `examples/`
- **RNA dinucleotide** — minimal RNA fixture in `tests/conftest.py`

### External dependencies

- **BioPython** — structure parsing (pip installed)
- **mkdssp** — DSSP secondary structure (system installed: `apt install dssp`)

## Current label builtins (17)

rama_category/rama6/rama5/rama4/rama3, dssp/dssp3,
is_cis/is_trans/peptide_bond, is_gly/is_pro/is_ileval/is_prepro,
has_all_mc/has_all_sc, is_left/is_right/chirality

## Planned work

- **property type** — numeric per-residue values (cbdev, mc_bfactor, sc_bfactor)
- **Input improvements** — ~~file lists, piping~~ (done), web fetch
- ~~**JSON output**~~ — JSONL format via `-o jsonl` (done)
- ~~**Parallelization**~~ — multiprocessing support (done: `-j`/`--jobs` flag)
- **More Java Dangle features** — see DIFFERENCES.md

## Style preferences

- Do not use bullet points or numbered lists for explanations in conversation.
- Always run `ruff check --fix` before `make all`.
- Update all documentation (README.md, examples/README.md, DIFFERENCES.md) with each feature addition.
- Test against `examples/1ubq.pdb` for smoke tests.
- Use conventional commit messages.
