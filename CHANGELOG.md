# Changelog

All notable changes to this project will be documented in this file.
This file is auto-generated from git history after each commit.

## Unreleased

- feat: add benchmark script comparing pydangle vs Java Dangle (3747e55)
- feat: add file-list, glob, and directory input methods (36833cd)
- docs: add DIFFERENCES.md and update documentation for dssp/rama variants (eff1a19)
- feat: add dssp and dssp3 labels via external mkdssp (bbd6a4d)
- feat: add rama3/rama4/rama5/rama6 Ramachandran category variants (3d0a2d9)
- feat: add DEVELOPMENT.md cheatsheet and make bump target (9ded5ce)
- fix: add --version flag, fix email typo in README (1df4ae7)
- docs: update README with residue classification labels (5e8632b)

## v0.3.0 — 2026-03-03

- refactor: replace synthetic hexapeptide with real ubiquitin coordinates (f8699ef)
- feat: add is_left/is_right chirality labels (83e4050)
- feat: add has_all_sc label for sidechain atom completeness (5cc303e)
- feat: add residue classification labels and examples documentation (952e058)

## v0.2.0 — 2026-03-03

- docs: update README with development story and Dangle differences (a042e27)
- fix: use Java Dangle residue-name whitelist for polymer filter (9b25c27)
- fix: use PPBuilder for chain gap detection, fix RNA fallback (07a88de)
- fix: correct omega definition to IUPAC convention (8213147)
- feat: accept multiple structure files on command line (e19c3dc)
- refactor: move tau to PROTEIN_BACKBONE, add pbNCAC alias tau is a bond angle (N-CA-C), not a dihedral. Moved to PROTEIN_BACKBONE with a pbNCAC alias following the naming convention of other backbone angle specs. The tau key is retained for Java Dangle compatibility. (0d3c5e0)
- feat: add PRO-specific chi3 definition (CB-CG-CD-N) (02ab062)
- fix: add type annotations for mypy strict mode (b0fa972)
- to correct lint errors ran ruff check --fix src/ tests/ (3323e13)
- feat: initial pydangle-biopython project (80e2978)

