---
name: TDA pipeline design
description: Planned TDA analysis approach - superlevel set cubical persistence via GUDHI on von Mises KDE density grids of Ramachandran data
type: project
---

**Goal:** Create a TDA-based parallel to the traditional adaptive nonparametric von Mises KDE
analysis of Ramachandran maps used for MolProbity validation (98%/99.5% contour boundaries).

**Approach — superlevel set cubical persistence:**
1. Pydangle extracts phi, psi, rama_category for ~2.5M residues from top2018 dataset
2. Stratify by rama category (General, Gly, IleVal, TransPro, CisPro, PrePro)
3. Compute von Mises KDE on periodic grid (360×360 or 720×720) per category
4. Run GUDHI CubicalComplex with periodic_dimensions=[True,True] for superlevel set persistence
5. Analyze persistence diagrams — H0 (density peaks/clusters), H1 (voids/forbidden regions)

**Why cubical persistence, not Ripser:** The density function on a grid is the natural
input. Ripser on 2.5M points would need ~50TB distance matrix. Cubical persistence
on a 720×720 grid takes seconds and ~50MB RAM.

**Why GUDHI:** Supports periodic boundary conditions (essential for toroidal phi/psi space).
Ripser doesn't do cubical complexes. GUDHI also provides persistence diagram analysis tools.

**Scientific value:** Persistence diagrams reveal at which density thresholds topological
transitions occur (region merging, void closing) — providing topologically-motivated
validation boundaries rather than arbitrary percentile choices.

**Estimated full pipeline time (genfive, parallel):** ~25 minutes total (20 min pydangle + seconds for TDA).

**Why:** Collaboration with Ezra Miller's group at Duke Math on TDA for Ramachandran analysis.

**How to apply:** When implementing the TDA pipeline, use GUDHI CubicalComplex with periodic
boundaries. GPU (Tesla M40) does not help this workload. Genfive is the recommended platform.
