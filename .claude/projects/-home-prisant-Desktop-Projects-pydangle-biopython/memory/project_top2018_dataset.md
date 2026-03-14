---
name: top2018 dataset locations and timing
description: Locations of top2018 70% homology filtered PDB datasets and empirical processing time estimates for pydangle-biopython
type: project
---

The top2018 70% homology mc+sc filtered datasets are at:

- **Full filtered (mc+sc, production):** `/home/prisant/Desktop/Backups/prisant/Desktop/FullFilteredResidues2025/top2018_pdbs_full_filtered_hom70`
  - 12,125 pruned files (primary production target, `*_pruned_all.pdb`)
  - 11,840 original PDBs
  - Average file size: ~488 KB, 10.9 GB total
  - Hierarchical structure: `[2-char]/[pdbcode]/[pdbcode].pdb`

- **MC filtered:** `/home/prisant/Desktop/Backups/prisant/Desktop/ClevelandMCFilteredResidues/top2018_pdbs_mc_filtered_hom70`
  - 13,308 original + 13,677 pruned MC = 26,985 total

**Empirical timing** for `phi; psi; omega; tau; chi1; dssp; rama_category`:
- 1.50 s/file in batch mode (100-file sample, 2m30s)
- Without DSSP: ~1.0 s/file
- Production estimate (12,125 pruned files): **~5.1 hours serial, ~40 min with 8-core parallelization**

**Why:** This dataset is the target for TDA/Ramachandran analysis work with Ezra Miller's group at Duke Math.

**How to apply:** When discussing processing the top2018 set, use these estimates. The pruned files are the production target.
