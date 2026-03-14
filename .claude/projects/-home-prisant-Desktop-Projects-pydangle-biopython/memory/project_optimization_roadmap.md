---
name: performance optimization roadmap
description: Agreed optimization strategy for pydangle - parallelize first, then conda, then gemmi sibling
type: project
---

Performance optimization roadmap (discussed 2026-03-12):

1. **Parallelize pydangle-biopython** (first priority, ~0.5 day)
   - Add `multiprocessing.Pool` to batch processing in `cli.py`
   - Expected ~6-7x speedup on 8-core machine
   - Low risk — embarrassingly parallel, each file independent

2. **Conda packaging** (second priority, ~1 day)
   - Move to conda for all pydangle siblings
   - Key benefit: mkdssp and cctbx installable without root via conda-forge
   - Maintain dual support (pip still works, DSSP degrades gracefully)
   - `conda install -c conda-forge dssp` provides mkdssp

3. **Create pydangle-gemmi** (medium term, 2-3 days)
   - Gemmi C++ parser is 5-20x faster than BioPython
   - Combined with parallelization: 12K files in ~5-10 minutes
   - Share measurement logic from pydangle-biopython where possible
   - Gemmi available via `conda install -c conda-forge gemmi`

4. **Optimize BioPython code** — skip (ceiling ~1.2x, bottleneck is BioPython's parser)

**Why:** Serial processing of 12K files takes ~5 hours. Parallelization alone cuts to ~40 min. Gemmi+parallel could reach 5-10 min.

**How to apply:** When the user wants to start optimization work, begin with parallelization in cli.py. Don't waste time optimizing BioPython parsing code.
