---
name: development machines
description: Three development platforms (nuc81, x99i7, genfive) with specs and relative performance for pydangle benchmarking
type: reference
---

| Machine | CPU | Cores/Threads | Boost | RAM | Storage | OS |
|---------|-----|---------------|-------|-----|---------|-----|
| nuc81 | i7-8559U (Coffee Lake 2018) | 4C/8T | 4.5 GHz | 64 GB | Samsung 860 EVO SATA | Mint 22.3 |
| x99i7 | i7-6850K (Broadwell-E 2016) | 6C/12T | 4.0 GHz | 64 GB | Timetec 1TB SATA | Mint 21.3 |
| genfive | i7-12700H (Alder Lake 2022) | 6P+8E/20T | 4.7 GHz | 64 GB | Samsung 980 PRO NVMe | Mint 22.3 |

**Relative single-thread performance** (nuc81 = 1.0): x99i7 ~0.80, genfive ~1.5

**Pydangle top2018 estimates (12,125 pruned files, phi/psi/omega/tau/chi1/dssp/rama):**

| Machine | Serial | Parallel |
|---------|--------|----------|
| nuc81 | 5.1 hrs | ~61 min |
| x99i7 | 6.3 hrs | ~47 min |
| genfive | 3.4 hrs | ~20 min |

Recommendation: use genfive for production runs. All benchmarks in BENCHMARK.md were measured on nuc81.
