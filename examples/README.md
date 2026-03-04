# pydangle-biopython examples

## Example structures

- `1ubq.pdb` — ubiquitin (76 residues, chain A)
- `1ubiH-1AtomRecs.pdb` — ubiquitin with hydrogens (no chain ID)
- `AGPQVS.pdb` — hexapeptide Ala-Gly-Pro-Gln-Val-Ser

## Quick start

```bash
# Default measurements (phi, psi, chi1–chi4)
pydangle-biopython examples/1ubq.pdb

# Backbone dihedrals with Ramachandran category
pydangle-biopython -c 'phi; psi; rama_category' examples/1ubq.pdb

# Find cis peptide bonds
pydangle-biopython -c 'omega; is_cis' examples/1ubq.pdb | grep True

# Pre-proline residues with backbone angles
pydangle-biopython -c 'phi; psi; is_prepro' examples/1ubq.pdb | grep True

# Custom Cα–Cα distance
pydangle-biopython -c 'distance: CaCa: i _CA_, i+1 _CA_' examples/AGPQVS.pdb

# Virtual Cα angle
pydangle-biopython -c 'angle: vCa: i-1 _CA_, i _CA_, i+1 _CA_' examples/AGPQVS.pdb

# Multiple files
pydangle-biopython -c 'phi; psi' examples/*.pdb
```

## Builtin reference

### Protein backbone dihedrals

| Name    | Type     | Definition                              |
|---------|----------|-----------------------------------------|
| phi     | dihedral | i-1 C, i N, i CA, i C                   |
| psi     | dihedral | i N, i CA, i C, i+1 N                   |
| omega   | dihedral | i-1 CA, i-1 C, i N, i CA               |
| tau     | angle    | i N, i CA, i C                          |

### Protein sidechain dihedrals

| Name | Type     | Definition                                      |
|------|----------|-------------------------------------------------|
| chi1 | dihedral | i N, i CA, i CB, i xG                           |
| chi2 | dihedral | i CA, i CB, i xG, i xD                          |
| chi3 | dihedral | i CB, i xG, i xD, i xE                          |
| chi4 | dihedral | i xG, i xD, i xE, i xZ                          |

(xG, xD, xE, xZ = gamma, delta, epsilon, zeta atoms matched by regex)

### Protein backbone distances and angles

| Name    | Type     | Definition                              |
|---------|----------|-----------------------------------------|
| pbCACAd | distance | i CA, i+1 CA                            |
| pbCACd  | distance | i CA, i C                               |
| pbCOd   | distance | i C, i O                                |
| pbCNd   | distance | i C, i+1 N                              |
| pbNCAd  | distance | i N, i CA                               |
| pbCACOa | angle    | i CA, i C, i O                          |
| pbCACNa | angle    | i CA, i C, i+1 N                        |
| pbOCNa  | angle    | i O, i C, i+1 N                         |
| pbCNCAa | angle    | i C, i+1 N, i+1 CA                      |
| vCAd    | distance | i-1 CA, i CA                            |
| vCAa    | angle    | i-1 CA, i CA, i+1 CA                    |
| vCAt    | dihedral | i-1 CA, i CA, i+1 CA, i+2 CA           |

### Nucleic acid backbone dihedrals

| Name    | Type     | Definition                              |
|---------|----------|-----------------------------------------|
| alpha   | dihedral | i-1 O3', i P, i O5', i C5'             |
| beta    | dihedral | i P, i O5', i C5', i C4'               |
| gamma   | dihedral | i O5', i C5', i C4', i C3'             |
| delta   | dihedral | i C5', i C4', i C3', i O3'             |
| epsilon | dihedral | i C4', i C3', i O3', i+1 P             |
| zeta    | dihedral | i C3', i O3', i+1 P, i+1 O5'           |
| eta     | dihedral | i-1 C4', i P, i C4', i+1 P             |
| theta   | dihedral | i P, i C4', i+1 P, i+1 C4'             |
| chi_na  | dihedral | i O4', i C1', i N9/N1, i C4/C2         |

### Nucleic acid sugar pucker

| Name | Type     | Definition                              |
|------|----------|-----------------------------------------|
| nu0  | dihedral | i C4', i O4', i C1', i C2'             |
| nu1  | dihedral | i O4', i C1', i C2', i C3'             |
| nu2  | dihedral | i C1', i C2', i C3', i C4'             |
| nu3  | dihedral | i C2', i C3', i C4', i O4'             |
| nu4  | dihedral | i C3', i C4', i O4', i C1'             |

### Residue classification labels

| Name          | Type  | Values                                           |
|---------------|-------|--------------------------------------------------|
| is_cis        | label | True / False / \_\_?\_\_                         |
| is_trans       | label | True / False / \_\_?\_\_                         |
| is_gly         | label | True / False                                     |
| is_pro         | label | True / False                                     |
| is_ileval      | label | True / False                                     |
| is_prepro      | label | True / False / \_\_?\_\_                         |
| has_all_mc     | label | True / False                                     |
| has_all_sc     | label | True / False / \_\_?\_\_                         |
| rama_category  | label | General / Gly / IleVal / TransPro / CisPro / PrePro |

Ramachandran categories follow the wwPDB Validation Task Force conventions
(Read et al., Structure 19:1395–1412, 2011).

## Custom measurements

Use the explicit syntax `function_type: label: atom_specs` for measurements
not covered by builtins:

```bash
# Distance between CB atoms
pydangle-biopython -c 'distance: CbCb: i _CB_, i+1 _CB_' structure.pdb

# Angle with regex atom matching
pydangle-biopython -c 'angle: test: i _N__, i _CA_, i /_C[BG]_/' structure.pdb

# Dihedral across two residues
pydangle-biopython -c 'dihedral: cross: i-1 _CA_, i _N__, i _CA_, i+1 _N__' structure.pdb
```
