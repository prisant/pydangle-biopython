"""Built-in measurement command definitions for proteins and nucleic acids.

Each entry maps a short name to its expanded command string in the form:
    "function_type: label: atom_spec1, atom_spec2 [, atom_spec3 [, atom_spec4]]"

Atom specifications use:
    - Residue offset notation: i-1, i, i+1, i+2 (default is i)
    - Underscore (_) in place of spaces in atom names
    - /regexp/ for regular expression matching of atom names

Label builtins use the form:
    "label: label_name"
with no atom specifications; the label function determines what
structural data it needs internally.
"""

# ---------------------------------------------------------------------------
# Protein backbone measurements
# ---------------------------------------------------------------------------
PROTEIN_BACKBONE = {
    # Backbone bond distances
    "pbCACAd": "distance: pbCACAd: i _CA_, i+1 _CA_",
    "pbCACd": "distance: pbCACd:  i _CA_, i _C__",
    "pbCOd": "distance: pbCOd:   i _C__, i _O__",
    "pbCNd": "distance: pbCNd:   i _C__, i+1 _N__",
    "pbNCAd": "distance: pbNCAd:  i _N__, i _CA_",
    # Backbone bond angles
    "pbCACOa": "angle: pbCACOa: i _CA_, i _C__, i _O__",
    "pbCACNa": "angle: pbCACNa: i _CA_, i _C__, i+1 _N__",
    "pbOCNa": "angle: pbOCNa:  i _O__, i _C__, i+1 _N__",
    "pbCNCAa": "angle: pbCNCAa: i _C__, i+1 _N__, i+1 _CA_",
    "pbNCAC": "angle: tau:     i _N__,   i _CA_, i _C__",
    "tau": "angle: tau:   i _N__,   i _CA_, i _C__",
    # Virtual Cα measurements
    "vCAd": "distance: vCAd: i-1 _CA_, i _CA_",
    "vCAa": "angle:    vCAa: i-1 _CA_, i _CA_, i+1 _CA_",
    "vCAt": "dihedral: vCAt: i-1 _CA_, i _CA_, i+1 _CA_, i+2 _CA_",
}

# ---------------------------------------------------------------------------
# Protein backbone dihedral angles
# ---------------------------------------------------------------------------
PROTEIN_DIHEDRALS = {
    "phi": "dihedral: phi:   i-1 _C__, i _N__, i _CA_, i _C__",
    "psi": "dihedral: psi:   i _N__,   i _CA_, i _C__, i+1 _N__",
    "omega": "dihedral: omega: i-1 _CA_, i-1 _C__, i _N__,  i _CA_",
}

# ---------------------------------------------------------------------------
# Protein sidechain measurements
# ---------------------------------------------------------------------------
PROTEIN_SIDECHAIN = {
    "chi1": "dihedral: chi1:  i _N__,   i _CA_, i _CB_, i /_[ACNOS]G[_1]/",
    "chi2": "dihedral: chi2:  i _CA_,   i _CB_, i /_[ACNOS]G[_1]/, i /_[ACNOS]D[_1]/",
    "chi3": (
        "dihedral: chi3:"
        " i _CB_, i /_[ACNOS]G[_1]/,"
        " i /_[ACNOS]D[_1]/, i /_[ACNOS]E[_1]/"
        " | i _CB_, i _CG_, i _CD_, i _N__"
    ),
    "chi4": (
        "dihedral: chi4:"
        " i /_[ACNOS]G[_1]/, i /_[ACNOS]D[_1]/,"
        " i /_[ACNOS]E[_1]/, i /_[ACNOS]Z[_1]/"
    ),
    "scABG": "angle:    scABG: i _CA_,   i _CB_, i /_[ACNOS]G[_1]/",
    "scBGD": "angle:    scBGD: i _CB_,   i /_[ACNOS]G[_1]/, i /_[ACNOS]D[_1]/",
    "scGDE": "angle:    scGDE: i /_[ACNOS]G[_1]/, i /_[ACNOS]D[_1]/, i /_[ACNOS]E[_1]/",
    "scDEZ": "angle:    scDEZ: i /_[ACNOS]D[_1]/, i /_[ACNOS]E[_1]/, i /_[ACNOS]Z[_1]/",
}

# ---------------------------------------------------------------------------
# RNA/DNA backbone dihedral angles
#
# Standard nucleic acid backbone: alpha through zeta
#   alpha:   O3'(i-1) - P    - O5'  - C5'
#   beta:    P        - O5'  - C5'  - C4'
#   gamma:   O5'      - C5'  - C4'  - C3'
#   delta:   C5'      - C4'  - C3'  - O3'
#   epsilon: C4'      - C3'  - O3'  - P(i+1)
#   zeta:    C3'      - O3'  - P(i+1) - O5'(i+1)
#
# Pseudotorsion angles for RNA:
#   eta:     C4'(i-1) - P    - C4'  - P(i+1)
#   theta:   P        - C4'  - P(i+1) - C4'(i+1)
#
# Glycosidic torsion:
#   chi (purines):   O4' - C1' - N9  - C4
#   chi (pyrimidines): O4' - C1' - N1  - C2
#
# Note: Atom names use * in place of ' for compatibility.
# Both conventions (* and ') are matched by the regex patterns.
# ---------------------------------------------------------------------------
NUCLEIC_ACID_BACKBONE = {
    "alpha": "dihedral: alpha:   i-1 _O3*, i _P__, i _O5*, i _C5*",
    "beta": "dihedral: beta:    i _P__,   i _O5*, i _C5*, i _C4*",
    "gamma": "dihedral: gamma:   i _O5*,   i _C5*, i _C4*, i _C3*",
    "delta": "dihedral: delta:   i _C5*,   i _C4*, i _C3*, i _O3*",
    "epsilon": "dihedral: epsilon: i _C4*,   i _C3*, i _O3*, i+1 _P__",
    "zeta": "dihedral: zeta:    i _C3*,   i _O3*, i+1 _P__, i+1 _O5*",
    # Pseudotorsion angles
    "eta": "dihedral: eta:     i-1 _C4*, i _P__, i _C4*, i+1 _P__",
    "theta": "dihedral: theta:   i _P__,   i _C4*, i+1 _P__, i+1 _C4*",
    # Glycosidic torsion: use | (or) to handle purines and pyrimidines
    "chi_na": (
        "dihedral: chi: i _O4*, i _C1*, i _N9_, i _C4_ | i _O4*, i _C1*, i _N1_, i _C2_"
    ),
}

# ---------------------------------------------------------------------------
# RNA/DNA sugar pucker and backbone angles
# ---------------------------------------------------------------------------
NUCLEIC_ACID_ANGLES = {
    # Sugar pucker pseudo-rotation related
    "nu0": "dihedral: nu0: i _C4*, i _O4*, i _C1*, i _C2*",
    "nu1": "dihedral: nu1: i _O4*, i _C1*, i _C2*, i _C3*",
    "nu2": "dihedral: nu2: i _C1*, i _C2*, i _C3*, i _C4*",
    "nu3": "dihedral: nu3: i _C2*, i _C3*, i _C4*, i _O4*",
    "nu4": "dihedral: nu4: i _C3*, i _C4*, i _O4*, i _C1*",
}

# ---------------------------------------------------------------------------
# Residue classification labels
#
# These use the ``label`` function type rather than geometric measurements.
# The label functions are defined in pydangle_biopython.labels and
# dispatched by name via the LABEL_REGISTRY.
# ---------------------------------------------------------------------------
RESIDUE_LABELS = {
    "is_cis": "label: is_cis",
    "is_trans": "label: is_trans",
    "peptide_bond": "label: peptide_bond",
    "is_gly": "label: is_gly",
    "is_pro": "label: is_pro",
    "is_ileval": "label: is_ileval",
    "is_prepro": "label: is_prepro",
    "has_all_mc": "label: has_all_mc",
    "has_all_sc": "label: has_all_sc",
    "is_left": "label: is_left",
    "is_right": "label: is_right",
    "chirality": "label: chirality",
    "rama_category": "label: rama_category",
    "rama6": "label: rama6",
    "rama5": "label: rama5",
    "rama4": "label: rama4",
    "rama3": "label: rama3",
    "dssp": "label: dssp",
    "dssp3": "label: dssp3",
}

# ---------------------------------------------------------------------------
# Combined dictionary of all builtins
# ---------------------------------------------------------------------------
BUILTIN_COMMANDS: dict[str, str] = {}
BUILTIN_COMMANDS.update(PROTEIN_BACKBONE)
BUILTIN_COMMANDS.update(PROTEIN_DIHEDRALS)
BUILTIN_COMMANDS.update(PROTEIN_SIDECHAIN)
BUILTIN_COMMANDS.update(NUCLEIC_ACID_BACKBONE)
BUILTIN_COMMANDS.update(NUCLEIC_ACID_ANGLES)
BUILTIN_COMMANDS.update(RESIDUE_LABELS)
