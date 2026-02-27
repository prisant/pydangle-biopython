"""Tests for pydangle builtin command definitions."""

import pytest

from pydangle_biopython.builtins import (
    BUILTIN_COMMANDS,
    NUCLEIC_ACID_ANGLES,
    NUCLEIC_ACID_BACKBONE,
    PROTEIN_BACKBONE,
    PROTEIN_DIHEDRALS,
    PROTEIN_SIDECHAIN,
)
from pydangle_biopython.parser import command_string_parser


class TestBuiltinDefinitions:
    """Verify that all builtin definitions are well-formed and parseable."""

    def test_no_duplicate_keys(self):
        """All sub-dictionaries should have unique keys."""
        all_keys = []
        for d in (PROTEIN_BACKBONE, PROTEIN_DIHEDRALS, PROTEIN_SIDECHAIN,
                  NUCLEIC_ACID_BACKBONE, NUCLEIC_ACID_ANGLES):
            all_keys.extend(d.keys())
        assert len(all_keys) == len(set(all_keys)), "Duplicate builtin keys found"

    def test_combined_dict_has_all_entries(self):
        """BUILTIN_COMMANDS should contain all individual dictionaries."""
        expected_count = (
            len(PROTEIN_BACKBONE) + len(PROTEIN_DIHEDRALS) +
            len(PROTEIN_SIDECHAIN) + len(NUCLEIC_ACID_BACKBONE) +
            len(NUCLEIC_ACID_ANGLES)
        )
        assert len(BUILTIN_COMMANDS) == expected_count

    @pytest.mark.parametrize("name", list(BUILTIN_COMMANDS.keys()))
    def test_builtin_parses(self, name):
        """Each builtin should parse without error."""
        result = command_string_parser(name)
        assert len(result) == 1
        fun_key, label, arg_lists = result[0]
        assert fun_key in ('distance', 'angle', 'dihedral')
        assert label  # label should be non-empty
        assert len(arg_lists) >= 1

    @pytest.mark.parametrize("name", list(BUILTIN_COMMANDS.keys()))
    def test_builtin_has_three_fields(self, name):
        """Each builtin command string should have exactly 3 colon-separated fields."""
        cmd_str = BUILTIN_COMMANDS[name]
        fields = [f.strip() for f in cmd_str.split(':')]
        assert len(fields) == 3, (
            f"Builtin {name!r} has {len(fields)} fields instead of 3: {cmd_str!r}"
        )


class TestProteinBuiltins:
    """Spot-check protein-specific builtins."""

    def test_phi_is_dihedral(self):
        result = command_string_parser("phi")
        assert result[0][0] == 'dihedral'
        assert len(result[0][2][0]) == 4

    def test_psi_is_dihedral(self):
        result = command_string_parser("psi")
        assert result[0][0] == 'dihedral'

    def test_omega_is_dihedral(self):
        result = command_string_parser("omega")
        assert result[0][0] == 'dihedral'

    def test_tau_is_angle(self):
        result = command_string_parser("tau")
        assert result[0][0] == 'angle'
        assert len(result[0][2][0]) == 3

    def test_chi1_through_chi4(self):
        for chi in ('chi1', 'chi2', 'chi3', 'chi4'):
            result = command_string_parser(chi)
            assert result[0][0] == 'dihedral'
            assert len(result[0][2][0]) == 4


class TestNucleicAcidBuiltins:
    """Spot-check nucleic acid builtins."""

    def test_backbone_dihedrals(self):
        for name in ('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta'):
            result = command_string_parser(name)
            assert result[0][0] == 'dihedral', f"{name} should be a dihedral"
            assert len(result[0][2][0]) == 4

    def test_pseudotorsions(self):
        for name in ('eta', 'theta'):
            result = command_string_parser(name)
            assert result[0][0] == 'dihedral'
            assert len(result[0][2][0]) == 4

    def test_chi_has_alternatives(self):
        """Nucleic acid chi should have two alternatives for purine/pyrimidine."""
        result = command_string_parser("chi_na")
        assert result[0][0] == 'dihedral'
        assert len(result[0][2]) == 2  # two alternatives

    def test_sugar_pucker_dihedrals(self):
        for name in ('nu0', 'nu1', 'nu2', 'nu3', 'nu4'):
            result = command_string_parser(name)
            assert result[0][0] == 'dihedral'
            assert len(result[0][2][0]) == 4
