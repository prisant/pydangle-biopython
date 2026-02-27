"""Tests for the pydangle command string parser."""

import pytest

from pydangle_biopython.builtins import BUILTIN_COMMANDS
from pydangle_biopython.parser import (
    _tokenize,
    command_string_parser,
    expand_command_list,
    parse_alternative_arg_lists,
    parse_atom_list,
    parse_atom_position,
    parse_command_fields,
)


# ---------------------------------------------------------------------------
# _tokenize
# ---------------------------------------------------------------------------

class TestTokenize:
    def test_semicolon_split(self):
        result = _tokenize("phi; psi; omega", ';')
        assert result == ['phi', 'psi', 'omega']

    def test_colon_split(self):
        result = _tokenize("dihedral: phi: i-1 _C__, i _N__", ':')
        assert result == ['dihedral', 'phi', 'i-1 _C__, i _N__']

    def test_strips_whitespace(self):
        result = _tokenize("  a  ;  b  ;  c  ", ';')
        assert result == ['a', 'b', 'c']

    def test_single_token(self):
        result = _tokenize("phi", ';')
        assert result == ['phi']


# ---------------------------------------------------------------------------
# expand_command_list
# ---------------------------------------------------------------------------

class TestExpandCommandList:
    def test_builtin_expansion(self):
        result = expand_command_list("phi")
        assert len(result) == 1
        assert 'dihedral' in result[0].lower()

    def test_multiple_builtins(self):
        result = expand_command_list("phi; psi; omega")
        assert len(result) == 3

    def test_passthrough_custom_command(self):
        custom = "distance: myDist: i _CA_, i+1 _CA_"
        result = expand_command_list(custom)
        assert result == [custom]

    def test_mixed_builtin_and_custom(self):
        result = expand_command_list("phi; distance: myDist: i _CA_, i+1 _CA_")
        assert len(result) == 2
        assert 'dihedral' in result[0].lower()
        assert 'myDist' in result[1]

    def test_empty_tokens_skipped(self):
        result = expand_command_list("phi; ; psi")
        assert len(result) == 2

    def test_all_protein_builtins_expand(self):
        """Every entry in BUILTIN_COMMANDS should expand without error."""
        for name in BUILTIN_COMMANDS:
            result = expand_command_list(name)
            assert len(result) == 1
            assert ':' in result[0], f"Builtin {name!r} did not expand properly"

    def test_nucleic_acid_builtins_present(self):
        """Key nucleic acid builtins should be present."""
        for name in ('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                     'eta', 'theta', 'chi_na'):
            assert name in BUILTIN_COMMANDS, f"Missing nucleic acid builtin: {name}"


# ---------------------------------------------------------------------------
# parse_command_fields
# ---------------------------------------------------------------------------

class TestParseCommandFields:
    def test_standard_command(self):
        key, label, args = parse_command_fields(
            "dihedral: phi: i-1 _C__, i _N__, i _CA_, i _C__"
        )
        assert key == 'dihedral'
        assert label == 'phi'
        assert '_C__' in args

    def test_synonym_keys(self):
        for syn in ('tor', 'tors', 'torsion', 'dihedral'):
            cmd = f"{syn}: test: i _N__, i _CA_, i _CB_, i _CG_"
            key, _, _ = parse_command_fields(cmd)
            assert key == 'dihedral'

        for syn in ('dis', 'dist', 'distance'):
            key, _, _ = parse_command_fields(f"{syn}: test: i _CA_, i+1 _CA_")
            assert key == 'distance'

        for syn in ('ang', 'angl', 'angle'):
            key, _, _ = parse_command_fields(f"{syn}: test: i _N__, i _CA_, i _C__")
            assert key == 'angle'

    def test_invalid_key_raises(self):
        with pytest.raises(ValueError, match="Unknown function key"):
            parse_command_fields("bogus: test: i _CA_, i+1 _CA_")

    def test_wrong_field_count_raises(self):
        with pytest.raises(ValueError, match="3 colon-separated fields"):
            parse_command_fields("dihedral: phi")


# ---------------------------------------------------------------------------
# parse_alternative_arg_lists
# ---------------------------------------------------------------------------

class TestParseAlternativeArgLists:
    def test_single_list(self):
        result = parse_alternative_arg_lists("i _N__, i _CA_, i _CB_, i _CG_")
        assert len(result) == 1

    def test_two_alternatives(self):
        result = parse_alternative_arg_lists(
            "i _O4*, i _C1*, i _N9_, i _C4_ | i _O4*, i _C1*, i _N1_, i _C2_"
        )
        assert len(result) == 2
        assert '_N9_' in result[0]
        assert '_N1_' in result[1]


# ---------------------------------------------------------------------------
# parse_atom_list
# ---------------------------------------------------------------------------

class TestParseAtomList:
    def test_distance_two_atoms(self):
        result = parse_atom_list('distance', "i _CA_, i+1 _CA_")
        assert len(result) == 2

    def test_angle_three_atoms(self):
        result = parse_atom_list('angle', "i _N__, i _CA_, i _C__")
        assert len(result) == 3

    def test_dihedral_four_atoms(self):
        result = parse_atom_list('dihedral', "i-1 _C__, i _N__, i _CA_, i _C__")
        assert len(result) == 4

    def test_wrong_count_raises(self):
        with pytest.raises(ValueError, match="requires 4 atoms"):
            parse_atom_list('dihedral', "i _N__, i _CA_")


# ---------------------------------------------------------------------------
# parse_atom_position
# ---------------------------------------------------------------------------

class TestParseAtomPosition:
    def test_explicit_offset_zero(self):
        offset, regex = parse_atom_position("i _CA_")
        assert offset == 0
        assert regex.match(" CA ")

    def test_implicit_offset(self):
        offset, regex = parse_atom_position("_CA_")
        assert offset == 0

    def test_positive_offset(self):
        offset, regex = parse_atom_position("i+1 _CA_")
        assert offset == 1

    def test_negative_offset(self):
        offset, regex = parse_atom_position("i-1 _C__")
        assert offset == -1

    def test_large_offset(self):
        offset, regex = parse_atom_position("i+2 _CA_")
        assert offset == 2

    def test_underscore_to_space(self):
        offset, regex = parse_atom_position("i _CA_")
        assert regex.match(" CA ")
        assert not regex.match("_CA_")

    def test_regex_pattern(self):
        offset, regex = parse_atom_position("i /_[ACNOS]G[_1]/")
        # Should match atoms like " CG " or " CG1" or " OG " etc.
        assert regex.match(" CG ")
        assert regex.match(" CG1")
        assert regex.match(" OG ")
        assert regex.match(" SG ")
        assert not regex.match(" CD ")

    def test_nucleic_acid_star(self):
        """The * in nucleic acid atom names should match both * and '."""
        offset, regex = parse_atom_position("i _O3*")
        assert regex.match(" O3*")
        assert regex.match(" O3'")

    def test_bare_i_offset(self):
        offset, regex = parse_atom_position("i _N__")
        assert offset == 0


# ---------------------------------------------------------------------------
# command_string_parser (integration)
# ---------------------------------------------------------------------------

class TestCommandStringParser:
    def test_single_builtin(self):
        result = command_string_parser("phi")
        assert len(result) == 1
        fun_key, label, arg_lists = result[0]
        assert fun_key == 'dihedral'
        assert label == 'phi'
        assert len(arg_lists) == 1      # one alternative
        assert len(arg_lists[0]) == 4   # four atoms for dihedral

    def test_multiple_builtins(self):
        result = command_string_parser("phi; psi; omega")
        assert len(result) == 3
        labels = [cmd[1] for cmd in result]
        assert labels == ['phi', 'psi', 'omega']

    def test_custom_distance(self):
        result = command_string_parser("distance: myDist: i-1 _CA_, i _CA_")
        assert len(result) == 1
        assert result[0][0] == 'distance'
        assert len(result[0][2][0]) == 2

    def test_alternative_arg_lists(self):
        """chi for nucleic acids has two alternatives (purine/pyrimidine)."""
        result = command_string_parser("chi_na")
        assert len(result) == 1
        fun_key, label, arg_lists = result[0]
        assert fun_key == 'dihedral'
        assert len(arg_lists) == 2  # two alternatives

    def test_all_builtins_parse(self):
        """Every builtin should parse without raising."""
        for name in BUILTIN_COMMANDS:
            result = command_string_parser(name)
            assert len(result) == 1, f"Builtin {name!r} failed to parse"
