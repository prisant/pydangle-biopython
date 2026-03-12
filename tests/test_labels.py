"""Tests for pydangle residue classification labels."""

from pydangle_biopython.labels import (
    LABEL_REGISTRY,
    _compute_omega,
    label_has_all_mc,
    label_has_all_sc,
    label_is_cis,
    label_is_gly,
    label_is_ileval,
    label_is_left,
    label_is_prepro,
    label_is_pro,
    label_is_right,
    label_is_trans,
    label_rama3,
    label_rama4,
    label_rama5,
    label_rama_category,
)
from pydangle_biopython.measure import process_measurement_commands
from pydangle_biopython.parser import command_string_parser

UNK = "__?__"

# Hexapeptide sequence (from conftest.py, ubiquitin 34-39):
#   index 0: GLU
#   index 1: GLY
#   index 2: ILE  (pre-proline)
#   index 3: PRO
#   index 4: PRO
#   index 5: ASP


class TestLabelRegistry:
    """Verify the label registry is well-formed."""

    def test_all_labels_registered(self):
        expected = {
            "is_cis",
            "is_trans",
            "is_gly",
            "is_pro",
            "is_ileval",
            "is_prepro",
            "has_all_mc",
            "has_all_sc",
            "is_left",
            "is_right",
            "rama_category",
            "rama6",
            "rama3",
            "rama4",
            "rama5",
        }
        assert set(LABEL_REGISTRY.keys()) == expected

    def test_all_labels_callable(self):
        for name, func in LABEL_REGISTRY.items():
            assert callable(func), f"{name} is not callable"


class TestLabelParsing:
    """Verify label builtins parse correctly."""

    def test_parse_is_cis(self):
        result = command_string_parser("is_cis")
        assert len(result) == 1
        fun_key, label, arg_lists = result[0]
        assert fun_key == "label"
        assert label == "is_cis"
        assert arg_lists == []

    def test_parse_rama_category(self):
        result = command_string_parser("rama_category")
        assert len(result) == 1
        assert result[0][0] == "label"
        assert result[0][1] == "rama_category"

    def test_parse_mixed_commands(self):
        """Labels can be mixed with geometric measurements."""
        result = command_string_parser("phi; psi; rama_category; is_cis")
        assert len(result) == 4
        assert result[0][0] == "dihedral"  # phi
        assert result[1][0] == "dihedral"  # psi
        assert result[2][0] == "label"  # rama_category
        assert result[3][0] == "label"  # is_cis

    def test_parse_explicit_label_syntax(self):
        """Explicit label: syntax should work."""
        result = command_string_parser("label: is_gly")
        assert len(result) == 1
        assert result[0][0] == "label"
        assert result[0][1] == "is_gly"

    def test_all_label_builtins_parse(self):
        for name in LABEL_REGISTRY:
            result = command_string_parser(name)
            assert len(result) == 1
            assert result[0][0] == "label", f"{name} did not parse as label"


class TestOmegaComputation:
    """Test the omega helper used by cis/trans labels."""

    def test_first_residue_returns_none(self, hexapeptide_chain):
        residue_list = list(hexapeptide_chain.get_residues())
        assert _compute_omega(residue_list, 0) is None

    def test_second_residue_computes(self, hexapeptide_chain):
        residue_list = list(hexapeptide_chain.get_residues())
        omega = _compute_omega(residue_list, 1)
        assert omega is not None
        # Omega should be roughly +/- 180 for trans peptide bonds
        assert abs(abs(omega) - 180.0) < 30.0


class TestPrimitiveLabels:
    """Test individual label functions against the hexapeptide."""

    def test_is_gly_true(self, hexapeptide_chain):
        """Index 1 is GLY."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_gly(residue_list, 1, UNK) == "True"

    def test_is_gly_false(self, hexapeptide_chain):
        """Index 0 is GLU, not GLY."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_gly(residue_list, 0, UNK) == "False"

    def test_is_pro_true(self, hexapeptide_chain):
        """Index 3 is PRO."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_pro(residue_list, 3, UNK) == "True"

    def test_is_pro_false(self, hexapeptide_chain):
        """Index 0 is GLU, not PRO."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_pro(residue_list, 0, UNK) == "False"

    def test_is_ileval_true(self, hexapeptide_chain):
        """Index 2 is ILE."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_ileval(residue_list, 2, UNK) == "True"

    def test_is_ileval_false(self, hexapeptide_chain):
        """Index 0 is GLU, not ILE/VAL."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_ileval(residue_list, 0, UNK) == "False"

    def test_is_prepro_true(self, hexapeptide_chain):
        """Index 2 (ILE) precedes PRO at index 3."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_prepro(residue_list, 2, UNK) == "True"

    def test_is_prepro_false(self, hexapeptide_chain):
        """Index 0 (GLU) precedes GLY, not PRO."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_prepro(residue_list, 0, UNK) == "False"

    def test_is_prepro_last_residue(self, hexapeptide_chain):
        """Last residue has no successor; should return unknown."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert (
            label_is_prepro(
                residue_list,
                len(residue_list) - 1,
                UNK,
            )
            == UNK
        )

    def test_is_cis_first_residue(self, hexapeptide_chain):
        """First residue has no omega; should return unknown."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_cis(residue_list, 0, UNK) == UNK

    def test_is_trans_normal_peptide(self, hexapeptide_chain):
        """Normal peptide bonds should be trans."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_trans(residue_list, 1, UNK) == "True"

    def test_is_cis_normal_peptide(self, hexapeptide_chain):
        """Normal peptide bonds should not be cis."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_cis(residue_list, 1, UNK) == "False"


class TestRamaCategory:
    """Test Ramachandran category assignment."""

    def test_gly_category(self, hexapeptide_chain):
        """GLY at index 1 should be 'Gly'."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama_category(residue_list, 1, UNK) == "Gly"

    def test_ileval_category(self, hexapeptide_chain):
        """ILE at index 2 should be 'IleVal'."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama_category(residue_list, 2, UNK) == "IleVal"

    def test_transpro_category(self, hexapeptide_chain):
        """PRO at index 3 should be 'TransPro'."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama_category(residue_list, 3, UNK) == "TransPro"

    def test_general_category(self, hexapeptide_chain):
        """ASP at index 5 (last, not followed by PRO) should be 'General'."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama_category(residue_list, 5, UNK) == "General"

    def test_general_glu(self, hexapeptide_chain):
        """GLU at index 0 (followed by GLY, not PRO) should be 'General'."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama_category(residue_list, 0, UNK) == "General"

    def test_prepro_ile(self, hexapeptide_chain):
        """ILE at index 2 is IleVal, not PrePro (IleVal takes priority)."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama_category(residue_list, 2, UNK) == "IleVal"


class TestRamaVariants:
    """Test reduced Ramachandran category schemes."""

    def test_rama6_is_alias(self, hexapeptide_chain):
        """rama6 should produce same result as rama_category."""
        residue_list = list(hexapeptide_chain.get_residues())
        for i in range(len(residue_list)):
            assert label_rama_category(residue_list, i, UNK) == (
                LABEL_REGISTRY['rama6'](residue_list, i, UNK)
            )

    def test_rama3_gly(self, hexapeptide_chain):
        """GLY (index 1) should be 'Gly' in rama3."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama3(residue_list, 1, UNK) == 'Gly'

    def test_rama3_pro(self, hexapeptide_chain):
        """PRO (index 3) should be 'Pro' in rama3."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama3(residue_list, 3, UNK) == 'Pro'

    def test_rama3_ileval_maps_to_general(self, hexapeptide_chain):
        """ILE (index 2) maps to 'General' in rama3."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama3(residue_list, 2, UNK) == 'General'

    def test_rama4_prepro(self, hexapeptide_chain):
        """ILE (index 2) is IleVal in rama6 but PrePro check:
        ILE before PRO maps to General in rama4 (IleVal takes priority
        over PrePro in rama6, so it maps to General not PrePro)."""
        residue_list = list(hexapeptide_chain.get_residues())
        # ILE at index 2 is IleVal in rama6 -> General in rama4
        assert label_rama4(residue_list, 2, UNK) == 'General'

    def test_rama5_ileval(self, hexapeptide_chain):
        """ILE (index 2) should remain 'IleVal' in rama5."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama5(residue_list, 2, UNK) == 'IleVal'

    def test_rama5_pro(self, hexapeptide_chain):
        """PRO (index 3) should be 'Pro' in rama5 (not TransPro)."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_rama5(residue_list, 3, UNK) == 'Pro'

    def test_rama3_only_three_values(self, hexapeptide_chain):
        """All rama3 values should be General, Gly, or Pro."""
        residue_list = list(hexapeptide_chain.get_residues())
        valid = {'General', 'Gly', 'Pro'}
        for i in range(len(residue_list)):
            result = label_rama3(residue_list, i, UNK)
            assert result in valid, f"Index {i}: {result!r} not in {valid}"


class TestHasAllMc:
    """Test mainchain atom completeness."""

    def test_complete_residue(self, hexapeptide_chain):
        """Normal residues should have all mainchain atoms."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_has_all_mc(residue_list, 0, UNK) == "True"

    def test_all_residues_complete(self, hexapeptide_chain):
        """All hexapeptide residues should have complete mainchain."""
        residue_list = list(hexapeptide_chain.get_residues())
        for i in range(len(residue_list)):
            assert label_has_all_mc(residue_list, i, UNK) == "True"

    def test_has_all_mc_in_output(self, hexapeptide_structure):
        """has_all_mc should produce True/False values in output."""
        lines = process_measurement_commands(
            "test",
            hexapeptide_structure,
            "has_all_mc",
        )
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) > 0
        for line in data_lines:
            value = line.rsplit(":", 1)[-1]
            assert value in {"True", "False"}, f"Invalid has_all_mc value: {value!r}"


class TestHasAllSc:
    """Test sidechain atom completeness."""

    def test_glu_complete(self, hexapeptide_chain):
        """GLU (index 0) should have all sidechain atoms."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_has_all_sc(residue_list, 0, UNK) == "True"

    def test_gly_always_true(self, hexapeptide_chain):
        """GLY (index 1) has no sidechain; always True."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_has_all_sc(residue_list, 1, UNK) == "True"

    def test_pro_complete(self, hexapeptide_chain):
        """PRO (index 3) should have CB, CG, CD."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_has_all_sc(residue_list, 3, UNK) == "True"

    def test_ile_complete(self, hexapeptide_chain):
        """ILE (index 2) should have CB, CG1, CG2, CD1."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_has_all_sc(residue_list, 2, UNK) == "True"

    def test_has_all_sc_in_output(self, hexapeptide_structure):
        """has_all_sc should produce True/False values in output."""
        lines = process_measurement_commands(
            "test",
            hexapeptide_structure,
            "has_all_sc",
        )
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) > 0
        for line in data_lines:
            value = line.rsplit(":", 1)[-1]
            assert value in {"True", "False", UNK}, (
                f"Invalid has_all_sc value: {value!r}"
            )


class TestChirality:
    """Test Calpha chirality labels."""

    def test_gly_unknown(self, hexapeptide_chain):
        """GLY (index 1) has no CB; should return unknown."""
        residue_list = list(hexapeptide_chain.get_residues())
        assert label_is_left(residue_list, 1, UNK) == UNK
        assert label_is_right(residue_list, 1, UNK) == UNK

    def test_non_gly_returns_bool(self, hexapeptide_chain):
        """Non-GLY residues should return True or False, not unknown."""
        residue_list = list(hexapeptide_chain.get_residues())
        for i in range(len(residue_list)):
            if residue_list[i].get_resname() == "GLY":
                continue
            assert label_is_left(residue_list, i, UNK) in {"True", "False"}
            assert label_is_right(residue_list, i, UNK) in {"True", "False"}

    def test_left_right_exclusive(self, hexapeptide_chain):
        """is_left and is_right should be mutually exclusive."""
        residue_list = list(hexapeptide_chain.get_residues())
        for i in range(len(residue_list)):
            left = label_is_left(residue_list, i, UNK)
            right = label_is_right(residue_list, i, UNK)
            if left == UNK:
                assert right == UNK
            else:
                assert left != right, f"Residue {i}: is_left={left} is_right={right}"

    def test_chirality_in_output(self, hexapeptide_structure):
        """is_left should produce valid values in output."""
        lines = process_measurement_commands(
            "test",
            hexapeptide_structure,
            "is_left; is_right",
        )
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) > 0
        valid = {"True", "False", UNK}
        for line in data_lines:
            parts = line.split(":")
            assert parts[-1] in valid
            assert parts[-2] in valid


class TestLabelIntegration:
    """Test labels via the full measurement pipeline."""

    def test_rama_in_output(self, hexapeptide_structure):
        """All rama variants should appear in output lines."""
        lines = process_measurement_commands(
            "test",
            hexapeptide_structure,
            "rama_category; rama6; rama5; rama4; rama3",
        )
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) > 0
        valid = {
            "General", "Gly", "IleVal",
            "TransPro", "CisPro", "PrePro", "Pro",
        }
        for line in data_lines:
            # Check all 5 rama fields (last 5 colon-separated values)
            parts = line.split(":")
            rama_values = parts[-5:]
            for val in rama_values:
                assert val in valid, f"Invalid category: {val!r} in {line}"

    def test_mixed_output(self, hexapeptide_structure):
        """Labels mixed with measurements should all produce output."""
        lines = process_measurement_commands(
            "test",
            hexapeptide_structure,
            "phi; rama_category",
        )
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) > 0
        valid = {
            "General",
            "Gly",
            "IleVal",
            "TransPro",
            "CisPro",
            "PrePro",
        }
        for line in data_lines:
            parts = line.split(":")
            assert parts[-1] in valid

    def test_is_cis_in_output(self, hexapeptide_structure):
        """is_cis should produce True/False or __?__ values."""
        lines = process_measurement_commands(
            "test",
            hexapeptide_structure,
            "is_cis",
        )
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) > 0
        valid = {"True", "False", "__?__"}
        for line in data_lines:
            value = line.rsplit(":", 1)[-1]
            assert value in valid, f"Invalid is_cis value: {value!r}"
