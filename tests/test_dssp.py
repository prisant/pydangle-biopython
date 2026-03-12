"""Tests for DSSP secondary structure assignment."""

import os

import pytest

from pydangle_biopython.dssp import (
    find_mkdssp,
    label_dssp,
    label_dssp3,
    parse_dssp_output,
    run_dssp,
    set_dssp_assignments,
)
from pydangle_biopython.measure import process_measurement_commands

# Skip all tests if mkdssp is not installed
pytestmark = pytest.mark.skipif(
    find_mkdssp() is None,
    reason="mkdssp not installed",
)

UBQ_PATH = os.path.join(
    os.path.dirname(__file__), "..", "examples", "1ubq.pdb",
)


class TestFindMkdssp:
    def test_finds_binary(self):
        assert find_mkdssp() is not None


class TestRunDssp:
    def test_returns_output(self):
        result = run_dssp(UBQ_PATH)
        assert result is not None
        assert "#  RESIDUE AA" in result

    def test_nonexistent_file_returns_none(self):
        result = run_dssp("/nonexistent/file.pdb")
        assert result is None


class TestParseDsspOutput:
    def test_parses_ubiquitin(self):
        text = run_dssp(UBQ_PATH)
        assert text is not None
        assignments = parse_dssp_output(text)
        # Ubiquitin has 76 residues
        assert len(assignments) == 76

    def test_known_helix_residue(self):
        """Residue 23-34 should be alpha helix (H)."""
        text = run_dssp(UBQ_PATH)
        assert text is not None
        assignments = parse_dssp_output(text)
        assert assignments[("A", 26, "")] == "H"

    def test_known_strand_residue(self):
        """Residue 2 should be beta strand (E)."""
        text = run_dssp(UBQ_PATH)
        assert text is not None
        assignments = parse_dssp_output(text)
        assert assignments[("A", 2, "")] == "E"

    def test_coil_residue(self):
        """Residue 17 should be coil (C)."""
        text = run_dssp(UBQ_PATH)
        assert text is not None
        assignments = parse_dssp_output(text)
        assert assignments[("A", 17, "")] == "C"

    def test_all_codes_valid(self):
        text = run_dssp(UBQ_PATH)
        assert text is not None
        assignments = parse_dssp_output(text)
        valid_codes = {"H", "B", "E", "G", "I", "T", "S", "C", "P"}
        for code in assignments.values():
            assert code in valid_codes, f"Invalid DSSP code: {code!r}"


class TestLabelDssp:
    @pytest.fixture(autouse=True)
    def _setup_assignments(self):
        """Set DSSP assignments for ubiquitin before each test."""
        text = run_dssp(UBQ_PATH)
        assert text is not None
        assignments = parse_dssp_output(text)
        set_dssp_assignments(assignments)
        yield
        set_dssp_assignments(None)

    @pytest.fixture()
    def ubq_chain(self):
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("ubq", UBQ_PATH)
        return list(structure[0]["A"].get_residues())

    def test_helix_residue(self, ubq_chain):
        # Residue 26 (index 25) should be H
        result = label_dssp(ubq_chain, 25, "__?__")
        assert result == "H"

    def test_strand_residue(self, ubq_chain):
        # Residue 2 (index 1) should be E
        result = label_dssp(ubq_chain, 1, "__?__")
        assert result == "E"

    def test_dssp3_helix(self, ubq_chain):
        result = label_dssp3(ubq_chain, 25, "__?__")
        assert result == "H"

    def test_dssp3_strand(self, ubq_chain):
        result = label_dssp3(ubq_chain, 1, "__?__")
        assert result == "E"

    def test_dssp3_only_three_values(self, ubq_chain):
        valid = {"H", "E", "C"}
        for i in range(min(len(ubq_chain), 76)):
            result = label_dssp3(ubq_chain, i, "__?__")
            assert result in valid, f"Index {i}: {result!r}"

    def test_no_assignments_returns_unknown(self, ubq_chain):
        set_dssp_assignments(None)
        result = label_dssp(ubq_chain, 0, "__?__")
        assert result == "__?__"


class TestDsspIntegration:
    def test_dssp_in_output(self):
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("ubq", UBQ_PATH)
        lines = process_measurement_commands(
            "1ubq.pdb", structure, "dssp; dssp3",
            filepath=UBQ_PATH,
        )
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) == 76
        valid8 = {"H", "B", "E", "G", "I", "T", "S", "C", "P"}
        valid3 = {"H", "E", "C"}
        for line in data_lines:
            parts = line.split(":")
            assert parts[-2] in valid8, f"Invalid dssp: {line}"
            assert parts[-1] in valid3, f"Invalid dssp3: {line}"

    def test_mixed_with_phi_psi(self):
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("ubq", UBQ_PATH)
        lines = process_measurement_commands(
            "1ubq.pdb", structure, "phi; psi; dssp; dssp3",
            filepath=UBQ_PATH,
        )
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) > 0
        # Each data line should have phi:psi:dssp:dssp3 after the header fields
        for line in data_lines:
            parts = line.split(":")
            assert len(parts) >= 10  # header(6) + phi + psi + dssp + dssp3
