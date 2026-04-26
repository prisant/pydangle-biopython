"""Tests for DSSP secondary structure assignment."""

import os
import tempfile

import pytest

from pydangle_biopython.dssp import (
    _clean_pdb_for_dssp,
    _needs_pdb_cleanup,
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


class TestPdbCleanup:
    """Cleanup of mkdssp-hostile records before invoking mkdssp 4.x.

    Regression coverage for the blank-chain-SEQRES bug discovered while
    propagating DSSP fixes from top100 to top500: legacy single-chain
    PDB depositions can have a space at column 12 of SEQRES records,
    which mkdssp 4.x rejects.  Build_ersatz_only.py only filled column
    22 (ATOM) chains, so the blank SEQRES propagated through to
    pydangle's per-file mkdssp call, silently nulling the entire
    structure's DSSP output.
    """

    def _write_pdb(self, lines):
        """Write *lines* to a temp file, return its path."""
        fd, path = tempfile.mkstemp(suffix=".pdb")
        with os.fdopen(fd, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return path

    def test_blank_seqres_chain_triggers_cleanup(self):
        path = self._write_pdb([
            "HEADER    HYDROLASE",
            "SEQRES   1     93  SER LYS ALA",
            "ATOM      1  N   SER A   1       0.000   0.000   0.000",
        ])
        try:
            assert _needs_pdb_cleanup(path) is True
        finally:
            os.unlink(path)

    def test_filled_seqres_chain_does_not_trigger_cleanup(self):
        path = self._write_pdb([
            "HEADER    HYDROLASE",
            "SEQRES   1 A   3  SER LYS ALA",
            "ATOM      1  N   SER A   1       0.000   0.000   0.000",
        ])
        try:
            assert _needs_pdb_cleanup(path) is False
        finally:
            os.unlink(path)

    def test_cleanup_fills_blank_seqres_chain(self):
        path = self._write_pdb([
            "HEADER    HYDROLASE",
            "SEQRES   1     93  SER LYS ALA",
            "SEQRES   2     93  VAL LYS TYR",
            "ATOM      1  N   SER A   1       0.000   0.000   0.000",
        ])
        try:
            cleaned = _clean_pdb_for_dssp(path)
            try:
                with open(cleaned, encoding="utf-8") as fh:
                    seqres_lines = [
                        line for line in fh
                        if line.startswith("SEQRES")
                    ]
                assert len(seqres_lines) == 2
                for line in seqres_lines:
                    assert line[11] == "A", (
                        f"SEQRES chain ID not filled: {line!r}"
                    )
            finally:
                os.unlink(cleaned)
        finally:
            os.unlink(path)


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

    def test_loop_residue(self):
        """Residue 17 should be loop (None) — DSSP space."""
        text = run_dssp(UBQ_PATH)
        assert text is not None
        assignments = parse_dssp_output(text)
        assert assignments[("A", 17, "")] is None

    def test_all_codes_valid(self):
        text = run_dssp(UBQ_PATH)
        assert text is not None
        assignments = parse_dssp_output(text)
        valid_codes = {"H", "B", "E", "G", "I", "T", "S", "P", None}
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
        # dssp: SS codes or __?__ for loop/unassigned
        valid8 = {"H", "B", "E", "G", "I", "T", "S", "P", "__?__"}
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
