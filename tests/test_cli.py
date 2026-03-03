"""Tests for the pydangle command-line interface."""

import pytest

from pydangle_biopython.cli import _guess_format, main
from tests.conftest import HEXAPEPTIDE_PDB, RNA_DINUCLEOTIDE_PDB

# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

class TestGuessFormat:
    def test_pdb_extension(self):
        assert _guess_format("structure.pdb") == 'pdb'

    def test_ent_extension(self):
        assert _guess_format("1abc.ent") == 'pdb'

    def test_cif_extension(self):
        assert _guess_format("structure.cif") == 'cif'

    def test_mmcif_extension(self):
        assert _guess_format("structure.mmcif") == 'cif'

    def test_unknown_defaults_to_pdb(self):
        assert _guess_format("structure.xyz") == 'pdb'

    def test_case_insensitive(self):
        assert _guess_format("STRUCTURE.PDB") == 'pdb'
        assert _guess_format("STRUCTURE.CIF") == 'cif'


# ---------------------------------------------------------------------------
# CLI integration tests
# ---------------------------------------------------------------------------

class TestCLI:
    @pytest.fixture
    def pdb_file(self, tmp_path):
        """Write the hexapeptide PDB to a temporary file."""
        f = tmp_path / "test.pdb"
        f.write_text(HEXAPEPTIDE_PDB)
        return str(f)

    @pytest.fixture
    def rna_file(self, tmp_path):
        """Write the RNA PDB to a temporary file."""
        f = tmp_path / "rna.pdb"
        f.write_text(RNA_DINUCLEOTIDE_PDB)
        return str(f)

    def test_default_commands(self, pdb_file, capsys):
        ret = main([pdb_file])
        assert ret == 0
        output = capsys.readouterr().out
        assert 'pydangle' in output  # header line
        assert len(output.strip().split('\n')) > 1

    def test_custom_commands(self, pdb_file, capsys):
        ret = main(['-c', 'tau; phi; psi', pdb_file])
        assert ret == 0
        output = capsys.readouterr().out
        assert len(output.strip().split('\n')) > 1

    def test_force_pdb_format(self, pdb_file, capsys):
        ret = main(['--pdb', pdb_file])
        assert ret == 0

    def test_nonexistent_file(self, capsys):
        ret = main(['nonexistent.pdb'])
        assert ret == 1
        stderr = capsys.readouterr().err
        assert 'cannot find' in stderr

    def test_rna_commands(self, rna_file, capsys):
        ret = main(['-c', 'delta', rna_file])
        assert ret == 0
        output = capsys.readouterr().out
        assert len(output.strip().split('\n')) >= 1

    def test_multiple_files(self, pdb_file, rna_file, capsys):
        ret = main(['-c', 'phi; psi', pdb_file, pdb_file])
        assert ret == 0
        output = capsys.readouterr().out
        lines = [
            ln for ln in output.strip().split('\n')
            if not ln.startswith('#')
        ]
        # Two copies of the same file should double the output
        assert len(lines) > 2

    def test_multiple_files_one_missing(self, pdb_file, capsys):
        ret = main([pdb_file, 'nonexistent.pdb'])
        assert ret == 1
        stderr = capsys.readouterr().err
        assert 'cannot find' in stderr
