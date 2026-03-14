"""Tests for the pydangle command-line interface."""

import json

import pytest

from pydangle_biopython import __version__
from pydangle_biopython.cli import _guess_format, main
from tests.conftest import HEXAPEPTIDE_PDB, RNA_DINUCLEOTIDE_PDB

# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------


class TestGuessFormat:
    def test_pdb_extension(self):
        assert _guess_format("structure.pdb") == "pdb"

    def test_ent_extension(self):
        assert _guess_format("1abc.ent") == "pdb"

    def test_cif_extension(self):
        assert _guess_format("structure.cif") == "cif"

    def test_mmcif_extension(self):
        assert _guess_format("structure.mmcif") == "cif"

    def test_unknown_defaults_to_pdb(self):
        assert _guess_format("structure.xyz") == "pdb"

    def test_case_insensitive(self):
        assert _guess_format("STRUCTURE.PDB") == "pdb"
        assert _guess_format("STRUCTURE.CIF") == "cif"


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
        assert "pydangle" in output  # header line
        assert len(output.strip().split("\n")) > 1

    def test_custom_commands(self, pdb_file, capsys):
        ret = main(["-c", "tau; phi; psi", pdb_file])
        assert ret == 0
        output = capsys.readouterr().out
        assert len(output.strip().split("\n")) > 1

    def test_force_pdb_format(self, pdb_file, capsys):
        ret = main(["--pdb", pdb_file])
        assert ret == 0

    def test_nonexistent_file(self, capsys):
        ret = main(["nonexistent.pdb"])
        assert ret == 1
        stderr = capsys.readouterr().err
        assert "cannot find" in stderr

    def test_rna_commands(self, rna_file, capsys):
        ret = main(["-c", "delta", rna_file])
        assert ret == 0
        output = capsys.readouterr().out
        assert len(output.strip().split("\n")) >= 1

    def test_multiple_files(self, pdb_file, rna_file, capsys):
        ret = main(["-c", "phi; psi", pdb_file, pdb_file])
        assert ret == 0
        output = capsys.readouterr().out
        lines = [ln for ln in output.strip().split("\n") if not ln.startswith("#")]
        # Two copies of the same file should double the output
        assert len(lines) > 2

    def test_multiple_files_one_missing(self, pdb_file, capsys):
        ret = main([pdb_file, "nonexistent.pdb"])
        assert ret == 1
        stderr = capsys.readouterr().err
        assert "cannot find" in stderr

    def test_no_files_error(self, capsys):
        ret = main(["-c", "phi"])
        assert ret == 1
        stderr = capsys.readouterr().err
        assert "no structure files" in stderr

    def test_file_list_flag(self, pdb_file, tmp_path, capsys):
        listing = tmp_path / "files.txt"
        listing.write_text(pdb_file + "\n")
        ret = main(["-c", "phi", "-f", str(listing)])
        assert ret == 0
        output = capsys.readouterr().out
        assert "pydangle" in output

    def test_file_list_stdin(self, pdb_file, monkeypatch, capsys):
        from io import StringIO

        monkeypatch.setattr("sys.stdin", StringIO(pdb_file + "\n"))
        ret = main(["-c", "phi", "-f", "-"])
        assert ret == 0
        output = capsys.readouterr().out
        assert "pydangle" in output

    def test_glob_flag(self, pdb_file, capsys):
        import os

        pattern = os.path.join(os.path.dirname(pdb_file), "*.pdb")
        ret = main(["-c", "phi", "-g", pattern])
        assert ret == 0
        output = capsys.readouterr().out
        assert "pydangle" in output

    def test_directory_flag(self, pdb_file, capsys):
        import os

        ret = main(["-c", "phi", "-d", os.path.dirname(pdb_file)])
        assert ret == 0
        output = capsys.readouterr().out
        assert "pydangle" in output

    def test_combined_sources(self, pdb_file, tmp_path, capsys):
        # Second file via file list
        f2 = tmp_path / "second.pdb"
        from tests.conftest import HEXAPEPTIDE_PDB

        f2.write_text(HEXAPEPTIDE_PDB)
        listing = tmp_path / "files.txt"
        listing.write_text(str(f2) + "\n")
        ret = main(["-c", "phi", "-f", str(listing), pdb_file])
        assert ret == 0
        output = capsys.readouterr().out
        lines = [ln for ln in output.strip().split("\n") if not ln.startswith("#")]
        assert len(lines) > 2  # output from both files

    def test_version(self, capsys):
        with pytest.raises(SystemExit, match="0"):
            main(["--version"])
        assert __version__ in capsys.readouterr().out

    # -----------------------------------------------------------------------
    # Multiprocessing (-j / --jobs)
    # -----------------------------------------------------------------------

    def test_jobs_flag_serial(self, pdb_file, capsys):
        """``-j 1`` produces the same output as the default."""
        main(["-c", "phi; psi", pdb_file])
        default_out = capsys.readouterr().out

        main(["-c", "phi; psi", "-j", "1", pdb_file])
        serial_out = capsys.readouterr().out

        assert default_out == serial_out

    def test_jobs_flag_parallel(self, pdb_file, capsys):
        """``-j 2`` on multiple files produces correct output."""
        main(["-c", "phi; psi", pdb_file, pdb_file])
        serial_out = capsys.readouterr().out

        main(["-c", "phi; psi", "-j", "2", pdb_file, pdb_file])
        parallel_out = capsys.readouterr().out

        assert serial_out == parallel_out

    # -----------------------------------------------------------------------
    # JSONL output format (-o jsonl)
    # -----------------------------------------------------------------------

    def test_jsonl_output(self, pdb_file, capsys):
        """``-o jsonl`` produces valid JSONL (each line parses as JSON)."""
        ret = main(["-o", "jsonl", "-c", "phi; psi", pdb_file])
        assert ret == 0
        output = capsys.readouterr().out.strip()
        lines = output.split("\n")
        assert len(lines) >= 1
        for line in lines:
            json.loads(line)  # raises on invalid JSON

    def test_jsonl_has_expected_keys(self, pdb_file, capsys):
        """JSONL records contain file, model, chain, resnum, ins, resname."""
        main(["-o", "jsonl", "-c", "phi; psi", pdb_file])
        output = capsys.readouterr().out.strip()
        record = json.loads(output.split("\n")[0])
        for key in ("file", "model", "chain", "resnum", "ins", "resname"):
            assert key in record

    def test_jsonl_numeric_values(self, pdb_file, capsys):
        """Geometric measurements are JSON numbers, not strings."""
        main(["-o", "jsonl", "-c", "phi; psi", pdb_file])
        output = capsys.readouterr().out.strip()
        for line in output.split("\n"):
            record = json.loads(line)
            if record["phi"] is not None:
                assert isinstance(record["phi"], float)
            if record["psi"] is not None:
                assert isinstance(record["psi"], float)

    def test_jsonl_null_for_unknown(self, pdb_file, capsys):
        """Unmeasurable values are JSON null."""
        main(["-o", "jsonl", "-c", "phi; psi", pdb_file])
        output = capsys.readouterr().out.strip()
        records = [json.loads(ln) for ln in output.split("\n")]
        has_null = any(r["phi"] is None or r["psi"] is None for r in records)
        assert has_null

    def test_jsonl_no_comment_lines(self, pdb_file, capsys):
        """JSONL output should not contain comment lines."""
        main(["-o", "jsonl", "-c", "phi; psi", pdb_file])
        output = capsys.readouterr().out.strip()
        for line in output.split("\n"):
            assert not line.startswith("#")

    def test_default_format_unchanged(self, pdb_file, capsys):
        """No ``-o`` flag gives CSV output (same as ``-o csv``)."""
        main(["-c", "phi; psi", pdb_file])
        default_out = capsys.readouterr().out

        main(["-o", "csv", "-c", "phi; psi", pdb_file])
        csv_out = capsys.readouterr().out

        assert default_out == csv_out

    def test_jobs_zero_auto(self, pdb_file, capsys):
        """``-j 0`` auto-detects cores and doesn't crash."""
        ret = main(["-c", "phi; psi", "-j", "0", pdb_file])
        assert ret == 0
        output = capsys.readouterr().out
        assert "pydangle" in output


    def test_jobs_negative_error(self, capsys):
        """Negative --jobs value returns error."""
        ret = main(["-j", "-1", "dummy.pdb"])
        assert ret == 1
        stderr = capsys.readouterr().err
        assert "--jobs" in stderr

    # -----------------------------------------------------------------------
    # Parallel regression: real structure, broad command set
    # -----------------------------------------------------------------------

    def test_parallel_regression_1ubq(self, capsys):
        """Parallel output on 1ubq.pdb must be identical to serial output.

        Uses a broad command set (geometric + labels + DSSP) across
        multiple copies to exercise ordering and data integrity.
        """
        import os

        ubq = os.path.join(
            os.path.dirname(__file__), os.pardir, "examples", "1ubq.pdb"
        )
        if not os.path.isfile(ubq):
            pytest.skip("examples/1ubq.pdb not found")

        cmds = "phi; psi; omega; tau; chi1; rama_category; is_cis; dssp"
        files = [ubq, ubq, ubq]

        main(["-c", cmds, "-j", "1"] + files)
        serial_out = capsys.readouterr().out

        main(["-c", cmds, "-j", "2"] + files)
        parallel_out = capsys.readouterr().out

        assert serial_out == parallel_out
