"""Tests for file input collection utilities."""

import os
from io import StringIO

import pytest

from pydangle_biopython.fileinput import (
    collect_files,
    files_from_directory,
    files_from_glob,
    files_from_list,
)

# ---------------------------------------------------------------------------
# files_from_list
# ---------------------------------------------------------------------------


class TestFilesFromList:
    def test_basic(self, tmp_path):
        listing = tmp_path / "files.txt"
        listing.write_text("a.pdb\nb.pdb\nc.pdb\n")
        assert files_from_list(str(listing)) == ["a.pdb", "b.pdb", "c.pdb"]

    def test_comments_and_blanks(self, tmp_path):
        listing = tmp_path / "files.txt"
        listing.write_text("# header\na.pdb\n\n  \n# comment\nb.pdb\n")
        assert files_from_list(str(listing)) == ["a.pdb", "b.pdb"]

    def test_strips_whitespace(self, tmp_path):
        listing = tmp_path / "files.txt"
        listing.write_text("  a.pdb  \n  b.pdb\n")
        assert files_from_list(str(listing)) == ["a.pdb", "b.pdb"]

    def test_stdin(self, monkeypatch):
        monkeypatch.setattr("sys.stdin", StringIO("x.pdb\ny.pdb\n"))
        assert files_from_list("-") == ["x.pdb", "y.pdb"]

    def test_missing_file(self):
        with pytest.raises(FileNotFoundError, match="file list not found"):
            files_from_list("/nonexistent/list.txt")


# ---------------------------------------------------------------------------
# files_from_glob
# ---------------------------------------------------------------------------


class TestFilesFromGlob:
    def test_recursive(self, tmp_path):
        sub = tmp_path / "sub"
        sub.mkdir()
        (tmp_path / "a.pdb").touch()
        (sub / "b.pdb").touch()
        result = files_from_glob(str(tmp_path / "**" / "*.pdb"))
        assert len(result) == 2
        basenames = [os.path.basename(p) for p in result]
        assert "a.pdb" in basenames
        assert "b.pdb" in basenames

    def test_no_match(self, tmp_path):
        result = files_from_glob(str(tmp_path / "*.xyz"))
        assert result == []


# ---------------------------------------------------------------------------
# files_from_directory
# ---------------------------------------------------------------------------


class TestFilesFromDirectory:
    def test_finds_structure_files(self, tmp_path):
        (tmp_path / "a.pdb").touch()
        (tmp_path / "b.cif").touch()
        (tmp_path / "c.txt").touch()  # not a structure file
        sub = tmp_path / "sub"
        sub.mkdir()
        (sub / "d.mmcif").touch()
        result = files_from_directory(str(tmp_path))
        basenames = [os.path.basename(p) for p in result]
        assert "a.pdb" in basenames
        assert "b.cif" in basenames
        assert "d.mmcif" in basenames
        assert "c.txt" not in basenames

    def test_missing_directory(self):
        with pytest.raises(FileNotFoundError, match="directory not found"):
            files_from_directory("/nonexistent/dir")


# ---------------------------------------------------------------------------
# collect_files
# ---------------------------------------------------------------------------


class TestCollectFiles:
    def test_positional_only(self):
        result = collect_files(["a.pdb", "b.pdb"], None, None, None)
        assert result == ["a.pdb", "b.pdb"]

    def test_all_none(self):
        result = collect_files(None, None, None, None)
        assert result == []

    def test_deduplication(self, tmp_path):
        f = tmp_path / "a.pdb"
        f.touch()
        abs_path = str(f)
        result = collect_files([abs_path, abs_path], None, None, None)
        assert result == [abs_path]

    def test_preserves_order(self):
        result = collect_files(["c.pdb", "a.pdb", "b.pdb"], None, None, None)
        assert result == ["c.pdb", "a.pdb", "b.pdb"]

    def test_combined_sources(self, tmp_path):
        # Positional
        f1 = tmp_path / "pos.pdb"
        f1.touch()
        # File list
        f2 = tmp_path / "listed.pdb"
        f2.touch()
        listing = tmp_path / "files.txt"
        listing.write_text(str(f2) + "\n")
        # Glob
        f3 = tmp_path / "globbed.pdb"
        f3.touch()
        # Directory
        sub = tmp_path / "subdir"
        sub.mkdir()
        f4 = sub / "found.pdb"
        f4.touch()

        result = collect_files(
            [str(f1)],
            [str(listing)],
            [str(tmp_path / "globbed.*")],
            [str(sub)],
        )
        basenames = [os.path.basename(p) for p in result]
        assert "pos.pdb" in basenames
        assert "listed.pdb" in basenames
        assert "globbed.pdb" in basenames
        assert "found.pdb" in basenames
        assert len(result) == 4
