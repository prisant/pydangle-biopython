"""File input collection utilities for pydangle-biopython."""

import glob as _glob
import os
import sys

# Extensions recognized for directory recursion
STRUCTURE_EXTENSIONS = {".pdb", ".ent", ".pdb1", ".cif", ".mmcif"}


def files_from_list(path: str) -> list[str]:
    """Read file paths from a text file, one per line.

    If *path* is ``"-"``, read from stdin.  Blank lines and lines
    starting with ``#`` are skipped.
    """
    if path == "-":
        if sys.stdin.isatty():
            print(
                "pydangle-biopython: warning: reading file list from "
                "terminal stdin (use Ctrl-D to end)",
                file=sys.stderr,
            )
        lines = sys.stdin.read().splitlines()
    else:
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"file list not found: {path}",
            )
        with open(path, encoding="utf-8") as fh:
            lines = fh.read().splitlines()

    result: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            result.append(stripped)
    return result


def files_from_glob(pattern: str) -> list[str]:
    """Expand a Python glob pattern (supports ``**/`` recursive)."""
    return sorted(_glob.glob(pattern, recursive=True))


def files_from_directory(directory: str) -> list[str]:
    """Recursively find structure files under *directory*."""
    if not os.path.isdir(directory):
        raise FileNotFoundError(
            f"directory not found: {directory}",
        )
    result: list[str] = []
    for dirpath, _dirnames, filenames in os.walk(directory):
        for fname in sorted(filenames):
            ext = os.path.splitext(fname)[1].lower()
            if ext in STRUCTURE_EXTENSIONS:
                result.append(os.path.join(dirpath, fname))
    return result


def collect_files(
    positional: list[str] | None,
    file_lists: list[str] | None,
    glob_patterns: list[str] | None,
    directories: list[str] | None,
) -> list[str]:
    """Combine all input sources into a deduplicated, ordered file list.

    Deduplication is based on ``os.path.realpath`` so that symlinks and
    relative/absolute variants of the same path count as one entry.
    Insertion order is preserved.
    """
    seen: dict[str, str] = {}  # realpath -> original path

    def _add(paths: list[str]) -> None:
        for p in paths:
            rp = os.path.realpath(p)
            if rp not in seen:
                seen[rp] = p

    if positional:
        _add(positional)
    if file_lists:
        for fl in file_lists:
            _add(files_from_list(fl))
    if glob_patterns:
        for pat in glob_patterns:
            _add(files_from_glob(pat))
    if directories:
        for d in directories:
            _add(files_from_directory(d))

    return list(seen.values())
