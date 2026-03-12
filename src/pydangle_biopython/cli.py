"""Command-line interface for pydangle-biopython."""

import argparse
import os
import sys
from typing import Any

from Bio.PDB import MMCIFParser, PDBParser  # type: ignore[attr-defined]

from pydangle_biopython import __version__
from pydangle_biopython.fileinput import collect_files
from pydangle_biopython.measure import process_measurement_commands

# ---------------------------------------------------------------------------
# File format detection
# ---------------------------------------------------------------------------

_PDB_EXTENSIONS = {".pdb", ".ent", ".pdb1"}
_CIF_EXTENSIONS = {".cif", ".mmcif"}


def _guess_format(filepath: str) -> str:
    """Guess file format from extension.

    Returns
    -------
    str
        ``'pdb'`` or ``'cif'``.
    """
    ext = os.path.splitext(filepath)[1].lower()
    if ext in _CIF_EXTENSIONS:
        return "cif"
    return "pdb"


# ---------------------------------------------------------------------------
# Default commands
# ---------------------------------------------------------------------------

DEFAULT_COMMANDS = "phi; psi; chi1; chi2; chi3; chi4"


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    """Main entry point for the pydangle-biopython command-line tool.

    Parameters
    ----------
    argv : list[str] or None
        Command-line arguments.  If ``None``, ``sys.argv[1:]`` is used.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    ap = argparse.ArgumentParser(
        prog="pydangle-biopython",
        description=(
            "Compute distances, angles, and dihedral angles in protein "
            "and nucleic acid structures (BioPython backend)."
        ),
        epilog=(
            "Examples:\n"
            "  pydangle-biopython 1abc.pdb\n"
            "  pydangle-biopython *.pdb\n"
            "  pydangle-biopython -c 'phi; psi; omega'"
            " structure.pdb\n"
            "  pydangle-biopython -c 'alpha; beta; gamma;"
            " delta; epsilon; zeta' rna.cif\n"
            "  pydangle-biopython -c 'distance: Ca_Ca:"
            " i-1 _CA_, i _CA_' structure.pdb\n"
            "\n"
            "Bulk input:\n"
            "  pydangle-biopython -f file_list.txt\n"
            "  find /data -name '*.pdb' | pydangle-biopython"
            " -f -\n"
            "  pydangle-biopython -g '**/*.pdb'\n"
            "  pydangle-biopython -d /data/structures/\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    ap.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    ap.add_argument(
        "structure_files",
        nargs="*",
        help="One or more PDB or mmCIF structure files.",
    )
    ap.add_argument(
        "-f",
        "--file-list",
        action="append",
        dest="file_lists",
        metavar="FILE",
        help=(
            "Read structure file paths from FILE (one per line). "
            "Use - for stdin."
        ),
    )
    ap.add_argument(
        "-g",
        "--glob",
        action="append",
        dest="glob_patterns",
        metavar="PATTERN",
        help="Expand a Python glob pattern (e.g., '**/*.pdb').",
    )
    ap.add_argument(
        "-d",
        "--directory",
        action="append",
        dest="directories",
        metavar="DIR",
        help="Recursively find .pdb/.cif/.mmcif/.ent files under DIR.",
    )
    ap.add_argument(
        "-c",
        "--commands",
        dest="command_string",
        default=DEFAULT_COMMANDS,
        help=(
            f"Measurement command string (default: '{DEFAULT_COMMANDS}'). "
            "Use semicolons to separate multiple measurements."
        ),
    )

    format_group = ap.add_mutually_exclusive_group()
    format_group.add_argument(
        "-p",
        "--pdb",
        action="store_const",
        const="pdb",
        dest="file_format",
        help="Force PDB format parsing.",
    )
    format_group.add_argument(
        "-m",
        "--mmcif",
        action="store_const",
        const="cif",
        dest="file_format",
        help="Force mmCIF format parsing.",
    )

    args = ap.parse_args(argv)

    # Collect files from all input sources
    try:
        all_files = collect_files(
            args.structure_files or None,
            args.file_lists,
            args.glob_patterns,
            args.directories,
        )
    except FileNotFoundError as exc:
        print(
            f"pydangle-biopython: error: {exc}",
            file=sys.stderr,
        )
        return 1

    if not all_files:
        print(
            "pydangle-biopython: error: no structure files specified",
            file=sys.stderr,
        )
        return 1

    # Validate all input files before processing any
    for filepath in all_files:
        if not os.path.isfile(filepath):
            print(
                f"pydangle-biopython: error: cannot find {filepath}",
                file=sys.stderr,
            )
            return 1

    # Process each file
    total = len(all_files)
    show_progress = total > 10 and sys.stderr.isatty()
    for idx, filepath in enumerate(all_files):
        if show_progress:
            print(
                f"\r[{idx + 1}/{total}] {os.path.basename(filepath)}",
                end="",
                file=sys.stderr,
            )
        file_format: str = args.file_format or _guess_format(filepath)
        basename = os.path.basename(filepath)

        struct_parser: Any
        if file_format == "cif":
            struct_parser = MMCIFParser(  # type: ignore[no-untyped-call]
                QUIET=True,
            )
        else:
            struct_parser = PDBParser(  # type: ignore[no-untyped-call]
                QUIET=True,
            )

        structure = struct_parser.get_structure(
            "X",
            filepath,
        )

        output_lines = process_measurement_commands(
            basename,
            structure,
            args.command_string,
            filepath=filepath,
        )

        for line in output_lines:
            print(line)

    if show_progress:
        print("", file=sys.stderr)  # newline after progress
    return 0


if __name__ == "__main__":
    sys.exit(main())
