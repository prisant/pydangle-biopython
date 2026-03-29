"""Command-line interface for pydangle-biopython."""

import argparse
import multiprocessing
import os
import sys
import tempfile
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
_COORD_PREFIXES = ("ATOM  ", "HETATM", "MODEL ", "ENDMDL", "TER   ", "TER\n")


def _parse_pdb_resilient(filepath: str) -> Any:
    """Parse a PDB file, falling back to coordinate-only parsing on error.

    BioPython's PDB header parser can crash on some files (e.g. Reduce
    hydrogen-added files with unusual REMARK records).  When that happens,
    retry parsing with only ATOM/HETATM/MODEL/ENDMDL/TER lines.
    """
    parser: Any = PDBParser(QUIET=True)  # type: ignore[no-untyped-call]
    try:
        return parser.get_structure("X", filepath)
    except Exception:
        pass

    # Fallback: strip non-coordinate lines and parse from a temp file
    print(
        f"pydangle-biopython: warning: header parse failed for "
        f"{os.path.basename(filepath)}, retrying without headers",
        file=sys.stderr,
    )
    with open(filepath, encoding="utf-8") as fh:
        coord_lines = [
            line
            for line in fh
            if line.startswith(_COORD_PREFIXES) or line.strip() == "END"
        ]
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".pdb", delete=False
    ) as tf:
        tf.writelines(coord_lines)
        tmp_path = tf.name
    try:
        return parser.get_structure("X", tmp_path)
    finally:
        os.unlink(tmp_path)


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
# Per-file processing (top-level for picklability)
# ---------------------------------------------------------------------------


def _process_one_file(
    args_tuple: tuple[str, str, str, str, tuple[str, ...] | None],
) -> list[str]:
    """Process a single structure file.

    Designed as a top-level function for use with :func:`multiprocessing.Pool`.

    Parameters
    ----------
    args_tuple : tuple
        ``(filepath, file_format, command_string, output_format,
        chain_filter)`` where *file_format* is ``'pdb'``, ``'cif'``,
        or ``''`` (auto-detect), *output_format* is ``'csv'`` or
        ``'jsonl'``, and *chain_filter* is a tuple of chain IDs to
        include or ``None`` for all chains.

    Returns
    -------
    list[str]
        Output lines for this file.
    """
    filepath, file_format, command_string, output_format, chain_filter = (
        args_tuple
    )
    if not file_format:
        file_format = _guess_format(filepath)
    basename = os.path.basename(filepath)

    try:
        if file_format == "cif":
            struct_parser: Any = MMCIFParser(QUIET=True)  # type: ignore[no-untyped-call]
            structure = struct_parser.get_structure("X", filepath)
        else:
            structure = _parse_pdb_resilient(filepath)

        return process_measurement_commands(
            basename,
            structure,
            command_string,
            filepath=filepath,
            output_format=output_format,
            chain_filter=set(chain_filter) if chain_filter else None,
        )
    except Exception as exc:
        print(
            f"pydangle-biopython: warning: failed to process "
            f"{basename}: {exc}",
            file=sys.stderr,
        )
        return []


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

    ap.add_argument(
        "-C",
        "--chain",
        action="append",
        dest="chain_filter",
        metavar="CHAIN",
        help=(
            "Limit output to the specified chain ID(s). May be "
            "repeated (e.g. -C A -C B). DSSP still runs on the "
            "full structure for correct assignments."
        ),
    )

    ap.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Number of parallel worker processes (default: 1). "
            "Use 0 for auto-detection (one per CPU core)."
        ),
    )

    ap.add_argument(
        "-o",
        "--output-format",
        choices=["csv", "jsonl"],
        default="csv",
        dest="output_format",
        help=(
            "Output format: csv (colon-separated, default) or "
            "jsonl (one JSON object per line)."
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

    # Resolve jobs count
    jobs: int = args.jobs
    if jobs < 0:
        print(
            "pydangle-biopython: error: --jobs must be >= 0",
            file=sys.stderr,
        )
        return 1
    if jobs == 0:
        jobs = os.cpu_count() or 1

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
    file_format: str = args.file_format or ""
    command_string: str = args.command_string
    output_format: str = args.output_format
    chain_filter: tuple[str, ...] | None = (
        tuple(args.chain_filter) if args.chain_filter else None
    )

    if jobs == 1:
        # Serial processing (original behaviour, no multiprocessing overhead)
        for idx, filepath in enumerate(all_files):
            if show_progress:
                print(
                    f"\r[{idx + 1}/{total}] {os.path.basename(filepath)}",
                    end="",
                    file=sys.stderr,
                )
            output_lines = _process_one_file(
                (filepath, file_format, command_string, output_format,
                 chain_filter),
            )
            for line in output_lines:
                print(line)
    else:
        # Parallel processing
        work = [
            (fp, file_format, command_string, output_format, chain_filter)
            for fp in all_files
        ]
        with multiprocessing.Pool(jobs) as pool:
            for idx, lines in enumerate(pool.imap(_process_one_file, work)):
                if show_progress:
                    print(
                        f"\r[{idx + 1}/{total}]",
                        end="",
                        file=sys.stderr,
                    )
                for line in lines:
                    print(line)

    if show_progress:
        print("", file=sys.stderr)  # newline after progress
    return 0


if __name__ == "__main__":
    sys.exit(main())
