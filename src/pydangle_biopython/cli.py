"""Command-line interface for pydangle-biopython."""

import argparse
import os
import sys
from typing import Any

from Bio.PDB import MMCIFParser, PDBParser  # type: ignore[attr-defined]

from pydangle_biopython.measure import process_measurement_commands

# ---------------------------------------------------------------------------
# File format detection
# ---------------------------------------------------------------------------

_PDB_EXTENSIONS = {'.pdb', '.ent', '.pdb1'}
_CIF_EXTENSIONS = {'.cif', '.mmcif'}


def _guess_format(filepath: str) -> str:
    """Guess file format from extension.

    Returns
    -------
    str
        ``'pdb'`` or ``'cif'``.
    """
    ext = os.path.splitext(filepath)[1].lower()
    if ext in _CIF_EXTENSIONS:
        return 'cif'
    return 'pdb'


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
        prog='pydangle-biopython',
        description=(
            "Compute distances, angles, and dihedral angles in protein "
            "and nucleic acid structures (BioPython backend)."
        ),
        epilog=(
            "Examples:\n"
            "  pydangle-biopython 1abc.pdb\n"
            "  pydangle-biopython -c 'phi; psi; omega'"
            " structure.pdb\n"
            "  pydangle-biopython -c 'alpha; beta; gamma;"
            " delta; epsilon; zeta' rna.cif\n"
            "  pydangle-biopython -c 'distance: Ca_Ca:"
            " i-1 _CA_, i _CA_' structure.pdb\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    ap.add_argument(
        'structure_file',
        help="PDB or mmCIF structure file to process.",
    )
    ap.add_argument(
        '-c', '--commands',
        dest='command_string',
        default=DEFAULT_COMMANDS,
        help=(
            f"Measurement command string (default: '{DEFAULT_COMMANDS}'). "
            "Use semicolons to separate multiple measurements."
        ),
    )

    format_group = ap.add_mutually_exclusive_group()
    format_group.add_argument(
        '-p', '--pdb',
        action='store_const',
        const='pdb',
        dest='file_format',
        help="Force PDB format parsing.",
    )
    format_group.add_argument(
        '-m', '--mmcif',
        action='store_const',
        const='cif',
        dest='file_format',
        help="Force mmCIF format parsing.",
    )

    args = ap.parse_args(argv)

    # Validate input file
    if not os.path.isfile(args.structure_file):
        print(
            f"pydangle-biopython: error: cannot find {args.structure_file}",
            file=sys.stderr,
        )
        return 1

    # Determine format
    file_format = args.file_format or _guess_format(args.structure_file)

    # Parse structure
    basename = os.path.basename(args.structure_file)
    struct_parser: Any
    if file_format == 'cif':
        struct_parser = MMCIFParser(QUIET=True)  # type: ignore[no-untyped-call]
    else:
        struct_parser = PDBParser(QUIET=True)  # type: ignore[no-untyped-call]

    structure = struct_parser.get_structure(
        'X', args.structure_file,
    )

    # Run measurements
    output_lines = process_measurement_commands(
        basename, structure, args.command_string
    )

    for line in output_lines:
        print(line)

    return 0


if __name__ == '__main__':
    sys.exit(main())
