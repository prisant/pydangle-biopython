"""Command string parsing for pydangle measurement specifications.

Translates user-supplied command strings into structured measurement
specifications that can be evaluated against macromolecular structures.

Command string grammar:
    command_list  ::= command { ';' command }
    command       ::= builtin_name | function_type ':' label ':' arg_lists
    arg_lists     ::= arg_list { '|' arg_list }
    arg_list      ::= atom_pos { ',' atom_pos }
    atom_pos      ::= [residue_offset] atom_name
    residue_offset ::= 'i' ['+' | '-'] digit
    atom_name     ::= 4_char_name | '/' regexp '/'
"""

import re

from pydangle_biopython.builtins import BUILTIN_COMMANDS

# ---------------------------------------------------------------------------
# Function key synonyms → canonical form
# ---------------------------------------------------------------------------
FUNCTION_KEY_MAP: dict[str, str] = {
    'dis':      'distance',
    'dist':     'distance',
    'distance': 'distance',
    'ang':      'angle',
    'angl':     'angle',
    'angle':    'angle',
    'tor':      'dihedral',
    'tors':     'dihedral',
    'torsion':  'dihedral',
    'dihedral': 'dihedral',
    'label':    'label',
    'lbl':      'label',
}

# Expected argument counts per function type
_EXPECTED_ARGS: dict[str, int] = {
    'distance': 2,
    'angle':    3,
    'dihedral': 4,
    # 'label' has no atom arguments — handled separately
}


# Type aliases for parsed command specifications
AtomPos = tuple[int, re.Pattern[str]]
ArgList = list[AtomPos]
ParsedCommand = tuple[str, str, list[ArgList]]


# ---------------------------------------------------------------------------
# Low-level tokenisation helpers
# ---------------------------------------------------------------------------

def _tokenize(s: str, delimiter: str) -> list[str]:
    """Split string *s* on *delimiter* and strip whitespace from each token."""
    return [tok.strip() for tok in s.split(delimiter)]


# ---------------------------------------------------------------------------
# Parsing stages
# ---------------------------------------------------------------------------

def expand_command_list(command_string: str) -> list[str]:
    """Stage 1 – split on semicolons and expand builtin names.

    Parameters
    ----------
    command_string : str
        Semicolon-delimited list of commands (e.g. ``"phi; psi; chi1"``).

    Returns
    -------
    list[str]
        List of expanded command strings.
    """
    commands = _tokenize(command_string, ';')
    expanded = []
    for cmd in commands:
        cmd = cmd.strip()
        if not cmd:
            continue
        if cmd in BUILTIN_COMMANDS:
            expanded.append(BUILTIN_COMMANDS[cmd])
        else:
            expanded.append(cmd)
    return expanded


def parse_command_fields(command_string: str) -> tuple[str, str, str]:
    """Stage 2 – split a single command into (function_key, label, arg_string).

    Parameters
    ----------
    command_string : str
        A single expanded command such as
        ``"dihedral: phi: i-1 _C__, i _N__, i _CA_, i _C__"``.

    Returns
    -------
    tuple[str, str, str]
        Canonical function key, label, and raw argument string.

    Raises
    ------
    ValueError
        If the command does not have exactly three colon-separated fields
        or if the function key is unrecognised.
    """
    fields = _tokenize(command_string, ':')
    if len(fields) == 2:
        # Label commands: "label: name" (no atom args)
        raw_key, label = fields
        raw_key_lower = raw_key.strip().lower()
        if raw_key_lower not in FUNCTION_KEY_MAP:
            raise ValueError(f"Unknown function key: {raw_key!r}")
        canonical = FUNCTION_KEY_MAP[raw_key_lower]
        if canonical != 'label':
            raise ValueError(
                f"Two-field command requires 'label' function type, "
                f"got {raw_key!r}: {command_string!r}"
            )
        return 'label', label.strip(), ''
    if len(fields) != 3:
        raise ValueError(
            f"Expected 2 or 3 colon-separated fields, got {len(fields)}: "
            f"{command_string!r}"
        )
    raw_key, label, arg_string = fields
    raw_key_lower = raw_key.strip().lower()
    if raw_key_lower not in FUNCTION_KEY_MAP:
        raise ValueError(f"Unknown function key: {raw_key!r}")
    return FUNCTION_KEY_MAP[raw_key_lower], label.strip(), arg_string.strip()


def parse_alternative_arg_lists(arg_string: str) -> list[str]:
    """Stage 3 – split on ``|`` to get alternative argument lists.

    The ``|`` operator allows specifying alternative atom sets (e.g. for
    the nucleic acid glycosidic chi which differs for purines vs pyrimidines).

    Parameters
    ----------
    arg_string : str

    Returns
    -------
    list[str]
    """
    return [s for s in _tokenize(arg_string, '|') if s]


def parse_atom_list(function_key: str, arg_list_string: str) -> list[str]:
    """Stage 4 – split on commas to get individual atom position strings.

    Parameters
    ----------
    function_key : str
        Canonical function key (``'distance'``, ``'angle'``, ``'dihedral'``).
    arg_list_string : str

    Returns
    -------
    list[str]

    Raises
    ------
    ValueError
        If the number of atoms does not match the function type.
    """
    atoms = [s for s in _tokenize(arg_list_string, ',') if s]
    expected = _EXPECTED_ARGS.get(function_key)
    if expected is not None and len(atoms) != expected:
        raise ValueError(
            f"{function_key} requires {expected} atoms, got {len(atoms)}: "
            f"{arg_list_string!r}"
        )
    return atoms


def parse_atom_position(atom_pos_string: str) -> AtomPos:
    """Stage 5 – parse a single atom position into (offset, name_regex).

    An atom position string has the form ``[residue_offset] atom_name``
    where *residue_offset* defaults to ``i`` (i.e. offset 0) and
    *atom_name* is either a 4-character padded name (underscores → spaces)
    or a ``/regexp/``.

    Parameters
    ----------
    atom_pos_string : str

    Returns
    -------
    tuple[int, re.Pattern]
        Integer residue offset and compiled regex for atom name matching.
    """
    parts = atom_pos_string.split()
    if len(parts) == 1:
        # No explicit offset → default to i (offset 0)
        offset_str = 'i+0'
        name_str = parts[0]
    elif len(parts) == 2:
        offset_str = parts[0]
        name_str = parts[1]
    else:
        raise ValueError(
            f"Expected 1 or 2 space-separated tokens in atom position, "
            f"got {len(parts)}: {atom_pos_string!r}"
        )

    # Parse offset: "i" → 0, "i+1" → 1, "i-2" → -2
    offset_str = offset_str.strip()
    if offset_str == 'i':
        offset = 0
    else:
        # Remove the leading 'i'
        offset_tail = offset_str[1:]  # e.g. "+1" or "-2"
        try:
            offset = int(offset_tail)
        except ValueError:
            raise ValueError(
                f"Cannot parse residue offset: {offset_str!r}"
            ) from None

    # Parse atom name: strip slashes for regex, replace _ with space
    name_str = name_str.strip('/')
    # Replace * with [*'] to match both conventions in nucleic acids
    name_str = name_str.replace('*', "[*']")
    name_str = name_str.replace('_', ' ')
    regex = re.compile('^' + name_str + '$')

    return offset, regex


# ---------------------------------------------------------------------------
# Top-level parser
# ---------------------------------------------------------------------------

def command_string_parser(command_string: str) -> list[ParsedCommand]:
    """Parse a full command string into structured measurement specifications.

    Parameters
    ----------
    command_string : str
        User-supplied command string, e.g.
        ``"phi; psi; chi1; chi2; chi3; chi4"`` or
        ``"distance Ca--Ca i-1 _CA_, i _CA_"``.

    Returns
    -------
    list[tuple]
        Each element is ``(function_key, label, arg_lists)`` where:
        - *function_key* is ``'distance'``, ``'angle'``, or ``'dihedral'``
        - *label* is a user-readable name
        - *arg_lists* is a ``list[list[tuple[int, re.Pattern]]]``
          (outer list = alternatives via ``|``; inner list = atom positions).
    """
    parsed_commands: list[ParsedCommand] = []
    for cmd_str in expand_command_list(command_string):
        fun_key, label, arg_str = parse_command_fields(cmd_str)
        if fun_key == 'label':
            # Label commands have no atom arguments
            parsed_commands.append((fun_key, label, []))
            continue
        alternatives = []
        for alt_str in parse_alternative_arg_lists(arg_str):
            atom_positions = []
            for atom_str in parse_atom_list(fun_key, alt_str):
                atom_positions.append(parse_atom_position(atom_str))
            alternatives.append(atom_positions)
        parsed_commands.append((fun_key, label, alternatives))
    return parsed_commands
