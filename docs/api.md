# API Reference

## pydangle_biopython.builtins

### `BUILTIN_COMMANDS`

```python
BUILTIN_COMMANDS: dict[str, str]
```

Dictionary mapping short measurement names to their expanded command strings.
Includes protein backbone, sidechain, and nucleic acid measurements.

**Protein builtins:** `phi`, `psi`, `omega`, `tau`, `chi1`–`chi4`,
`vCAd`, `vCAa`, `vCAt`, `pbCACAd`, `pbCACd`, `pbCOd`, `pbCNd`, `pbNCAd`,
`pbCACOa`, `pbCACNa`, `pbOCNa`, `pbCNCAa`, `scABG`, `scBGD`, `scGDE`, `scDEZ`

**Nucleic acid builtins:** `alpha`, `beta`, `gamma`, `delta`, `epsilon`,
`zeta`, `eta`, `theta`, `chi_na`, `nu0`–`nu4`

---

## pydangle_biopython.parser

### `command_string_parser(command_string)`

```python
def command_string_parser(command_string: str) -> list[tuple]
```

Parse a measurement command string into structured specifications.

**Parameters:**

- `command_string` — Semicolon-delimited commands, e.g. `"phi; psi; chi1"`

**Returns:** List of `(function_key, label, arg_lists)` tuples where:

- `function_key` is `'distance'`, `'angle'`, or `'dihedral'`
- `label` is a human-readable name
- `arg_lists` is a list of alternative atom-position lists

**Example:**

```python
from pydangle_biopython.parser import command_string_parser

commands = command_string_parser("phi; psi")
for fun_key, label, arg_lists in commands:
    print(f"{label}: {fun_key} with {len(arg_lists[0])} atoms")
# phi: dihedral with 4 atoms
# psi: dihedral with 4 atoms
```

---

## pydangle_biopython.measure

### `process_measurement_commands(label, structure, commands, unknown_str="__?__")`

```python
def process_measurement_commands(
    label: str,
    structure: Bio.PDB.Structure.Structure,
    commands: str,
    unknown_str: str = "__?__",
) -> list[str]
```

Process all measurement commands on an entire structure.

**Parameters:**

- `label` — Identifier for the structure (typically the filename)
- `structure` — Parsed BioPython structure object
- `commands` — Measurement command string
- `unknown_str` — String used when a measurement cannot be computed

**Returns:** List of formatted output lines.

**Example:**

```python
from Bio.PDB import PDBParser
from pydangle_biopython.measure import process_measurement_commands

parser = PDBParser(QUIET=True)
structure = parser.get_structure('X', '1abc.pdb')
lines = process_measurement_commands('1abc.pdb', structure, 'tau; phi; psi')
```

### `compute_measurement(command, chain, residue_index, unknown_str)`

```python
def compute_measurement(command, chain, residue_index, unknown_str) -> str
```

Compute a single measurement for one residue.

**Parameters:**

- `command` — Parsed command tuple from `command_string_parser`
- `chain` — BioPython Chain object
- `residue_index` — Index into `chain.child_list`
- `unknown_str` — Fallback string

**Returns:** Formatted measurement or `unknown_str`.

### `process_measurement_for_residue(label, chain, residue_index, command_list, unknown_str="__?__")`

```python
def process_measurement_for_residue(label, chain, residue_index, command_list, unknown_str) -> str | None
```

Compute all measurements for a single residue.  Returns a formatted output
line, or `None` if no measurements were computable.

---

## pydangle_biopython.cli

### `main(argv=None)`

```python
def main(argv: list[str] | None = None) -> int
```

CLI entry point.  Parses command-line arguments and runs measurements.

---

## pydangle_biopython.exceptions

### `PydangleError`

Base exception for all pydangle-biopython errors.

### `ParseError`

Raised when a measurement command string cannot be parsed.

### `MeasurementError`

Raised when a measurement cannot be computed.

### `FileFormatError`

Raised when a structure file format is unrecognised.

## Measurement command syntax

### Built-in shortcuts

Simply name the measurement: `phi`, `psi`, `omega`, `alpha`, `delta`, etc.

### Custom measurements

```
distance: <label>: [i±n] <atom>, [i±n] <atom>
angle:    <label>: [i±n] <atom>, [i±n] <atom>, [i±n] <atom>
dihedral: <label>: [i±n] <atom>, [i±n] <atom>, [i±n] <atom>, [i±n] <atom>
```

- Use underscores (`_`) for spaces in 4-character atom names
- Use `/regexp/` for regex matching
- Use `*` in nucleic acid names to match both `*` and `'` conventions
- Use `|` to separate alternative atom sets

### Output format

```
filename:model:chain:resnum:icode:resname:measurement1:measurement2:...
```
