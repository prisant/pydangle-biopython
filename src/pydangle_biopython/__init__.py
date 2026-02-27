"""pydangle-biopython: distances, angles, and dihedrals in structures.

A Python translation of the Java Dangle program from the Richardson
Laboratory at Duke University, using BioPython as the structure-parsing
backend.  Measures backbone and sidechain geometry in protein and nucleic
acid structures from PDB or mmCIF files.
"""

__version__ = "0.1.0"

from pydangle_biopython.builtins import (
    BUILTIN_COMMANDS as BUILTIN_COMMANDS,
)
from pydangle_biopython.measure import (
    process_measurement_commands as process_measurement_commands,
)
from pydangle_biopython.parser import (
    command_string_parser as command_string_parser,
)
