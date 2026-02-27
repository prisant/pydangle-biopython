"""Custom exception classes for pydangle-biopython."""


class PydangleError(Exception):
    """Base exception for all pydangle-biopython errors."""


class ParseError(PydangleError):
    """Raised when a measurement command string cannot be parsed."""


class MeasurementError(PydangleError):
    """Raised when a measurement cannot be computed due to missing or
    degenerate geometry."""


class FileFormatError(PydangleError):
    """Raised when a structure file cannot be read or its format is
    unrecognised."""
