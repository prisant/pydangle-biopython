"""Tests for the --chain / -C chain filter feature."""

from pydangle_biopython.cli import main


def test_chain_filter_matching(capsys: object) -> None:
    """Chain A filter on single-chain-A file gives same output as no filter."""

    # Run without filter
    main(["-c", "phi", "-o", "jsonl", "examples/1ubq.pdb"])
    captured = getattr(capsys, "readouterr")()
    no_filter = captured.out.strip().splitlines()

    # Run with -C A (1ubq is chain A)
    main(["-c", "phi", "-o", "jsonl", "-C", "A", "examples/1ubq.pdb"])
    captured = getattr(capsys, "readouterr")()
    with_filter = captured.out.strip().splitlines()

    assert len(no_filter) > 0
    assert no_filter == with_filter


def test_chain_filter_nonexistent(capsys: object) -> None:
    """Chain B filter on single-chain-A file gives empty output."""
    main(["-c", "phi", "-o", "jsonl", "-C", "B", "examples/1ubq.pdb"])
    captured = getattr(capsys, "readouterr")()
    lines = [
        line for line in captured.out.strip().splitlines()
        if line and not line.startswith("#")
    ]
    assert len(lines) == 0


def test_chain_filter_multiple(capsys: object) -> None:
    """Multiple -C flags work together."""
    main(["-c", "phi", "-o", "jsonl", "-C", "A", "-C", "B",
          "examples/1ubq.pdb"])
    captured = getattr(capsys, "readouterr")()
    # Should get chain A output (B doesn't exist but A does)
    lines = [
        line for line in captured.out.strip().splitlines()
        if line and not line.startswith("#")
    ]
    assert len(lines) > 0
