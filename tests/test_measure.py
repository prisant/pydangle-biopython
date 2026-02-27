"""Tests for the pydangle measurement computation module."""

import math

import pytest
from Bio.PDB.vectors import Vector

from pydangle_biopython.measure import (
    _add_jitter,
    _is_origin,
    angle_to_string,
    calc_dist,
    calc_wrapper,
    compute_measurement,
    number_to_string,
    process_measurement_commands,
    process_measurement_for_residue,
)
from pydangle_biopython.parser import command_string_parser


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

class TestAngleToString:
    def test_zero(self):
        assert angle_to_string(0.0) == '0'

    def test_ninety_degrees(self):
        result = angle_to_string(math.radians(90.0))
        assert abs(float(result) - 90.0) < 0.01

    def test_negative(self):
        result = angle_to_string(math.radians(-120.5))
        assert abs(float(result) - (-120.5)) < 0.01

    def test_trailing_zeros_stripped(self):
        # 180 degrees = pi radians
        result = angle_to_string(math.pi)
        assert '000' not in result


class TestNumberToString:
    def test_integer_value(self):
        result = number_to_string(3.0)
        assert result == '3'

    def test_decimal_value(self):
        result = number_to_string(1.523)
        assert result == '1.523'

    def test_trailing_zeros_stripped(self):
        result = number_to_string(2.100)
        assert result == '2.1'


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

class TestCalcDist:
    def test_zero_distance(self):
        v = Vector(1.0, 2.0, 3.0)
        assert calc_dist(v, v) == pytest.approx(0.0, abs=1e-10)

    def test_unit_distance(self):
        v1 = Vector(0.0, 0.0, 0.0)
        v2 = Vector(1.0, 0.0, 0.0)
        assert calc_dist(v1, v2) == pytest.approx(1.0, abs=1e-10)

    def test_diagonal_distance(self):
        v1 = Vector(0.0, 0.0, 0.0)
        v2 = Vector(1.0, 1.0, 1.0)
        assert calc_dist(v1, v2) == pytest.approx(math.sqrt(3.0), abs=1e-10)


class TestIsOrigin:
    def test_origin(self):
        assert _is_origin(Vector(0.0, 0.0, 0.0))

    def test_not_origin(self):
        assert not _is_origin(Vector(1.0, 0.0, 0.0))

    def test_near_origin(self):
        # Below the threshold
        assert _is_origin(Vector(1e-15, 1e-15, 1e-15))


class TestAddJitter:
    def test_modifies_vector(self):
        v = Vector(0.0, 0.0, 0.0)
        jittered = _add_jitter(v)
        # Should not be exactly zero anymore
        assert not _is_origin(jittered)


# ---------------------------------------------------------------------------
# calc_wrapper
# ---------------------------------------------------------------------------

class TestCalcWrapper:
    def test_distance(self):
        v1 = Vector(0.0, 0.0, 0.0)
        v2 = Vector(3.0, 4.0, 0.0)
        result = calc_wrapper('distance', [v1, v2], '__?__')
        assert abs(float(result) - 5.0) < 0.01

    def test_angle(self):
        v1 = Vector(1.0, 0.0, 0.0)
        v2 = Vector(0.0, 0.0, 0.0)
        v3 = Vector(0.0, 1.0, 0.0)
        result = calc_wrapper('angle', [v1, v2, v3], '__?__')
        assert abs(float(result) - 90.0) < 0.01

    def test_dihedral(self):
        v1 = Vector(1.0, 0.0, 0.0)
        v2 = Vector(0.0, 0.0, 0.0)
        v3 = Vector(0.0, 1.0, 0.0)
        v4 = Vector(0.0, 1.0, 1.0)
        result = calc_wrapper('dihedral', [v1, v2, v3, v4], '__?__')
        assert result != '__?__'

    def test_wrong_arg_count_returns_unknown(self):
        v1 = Vector(1.0, 0.0, 0.0)
        v2 = Vector(0.0, 0.0, 0.0)
        assert calc_wrapper('dihedral', [v1, v2], '__?__') == '__?__'
        assert calc_wrapper('angle', [v1], '__?__') == '__?__'
        assert calc_wrapper('distance', [v1], '__?__') == '__?__'


# ---------------------------------------------------------------------------
# compute_measurement (requires structure fixtures from conftest.py)
# ---------------------------------------------------------------------------

class TestComputeMeasurement:
    def test_phi_middle_residue(self, hexapeptide_chain):
        """Phi should be computable for residue index 2 (PRO)."""
        commands = command_string_parser("phi")
        result = compute_measurement(
            commands[0], hexapeptide_chain, 2, '__?__'
        )
        assert result != '__?__'
        # Should be a valid number
        float(result)

    def test_phi_first_residue_is_unknown(self, hexapeptide_chain):
        """Phi requires i-1 C, which doesn't exist for the first residue."""
        commands = command_string_parser("phi")
        result = compute_measurement(
            commands[0], hexapeptide_chain, 0, '__?__'
        )
        assert result == '__?__'

    def test_psi_last_residue_is_unknown(self, hexapeptide_chain):
        """Psi requires i+1 N, which doesn't exist for the last residue."""
        commands = command_string_parser("psi")
        result = compute_measurement(
            commands[0], hexapeptide_chain, 5, '__?__'
        )
        assert result == '__?__'

    def test_tau_angle(self, hexapeptide_chain):
        """Tau (N-CA-C angle) should work for any middle residue."""
        commands = command_string_parser("tau")
        result = compute_measurement(
            commands[0], hexapeptide_chain, 2, '__?__'
        )
        assert result != '__?__'
        angle = float(result)
        # Tau should be roughly 100-120 degrees for standard protein geometry
        assert 80.0 < angle < 140.0

    def test_chi1_alanine_is_unknown(self, hexapeptide_chain):
        """ALA has CB but no CG/OG/SG → chi1 should be unknown."""
        commands = command_string_parser("chi1")
        result = compute_measurement(
            commands[0], hexapeptide_chain, 0, '__?__'
        )
        assert result == '__?__'

    def test_chi1_glutamine(self, hexapeptide_chain):
        """GLN (residue index 3) has CG → chi1 should be computable."""
        commands = command_string_parser("chi1")
        result = compute_measurement(
            commands[0], hexapeptide_chain, 3, '__?__'
        )
        assert result != '__?__'
        float(result)

    def test_distance_ca_ca(self, hexapeptide_chain):
        """Cα–Cα distance between adjacent residues."""
        commands = command_string_parser("vCAd")
        # Residue index 1 needs i-1 CA
        result = compute_measurement(
            commands[0], hexapeptide_chain, 1, '__?__'
        )
        assert result != '__?__'
        dist = float(result)
        # Virtual CA distance should be roughly 3.0-4.5 Å
        assert 2.0 < dist < 6.0

    def test_custom_distance(self, hexapeptide_chain):
        """Custom distance measurement."""
        commands = command_string_parser(
            "distance: test: i _CA_, i _C__"
        )
        result = compute_measurement(
            commands[0], hexapeptide_chain, 0, '__?__'
        )
        assert result != '__?__'
        dist = float(result)
        # CA-C bond length should be roughly 1.5 Å
        assert 1.0 < dist < 2.0


# ---------------------------------------------------------------------------
# process_measurement_for_residue
# ---------------------------------------------------------------------------

class TestProcessMeasurementForResidue:
    def test_returns_string_for_valid_residue(self, hexapeptide_chain):
        commands = command_string_parser("tau; phi; psi")
        result = process_measurement_for_residue(
            "test:1:A:", hexapeptide_chain, 2, commands
        )
        assert result is not None
        assert 'test:1:A:' in result

    def test_returns_none_when_all_unknown(self, hexapeptide_chain):
        """First residue has no phi; if we only measure phi, all is unknown."""
        commands = command_string_parser("phi")
        result = process_measurement_for_residue(
            "test:1:A:", hexapeptide_chain, 0, commands
        )
        assert result is None


# ---------------------------------------------------------------------------
# process_measurement_commands (end-to-end)
# ---------------------------------------------------------------------------

class TestProcessMeasurementCommands:
    def test_protein_phi_psi(self, hexapeptide_structure):
        lines = process_measurement_commands(
            "AGPQVS", hexapeptide_structure, "phi; psi"
        )
        # First line is the header comment
        assert lines[0].startswith('#')
        # Should have output for middle residues
        assert len(lines) > 1

    def test_protein_tau_phi_psi(self, hexapeptide_structure):
        lines = process_measurement_commands(
            "AGPQVS", hexapeptide_structure, "tau; phi; psi"
        )
        assert len(lines) > 1
        # Each data line should have colon-separated fields
        for line in lines[1:]:
            fields = line.split(':')
            assert len(fields) >= 5  # file:model:chain:resnum:icode:resname:...

    def test_protein_chi_values(self, hexapeptide_structure):
        lines = process_measurement_commands(
            "AGPQVS", hexapeptide_structure, "chi1; chi2"
        )
        # Should have at least some valid chi values (e.g., for GLN, VAL, SER)
        data_lines = [line for line in lines if not line.startswith('#')]
        assert len(data_lines) > 0

    def test_protein_distances(self, hexapeptide_structure):
        lines = process_measurement_commands(
            "AGPQVS", hexapeptide_structure, "vCAd; pbCACd"
        )
        data_lines = [line for line in lines if not line.startswith('#')]
        assert len(data_lines) > 0

    def test_output_format_consistency(self, hexapeptide_structure):
        """All data lines should have the same number of fields."""
        lines = process_measurement_commands(
            "AGPQVS", hexapeptide_structure, "tau; phi; psi"
        )
        data_lines = [line for line in lines if not line.startswith('#')]
        if data_lines:
            n_fields = len(data_lines[0].split(':'))
            for line in data_lines:
                assert len(line.split(':')) == n_fields

    def test_rna_backbone_angles(self, rna_structure):
        """RNA backbone angles should produce output for nucleotide residues."""
        lines = process_measurement_commands(
            "RNA", rna_structure, "delta"
        )
        data_lines = [line for line in lines if not line.startswith('#')]
        # delta is within a single residue so should work
        assert len(data_lines) >= 1

    def test_unknown_string_customizable(self, hexapeptide_structure):
        lines = process_measurement_commands(
            "AGPQVS", hexapeptide_structure, "phi; psi",
            unknown_str="N/A"
        )
        # Should still produce output
        assert len(lines) > 1
