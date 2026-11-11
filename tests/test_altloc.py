"""Tests for alternate-conformation selection (altloc.py).

Fixtures are real coordinates from the ersatz files of the Richardson
datasets.  Expected angles are computed here from the atoms of the
conformer that should be used, except for 6b6w, where they are the values
in C. J. Williams' Top2018 database (cctbx): 6b6w, 3zeu, 5k0a.
"""

import io
import json
import math

import pytest
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_angle, calc_dihedral

from pydangle_biopython.measure import process_measurement_commands
from tests.conftest import ALTLOC_NTERM_PDB

# 1din A122-124 (top500): CSD altloc A written as HETATM, CYS altloc B as
# ATOM -- BioPython makes them two residues, ('H_CSD', 123) and (' ', 123).
SPLIT_MICROHET_PDB = """\
ATOM    927  N   TYR A 122      17.351  49.543   9.416  1.00  7.19           N
ATOM    928  CA  TYR A 122      17.068  50.014   8.065  1.00  9.39           C
ATOM    929  C   TYR A 122      15.617  49.640   7.780  1.00 10.73           C
ATOM    930  O   TYR A 122      15.128  48.520   7.971  1.00 11.30           O
ATOM    931  CB  TYR A 122      17.942  49.245   7.034  1.00 12.82           C
HETATM  939  N  ACSD A 123      14.951  50.672   7.287  0.60 10.41           N
HETATM  940  CA ACSD A 123      13.530  50.632   6.921  0.60 22.56           C
HETATM  941  CB ACSD A 123      13.500  49.608   5.835  0.60 16.52           C
HETATM  942  SG ACSD A 123      12.125  49.265   4.790  0.60 35.39           S
HETATM  943  C  ACSD A 123      12.602  50.101   8.026  0.60 14.09           C
HETATM  944  O  ACSD A 123      12.546  51.002   8.873  0.60 16.11           O
HETATM  945  OD1ACSD A 123      12.471  48.128   4.979  0.60 25.92           O
HETATM  946  OD2ACSD A 123      11.149  50.292   4.778  0.60 32.06           O
ATOM    947  N  BCYS A 123      14.951  50.672   7.287  0.40 10.41           N
ATOM    948  CA BCYS A 123      13.530  50.632   6.921  0.40 22.56           C
ATOM    949  C  BCYS A 123      12.602  50.101   8.026  0.40 14.09           C
ATOM    950  O  BCYS A 123      12.546  51.002   8.873  0.40 16.11           O
ATOM    951  CB BCYS A 123      13.500  49.608   5.835  0.40 16.52           C
ATOM    952  SG BCYS A 123      14.425  49.910   4.378  0.40 21.48           S
ATOM    953  N   LEU A 124      11.990  48.849   8.089  1.00 13.58           N
ATOM    954  CA  LEU A 124      11.158  48.493   9.201  1.00 11.40           C
ATOM    955  C   LEU A 124      12.131  48.690  10.371  1.00 12.06           C
ATOM    956  O   LEU A 124      11.716  49.254  11.380  1.00 12.09           O
ATOM    957  CB  LEU A 124      10.891  47.001   9.058  1.00 11.85           C
END
"""

# 2zs0 D30-32 (top2018_mc): ILE altloc A / VAL altloc B, occupancy 0.50
# each, all ATOM -- a BioPython DisorderedResidue.  Java Dangle 1.07
# measured the Ile; pydangle before v0.6 the Val.
DISORDERED_RESIDUE_PDB = """\
ATOM   3613  N   ALA D  30      30.091  21.296  28.133  1.00 19.88           N
ATOM   3614  CA  ALA D  30      30.004  20.470  29.328  1.00 19.24           C
ATOM   3615  C   ALA D  30      30.512  19.052  29.063  1.00 20.04           C
ATOM   3616  O   ALA D  30      31.224  18.471  29.888  1.00 21.54           O
ATOM   3617  CB  ALA D  30      28.554  20.435  29.827  1.00 19.89           C
ATOM   3618  N  AILE D  31      30.110  18.503  27.933  0.50 21.22           N
ATOM   3619  CA AILE D  31      30.583  17.150  27.543  0.50 21.88           C
ATOM   3620  C  AILE D  31      32.101  17.171  27.404  0.50 21.55           C
ATOM   3621  O  AILE D  31      32.791  16.304  27.930  0.50 21.75           O
ATOM   3622  CB AILE D  31      29.974  16.708  26.183  0.50 21.66           C
ATOM   3623  CG1AILE D  31      28.454  16.528  26.310  0.50 26.29           C
ATOM   3624  CG2AILE D  31      30.636  15.412  25.711  0.50 25.18           C
ATOM   3625  CD1AILE D  31      28.037  15.386  27.207  0.50 28.03           C
ATOM   3626  N  BVAL D  31      30.132  18.493  27.922  0.50 19.72           N
ATOM   3627  CA BVAL D  31      30.572  17.154  27.560  0.50 19.82           C
ATOM   3628  C  BVAL D  31      32.089  17.161  27.399  0.50 19.74           C
ATOM   3629  O  BVAL D  31      32.787  16.305  27.926  0.50 20.32           O
ATOM   3630  CB BVAL D  31      29.917  16.689  26.234  0.50 18.53           C
ATOM   3631  CG1BVAL D  31      30.530  15.366  25.776  0.50 23.03           C
ATOM   3632  CG2BVAL D  31      28.404  16.526  26.437  0.50 22.61           C
ATOM   3633  N   GLY D  32      32.582  18.157  26.675  1.00 20.49           N
ATOM   3634  CA  GLY D  32      34.009  18.266  26.425  1.00 18.19           C
ATOM   3635  C   GLY D  32      34.873  18.469  27.654  1.00 18.20           C
ATOM   3636  O   GLY D  32      35.960  17.877  27.765  1.00 18.67           O
END
"""

# 6b6w B300-303 (top2018_full): CYS B301 exists only as altloc B; CYS B302
# has A (0.20) and B (0.80); THR B303 has none.
B_ONLY_RESIDUE_PDB = """\
ATOM   7094  N   LEU B 300      14.502  46.480  17.817  1.00 12.59           N
ATOM   7095  CA  LEU B 300      14.179  45.764  16.594  1.00 11.38           C
ATOM   7096  C   LEU B 300      13.703  44.365  16.961  1.00 11.19           C
ATOM   7097  O   LEU B 300      14.103  43.810  17.989  1.00 15.22           O
ATOM   7098  CB  LEU B 300      15.391  45.648  15.649  1.00 12.06           C
ATOM   7102  N  BCYS B 301      12.834  43.814  16.120  1.00 16.49           N
ATOM   7103  CA BCYS B 301      12.459  42.393  16.168  1.00 16.12           C
ATOM   7104  C  BCYS B 301      11.823  42.070  17.526  1.00 16.30           C
ATOM   7105  O  BCYS B 301      11.135  42.904  18.134  1.00 15.44           O
ATOM   7106  CB BCYS B 301      13.683  41.543  15.820  1.00 15.03           C
ATOM   7107  SG BCYS B 301      13.376  39.752  15.513  1.00 23.40           S
ATOM   7108  N  ACYS B 302      12.296  40.961  17.963  0.20 12.79           N
ATOM   7109  N  BCYS B 302      12.044  40.851  18.012  0.80 12.40           N
ATOM   7110  CA ACYS B 302      11.411  40.595  19.064  0.20 13.23           C
ATOM   7111  CA BCYS B 302      11.359  40.392  19.214  0.80 12.21           C
ATOM   7112  C  ACYS B 302      11.956  41.064  20.410  0.20 12.66           C
ATOM   7113  C  BCYS B 302      11.960  40.983  20.481  0.80 13.02           C
ATOM   7114  O  ACYS B 302      11.187  41.356  21.324  0.20 14.03           O
ATOM   7115  O  BCYS B 302      11.228  41.292  21.428  0.80 14.19           O
ATOM   7116  CB ACYS B 302      11.186  39.084  19.086  0.20 16.22           C
ATOM   7117  CB BCYS B 302      11.388  38.867  19.272  0.80 14.15           C
ATOM   7120  N   THR B 303      13.284  41.123  20.537  1.00 13.48           N
ATOM   7121  CA  THR B 303      13.874  41.730  21.721  1.00 15.16           C
ATOM   7122  C   THR B 303      13.419  43.178  21.851  1.00 10.76           C
ATOM   7123  O   THR B 303      13.223  43.671  22.968  1.00 13.81           O
END
"""


# 3zeu A1 and B167-169 (top2018_mc): the file's first alternates are MET A1
# conformers B and C (no A); LYS B168 has A (0.27, written first) and B
# (0.73).  The Richardson Lab measured B168 in conformer A.
LOCAL_ORDER_PDB = """\
ATOM      1  N  BMET A   1       9.680   4.996  24.524  0.55 38.87           N
ATOM      2  N  CMET A   1       9.700   5.567  24.971  0.45 34.77           N
ATOM      3  CA BMET A   1       8.692   5.981  24.103  0.55 35.62           C
ATOM      4  CA CMET A   1       8.688   5.735  23.936  0.45 35.58           C
ATOM      5  C  BMET A   1       7.299   5.574  24.574  0.55 26.00           C
ATOM      6  C  CMET A   1       7.300   5.420  24.469  0.45 25.87           C
ATOM      7  O  BMET A   1       6.288   6.001  24.015  0.55 32.55           O
ATOM      8  O  CMET A   1       6.294   5.796  23.870  0.45 31.59           O
TER
ATOM   3049  N   ASP B 167      -4.360  26.605  73.362  1.00 11.22           N
ATOM   3050  CA  ASP B 167      -4.059  27.535  72.269  1.00 15.16           C
ATOM   3051  C   ASP B 167      -4.880  28.821  72.371  1.00 11.41           C
ATOM   3052  O   ASP B 167      -4.388  29.920  72.061  1.00 12.79           O
ATOM   3057  N  ALYS B 168      -6.130  28.688  72.811  0.27 11.78           N
ATOM   3058  N  BLYS B 168      -6.126  28.698  72.824  0.73 11.65           N
ATOM   3059  CA ALYS B 168      -7.016  29.843  72.935  0.27 13.64           C
ATOM   3060  CA BLYS B 168      -6.987  29.879  72.913  0.73 13.55           C
ATOM   3061  C  ALYS B 168      -6.547  30.793  74.034  0.27 14.88           C
ATOM   3062  C  BLYS B 168      -6.545  30.807  74.042  0.73 14.93           C
ATOM   3063  O  ALYS B 168      -6.535  32.009  73.847  0.27 15.60           O
ATOM   3064  O  BLYS B 168      -6.546  32.026  73.885  0.73 15.59           O
ATOM   3075  N   THR B 169      -6.162  30.239  75.180  1.00 11.52           N
ATOM   3076  CA  THR B 169      -5.648  31.066  76.267  1.00 14.31           C
ATOM   3077  C   THR B 169      -4.281  31.657  75.917  1.00 13.71           C
ATOM   3078  O   THR B 169      -3.995  32.810  76.244  1.00 13.89           O
END
"""

# 5k0a A155-159 (top2018_mc): CSO altloc A written as HETATM (first), CYS
# altloc B as ATOM.  The Richardson Lab's phi/psi at A156-158 are those of
# the CYS conformer.
MODIFIED_FIRST_PDB = """\
ATOM   2146  N   HIS A 155       9.625  30.130  -5.182  1.00 22.88           N
ATOM   2147  CA  HIS A 155       8.873  29.425  -6.203  1.00 21.25           C
ATOM   2148  C   HIS A 155       8.345  28.096  -5.686  1.00 29.56           C
ATOM   2149  O   HIS A 155       8.808  27.562  -4.674  1.00 28.57           O
ATOM   2150  CB  HIS A 155       9.749  29.210  -7.432  1.00 26.38           C
ATOM   2163  N   VAL A 156       7.342  27.574  -6.386  1.00 25.03           N
ATOM   2164  CA  VAL A 156       6.795  26.261  -6.062  1.00 24.83           C
ATOM   2165  C   VAL A 156       6.989  25.253  -7.176  1.00 28.44           C
ATOM   2166  O   VAL A 156       6.722  24.062  -6.961  1.00 27.55           O
ATOM   2167  CB  VAL A 156       5.303  26.340  -5.669  1.00 25.09           C
HETATM 2179  N  ACSO A 157       7.358  25.745  -8.354  0.49 31.03           N
HETATM 2180  CA ACSO A 157       7.580  24.893  -9.514  0.49 24.25           C
HETATM 2181  CB ACSO A 157       6.569  25.244 -10.610  0.49 28.08           C
HETATM 2182  SG ACSO A 157       6.852  24.264 -12.046  0.49 34.17           S
HETATM 2183  C  ACSO A 157       8.980  25.066 -10.036  0.49 28.57           C
HETATM 2184  O  ACSO A 157       9.241  26.137 -10.621  0.49 30.52           O
HETATM 2185  OD ACSO A 157       6.749  22.641 -11.448  0.49 35.77           O
ATOM   2189  N  BCYS A 157       7.521  25.676  -8.326  0.51 21.86           N
ATOM   2190  CA BCYS A 157       7.525  24.886  -9.556  0.51 27.11           C
ATOM   2191  C  BCYS A 157       8.895  25.029 -10.214  0.51 23.57           C
ATOM   2192  O  BCYS A 157       9.231  26.105 -10.716  0.51 25.17           O
ATOM   2193  CB BCYS A 157       6.410  25.362 -10.496  0.51 22.73           C
ATOM   2194  SG BCYS A 157       6.526  24.778 -12.204  0.51 27.33           S
ATOM   2200  N   MET A 158       9.683  23.953 -10.230  1.00 25.15           N
ATOM   2201  CA  MET A 158      11.035  24.068 -10.743  1.00 38.69           C
ATOM   2202  C   MET A 158      11.048  24.303 -12.249  1.00 27.46           C
ATOM   2203  O   MET A 158      11.840  25.104 -12.746  1.00 31.12           O
ATOM   2204  CB  MET A 158      11.812  22.820 -10.390  1.00 42.58           C
ATOM   2217  N   ILE A 159      10.186  23.618 -13.000  1.00 24.19           N
ATOM   2218  CA  ILE A 159      10.262  23.803 -14.448  1.00 26.34           C
ATOM   2219  C   ILE A 159       9.766  25.187 -14.836  1.00 30.97           C
ATOM   2220  O   ILE A 159      10.254  25.773 -15.816  1.00 29.58           O
ATOM   2221  CB  ILE A 159       9.519  22.686 -15.200  1.00 30.47           C
END
"""


def _structure(text, name="x"):
    return PDBParser(QUIET=True).get_structure(name, io.StringIO(text))


def _records(structure, commands, mode):
    lines = process_measurement_commands(
        "x", structure, commands, output_format="jsonl", altloc_mode=mode
    )
    return [json.loads(line) for line in lines]


def _by_key(records):
    return {(r["chain"], r["resnum"]): r for r in records}


def _deg(radians):
    return math.degrees(radians)


class TestFirstConformer:
    """Default mode: the first conformer, whole-residue."""

    def test_tau_from_conformer_a(self, altloc_nterm_structure):
        res = altloc_nterm_structure[0]["A"][2]
        expected = _deg(calc_angle(
            res["N"].disordered_get("A").get_vector(),
            res["CA"].disordered_get("A").get_vector(),
            res["C"].disordered_get("A").get_vector(),
        ))
        records = _by_key(_records(altloc_nterm_structure, "tau", "first"))
        assert records[("A", 2)]["tau"] == pytest.approx(expected, abs=0.001)

    def test_occupancy_mode_unchanged(self, altloc_nterm_structure):
        records = _by_key(_records(altloc_nterm_structure, "tau", "occupancy"))
        # Conformer B (occupancy 0.51), as before v0.6
        assert records[("A", 2)]["tau"] == pytest.approx(106.719, abs=0.001)

    def test_every_disordered_atom_selects_first(self, altloc_nterm_structure):
        _records(altloc_nterm_structure, "tau", "first")
        selected = {
            atom.get_altloc()
            for atom in altloc_nterm_structure.get_atoms()
            if atom.is_disordered()
        }
        assert selected == {"A"}

    def test_first_is_per_residue_not_whole_file(self):
        # 3zeu: B and C appear first in the file (A1), but B168 lists A first.
        records = _by_key(_records(_structure(LOCAL_ORDER_PDB), "phi; psi", "first"))
        # Williams Top2018 SQLite (cctbx), 3zeu B168
        assert records[("B", 168)]["phi"] == pytest.approx(-67.634, abs=0.002)
        assert records[("B", 168)]["psi"] == pytest.approx(-44.071, abs=0.002)

    def test_order_is_file_order_not_alphabetical(self):
        # Swap the labels: B now comes first in the file and carries the
        # coordinates that A had, so the result must equal the A result.
        swapped = "".join(
            line[:16] + {"A": "B", "B": "A"}.get(line[16], line[16]) + line[17:]
            if line.startswith("ATOM") else line
            for line in ALTLOC_NTERM_PDB.splitlines(keepends=True)
        )
        original = _by_key(_records(_structure(ALTLOC_NTERM_PDB), "tau", "first"))
        relabelled = _by_key(_records(_structure(swapped), "tau", "first"))
        assert relabelled[("A", 2)]["tau"] == pytest.approx(
            original[("A", 2)]["tau"], abs=1e-6
        )


class TestMicroheterogeneity:
    """One residue type per position in first mode: the first standard one."""

    def test_split_hetatm_atom_pair_one_record(self):
        # CSD (altloc A) is first, but CYS is the standard residue type.
        records = _records(_structure(SPLIT_MICROHET_PDB), "phi; psi", "first")
        at_123 = [r for r in records if r["resnum"] == 123]
        assert [r["resname"] for r in at_123] == ["CYS"]

    def test_split_pair_neighbours_use_kept_type(self):
        s = _structure(SPLIT_MICROHET_PDB)
        cys = s[0]["A"][123]
        leu = s[0]["A"][124]
        expected_phi = _deg(calc_dihedral(
            cys["C"].get_vector(), leu["N"].get_vector(),
            leu["CA"].get_vector(), leu["C"].get_vector(),
        ))
        records = _by_key(_records(s, "phi; psi", "first"))
        assert records[("A", 124)]["phi"] == pytest.approx(expected_phi, abs=0.001)
        assert records[("A", 123)]["psi"] is not None

    def test_standard_type_preferred_matches_lab_database(self):
        # 5k0a A157: CSO is conformer A, CYS conformer B; the lab's values
        # at A156-158 are those of the CYS conformer.
        records = _by_key(
            _records(_structure(MODIFIED_FIRST_PDB), "phi; psi", "first")
        )
        # Williams Top2018 SQLite (cctbx), 5k0a A156-158
        assert records[("A", 156)]["psi"] == pytest.approx(-4.937, abs=0.002)
        assert records[("A", 157)]["resname"] == "CYS"
        assert records[("A", 157)]["phi"] == pytest.approx(-137.641, abs=0.002)
        assert records[("A", 157)]["psi"] == pytest.approx(111.071, abs=0.002)
        assert records[("A", 158)]["phi"] == pytest.approx(-69.090, abs=0.002)

    def test_split_pair_occupancy_mode_two_records(self):
        records = _records(_structure(SPLIT_MICROHET_PDB), "phi; psi", "occupancy")
        assert sorted(r["resname"] for r in records if r["resnum"] == 123) == [
            "CSD", "CYS"
        ]

    def test_disordered_residue_first_type(self):
        s = _structure(DISORDERED_RESIDUE_PDB)
        ile = s[0]["D"][31].disordered_get("ILE")
        expected_tau = _deg(calc_angle(
            ile["N"].get_vector(), ile["CA"].get_vector(), ile["C"].get_vector()
        ))
        records = _by_key(_records(s, "tau", "first"))
        assert records[("D", 31)]["resname"] == "ILE"
        assert records[("D", 31)]["tau"] == pytest.approx(expected_tau, abs=0.001)

    def test_disordered_residue_occupancy_mode(self):
        records = _by_key(
            _records(_structure(DISORDERED_RESIDUE_PDB), "tau", "occupancy")
        )
        assert records[("D", 31)]["resname"] == "VAL"


class TestConformerConsistentAcrossNeighbours:
    """A residue present only as B is measured with its neighbours' B atoms."""

    def test_b_only_residue_and_next_match_lab_database(self):
        records = _by_key(_records(_structure(B_ONLY_RESIDUE_PDB), "phi; psi", "first"))
        # Williams Top2018 SQLite (cctbx), 6b6w B301 and B302
        assert records[("B", 301)]["phi"] == pytest.approx(60.608, abs=0.002)
        assert records[("B", 301)]["psi"] == pytest.approx(-146.704, abs=0.002)
        assert records[("B", 302)]["phi"] == pytest.approx(-76.932, abs=0.002)
        assert records[("B", 302)]["psi"] == pytest.approx(-38.206, abs=0.002)
