"""Shared test fixtures and sample structure data for pydangle tests."""

import io

import pytest
from Bio.PDB import PDBParser

# ---------------------------------------------------------------------------
# Sample PDB: hexapeptide from ubiquitin residues 34-39 (chain A)
# Sequence: GLU GLY ILE PRO PRO ASP
# Source: 1ubq.pdb (Vijay-Kumar et al., 1987)
# ---------------------------------------------------------------------------
HEXAPEPTIDE_PDB = """\
HEADER    HEXAPEPTIDE FROM UBIQUITIN 34-39
ATOM      1  N   GLU A   1      39.655  34.335  11.285  1.00 10.11           N
ATOM      2  CA  GLU A   1      39.676  35.547  12.072  1.00 10.07           C
ATOM      3  C   GLU A   1      40.675  35.527  13.200  1.00  9.32           C
ATOM      4  O   GLU A   1      40.814  36.528  13.911  1.00 11.61           O
ATOM      5  CB  GLU A   1      38.290  35.814  12.698  1.00 14.77           C
ATOM      6  CG  GLU A   1      37.156  35.985  11.688  1.00 18.75           C
ATOM      7  CD  GLU A   1      37.192  37.361  11.033  1.00 22.28           C
ATOM      8  OE1 GLU A   1      37.519  38.360  11.645  1.00 21.95           O
ATOM      9  OE2 GLU A   1      36.861  37.320   9.822  1.00 25.19           O
ATOM     10  N   GLY A   2      41.317  34.393  13.432  1.00  7.22           N
ATOM     11  CA  GLY A   2      42.345  34.269  14.431  1.00  6.29           C
ATOM     12  C   GLY A   2      41.949  34.076  15.842  1.00  6.93           C
ATOM     13  O   GLY A   2      42.829  34.000  16.739  1.00  7.41           O
ATOM     14  N   ILE A   3      40.642  33.916  16.112  1.00  5.86           N
ATOM     15  CA  ILE A   3      40.226  33.716  17.509  1.00  6.07           C
ATOM     16  C   ILE A   3      40.449  32.278  17.945  1.00  6.36           C
ATOM     17  O   ILE A   3      39.936  31.336  17.315  1.00  6.18           O
ATOM     18  CB  ILE A   3      38.693  34.106  17.595  1.00  7.47           C
ATOM     19  CG1 ILE A   3      38.471  35.546  17.045  1.00  8.52           C
ATOM     20  CG2 ILE A   3      38.146  33.932  19.027  1.00  7.36           C
ATOM     21  CD1 ILE A   3      36.958  35.746  16.680  1.00  9.49           C
ATOM     22  N   PRO A   4      41.189  32.085  19.031  1.00  8.65           N
ATOM     23  CA  PRO A   4      41.461  30.751  19.594  1.00  9.18           C
ATOM     24  C   PRO A   4      40.168  30.026  19.918  1.00  9.85           C
ATOM     25  O   PRO A   4      39.264  30.662  20.521  1.00  8.51           O
ATOM     26  CB  PRO A   4      42.195  31.142  20.913  1.00 11.42           C
ATOM     27  CG  PRO A   4      42.904  32.414  20.553  1.00  9.27           C
ATOM     28  CD  PRO A   4      41.822  33.188  19.813  1.00  8.33           C
ATOM     29  N   PRO A   5      40.059  28.758  19.607  1.00  8.71           N
ATOM     30  CA  PRO A   5      38.817  28.020  19.889  1.00  9.08           C
ATOM     31  C   PRO A   5      38.421  28.048  21.341  1.00  9.28           C
ATOM     32  O   PRO A   5      37.213  28.036  21.704  1.00  6.50           O
ATOM     33  CB  PRO A   5      39.090  26.629  19.325  1.00 10.31           C
ATOM     34  CG  PRO A   5      40.082  26.904  18.198  1.00 10.81           C
ATOM     35  CD  PRO A   5      41.035  27.909  18.879  1.00 12.00           C
ATOM     36  N   ASP A   6      39.374  28.090  22.240  1.00 11.20           N
ATOM     37  CA  ASP A   6      39.063  28.063  23.695  1.00 14.96           C
ATOM     38  C   ASP A   6      38.365  29.335  24.159  1.00 13.99           C
ATOM     39  O   ASP A   6      37.684  29.390  25.221  1.00 13.75           O
ATOM     40  CB  ASP A   6      40.340  27.692  24.468  1.00 24.16           C
ATOM     41  CG  ASP A   6      40.559  28.585  25.675  1.00 31.06           C
ATOM     42  OD1 ASP A   6      40.716  29.809  25.456  1.00 35.55           O
ATOM     43  OD2 ASP A   6      40.549  28.090  26.840  1.00 34.22           O
END
"""

# ---------------------------------------------------------------------------
# Sample PDB: short RNA duplex (minimal)
# ---------------------------------------------------------------------------
RNA_DINUCLEOTIDE_PDB = """\
HEADER    RNA DINUCLEOTIDE
ATOM      1  P     G A   1      27.240  25.930  17.060  1.00 30.00           P
ATOM      2  OP1   G A   1      27.980  27.140  16.680  1.00 30.00           O
ATOM      3  OP2   G A   1      27.720  24.630  16.560  1.00 30.00           O
ATOM      4  O5'   G A   1      25.750  26.120  16.620  1.00 30.00           O
ATOM      5  C5'   G A   1      24.780  25.150  16.990  1.00 30.00           C
ATOM      6  C4'   G A   1      23.420  25.690  16.630  1.00 30.00           C
ATOM      7  O4'   G A   1      23.050  26.780  17.510  1.00 30.00           O
ATOM      8  C3'   G A   1      23.370  26.310  15.240  1.00 30.00           C
ATOM      9  O3'   G A   1      22.470  25.570  14.430  1.00 30.00           O
ATOM     10  C2'   G A   1      22.870  27.710  15.570  1.00 30.00           C
ATOM     11  O2'   G A   1      21.480  27.710  15.800  1.00 30.00           O
ATOM     12  C1'   G A   1      23.100  28.390  14.230  1.00 30.00           C
ATOM     13  N9    G A   1      24.450  28.910  14.100  1.00 30.00           N
ATOM     14  C4    G A   1      24.850  30.070  13.490  1.00 30.00           C
ATOM     15  P     C A   2      22.670  24.070  13.860  1.00 30.00           P
ATOM     16  OP1   C A   2      24.100  23.780  13.690  1.00 30.00           O
ATOM     17  OP2   C A   2      21.780  23.170  14.650  1.00 30.00           O
ATOM     18  O5'   C A   2      22.100  24.050  12.380  1.00 30.00           O
ATOM     19  C5'   C A   2      22.830  24.670  11.330  1.00 30.00           C
ATOM     20  C4'   C A   2      22.000  24.650  10.070  1.00 30.00           C
ATOM     21  O4'   C A   2      21.920  23.280   9.580  1.00 30.00           O
ATOM     22  C3'   C A   2      20.550  25.100  10.230  1.00 30.00           C
ATOM     23  O3'   C A   2      20.370  26.500  10.080  1.00 30.00           O
ATOM     24  C2'   C A   2      19.850  24.290   9.150  1.00 30.00           C
ATOM     25  O2'   C A   2      19.970  24.880   7.870  1.00 30.00           O
ATOM     26  C1'   C A   2      20.570  22.950   9.230  1.00 30.00           C
ATOM     27  N1    C A   2      20.590  22.190  10.500  1.00 30.00           N
ATOM     28  C2    C A   2      19.420  21.570  10.920  1.00 30.00           C
END
"""


# ---------------------------------------------------------------------------
# Pytest fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def hexapeptide_structure():
    """Return a BioPython Structure for the ubiquitin hexapeptide."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure("EGIPpd", io.StringIO(HEXAPEPTIDE_PDB))


@pytest.fixture
def hexapeptide_chain(hexapeptide_structure):
    """Return chain A from the hexapeptide structure."""
    return hexapeptide_structure[0]["A"]


@pytest.fixture
def rna_structure():
    """Return a BioPython Structure for the RNA dinucleotide."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure("RNA", io.StringIO(RNA_DINUCLEOTIDE_PDB))


@pytest.fixture
def rna_chain(rna_structure):
    """Return chain A from the RNA structure."""
    return rna_structure[0]["A"]
