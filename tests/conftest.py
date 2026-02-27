"""Shared test fixtures and sample structure data for pydangle tests."""

import io

import pytest
from Bio.PDB import PDBParser

# ---------------------------------------------------------------------------
# Sample PDB: hexapeptide AGPQVS
# ---------------------------------------------------------------------------
HEXAPEPTIDE_PDB = """\
HEADER    HEXAPEPTIDE AGPQVS
TITLE     MODELED STRUCTURE OF AGPQVS HEXAPEPTIDE
ATOM      1  N   ALA A   1      -0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.010   1.410   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.270   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       2.000  -0.780   1.200  1.00  0.00           C
ATOM      6  N   GLY A   2       3.330   1.520   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       4.010   2.800   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2       3.590   3.640  -1.200  1.00  0.00           C
ATOM      9  O   GLY A   2       3.920   4.830  -1.240  1.00  0.00           O
ATOM     10  N   PRO A   3       2.870   2.990  -2.140  1.00  0.00           N
ATOM     11  CA  PRO A   3       2.430   3.720  -3.330  1.00  0.00           C
ATOM     12  C   PRO A   3       3.350   3.440  -4.520  1.00  0.00           C
ATOM     13  O   PRO A   3       3.340   4.170  -5.520  1.00  0.00           O
ATOM     14  CB  PRO A   3       0.990   3.210  -3.480  1.00  0.00           C
ATOM     15  CG  PRO A   3       0.790   2.310  -2.290  1.00  0.00           C
ATOM     16  CD  PRO A   3       1.750   2.660  -1.180  1.00  0.00           C
ATOM     17  N   GLN A   4       4.150   2.390  -4.400  1.00  0.00           N
ATOM     18  CA  GLN A   4       5.080   2.030  -5.460  1.00  0.00           C
ATOM     19  C   GLN A   4       4.440   0.940  -6.310  1.00  0.00           C
ATOM     20  O   GLN A   4       4.990   0.670  -7.390  1.00  0.00           O
ATOM     21  CB  GLN A   4       6.440   1.580  -4.900  1.00  0.00           C
ATOM     22  CG  GLN A   4       7.310   2.730  -4.360  1.00  0.00           C
ATOM     23  CD  GLN A   4       8.730   2.250  -4.140  1.00  0.00           C
ATOM     24  OE1 GLN A   4       9.140   1.170  -4.590  1.00  0.00           O
ATOM     25  NE2 GLN A   4       9.470   3.070  -3.410  1.00  0.00           N
ATOM     26  N   VAL A   5       3.260   0.310  -5.830  1.00  0.00           N
ATOM     27  CA  VAL A   5       2.560  -0.740  -6.550  1.00  0.00           C
ATOM     28  C   VAL A   5       1.150  -0.340  -6.930  1.00  0.00           C
ATOM     29  O   VAL A   5       0.460  -1.070  -7.650  1.00  0.00           O
ATOM     30  CB  VAL A   5       2.520  -2.070  -5.750  1.00  0.00           C
ATOM     31  CG1 VAL A   5       1.850  -3.150  -6.580  1.00  0.00           C
ATOM     32  CG2 VAL A   5       3.920  -2.500  -5.380  1.00  0.00           C
ATOM     33  N   SER A   6       0.730   0.810  -6.420  1.00  0.00           N
ATOM     34  CA  SER A   6      -0.590   1.300  -6.720  1.00  0.00           C
ATOM     35  C   SER A   6      -0.650   2.780  -6.370  1.00  0.00           C
ATOM     36  O   SER A   6      -1.680   3.430  -6.580  1.00  0.00           O
ATOM     37  CB  SER A   6      -1.010   1.090  -8.180  1.00  0.00           C
ATOM     38  OG  SER A   6      -0.230   1.870  -9.050  1.00  0.00           O
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
    """Return a BioPython Structure for the AGPQVS hexapeptide."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure('AGPQVS', io.StringIO(HEXAPEPTIDE_PDB))


@pytest.fixture
def hexapeptide_chain(hexapeptide_structure):
    """Return chain A from the hexapeptide structure."""
    return hexapeptide_structure[0]['A']


@pytest.fixture
def rna_structure():
    """Return a BioPython Structure for the RNA dinucleotide."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure('RNA', io.StringIO(RNA_DINUCLEOTIDE_PDB))


@pytest.fixture
def rna_chain(rna_structure):
    """Return chain A from the RNA structure."""
    return rna_structure[0]['A']
