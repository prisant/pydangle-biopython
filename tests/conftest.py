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


# ---------------------------------------------------------------------------
# Sample PDB: 1gk8 chain C residues 254-258 (ALA VAL SMC ALA LYS), heavy atoms
# SMC 256 (S-methylcysteine) is a HETATM modified residue inside a helix,
# not on the Java Dangle whitelist.  Regression fixture for the continuity
# bug reported by C. Williams 2026-11-07.  Coordinates from the top2018
# ersatz file (Reduce hydrogens removed).  Reference values, Java Dangle 1.07
# on the full 1gk8_ersatz.pdb:
#   VAL 255  psi -42.763
#   ALA 257  phi -56.602  psi -47.184  omega 173.232
# Source: 1gk8.pdb (Taylor et al., 2001)
# ---------------------------------------------------------------------------
SMC_HELIX_PDB = """\
HEADER    1GK8 CHAIN C 254-258 WITH SMC 256
ATOM   5631  N   ALA C 254      19.521  30.626   2.278  1.00  6.96           N
ATOM   5632  CA  ALA C 254      19.580  31.938   2.895  1.00  7.16           C
ATOM   5633  C   ALA C 254      18.535  32.058   3.994  1.00  7.26           C
ATOM   5634  O   ALA C 254      18.816  32.627   5.068  1.00  7.13           O
ATOM   5635  CB  ALA C 254      19.372  33.008   1.835  1.00  7.07           C
ATOM   5636  N   VAL C 255      17.314  31.582   3.742  1.00  7.12           N
ATOM   5637  CA  VAL C 255      16.258  31.637   4.748  1.00  7.55           C
ATOM   5638  C   VAL C 255      16.646  30.874   6.006  1.00  7.07           C
ATOM   5639  O   VAL C 255      16.413  31.333   7.130  1.00  7.71           O
ATOM   5640  CB  VAL C 255      14.929  31.100   4.175  1.00  7.63           C
ATOM   5641  CG1 VAL C 255      13.897  30.852   5.261  1.00  7.75           C
ATOM   5642  CG2 VAL C 255      14.384  32.077   3.139  1.00  7.38           C
HETATM 5643  N   SMC C 256      17.238  29.695   5.838  1.00  6.89           N
HETATM 5644  CA  SMC C 256      17.659  28.913   7.015  1.00  7.83           C
HETATM 5645  CB  SMC C 256      18.177  27.558   6.551  1.00  8.04           C
HETATM 5646  SG  SMC C 256      18.961  26.581   7.847  1.00  6.85           S
HETATM 5647  CS  SMC C 256      17.632  26.423   9.014  1.00  8.65           C
HETATM 5648  C   SMC C 256      18.728  29.660   7.795  1.00  7.77           C
HETATM 5649  O   SMC C 256      18.619  29.726   9.029  1.00  6.93           O
ATOM   5650  N   ALA C 257      19.712  30.272   7.149  1.00  6.80           N
ATOM   5651  CA  ALA C 257      20.669  31.118   7.858  1.00  7.56           C
ATOM   5652  C   ALA C 257      19.977  32.228   8.623  1.00  8.15           C
ATOM   5653  O   ALA C 257      20.306  32.497   9.792  1.00  8.05           O
ATOM   5654  CB  ALA C 257      21.673  31.696   6.870  1.00  8.00           C
ATOM   5655  N   LYS C 258      19.041  32.948   8.020  1.00  7.94           N
ATOM   5656  CA  LYS C 258      18.232  33.954   8.685  1.00  7.97           C
ATOM   5657  C   LYS C 258      17.586  33.392   9.943  1.00  7.79           C
ATOM   5658  O   LYS C 258      17.638  33.993  11.022  1.00  7.92           O
ATOM   5659  CB  LYS C 258      17.146  34.490   7.735  1.00  8.62           C
ATOM   5660  CG  LYS C 258      16.130  35.415   8.381  1.00  9.44           C
ATOM   5661  CD  LYS C 258      14.944  35.736   7.485  1.00 10.13           C
ATOM   5662  CE  LYS C 258      13.855  36.471   8.234  1.00  9.89           C
ATOM   5663  NZ  LYS C 258      14.173  37.896   8.494  1.00 10.16           N
END
"""


@pytest.fixture
def smc_helix_structure():
    """Return 1gk8 C 254-258, with the SMC 256 HETATM residue."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure("1gk8", io.StringIO(SMC_HELIX_PDB))


@pytest.fixture
def smc_deleted_structure():
    """Return 1gk8 C 254-258 with SMC 256 removed (a real chain break)."""
    text = "".join(
        line + "\n"
        for line in SMC_HELIX_PDB.splitlines()
        if not line.startswith("HETATM")
    )
    parser = PDBParser(QUIET=True)
    return parser.get_structure("1gk8gap", io.StringIO(text))


# ---------------------------------------------------------------------------
# Sample PDB: 1ctj chain A residues 1-3 (GLU ALA ASP), heavy atoms.
# GLU 1 and ALA 2 have whole-residue alternates A (occ 0.49) and B (0.51);
# C(1)--N(2) is bonded only within altloc A.  BioPython's PPBuilder would
# switch just those two atoms to altloc A, mixing conformers (ALA 2 tau
# 83.8 deg).  Fragment building must leave the default (B) selection
# alone: ALA 2 tau from conformer B is 106.719 deg.
# Source: 1ctj.pdb, from the top100 ersatz file (Reduce hydrogens removed)
# ---------------------------------------------------------------------------
ALTLOC_NTERM_PDB = """\
HEADER    1CTJ CHAIN A 1-3 WITH WHOLE-RESIDUE ALTLOCS
ATOM      1  N  AGLU A   1       4.127  26.179  -7.903  0.49 57.53           N
ATOM      2  N  BGLU A   1       3.535  25.488 -12.889  0.51 54.52           N
ATOM      3  CA AGLU A   1       5.490  26.607  -8.207  0.49 52.50           C
ATOM      4  CA BGLU A   1       2.754  26.395 -12.051  0.51 51.27           C
ATOM      5  C  AGLU A   1       5.550  27.734  -9.233  0.49 47.55           C
ATOM      6  C  BGLU A   1       3.471  27.729 -11.891  0.51 49.70           C
ATOM      7  O  AGLU A   1       4.982  28.809  -8.942  0.49 50.34           O
ATOM      8  O  BGLU A   1       2.870  28.781 -12.098  0.51 49.78           O
ATOM      9  CB AGLU A   1       6.372  25.388  -8.512  0.49 54.41           C
ATOM     10  CB BGLU A   1       2.416  25.770 -10.696  0.51 52.50           C
ATOM     11  CG AGLU A   1       7.853  25.736  -8.366  0.49 54.45           C
ATOM     12  CG BGLU A   1       1.268  26.491 -10.006  0.51 54.28           C
ATOM     13  CD AGLU A   1       8.764  24.531  -8.316  0.49 56.31           C
ATOM     14  CD BGLU A   1       1.720  27.296  -8.799  0.51 55.22           C
ATOM     15  OE1AGLU A   1       8.298  23.411  -8.022  0.49 60.72           O
ATOM     16  OE1BGLU A   1       2.727  28.028  -8.909  0.51 55.30           O
ATOM     17  OE2AGLU A   1       9.986  24.743  -8.482  0.49 59.28           O
ATOM     18  OE2BGLU A   1       1.068  27.166  -7.733  0.51 55.05           O
ATOM     19  N  AALA A   2       6.254  27.610 -10.341  0.49 39.82           N
ATOM     20  N  BALA A   2       4.759  27.656 -11.589  0.51 47.24           N
ATOM     21  CA AALA A   2       6.531  28.549 -11.391  0.49 40.63           C
ATOM     22  CA BALA A   2       5.670  28.786 -11.571  0.51 44.89           C
ATOM     23  C  AALA A   2       6.980  27.926 -12.710  0.49 38.03           C
ATOM     24  C  BALA A   2       5.766  29.304 -13.027  0.51 42.80           C
ATOM     25  O  AALA A   2       7.831  27.034 -12.721  0.49 32.72           O
ATOM     26  O  BALA A   2       5.165  30.351 -13.297  0.51 41.12           O
ATOM     27  CB AALA A   2       7.686  29.468 -10.955  0.49 40.12           C
ATOM     28  CB BALA A   2       7.039  28.354 -11.070  0.51 43.42           C
ATOM     29  N   ASP A   3       6.482  28.531 -13.792  1.00 36.71           N
ATOM     30  CA  ASP A   3       6.844  28.172 -15.120  1.00 32.51           C
ATOM     31  C   ASP A   3       8.056  28.945 -15.623  1.00 30.43           C
ATOM     32  O   ASP A   3       7.982  29.944 -16.316  1.00 26.86           O
ATOM     33  CB  ASP A   3       5.712  28.262 -16.130  1.00 31.89           C
ATOM     34  CG  ASP A   3       5.940  27.498 -17.420  1.00 34.79           C
ATOM     35  OD1 ASP A   3       7.066  27.105 -17.813  1.00 31.55           O
ATOM     36  OD2 ASP A   3       4.871  27.304 -18.059  1.00 41.40           O
END
"""


@pytest.fixture
def altloc_nterm_structure():
    """Return 1ctj A 1-3, whose first two residues have A/B alternates."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure("1ctj", io.StringIO(ALTLOC_NTERM_PDB))
