"""
Tests for Z-matrix (internal coordinates) functionality.

Covers InitStdZMat, InitAllXYZ, LoadFromString, SaveToString,
and coordinate manipulation.

Run with: pytest tests/test_zmat.py -v
"""
import os
import pytest

import molset

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_HOH2_HLM = os.path.join(REPO_ROOT, "tests", "data", "HOH2.hlm")
DATA_MG_HOH7_HLM = os.path.join(REPO_ROOT, "tests", "data", "MG_HOH7.hlm")

# Z-matrix string with numeric values for Mg(H2O)7 (22 atoms)
ZMAT_MG_HOH7_NUMERIC = """\
Mg
O     1     2.11247309
H     2     0.95336722        1     126.4789322
H     2      0.95341086       3     107.4375181        1     -176.1937715
O     1     2.11993303        2     86.53210581        3     -175.2259583
H     5     0.95291499        1     126.0634331        2     90.3834173
H     5     0.95331212        6     106.9683016       1     -167.3328815
O     1      2.09078869       2     93.03391603        3     -85.98541282
H     8     0.95248506        1     128.8425032        2     3.92765029
H     8      0.96159493       9     108.2715025        1     176.3897237
O     1     2.09467801        5     94.16458312       6     -89.65408357
H     11     0.95212312       1     128.412064        5      -134.3273927
H     11    0.96103362       12     107.0558795       1     157.0655298
O     1      2.12255818       2    85.86304319        3     5.80606768
H     14     0.95345963       1     124.8395297        2     -103.4893761
H     14     0.95341795      15     106.9748438        1     168.8434896
O     1      2.2              2     89.62591725        3     95.78766795
H     17     0.95391008        1     123.9207953        2    -2.2583911
H     17     0.95369514        18    106.8334524        1     -174.5228904
O     1      3.8               2     135.9682558       3     -95.28277733
H     20      0.9538636        1     -95.28277733        2    99.61542168
H     20     0.953690541       21    105.2902395        1    -178.8382717
"""

# Z-matrix string with symbolic variable names for Mg(H2O)7
ZMAT_MG_HOH7_SYMBOLIC = """\
Mg
O     1     rMO1
H     2     0.95          1     aMOH1
H     2     rOH12        3     aHOH1        1     dMHOH1
O     1     rMO2         2     aOMO2        3     dOMOH2
H     5     rOH21        1     aMOH2        2     dHOMO2
H     5     rOH22        6     aHOH2        1     dMHOH2
O     1     rMO3         2     aOMO3        3     dOMOH3
H     8     rOH31        1     aMOH3        2     dHOMO3
H     8     rOH32        9     aHOH3        1     dMHOH3
O     1     rMO4         5     aOMO4        6     dOMOH4
H     11     rOH41        1     aMOH4        5     dHOMO4
H     11     rOH42        12     aHOH4        1     dMHOH4
O     1     rMO5         2     aOMO5        3     dOMOH5
H     14     rOH51        1     aMOH5        2     dHOMO5
H     14     rOH52        15     aHOH5        1     dMHOH5
O     1     rMO6         2     aOMO6        3     dOMOH6
H     17     rOH61        1     aMOH6        2     dHOMO6
H     17     rOH62        18     aHOH6        1     dMHOH6
O     1     rMO7     2     aOMO7        3     dOMOH7
H     20     rOH71        1     aMOH7        2     dHOMO7
H     20     rOH72        21     aHOH7        1     dMHOH7

rMO1   2.11247309
rOH12    0.95341086
aMOH1   126.4789322
aHOH1   107.4375181
dMHOH1    -176.1937715
rMO2    2.11993303
rOH21    0.95291499
rOH22    0.95331212
aMOH2    126.0634331
aHOH2    106.9683016
dMHOH2   -167.3328815
aOMO2    86.53210581
dOMOH2   -175.2259583
dHOMO2    90.3834173
rMO3    2.09078869
rOH31    0.95248506
rOH32    0.96159493
aMOH3    128.8425032
aHOH3    108.2715025
dMHOH3   176.3897237
aOMO3    93.03391603
dOMOH3   -85.98541282
dHOMO3   3.92765029
rMO4   2.09467801
rOH41    0.95212312
rOH42    0.96103362
aMOH4   128.412064
aHOH4   107.0558795
dMHOH4   157.0655298
aOMO4    94.16458312
dOMOH4   -89.65408357
dHOMO4   -134.3273927
rMO5   2.12255818
rOH51    0.95345963
rOH52    0.95341795
aMOH5   124.8395297
aHOH5   106.9748438
dMHOH5   168.8434896
aOMO5    85.86304319
dOMOH5   5.80606768
dHOMO5   -103.4893761
rOH61   0.95391008
rOH62    0.95369514
aMOH6   123.9207953
aHOH6    106.8334524
dMHOH6   -174.5228904
aOMO6    89.62591725
dOMOH6   95.78766795
dHOMO6   -2.2583911
rOH71    0.9538636
rOH72    0.95369054
aMOH7    129.5063823
aHOH7    105.2902395
dMHOH7    -178.8382717
aOMO7    135.9682558
dOMOH7    -95.28277733
dHOMO7    99.61542168
rMO7   3.8

rMO6    2.2
"""


def new_molset_from_hlm(path):
    """Load an HLM file into a fresh MolSet and return it."""
    ms = molset.MolSet()
    result = ms.FetchFile(molset.FormatHarlem, path)
    assert result, f"FetchFile failed for {path}"
    return ms


class TestInitStdZMat:
    """Test InitStdZMat – build a standard Z-matrix from Cartesian coords."""

    def test_init_std_zmat_produces_output(self):
        """InitStdZMat should produce a non-empty Z-matrix string."""
        ms = new_molset_from_hlm(DATA_HOH2_HLM)
        zm = ms.GetZMat()
        zm.InitStdZMat()

        text = zm.SaveToString()
        assert text is not None
        assert len(text.strip()) > 0

    def test_init_std_zmat_line_count(self):
        """SaveToString should have one line per atom (6 for HOH2)."""
        ms = new_molset_from_hlm(DATA_HOH2_HLM)
        zm = ms.GetZMat()
        zm.InitStdZMat()

        lines = [l for l in zm.SaveToString().splitlines() if l.strip()]
        assert len(lines) == ms.GetNAtoms()

    def test_init_std_zmat_preserves_nz(self):
        """GetNZ should equal NAtoms after InitStdZMat."""
        ms = new_molset_from_hlm(DATA_HOH2_HLM)
        zm = ms.GetZMat()
        zm.InitStdZMat()
        assert zm.GetNZ() == ms.GetNAtoms()

    def test_init_std_zmat_ncrd(self):
        """GetNCrd should be 3*N - 6 for a non-linear molecule (N>=3)."""
        ms = new_molset_from_hlm(DATA_HOH2_HLM)
        zm = ms.GetZMat()
        zm.InitStdZMat()
        n = ms.GetNAtoms()
        # Internal coords: 3N - 6 for non-linear molecules
        assert zm.GetNCrd() == 3 * n - 6


class TestInitAllXYZ:
    """Test InitAllXYZ – convert Z-matrix to Cartesian coordinates."""

    def test_init_all_xyz_produces_output(self):
        """InitAllXYZ on a pre-existing zmat should produce XYZ output."""
        ms = new_molset_from_hlm(DATA_HOH2_HLM)
        zm = ms.GetZMat()
        zm.InitAllXYZ()

        text = zm.SaveToString()
        assert text is not None
        assert len(text.strip()) > 0

    def test_init_all_xyz_line_count(self):
        """After InitAllXYZ, SaveToString should have one line per atom."""
        ms = new_molset_from_hlm(DATA_HOH2_HLM)
        zm = ms.GetZMat()
        zm.InitAllXYZ()

        lines = [l for l in zm.SaveToString().splitlines() if l.strip()]
        assert len(lines) == ms.GetNAtoms()

    def test_init_all_xyz_ncrd(self):
        """After InitAllXYZ, NCrd should be 3*N (all Cartesian coords)."""
        ms = new_molset_from_hlm(DATA_HOH2_HLM)
        zm = ms.GetZMat()
        zm.InitAllXYZ()
        assert zm.GetNCrd() == 3 * ms.GetNAtoms()


class TestLoadFromStringNumeric:
    """Test LoadFromString with numeric Z-matrix values (Mg(H2O)7)."""

    def test_load_numeric_zmat_nz(self):
        """After loading numeric zmat, NZ should match atom count."""
        ms = new_molset_from_hlm(DATA_MG_HOH7_HLM)
        zm = ms.GetZMat()
        zm.LoadFromString(ZMAT_MG_HOH7_NUMERIC)
        assert zm.GetNZ() == 22

    def test_load_numeric_zmat_roundtrip(self):
        """SaveToString after LoadFromString should produce 22 non-empty lines."""
        ms = new_molset_from_hlm(DATA_MG_HOH7_HLM)
        zm = ms.GetZMat()
        zm.LoadFromString(ZMAT_MG_HOH7_NUMERIC)

        lines = [l for l in zm.SaveToString().splitlines() if l.strip()]
        assert len(lines) == 22

    def test_load_numeric_zmat_ncrd(self):
        """NCrd should be 3*22 - 6 = 60 internal coordinates."""
        ms = new_molset_from_hlm(DATA_MG_HOH7_HLM)
        zm = ms.GetZMat()
        zm.LoadFromString(ZMAT_MG_HOH7_NUMERIC)
        assert zm.GetNCrd() == 60


class TestLoadFromStringSymbolic:
    """Test LoadFromString with symbolic variables and coordinate manipulation."""

    def test_load_symbolic_zmat(self):
        """Symbolic zmat with variable definitions should load correctly."""
        ms = new_molset_from_hlm(DATA_MG_HOH7_HLM)
        zm = ms.GetZMat()
        zm.LoadFromString(ZMAT_MG_HOH7_SYMBOLIC)
        assert zm.GetNZ() == 22
        assert zm.GetNCrd() == 60

    def test_get_r_crd_for_atom(self):
        """GetRCrdForAtom should return the R coordinate for a named atom."""
        ms = new_molset_from_hlm(DATA_MG_HOH7_HLM)
        zm = ms.GetZMat()
        zm.LoadFromString(ZMAT_MG_HOH7_SYMBOLIC)

        aptr = ms.GetAtomByRef("HOH2.O")
        elc = zm.GetRCrdForAtom(aptr)
        assert abs(elc.GetValue() - 2.11247309) < 1e-5

    def test_set_value(self):
        """SetValue should update the coordinate value."""
        ms = new_molset_from_hlm(DATA_MG_HOH7_HLM)
        zm = ms.GetZMat()
        zm.LoadFromString(ZMAT_MG_HOH7_SYMBOLIC)

        aptr = ms.GetAtomByRef("HOH2.O")
        elc = zm.GetRCrdForAtom(aptr)
        elc.SetValue(2.5)
        assert abs(elc.GetValue() - 2.5) < 1e-10

    def test_set_frozen(self):
        """SetFrozen should freeze the coordinate."""
        ms = new_molset_from_hlm(DATA_MG_HOH7_HLM)
        zm = ms.GetZMat()
        zm.LoadFromString(ZMAT_MG_HOH7_SYMBOLIC)

        aptr = ms.GetAtomByRef("HOH3.O")
        elc = zm.GetRCrdForAtom(aptr)
        assert not elc.IsFrozen()
        elc.SetFrozen()
        assert elc.IsFrozen()

    def test_modified_values_in_output(self):
        """Modified coordinate values should appear in SaveToString output."""
        ms = new_molset_from_hlm(DATA_MG_HOH7_HLM)
        zm = ms.GetZMat()
        zm.LoadFromString(ZMAT_MG_HOH7_SYMBOLIC)

        aptr = ms.GetAtomByRef("HOH2.O")
        elc = zm.GetRCrdForAtom(aptr)
        elc.SetFrozen()
        elc.SetValue(2.5)

        aptr3 = ms.GetAtomByRef("HOH3.O")
        elc3 = zm.GetRCrdForAtom(aptr3)
        elc3.SetFrozen()
        elc3.SetValue(3.0)

        output = zm.SaveToString()
        assert "2.5" in output
        assert "3.0" in output
