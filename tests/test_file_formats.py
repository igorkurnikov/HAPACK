"""
Tests for loading and saving molecular files in different formats.

Uses small test molecules from the examples directory.
Run with: pytest tests/test_file_formats.py -v
"""
import os
import sys
import pytest

import molset

# Paths to test data files (relative to repo root)
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_GLY_PDB = os.path.join(REPO_ROOT, "examples", "MD", "TEST_MD_1_GLY", "gly1.pdb")
DATA_TYR_PDB = os.path.join(REPO_ROOT, "examples", "VRML", "tyr.pdb")
DATA_CROWN_SDF = os.path.join(REPO_ROOT, "MORT_LIB", "python", "crown.sdf")
DATA_PEPT_XYZ = os.path.join(REPO_ROOT, "examples", "MD", "PEPT3", "pept_3.xyz")
DATA_LIG_HIN = os.path.join(REPO_ROOT, "examples", "MD", "EDIT_MUT_MAP", "1b.hin")
DATA_MULTI_XYZ = os.path.join(REPO_ROOT, "tests", "data", "multi_mol.xyz")
DATA_TYK2_SDF = os.path.join(REPO_ROOT, "tests", "data", "tyk2_ligands.sdf")


# ============================================================
#  Helper
# ============================================================

def new_molset():
    """Create a fresh MolSet instance."""
    return molset.MolSet()


# ============================================================
#  Loading tests
# ============================================================

class TestLoadPDB:
    """Test loading PDB files."""

    def test_load_gly_pdb(self):
        mset = new_molset()
        result = mset.LoadPDBFile(DATA_GLY_PDB)
        assert result, f"LoadPDBFile failed for {DATA_GLY_PDB}"
        assert mset.GetNAtoms() == 8
        assert mset.GetNMol() == 1

    def test_load_tyr_pdb(self):
        mset = new_molset()
        result = mset.LoadPDBFile(DATA_TYR_PDB)
        assert result, f"LoadPDBFile failed for {DATA_TYR_PDB}"
        assert mset.GetNAtoms() == 12
        assert mset.GetNMol() == 1


class TestLoadSDF:
    """Test loading SDF files."""

    def test_load_crown_sdf(self):
        """Single-molecule SDF file."""
        mset = new_molset()
        result = mset.LoadSDFFile(DATA_CROWN_SDF)
        assert result, f"LoadSDFFile failed for {DATA_CROWN_SDF}"
        assert mset.GetNAtoms() == 54
        assert mset.GetNMol() == 1
        assert mset.GetNBonds() > 0

    def test_load_tyk2_multi_mol_sdf(self):
        """Multi-molecule SDF file (10 TYK2 ligands, 366 atoms total)."""
        mset = new_molset()
        result = mset.LoadSDFFile(DATA_TYK2_SDF)
        assert result, f"LoadSDFFile failed for {DATA_TYK2_SDF}"
        assert mset.GetNMol() == 10
        assert mset.GetNAtoms() == 366
        assert mset.GetNBonds() > 0


class TestLoadHIN:
    """Test loading HIN (Arbalest) files."""

    def test_load_lig_hin(self):
        mset = new_molset()
        result = mset.LoadHINFile(DATA_LIG_HIN)
        assert result, f"LoadHINFile failed for {DATA_LIG_HIN}"
        assert mset.GetNAtoms() > 0
        assert mset.GetNMol() == 1


class TestLoadXYZ:
    """Test loading TINKER XYZ files."""

    def test_load_pept_xyz(self):
        mset = new_molset()
        result = mset.LoadXYZFile(DATA_PEPT_XYZ)
        assert result, f"LoadXYZFile failed for {DATA_PEPT_XYZ}"
        assert mset.GetNAtoms() == 43
        assert mset.GetNMol() == 1

    def test_load_multi_mol_xyz(self):
        mset = new_molset()
        result = mset.LoadXYZFile(DATA_MULTI_XYZ)
        assert result, f"LoadXYZFile failed for {DATA_MULTI_XYZ}"
        assert mset.GetNMol() == 3
        mol1 = mset.GetMolByIdx(0)
        mol2 = mset.GetMolByIdx(1)
        mol3 = mset.GetMolByIdx(2)
        assert mol1.GetNAtom() == 3   # water
        assert mol2.GetNAtom() == 5   # methane
        assert mol3.GetNAtom() == 3   # water


class TestFetchFile:
    """Test the FetchFile dispatch with explicit format codes."""

    def test_fetch_pdb(self):
        mset = new_molset()
        result = mset.FetchFile(molset.FormatPDB, DATA_GLY_PDB)
        assert result
        assert mset.GetNAtoms() == 8

    def test_fetch_sdf(self):
        mset = new_molset()
        result = mset.FetchFile(molset.FormatSDF, DATA_CROWN_SDF)
        assert result
        assert mset.GetNAtoms() == 54

    def test_fetch_hin(self):
        mset = new_molset()
        result = mset.FetchFile(molset.FormatHIN, DATA_LIG_HIN)
        assert result
        assert mset.GetNAtoms() > 0

    def test_fetch_xyz(self):
        mset = new_molset()
        result = mset.FetchFile(molset.FormatXYZ, DATA_PEPT_XYZ)
        assert result
        assert mset.GetNAtoms() == 43


# ============================================================
#  Save and reload (round-trip) tests
# ============================================================

class TestSaveLoadPDB:
    """Test saving to PDB and reloading."""

    def test_roundtrip_pdb(self, tmp_path):
        mset = new_molset()
        mset.LoadPDBFile(DATA_GLY_PDB)
        natoms_orig = mset.GetNAtoms()

        out_file = str(tmp_path / "out.pdb")
        mset.SavePDBFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        result = mset2.LoadPDBFile(out_file)
        assert result
        assert mset2.GetNAtoms() == natoms_orig


class TestSaveLoadHIN:
    """Test saving to HIN and reloading."""

    def test_roundtrip_hin_from_pdb(self, tmp_path):
        mset = new_molset()
        mset.LoadPDBFile(DATA_GLY_PDB)
        natoms_orig = mset.GetNAtoms()

        out_file = str(tmp_path / "out.hin")
        mset.SaveHINFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        result = mset2.LoadHINFile(out_file)
        assert result
        assert mset2.GetNAtoms() == natoms_orig

    def test_roundtrip_hin(self, tmp_path):
        mset = new_molset()
        mset.LoadHINFile(DATA_LIG_HIN)
        natoms_orig = mset.GetNAtoms()
        nbonds_orig = mset.GetNBonds()

        out_file = str(tmp_path / "out.hin")
        mset.SaveHINFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        result = mset2.LoadHINFile(out_file)
        assert result
        assert mset2.GetNAtoms() == natoms_orig
        assert mset2.GetNBonds() == nbonds_orig


class TestSaveFile:
    """Test the SaveFile method that auto-detects format from extension."""

    def test_savefile_pdb(self, tmp_path):
        mset = new_molset()
        mset.LoadPDBFile(DATA_GLY_PDB)

        out_file = str(tmp_path / "out.pdb")
        mset.SaveFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        assert mset2.LoadPDBFile(out_file)
        assert mset2.GetNAtoms() == mset.GetNAtoms()

    def test_savefile_hin(self, tmp_path):
        mset = new_molset()
        mset.LoadPDBFile(DATA_GLY_PDB)

        out_file = str(tmp_path / "out.hin")
        mset.SaveFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        assert mset2.LoadHINFile(out_file)
        assert mset2.GetNAtoms() == mset.GetNAtoms()

    def test_savefile_sdf(self, tmp_path):
        mset = new_molset()
        mset.LoadSDFFile(DATA_CROWN_SDF)

        out_file = str(tmp_path / "out.sdf")
        mset.SaveFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        assert mset2.LoadSDFFile(out_file)
        assert mset2.GetNAtoms() == mset.GetNAtoms()


class TestSaveLoadSDF:
    """Test saving to SDF and reloading."""

    def test_roundtrip_sdf(self, tmp_path):
        """Load SDF, save SDF, reload and check atoms/bonds."""
        mset = new_molset()
        mset.LoadSDFFile(DATA_CROWN_SDF)
        natoms_orig = mset.GetNAtoms()
        nbonds_orig = mset.GetNBonds()

        out_file = str(tmp_path / "out.sdf")
        mset.SaveSDFFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        result = mset2.LoadSDFFile(out_file)
        assert result
        assert mset2.GetNAtoms() == natoms_orig
        assert mset2.GetNBonds() == nbonds_orig

    def test_roundtrip_multi_mol_sdf(self, tmp_path):
        """Multi-molecule SDF round trip."""
        mset = new_molset()
        mset.LoadSDFFile(DATA_TYK2_SDF)
        natoms_orig = mset.GetNAtoms()
        nmol_orig = mset.GetNMol()

        out_file = str(tmp_path / "out.sdf")
        mset.SaveSDFFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        result = mset2.LoadSDFFile(out_file)
        assert result
        assert mset2.GetNMol() == nmol_orig
        assert mset2.GetNAtoms() == natoms_orig


# ============================================================
#  Cross-format conversion tests
# ============================================================

class TestCrossFormat:
    """Test loading in one format and saving in another."""

    def test_pdb_to_hin_to_pdb(self, tmp_path):
        """PDB -> HIN -> PDB round trip."""
        mset1 = new_molset()
        mset1.LoadPDBFile(DATA_GLY_PDB)
        natoms = mset1.GetNAtoms()

        hin_file = str(tmp_path / "step1.hin")
        mset1.SaveHINFile(hin_file)

        mset2 = new_molset()
        mset2.LoadHINFile(hin_file)
        assert mset2.GetNAtoms() == natoms

        pdb_file = str(tmp_path / "step2.pdb")
        mset2.SavePDBFile(pdb_file)

        mset3 = new_molset()
        mset3.LoadPDBFile(pdb_file)
        assert mset3.GetNAtoms() == natoms

    def test_sdf_to_pdb(self, tmp_path):
        """SDF -> PDB save."""
        mset = new_molset()
        mset.LoadSDFFile(DATA_CROWN_SDF)
        natoms = mset.GetNAtoms()

        out_file = str(tmp_path / "crown.pdb")
        mset.SavePDBFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        mset2.LoadPDBFile(out_file)
        assert mset2.GetNAtoms() == natoms

    def test_sdf_to_hin(self, tmp_path):
        """SDF -> HIN save."""
        mset = new_molset()
        mset.LoadSDFFile(DATA_CROWN_SDF)
        natoms = mset.GetNAtoms()

        out_file = str(tmp_path / "crown.hin")
        mset.SaveHINFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        mset2.LoadHINFile(out_file)
        assert mset2.GetNAtoms() == natoms

    def test_pdb_to_sdf(self, tmp_path):
        """PDB -> SDF save."""
        mset = new_molset()
        mset.LoadPDBFile(DATA_GLY_PDB)
        natoms = mset.GetNAtoms()

        out_file = str(tmp_path / "gly.sdf")
        mset.SaveSDFFile(out_file)
        assert os.path.isfile(out_file)

        mset2 = new_molset()
        mset2.LoadSDFFile(out_file)
        assert mset2.GetNAtoms() == natoms


# ============================================================
#  Error handling tests
# ============================================================

class TestErrorHandling:
    """Test that loading nonexistent files returns failure."""

    def test_load_nonexistent_pdb(self):
        mset = new_molset()
        result = mset.LoadPDBFile("nonexistent_file.pdb")
        assert not result

    def test_load_nonexistent_sdf(self):
        mset = new_molset()
        result = mset.LoadSDFFile("nonexistent_file.sdf")
        assert not result

    def test_load_nonexistent_hin(self):
        mset = new_molset()
        result = mset.LoadHINFile("nonexistent_file.hin")
        assert not result

    def test_load_nonexistent_xyz(self):
        mset = new_molset()
        result = mset.LoadXYZFile("nonexistent_file.xyz")
        assert not result
