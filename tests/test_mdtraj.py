"""
Smoke tests: verify that molset and its optional scientific dependencies
can be imported.

Run with: pytest tests/test_mdtraj.py -v
"""
import pytest


def test_import_molset():
    import molset  # noqa: F401


def test_import_mdtraj():
    mdtraj = pytest.importorskip("mdtraj")


def test_import_rdkit():
    rdkit = pytest.importorskip("rdkit")


def test_import_openbabel():
    openbabel = pytest.importorskip("openbabel")


def test_import_pybel():
    openbabel = pytest.importorskip("openbabel")
    from openbabel import pybel  # noqa: F401
