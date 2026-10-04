"""
Shared helpers for the unit tests.

analyse() mirrors the IC generation in nomodeco.main() (nomodeco.py, "IC Generation" section)
for single molecules. It is replaced by pipeline.generate() in step 1.3 of the restructure plan.
"""

import contextlib
import io
import logging
from types import SimpleNamespace

import numpy as np
import pymatgen.core as mg
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

from nomodeco.libraries import icsel, specifications
from nomodeco.nomodeco import strip_numbers


def null_logger(name="nomodeco-tests"):
    log = logging.getLogger(name)
    log.addHandler(logging.NullHandler())
    log.propagate = False
    return log


def coordinates(mol) -> np.ndarray:
    return np.array([a.coordinates for a in mol], dtype=float)


def analyse(mol) -> SimpleNamespace:
    degofc = mol.degree_of_covalance()
    cov_bonds = sorted(mol.covalent_bonds(degofc))
    _, _, sub_symbols = mol.detect_submolecules(degofc)
    h_bonds = sorted(mol.intermolecular_h_bond(degofc, sub_symbols))
    bonds = sorted(set(cov_bonds) | set(h_bonds))

    angles, linear_angles = mol.generate_angles(bonds)
    dihedrals = mol.generate_dihedrals(bonds)

    pg = PointGroupAnalyzer(
        mg.Molecule([strip_numbers(a.symbol) for a in mol], [a.coordinates for a in mol])
    )
    spec = specifications.calculation_specification({}, mol, pg, bonds, angles, linear_angles)

    if spec["planar"] == "yes":
        out_of_plane = mol.generate_out_of_plane(bonds)
    elif spec["planar submolecule(s)"]:
        out_of_plane = mol.generate_oop_planar_subunits(bonds, spec["planar submolecule(s)"])
    else:
        out_of_plane = []

    idof = mol.idof_linear() if spec["linearity"] == "fully linear" else mol.idof_general()

    return SimpleNamespace(
        cov_bonds=cov_bonds, h_bonds=h_bonds, bonds=bonds, angles=angles,
        linear_angles=linear_angles, dihedrals=dihedrals, out_of_plane=out_of_plane,
        spec=spec, idof=idof, point_group=pg.sch_symbol,
    )


def generate_sets(mol, a=None):
    """Run icsel.get_sets for mol (stdout silenced). Returns (ic_dict, analysis)."""
    a = a or analyse(mol)
    with contextlib.redirect_stdout(io.StringIO()):
        ic_dict = icsel.get_sets(
            a.idof, null_logger(), mol, a.bonds, a.angles, a.linear_angles,
            a.out_of_plane, a.dihedrals, a.spec,
        )
    return ic_dict, a


def n_ics(ic_set) -> int:
    return sum(len(v) for v in ic_set.values())
