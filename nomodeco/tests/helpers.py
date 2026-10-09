"""
Shared helpers for the unit tests.

analyse() mirrors the IC generation in nomodeco.main() (nomodeco.py, "IC Generation" section)
for single molecules. It is replaced by pipeline.generate() in step 1.3 of the restructure plan.
"""

import contextlib
import io
import logging
import sys
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pymatgen.core as mg
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

from nomodeco.libraries import icsel, specifications
from nomodeco.libraries import topology as tp
from nomodeco.libraries.ic_class import InternalCoordinates
from nomodeco.nomodeco import strip_numbers


def null_logger(name="nomodeco-tests"):
    log = logging.getLogger(name)
    log.addHandler(logging.NullHandler())
    log.propagate = False
    return log


def coordinates(mol) -> np.ndarray:
    return np.array([a.coordinates for a in mol], dtype=float)


def analyse(mol) -> SimpleNamespace:
    G = mol.graph()
    cov_bonds = sorted(G.graph["cov_bonds"])
    h_bonds = sorted(G.graph["h_bonds"])
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


def total_ic_dict(mol, a) -> InternalCoordinates:
    """The cov_/h_bond_/acc_don_ IC dict that main() builds for the intermolecular topology variants."""
    acc_don = list(mol.graph().graph["acc_don_bonds"])
    bond_acc_don = sorted(set(a.cov_bonds) | set(acc_don))

    d = InternalCoordinates()
    d.add_coordinate("cov_bond", a.cov_bonds)
    d.add_coordinate("h_bond", a.h_bonds)
    d.add_coordinate("acc_don", acc_don)
    cov_angles, cov_linear_angles = mol.generate_angles(a.cov_bonds)
    d.add_coordinate("cov_angles", cov_angles)
    d.add_coordinate("cov_linear_angles", cov_linear_angles)
    d.add_coordinate("cov_dihedrals", mol.generate_dihedrals(a.cov_bonds))
    for prefix, bonds in (("h_bond", a.bonds), ("acc_don", bond_acc_don)):
        angles, linear_angles = mol.generate_angles(bonds)
        d.add_coord_diff(f"{prefix}_angles", angles, d["cov_angles"])
        d.add_coord_diff_linear(f"{prefix}_linear_angles", linear_angles, d["cov_linear_angles"])
        d.add_coord_diff(f"{prefix}_dihedrals", mol.generate_dihedrals(bonds), d["cov_dihedrals"])

    for key in ("cov_oop", "h_bond_oop", "acc_don_oop"):
        d.add_coordinate(key, [])
    if a.spec["planar"] == "yes":
        d.add_coordinate("cov_oop", mol.generate_out_of_plane(a.cov_bonds))
        d.add_coord_diff("h_bond_oop", mol.generate_out_of_plane(a.bonds), d["cov_oop"])
        d.add_coord_diff("acc_don_oop", mol.generate_out_of_plane(bond_acc_don), d["cov_oop"])
    return d


def generate_sets(mol, a=None):
    """Run icsel.get_sets for mol (stdout silenced). Returns (ic_dict, analysis)."""
    a = a or analyse(mol)
    # module globals read by the intermolecular variants; passed explicitly after step 1.3
    tp.atoms_list = mol
    tp.Total_IC_dict = total_ic_dict(mol, a)
    # the intermolecular variants call arguments.get_args() for --comb: give them the CLI
    # defaults, not pytest's argv (until step 0.3 passes comb as a parameter)
    with contextlib.redirect_stdout(io.StringIO()), mock.patch.object(sys, "argv", ["nomodeco"]):
        ic_dict = icsel.get_sets(
            a.idof, null_logger(), mol, a.bonds, a.angles, a.linear_angles,
            a.out_of_plane, a.dihedrals, a.spec,
        )
    return ic_dict, a


def n_ics(ic_set) -> int:
    return sum(len(v) for v in ic_set.values())
