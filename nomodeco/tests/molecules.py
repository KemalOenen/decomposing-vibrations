"""
Small test molecules built in code (no file I/O), one per branch of the topology decision tree.

Geometries are approximate equilibrium structures in Angstrom; they only need to be good enough
for bond detection, planarity/linearity classification and point-group symmetry.
Atom labels follow the parsers' convention (orca_parser.numerate_strings): an element is numbered
only if it occurs more than once, e.g. CO2 -> C, O1, O2.
"""

import numpy as np

from nomodeco.libraries.molecule_class import Molecule
from nomodeco.libraries.orca_parser import numerate_strings


def build(symbols, coords) -> Molecule:
    labels = numerate_strings(symbols)
    return Molecule([Molecule.Atom(l, tuple(float(x) for x in c)) for l, c in zip(labels, coords)])


def h2o():
    """bent, planar, acyclic, no dihedrals"""
    return build(["O", "H", "H"], [(0, 0, 0.1173), (0, 0.7572, -0.4692), (0, -0.7572, -0.4692)])


def co2():
    """fully linear"""
    return build(["C", "O", "O"], [(0, 0, 0), (0, 0, 1.16), (0, 0, -1.16)])


def nh3():
    """general (non-planar), acyclic"""
    return build(
        ["N", "H", "H", "H"],
        [(0, 0, 0.1162), (0, 0.9377, -0.2711), (0.8121, -0.4689, -0.2711), (-0.8121, -0.4689, -0.2711)],
    )


def ethylene():
    """planar, acyclic, has dihedrals and out-of-plane coordinates"""
    return build(
        ["C", "C", "H", "H", "H", "H"],
        [(0, 0, 0.6695), (0, 0, -0.6695),
         (0, 0.9289, 1.2321), (0, -0.9289, 1.2321), (0, 0.9289, -1.2321), (0, -0.9289, -1.2321)],
    )


def propyne():
    """general, acyclic, with a linear submolecule (C-C#C-H)"""
    methyl_h = [(1.03 * np.cos(a), 1.03 * np.sin(a), -0.39) for a in np.radians([0, 120, 240])]
    return build(
        ["C", "C", "C", "H", "H", "H", "H"],
        [(0, 0, 0), (0, 0, 1.46), (0, 0, 2.67), (0, 0, 3.73), *methyl_h],
    )


def hcocn():
    """planar, acyclic, with a linear submolecule (C-C#N) ending at an sp2 center that needs an oop"""
    return build(
        ["C", "O", "H", "C", "N"],
        [(0, 0, 0), (-0.60, 1.04, 0), (-0.55, -0.95, 0), (1.47, 0, 0), (2.63, 0, 0)],
    )


def acetyl_cyanide():
    """general (methyl), acyclic, with a linear submolecule (C-C#N) ending at a planar sp2 center"""
    methyl_h = [(-0.25, -2.25, 0.0), (-1.38, -1.20, 0.89), (-1.38, -1.20, -0.89)]
    return build(
        ["C", "O", "C", "C", "N", "H", "H", "H"],
        [(0, 0, 0), (-0.60, 1.04, 0), (-0.75, -1.30, 0), (1.47, 0, 0), (2.63, 0, 0), *methyl_h],
    )


def benzene():
    """planar, cyclic (mu = 1)"""
    angles = np.radians(np.arange(0, 360, 60))
    carbons = [(1.39 * np.cos(a), 1.39 * np.sin(a), 0) for a in angles]
    hydrogens = [(2.48 * np.cos(a), 2.48 * np.sin(a), 0) for a in angles]
    return build(["C"] * 6 + ["H"] * 6, carbons + hydrogens)


def water_dimer():
    """intermolecular: two covalent submolecules joined by one hydrogen bond"""
    return build(
        ["O", "H", "H", "O", "H", "H"],
        [(-1.5510, -0.1145, 0.0), (-1.9343, 0.7625, 0.0), (-0.5997, 0.0407, 0.0),
         (1.3505, 0.1116, 0.0), (1.6803, -0.3733, -0.7615), (1.6803, -0.3733, 0.7615)],
    )


def hcn_h2o():
    """intermolecular, planar: N#C-H...OH2, the acceptor O (oop center) ends the linear unit C-H...O"""
    return build(
        ["N", "C", "H", "O", "H", "H"],
        [(-3.27, 0, 0), (-2.12, 0, 0), (-1.05, 0, 0), (1.10, 0, 0), (1.69, 0.76, 0), (1.69, -0.76, 0)],
    )


def nh3_dimer():
    """intermolecular, cyclic (C2h): two bent N-H...N contacts at 2.63 A, degree of covalence 0.21"""
    return build(
        ["N", "H", "H", "H", "N", "H", "H", "H"],
        [(0.0472, -1.6438, 0.0), (-0.1094, -2.2155, -0.8067), (-0.1094, -2.2155, 0.8067),
         (-0.6459, -0.9203, 0.0), (-0.0472, 1.6438, 0.0), (0.1094, 2.2155, -0.8067),
         (0.6459, 0.9203, 0.0), (0.1094, 2.2155, 0.8067)],
    )


def cyclopropanol_water():
    """intermolecular, cyclic (cyclopropane ring) with a linear H-bond unit O1-H6...O2 (175.1 deg)"""
    return build(
        ["C", "C", "C", "H", "H", "H", "H", "H", "O", "H", "O", "H", "H"],
        [(-0.0570, -0.2206, -2.1091), (-0.1589, 1.0528, -1.3026), (0.4609, -0.1770, -0.7189),
         (-0.9631, -0.7787, -2.2607), (0.6441, -0.2579, -2.9242), (-1.1341, 1.3247, -0.9394),
         (0.4764, 1.8802, -1.5670), (1.5318, -0.1758, -0.5864), (-0.2388, -0.9014, 0.2229),
         (-0.1634, -0.4719, 1.0616), (0.0504, 0.2867, 2.8999), (-0.6586, 0.7673, 3.2910),
         (0.3380, -0.3339, 3.5472)],
    )


ALL = {
    "h2o": h2o,
    "co2": co2,
    "nh3": nh3,
    "ethylene": ethylene,
    "propyne": propyne,
    "hcocn": hcocn,
    "acetyl_cyanide": acetyl_cyanide,
    "benzene": benzene,
    "water_dimer": water_dimer,
    "hcn_h2o": hcn_h2o,
    "nh3_dimer": nh3_dimer,
    "cyclopropanol_water": cyclopropanol_water,
}