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


ALL = {
    "h2o": h2o,
    "co2": co2,
    "nh3": nh3,
    "ethylene": ethylene,
    "propyne": propyne,
    "benzene": benzene,
    "water_dimer": water_dimer,
}