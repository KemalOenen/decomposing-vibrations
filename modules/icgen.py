from __future__ import annotations

import itertools
from typing import NamedTuple
import pprint
import numpy as np
from scipy import constants
import argparse

from . import bond
from . import angle
from . import oop
from . import dihedral
from . import molpro_parser

class Atom(NamedTuple):
    symbol: str
    coordinates: tuple

#TODO: adaption for other output formats like Gaussian, Orca
def parse_xyz_from_outfile(outfile):
    bohr = constants.value(u'Bohr radius')
    angstrom = constants.value(u'Angstrom star')
    BOHR_PER_ANGSTROM = angstrom/bohr 
    for line in outfile:
        if line.strip().startswith('ATOMIC COORDINATES'):
            break
    next(outfile)
    next(outfile)
    next(outfile)
    names = []
    for line in outfile:
        entries = line.strip().split()
        if len(entries) == 0:
            break
        names.append(entries[1])
    for line in outfile:
        if line.strip().startswith('FREQUENCIES * CALCULATION OF NORMAL MODES'):
            break
    for _ in range(6):
        next(outfile)
    coordinates = []
    for line in outfile:
        entries = line.strip().split()
        if len(entries) == 0:
            break
        xyz = [float(f) / BOHR_PER_ANGSTROM for f in entries[3:]]
        coordinates.append(xyz)

    return [Atom(name, tuple(coordinate)) for name, coordinate in zip(names, coordinates)]

def initialize_bonds(atoms):
    hits = []
    for atom_a, atom_b in itertools.combinations(atoms, 2):
        if bond.is_valid(atom_a, atom_b):
            hits.append((atom_a.symbol, atom_b.symbol))
    return hits

def initialize_angles(atoms):
    hits = []
    linear_hits = []
    for a,b,c in itertools.permutations(atoms, 3):
        if angle.is_valid(a, b, c) and not (c.symbol, b.symbol, a.symbol) in hits and not (c.symbol, b.symbol, a.symbol) in linear_hits:
            if np.isclose(angle.bond_angle(a.coordinates, b.coordinates, c.coordinates) / np.pi, 0.99, atol=0.15):
                linear_hits.append((a.symbol, b.symbol, c.symbol))
                linear_hits.append((a.symbol, b.symbol, c.symbol))
            else:
                hits.append((a.symbol, b.symbol, c.symbol))
    # we need both DOF!
    #if len(hits) != 0:
    #    linear_hits = [*set(linear_hits)]
    return hits, linear_hits

def initialize_oop(atoms):
    hits = []
    for a,b,c,d in itertools.permutations(atoms, 4):
        if oop.is_valid(a, b, c, d) and not (a.symbol, b.symbol, d.symbol, c.symbol) in hits:
            hits.append((a.symbol, b.symbol, c.symbol, d.symbol))
    return hits

# if the molecule is general but contains a planar submolecule, only sample these coordinates
def initialize_oop_planar_subunits(atoms, central_atoms_list):
    hits = []
    for central_atom in central_atoms_list:
        for a,b,c,d in itertools.permutations(atoms, 4):
            if a.symbol == central_atom[0] and oop.is_valid(a, b, c, d) and not (a.symbol, b.symbol, d.symbol, c.symbol) in hits:
                hits.append((a.symbol, b.symbol, c.symbol, d.symbol))
    return hits

def initialize_dihedrals(atoms):
    hits = []
    for a,b,c,d in itertools.permutations(atoms, 4):
        if dihedral.is_valid(a, b, c, d) and not (d.symbol, c.symbol, b.symbol, a.symbol) in hits:
            hits.append((a.symbol, b.symbol, c.symbol, d.symbol))
    return hits

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    with open(args.output) as outfile:
        xyz= molpro_parser.parse_xyz_from_inputfile(outfile)
        n_atoms = len(xyz)
    bonds = initialize_bonds(xyz)
    print("The following", len(bonds), "bonds have been found:")
    pprint.pprint(bonds)
    angles = initialize_angles(xyz)
    print("The following", len(angles), "angles have been found:") 
    pprint.pprint(angles)
#    out_of_plane = initialize_oop(xyz)
#    print("The following possible", len(out_of_plane), "out of plane angle definitions have been found:") 
#    pprint.pprint(out_of_plane)
    dihedrals = initialize_dihedrals(xyz) 
    print("The following", len(dihedrals), "dihedrals have been found:") 
    pprint.pprint(dihedrals)


if __name__ == '__main__':
    main()
