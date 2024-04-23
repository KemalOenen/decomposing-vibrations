from __future__ import annotations

from itertools import combinations
import bond
import numpy as np

#TODO: clean bond_angle function
def bond_angle(coord_a: Iterable, coord_b: Iterable, coord_c: Iterable) -> float: #RAD
    coord_a = np.array(coord_a)
    coord_b = np.array(coord_b)
    coord_c = np.array(coord_c)
    cosine_angle = np.clip((np.inner((coord_a - coord_b), (coord_c - coord_b))) / (
            bond.actual_length(coord_a, coord_b) * (bond.actual_length(coord_c, coord_b))), -1.0, 1.0)
    return np.arccos(cosine_angle)

def is_valid(atom_a: Atom, atom_b: Atom, atom_c: Atom) -> bool:
    """returns True if any order of atoms is linked from 1->2->3 and if the the atoms do not form a
    pseudo angle."""
    if is_pseudo(atom_a, atom_b, atom_c):
        return False
    if bond.is_valid(atom_a, atom_b) and bond.is_valid(atom_b, atom_c):
        return True
#    for x, y, z in permutations((atom_a, atom_b, atom_c)):
#        if bond.is_valid(x, y) and bond.is_valid(y, z):
#            return True
    return False


def is_pseudo(atom_a: Atom, atom_b: Atom, atom_c: Atom) -> bool:
    for x, y in combinations([atom_a, atom_b, atom_c], 2):
        if bond.actual_length(x.coordinates, y.coordinates) == 0:
            #TODO: current check could lead to bugs, maybe np.isclose instead?
            return True
    return False

