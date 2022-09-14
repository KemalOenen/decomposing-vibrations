from __future__ import annotations

from itertools import permutations, combinations
import bond


def is_valid(atom_a: Atom, atom_b: Atom, atom_c: Atom, atom_d: Atom) -> bool:
    """returns True if any order of atoms 2,3,4 is linked to a central atom 1 and if the the atoms do not form a
    pseudo out-of-plane angle."""
    if is_pseudo(atom_a, atom_b, atom_c, atom_d):
        return False
    if bond.is_valid(atom_a, atom_b) and bond.is_valid(atom_a, atom_c) and bond.is_valid(atom_a, atom_d):
        return True
    return False


def is_pseudo(atom_a: Atom, atom_b: Atom, atom_c: Atom, atom_d: Atom) -> bool:
    for x, y in combinations([atom_a, atom_b, atom_c, atom_d], 2):
        if bond.actual_length(x.coordinates, y.coordinates) == 0:
            # maybe np.isclose instead?
            return True
    return False

