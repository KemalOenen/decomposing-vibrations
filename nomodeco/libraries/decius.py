"""
Decius counts: how many internal coordinates of each kind a complete, non-redundant set needs.

J. C. Decius, J. Chem. Phys. 17, 1315 (1949). With
    b    bonds (for a cyclic molecule: after opening the mu rings, b = b_total - mu)
    a    atoms
    a_1  terminal atoms (one bond)
    l    atoms in the linear submolecule(s) (specification["length of linear submolecule(s) l"])
the counts are

                n_r   n_phi                 n_phi' (linear bends)  n_gamma          n_tau
    planar      b     2b - a - (l - 1)      2 (l - 1)               2(b - a) + a_1   b - a_1
    general     b     4b - 3a + a_1 - (l-1) 2 (l - 1)               0                b - a_1 - (l - 1)

with the (l - 1) terms only for molecules with linear submolecules. Without linear units both
rows add up to 6b - 3a, which is 3a - 6 = idof for an acyclic molecule (b = a - 1).

Planar molecules keep n_tau = b - a_1 here; the topology functions subtract (l - 1) once per
terminal linear bond and remove one oop per linear center afterwards, because those corrections
depend on the generated ICs.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class DeciusCounts:
    n_r: int
    n_phi: int
    n_gamma: int
    n_tau: int
    n_phi_prime: int = 0

    def __iter__(self):
        """n_r, n_phi, n_gamma, n_tau, n_phi_prime = decius.planar(...)"""
        return iter((self.n_r, self.n_phi, self.n_gamma, self.n_tau, self.n_phi_prime))

    def total(self) -> int:
        return self.n_r + self.n_phi + self.n_phi_prime + self.n_gamma + self.n_tau


def _linear_dof(l) -> int:
    """l - 1 for a linear unit of l atoms, 0 without linear units (l=None)."""
    return 0 if l is None else l - 1


def planar(b, a, a_1, l=None) -> DeciusCounts:
    lin = _linear_dof(l)
    return DeciusCounts(
        n_r=b,
        n_phi=2 * b - a - lin,
        n_gamma=2 * (b - a) + a_1,
        n_tau=b - a_1,
        n_phi_prime=2 * lin,
    )


def general(b, a, a_1, l=None) -> DeciusCounts:
    lin = _linear_dof(l)
    return DeciusCounts(
        n_r=b,
        n_phi=4 * b - 3 * a + a_1 - lin,
        n_gamma=0,
        n_tau=b - a_1 - lin,
        n_phi_prime=2 * lin,
    )


def planar_submolecule_counts(planar_atoms, multiplicity_list, bonds) -> tuple[int, int]:
    """
    (n_phi, n_gamma) of a non-planar molecule with planar submolecules, summed per atom with
    m neighbours (m > 1): a planar center needs m - 1 angles and m - 2 oops, any other center
    2m - 3 angles. Neighbours are counted from bonds, which may be ring-opened;
    multiplicity_list (from the intact molecule) only gives the atoms.
    """
    n_neighbours = Counter(atom for bond in bonds for atom in bond)
    planar_atoms = set(planar_atoms)
    n_phi = 0
    n_gamma = 0
    for atom, _ in multiplicity_list:
        m = n_neighbours[atom]
        if m > 1:
            if atom in planar_atoms:
                n_phi += m - 1
                n_gamma += m - 2
            else:
                n_phi += 2 * m - 3
    return n_phi, n_gamma
