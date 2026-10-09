"""
decius: number of ICs of each kind in a complete set (J. C. Decius, J. Chem. Phys. 17, 1315 (1949)).
"""

import pytest

from nomodeco.libraries import decius
from nomodeco.tests import molecules
from nomodeco.tests.helpers import analyse


def terminal_atoms(bonds):
    neighbours = {}
    for a, b in bonds:
        neighbours[a] = neighbours.get(a, 0) + 1
        neighbours[b] = neighbours.get(b, 0) + 1
    return sum(1 for n in neighbours.values() if n == 1)


@pytest.mark.parametrize("b, a, a_1", [(2, 3, 2), (5, 6, 4), (11, 12, 6), (7, 8, 5)])
def test_totals_without_linear_units(b, a, a_1):
    assert decius.planar(b, a, a_1).total() == 6 * b - 3 * a
    assert decius.general(b, a, a_1).total() == 6 * b - 3 * a


@pytest.mark.parametrize("name, kind", [("h2o", "planar"), ("ethylene", "planar"), ("nh3", "general")])
def test_acyclic_molecules_get_idof_coordinates(name, kind):
    mol = molecules.ALL[name]()
    a = analyse(mol)
    counts = getattr(decius, kind)(len(a.bonds), len(mol), terminal_atoms(a.bonds))
    assert counts.total() == a.idof


def test_ethylene_counts():
    """C2H4: 5 bonds, 4 angles, 2 oops, 1 torsion."""
    assert tuple(decius.planar(5, 6, 4)) == (5, 4, 2, 1, 0)


def test_ammonia_counts():
    """NH3: 3 bonds, 3 angles, no torsions."""
    assert tuple(decius.general(3, 4, 3)) == (3, 3, 0, 0, 0)


def test_linear_unit_terms():
    """(l - 1) fewer angles, 2 (l - 1) linear bends; general also (l - 1) fewer torsions."""
    b, a, a_1, l = 6, 7, 4, 4
    plain_p, lin_p = decius.planar(b, a, a_1), decius.planar(b, a, a_1, l)
    assert lin_p.n_phi == plain_p.n_phi - (l - 1)
    assert lin_p.n_phi_prime == 2 * (l - 1)
    assert lin_p.n_tau == plain_p.n_tau          # planar: corrected later per terminal linear bond
    assert lin_p.total() == 6 * b - 3 * a + (l - 1)

    plain_g, lin_g = decius.general(b, a, a_1), decius.general(b, a, a_1, l)
    assert lin_g.n_tau == plain_g.n_tau - (l - 1)
    assert lin_g.total() == 6 * b - 3 * a


def test_no_linear_unit_is_none_not_zero():
    """l = 0 is a (degenerate) value, not 'no linear unit': it gives l - 1 = -1."""
    assert decius.planar(5, 6, 4, 0).n_phi == decius.planar(5, 6, 4).n_phi + 1


def test_planar_submolecule_counts():
    """C1 planar with 3 neighbours: 2 angles + 1 oop; C2 non-planar with 4 neighbours: 5 angles."""
    bonds = [("C1", "O"), ("C1", "H1"), ("C1", "C2"), ("C2", "H2"), ("C2", "H3"), ("C2", "H4")]
    mult = [("C1", 3), ("O", 1), ("H1", 1), ("C2", 4), ("H2", 1), ("H3", 1), ("H4", 1)]
    assert decius.planar_submolecule_counts({"C1"}, mult, bonds) == (2 + 5, 1)
    assert decius.planar_submolecule_counts(set(), mult, bonds) == (3 + 5, 0)


def test_counts_unpack_in_field_order():
    n_r, n_phi, n_gamma, n_tau, n_phi_prime = decius.planar(6, 7, 4, 4)
    c = decius.planar(6, 7, 4, 4)
    assert (n_r, n_phi, n_gamma, n_tau, n_phi_prime) == (c.n_r, c.n_phi, c.n_gamma, c.n_tau, c.n_phi_prime)
