"""
Normal-mode aligned linear bends (linear_bends.linear_bend_references).

The Hessian is a model one, F = B^T K B from a diagonal valence force field. Both bends of a
linear triple get the same force constant, so F is isotropic about the linear axis and has the
full symmetry of the molecule. The bend content of mode k in bend i is reported as
D[i, k]^2 / sum_j D[j, k]^2 over the two bends of the triple (equal to the PED share for a
diagonal internal F).
"""

import numpy as np
import pytest

from nomodeco.libraries import bmatrix, linear_bends
from nomodeco.libraries.molecule_class import Molecule
from nomodeco.nomodeco import reciprocal_square_massvector
from nomodeco.tests import molecules
from nomodeco.tests.helpers import analyse

TOL = 1e-8


def rotation(theta_deg, phi_deg):
    """Tilt by theta about x, then turn by phi about z."""
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    rx = np.array([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])
    rz = np.array([[np.cos(p), -np.sin(p), 0], [np.sin(p), np.cos(p), 0], [0, 0, 1]])
    return rz @ rx


def rotated(mol, R):
    return Molecule([Molecule.Atom(a.symbol, tuple(R @ np.array(a.coordinates))) for a in mol])


def both_bends(linear_angles):
    return [t for t in dict.fromkeys(map(tuple, linear_angles)) for _ in range(2)]


def model_normal_modes(mol, a):
    """(eigenvalues, L, diag_reciprocal_square) of the mass-weighted model Hessian."""
    bends = both_bends(a.linear_angles)
    B = bmatrix.b_matrix(mol, a.bonds, a.angles, bends, a.out_of_plane, a.dihedrals, 0)
    k = np.concatenate([
        np.linspace(0.5, 0.9, len(a.bonds)),
        np.linspace(0.10, 0.20, len(a.angles)),
        np.repeat(np.linspace(0.05, 0.08, len(bends) // 2), 2),
        np.linspace(0.02, 0.04, len(a.out_of_plane)),
        np.linspace(0.010, 0.015, len(a.dihedrals)),
    ])
    F = B.T @ np.diag(k) @ B
    m = reciprocal_square_massvector(mol)
    eigenvalues, L = np.linalg.eigh(m[:, None] * F * m[None, :])
    return eigenvalues, L, m


def bend_shares(mol, triple, L, m, bend_refs):
    """2 x 3N array: share of each mode in the two bends of triple."""
    B = bmatrix.b_matrix(mol, [], [], [triple, triple], [], [], 0, bend_refs)
    D = B @ (m[:, None] * L)
    norm = (D ** 2).sum(axis=0)
    return np.divide(D ** 2, norm, out=np.zeros_like(D), where=norm > 1e-12)


def aligned(mol, scramble_deg=0.0):
    """Run linear_bend_references; scramble_deg first rotates every degenerate pair of L."""
    a = analyse(mol)
    eigenvalues, L, m = model_normal_modes(mol, a)
    if scramble_deg:
        c, s = np.cos(np.radians(scramble_deg)), np.sin(np.radians(scramble_deg))
        for i in range(L.shape[1] - a.idof, L.shape[1] - 1):
            if abs(eigenvalues[i] - eigenvalues[i + 1]) < 1e-12:
                L[:, [i, i + 1]] = L[:, [i, i + 1]] @ np.array([[c, -s], [s, c]])
    L, refs = linear_bends.linear_bend_references(mol, a.linear_angles, L, eigenvalues, m, a.idof)
    return a, L, m, refs


@pytest.mark.parametrize("theta, phi", [(0, 0), (15, 0), (30, 20), (45, 45), (60, 70), (90, 30)])
@pytest.mark.parametrize("scramble", [0.0, 37.0])
def test_co2_degenerate_bends_are_100_0_at_every_tilt(theta, phi, scramble):
    """
    The two bending modes of CO2 are degenerate, so eigh returns them in an arbitrary rotation.
    After alignment each one must be a pure bend: 100/0 and 0/100, for any orientation of the
    molecule and any rotation of the degenerate pair.
    """
    mol = rotated(molecules.co2(), rotation(theta, phi))
    a, L, m, refs = aligned(mol, scramble)
    (triple,) = refs
    shares = bend_shares(mol, triple, L, m, refs)
    bending = [k for k in range(L.shape[1] - a.idof, L.shape[1]) if shares[:, k].any()]
    assert len(bending) == 2
    first, second = bending
    assert shares[:, first] == pytest.approx([1, 0], abs=TOL)
    assert shares[:, second] == pytest.approx([0, 1], abs=TOL)


def test_reference_vectors_are_orthonormal_and_perpendicular_to_the_axis():
    mol = rotated(molecules.co2(), rotation(30, 20))
    _, _, _, refs = aligned(mol)
    xyz = {atom.symbol: np.array(atom.coordinates) for atom in mol}
    for (M, O, _), (w1, w2) in refs.items():
        u = (xyz[M] - xyz[O]) / np.linalg.norm(xyz[M] - xyz[O])
        W = np.array([w1, w2, u])
        assert W @ W.T == pytest.approx(np.eye(3), abs=TOL)


@pytest.mark.parametrize("theta, phi", [(0, 0), (35, 50), (80, 10)])
def test_hcocn_bends_follow_the_cs_mirror_plane(theta, phi):
    """
    HCOCN is planar (Cs). Its modes are A' (in plane) or A'' (out of plane), so one bend
    reference must be the plane normal and the other must lie in the plane, and no mode may mix
    the two bends.
    """
    R = rotation(theta, phi)
    mol = rotated(molecules.hcocn(), R)
    a, L, m, refs = aligned(mol)
    normal = R @ np.array([0.0, 0.0, 1.0])
    (triple,) = refs
    w1, w2 = refs[triple]
    assert sorted([abs(w1 @ normal), abs(w2 @ normal)]) == pytest.approx([0, 1], abs=TOL)

    shares = bend_shares(mol, triple, L, m, refs)
    vib = range(L.shape[1] - a.idof, L.shape[1])
    assert all(shares[:, k].min() < TOL for k in vib)


def test_no_linear_angles_leaves_modes_unchanged():
    mol = molecules.h2o()
    a = analyse(mol)
    eigenvalues, L, m = model_normal_modes(mol, a)
    L_out, refs = linear_bends.linear_bend_references(mol, a.linear_angles, L, eigenvalues, m, a.idof)
    assert refs == {}
    assert np.array_equal(L_out, L)
