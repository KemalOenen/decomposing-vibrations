"""
icset_opt.find_optimal_coordinate_set: the vectorised per-set PED (precomputed G/H/D blocks,
indexed per set) must give the same metric as the plain per-set calculation that main() does
(B_aug -> G -> G_inv -> B_inv -> F_int -> D -> PED). The plain version below is the oracle for
the ped_calc extraction (step 1.1 of the restructure plan).

Hessian: model valence force field F = B^T K B over all generated ICs.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from nomodeco.libraries import bmatrix, icsel, icset_opt
from nomodeco.nomodeco import reciprocal_square_massvector
from nomodeco.tests import molecules
from nomodeco.tests.helpers import generate_sets

ARGS = SimpleNamespace(log=False, matrix_opt="contr")
KINDS = ("bonds", "angles", "linear valence angles", "out of plane angles", "dihedrals")


def model_hessian(mol, a):
    ics = [a.bonds, a.angles, a.linear_angles, a.out_of_plane, a.dihedrals]
    B = bmatrix.b_matrix(mol, *ics, 0)
    k = np.concatenate([
        np.linspace(0.5, 0.9, len(a.bonds)),
        np.linspace(0.10, 0.20, len(a.angles)),
        np.full(len(a.linear_angles), 0.06),
        np.linspace(0.02, 0.04, len(a.out_of_plane)),
        np.linspace(0.010, 0.015, len(a.dihedrals)),
    ])
    return B.T @ np.diag(k) @ B


@pytest.fixture(scope="module", params=["h2o", "nh3", "ethylene"])
def case(request):
    mol = molecules.ALL[request.param]()
    ic_dict, a = generate_sets(mol)
    F = model_hessian(mol, a)
    m = reciprocal_square_massvector(mol)
    eigenvalues, L = np.linalg.eigh(m[:, None] * F * m[None, :])
    rottra = L[:, : 3 * len(mol) - a.idof]
    return SimpleNamespace(mol=mol, ic_dict=ic_dict, a=a, F=F, m=m, L=L, rottra=rottra)


def find_optimal(case, ic_dict):
    return icset_opt.find_optimal_coordinate_set(
        ic_dict, ARGS, case.a.idof, np.diag(case.m ** 2), np.diag(case.m), case.rottra,
        case.F, case.mol, {}, case.L, 0.0, 0.0,
    )


def reference_metric(case, ic_set):
    """Per-set PED as in nomodeco.main(); None if the set is not complete."""
    mol, idof, m = case.mol, case.a.idof, case.m
    ics = [ic_set[k] for k in KINDS]
    n_internals = sum(map(len, ics))
    red = n_internals - idof
    B = np.concatenate([bmatrix.b_matrix(mol, *ics, 0), case.rottra.T])
    G = B @ np.diag(m ** 2) @ B.T
    e, K = np.linalg.eigh(G)
    order = e.argsort()[::-1]
    e, K = e[order], K[:, order]
    if red > 0:
        e, K = e[:-red], K[:, :-red]
    G_inv = K @ np.diag(1 / e) @ K.T
    B_inv = np.diag(m ** 2) @ B.T @ G_inv
    F_int = B_inv.T @ case.F @ B_inv
    if not np.allclose(B.T @ F_int @ B, case.F):
        return None
    D = B @ (np.diag(m) @ case.L)
    n_rottra = case.rottra.shape[1]
    vib = range(n_rottra, n_rottra + idof)
    eigenvalues = np.diag(D.T @ F_int @ D)
    diag = np.array([[D[i, k] ** 2 * F_int[i, i] / eigenvalues[k] for k in vib]
                     for i in range(n_internals)])
    contribution = diag / diag.sum(axis=0) * 100
    return icsel.Kemalian_metric(contribution, diag, 0, 0.0, 0.0, None)


def test_every_set_matches_the_plain_calculation(case):
    for k in case.ic_dict.keys():
        single = {0: case.ic_dict[k]}
        expected = reference_metric(case, case.ic_dict[k])
        result = find_optimal(case, single)
        if expected is None:
            assert result["set"] is None
        else:
            assert result["metric"] == pytest.approx(expected, rel=1e-8)


def test_selects_the_set_with_the_highest_metric(case):
    metrics = {k: reference_metric(case, case.ic_dict[k]) for k in case.ic_dict.keys()}
    metrics = {k: v for k, v in metrics.items() if v is not None}
    result = find_optimal(case, case.ic_dict)
    assert result["metric"] == pytest.approx(max(metrics.values()), rel=1e-8)
    assert metrics[result["best_key"]] == pytest.approx(result["metric"], rel=1e-8)
    assert result["set"] == case.ic_dict[result["best_key"]]


def test_empty_ic_dict(case):
    assert find_optimal(case, {}) == {"best_key": None, "metric": 0, "set": None}


def test_incomplete_set_is_never_selected():
    """H2O with the angle missing spans only 2 of 3 vibrations."""
    mol = molecules.h2o()
    ic_dict, a = generate_sets(mol)
    F = model_hessian(mol, a)
    m = reciprocal_square_massvector(mol)
    _, L = np.linalg.eigh(m[:, None] * F * m[None, :])
    case = SimpleNamespace(mol=mol, a=a, F=F, m=m, L=L, rottra=L[:, : 3 * len(mol) - a.idof])
    incomplete = dict(ic_dict[0], angles=[])
    assert reference_metric(case, incomplete) is None
    result = find_optimal(case, {0: incomplete, 1: ic_dict[0]})
    assert result["best_key"] == 1
