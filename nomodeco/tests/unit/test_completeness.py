"""
completeness: every method must give the same complete sets as the reference rank filter
icset_opt.complete_set_keys (batched SVD of the full B matrix).
"""

import numpy as np
import pytest

from nomodeco.libraries import completeness, icset_opt
from nomodeco.tests import molecules
from nomodeco.tests.helpers import generate_sets

METHODS = ["svd", "gram", "product"]


@pytest.fixture(scope="module", params=sorted(molecules.ALL))
def case(request):
    mol = molecules.ALL[request.param]()
    ic_dict, a = generate_sets(mol)
    return mol, ic_dict, a


@pytest.mark.parametrize("method", METHODS)
def test_agrees_with_the_reference_rank_filter(case, method):
    mol, ic_dict, a = case
    expected = icset_opt.complete_set_keys(ic_dict, mol, a.idof)
    assert list(completeness.complete_keys(ic_dict, mol, a.idof, method=method)) == expected


@pytest.mark.parametrize("method", METHODS)
def test_rejects_missing_and_dependent_coordinates(method):
    """NH3 needs 6 ICs: dropping an angle loses a motion, three copies of one angle too."""
    mol = molecules.nh3()
    ic_dict, a = generate_sets(mol)
    full = ic_dict[0]
    missing = dict(full, angles=full["angles"][:-1])
    dependent = dict(full, angles=[full["angles"][0]] * 3)
    sets = {0: missing, 1: full, 2: dependent}
    assert list(completeness.complete_keys(sets, mol, a.idof, method=method)) == [1]


def test_off_keeps_every_set():
    mol = molecules.ethylene()
    ic_dict, a = generate_sets(mol)
    assert list(completeness.complete_keys(ic_dict, mol, a.idof, method="off")) == list(ic_dict.keys())


def test_product_mask_is_in_lazy_ic_dict_order():
    """mask[k] must belong to ic_dict[k] (index = (angle, oop, dihedral) in C order)."""
    mol = molecules.acetyl_cyanide()          # 24 sets, 12 complete: the order matters
    ic_dict, a = generate_sets(mol)
    mask = completeness.complete_mask_product(ic_dict, mol, a.idof)
    for k in ic_dict.keys():
        B = completeness.set_b_matrix(mol, ic_dict[k])
        assert mask[k] == completeness.is_complete_svd(B, a.idof)


def test_dependent_dihedral_adds_no_dimension():
    """
    A dihedral that is a combination of the prefix rows leaves only rounding noise after the
    projection; that noise must not count as a new direction.
    """
    rng = np.random.default_rng(0)
    idof = 6
    prefix = rng.normal(size=(5, idof))
    _, s, Vt = np.linalg.svd(prefix)
    N = Vt[5:].T
    row = rng.normal(size=5) @ prefix
    groups = completeness.subsets_by_size([[("d",)]], {("d",): 0})
    assert not completeness.dihedrals_complete(row[None, :] @ N, groups, 1, 1, s[0])[0]


@pytest.mark.parametrize("name, n_rigid", [("h2o", 6), ("nh3", 6), ("co2", 5)])
def test_shape_space_basis(name, n_rigid):
    """Q is orthonormal, has 3N - 6 (3N - 5 linear) columns and B Q keeps the rank of B."""
    mol = molecules.ALL[name]()
    ic_dict, a = generate_sets(mol)
    Q = completeness.shape_space_basis(mol)
    assert Q.shape == (3 * len(mol), 3 * len(mol) - n_rigid)
    assert np.allclose(Q.T @ Q, np.eye(Q.shape[1]))
    B = completeness.set_b_matrix(mol, ic_dict[0])
    assert np.allclose(B @ Q @ Q.T, B)
