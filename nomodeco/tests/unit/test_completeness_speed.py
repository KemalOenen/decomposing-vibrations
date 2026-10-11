"""
Speed of the completeness check: the product method against the reference rank filter
(icset_opt.complete_set_keys). Run with `pytest -m benchmark`.

The test molecules give only a few sets, so the LazyIcDict is built like an unsymmetric search:
all angle subsets x all dihedral subsets of the sizes the topology picks.
"""

import itertools
import time

import numpy as np
import pytest

from nomodeco.libraries import completeness, icset_opt, orca_parser
from nomodeco.libraries.molecule_class import Molecule
from nomodeco.libraries.topology import LazyIcDict
from nomodeco.tests.conftest import ORCA_DATA
from nomodeco.tests.helpers import generate_sets

pytestmark = pytest.mark.benchmark

CASES = {
    "ethanol": ORCA_DATA / "C2H6O" / "ethanol.property.txt",
    "acetic_acid": ORCA_DATA / "CH3COOH" / "CH3COOH.property.txt",
}
MIN_SPEEDUP = 2          # measured 4x (acetic acid) and 7x (ethanol), best of 3


def unsymmetric_search(mol, max_subsets=3000):
    """One product block: every angle subset x every dihedral subset (capped)."""
    ic_dict, a = generate_sets(mol)
    s = ic_dict[0]
    angle_subsets = itertools.combinations(a.angles, len(s["angles"]))
    dihedral_subsets = itertools.combinations(a.dihedrals, len(s["dihedrals"]))
    lazy = LazyIcDict()
    lazy.add_product(
        s["bonds"], s["linear valence angles"],
        [list(x) for x in itertools.islice(angle_subsets, 100)],
        [s["out of plane angles"]],
        [list(x) for x in itertools.islice(dihedral_subsets, max_subsets)],
    )
    return lazy, a.idof


def best_time(f, repeat=3):
    times = []
    for _ in range(repeat):
        start = time.perf_counter()
        result = f()
        times.append(time.perf_counter() - start)
    return min(times), result


@pytest.mark.parametrize("name", sorted(CASES))
def test_product_is_faster_than_the_reference(name):
    with open(CASES[name]) as f:
        mol = Molecule(orca_parser.parse_xyz_from_inputfile(f))
    ic_dict, idof = unsymmetric_search(mol)

    t_ref, expected = best_time(lambda: icset_opt.complete_set_keys(ic_dict, mol, idof))
    t_product, keys = best_time(lambda: completeness.complete_keys(ic_dict, mol, idof))
    print(f"\n{name}: {len(ic_dict)} sets, {len(expected)} complete, "
          f"reference {t_ref:.3f}s, product {t_product:.3f}s, {t_ref / t_product:.1f}x")

    assert list(keys) == expected
    assert 0 < len(expected) < len(ic_dict)      # the check has to separate something
    assert t_ref / t_product > MIN_SPEEDUP
