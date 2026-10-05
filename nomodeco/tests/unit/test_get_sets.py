"""
IC-set generation (icsel.get_sets and the topology variants): every generated set must have
exactly idof coordinates of the right kinds.
"""

from nomodeco.tests import molecules
from nomodeco.tests.helpers import generate_sets, n_ics


def test_planar_linunit_keeps_oop_at_end_of_linear_unit():
    """
    H-C1(=O)-C2#N: C1 is the end of the linear unit C1-C2-N, not its center, so its oop is
    needed (4 bonds + 2 angles + 2 linear bends + 1 oop = idof 9). It used to be removed
    because it contains the linear bond C1-C2, which left 8 ICs.
    """
    ic_dict, a = generate_sets(molecules.hcocn())
    assert a.spec["linearity"] == "linear submolecules found" and a.spec["planar"] == "yes"
    assert len(ic_dict) > 0
    for k in range(len(ic_dict)):
        ic_set = ic_dict[k]
        assert n_ics(ic_set) == a.idof
        assert len(ic_set["out of plane angles"]) == 1
        assert ic_set["out of plane angles"][0][0] == "C1"
