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


def test_general_linunit_keeps_oop_at_end_of_linear_unit():
    """
    CH3-C1(=O)-C3#N: non-planar (methyl), C1 is a planar submolecule at the end of the linear
    unit C1-C3-N, so one oop at C1 is needed (idof 18). Removing every oop that contains the
    linear bond C1-C3 left no valid set at all.
    """
    ic_dict, a = generate_sets(molecules.acetyl_cyanide())
    assert a.spec["linearity"] == "linear submolecules found" and a.spec["planar"] == "no"
    assert len(ic_dict) > 0
    for k in range(len(ic_dict)):
        ic_set = ic_dict[k]
        assert n_ics(ic_set) == a.idof
        assert len(ic_set["out of plane angles"]) == 1
        assert ic_set["out of plane angles"][0][0] == "C1"


def test_intermolecular_cyclic_linunit_uses_the_linear_bends_of_the_h_bond():
    """
    Cyclopropanol...water: O1-H6...O2 is linear (175.1 deg), so every set needs its 2 linear
    bends. The linear-angle pool used to be built from the ordinary intermolecular angles
    (typo), which put acceptor-donor angles such as C3-O1...O2 into the linear-bend slot.
    """
    ic_dict, a = generate_sets(molecules.cyclopropanol_water())
    assert a.spec["intermolecular"] == "yes" and a.spec["linearity"] == "linear submolecules found"
    assert len(ic_dict) > 0
    for ic_set in ic_dict.values():
        assert n_ics(ic_set) == a.idof
        bends = [tuple(t) for t in ic_set["linear valence angles"]]
        assert len(bends) == 2
        assert all(t in (("O2", "H6", "O1"), ("O1", "H6", "O2")) for t in bends)


def test_intermolecular_planar_linunit_keeps_oop_at_end_of_linear_unit():
    """
    N#C-H...OH2: the h-bond oop at the acceptor O (wings H1, H2, H3) is at the end of the linear
    unit C-H1...O, not its center, so it is needed (idof 12). It used to be removed because it
    contains the linear bond H1...O, which left 11 ICs.
    """
    ic_dict, a = generate_sets(molecules.hcn_h2o())
    assert a.spec["intermolecular"] == "yes" and a.spec["planar"] == "yes"
    assert a.spec["linearity"] == "linear submolecules found"
    assert len(ic_dict) > 0
    for k in range(len(ic_dict)):
        ic_set = ic_dict[k]
        assert n_ics(ic_set) == a.idof
        assert len(ic_set["out of plane angles"]) == 1
        assert ic_set["out of plane angles"][0][0] == "O"
