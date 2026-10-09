"""
Helper functions of the topology module: ring cutting, IC removal after a cut and the
planar-submolecule parameters.
"""

import logging

from nomodeco.libraries import topology as tp

# cyclobutane-like ring C1-C2-C3-C4 with one hydrogen on C1
RING = [("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C4", "C1"), ("C1", "H1")]
RING_MULT = [("C1", 3), ("C2", 2), ("C3", 2), ("C4", 2), ("H1", 1)]


def test_molecule_is_not_split():
    assert tp.molecule_is_not_split(RING)
    assert tp.molecule_is_not_split(RING[1:])
    assert not tp.molecule_is_not_split([("C1", "C2"), ("C3", "C4")])


def test_eliminate_symmetric_tuples_keeps_first_orientation():
    assert tp.eliminate_symmetric_tuples([("A", "B"), ("B", "A"), ("A", "B"), ("B", "C")]) == [
        ("A", "B"), ("B", "C")]


def test_valide_atoms_to_cut_excludes_terminal_atoms():
    assert tp.valide_atoms_to_cut(RING, RING_MULT) == {"C1", "C2", "C3", "C4"}


# --- ring cutting ---------------------------------------------------------------------


def test_delete_bonds_symmetry_cuts_one_ring_bond():
    valid = tp.valide_atoms_to_cut(RING, RING_MULT)
    group = [("C1", "C2"), ("C3", "C4")]
    removed, bonds = tp.delete_bonds_symmetry(group, RING, 1, valid)
    assert removed == [("C1", "C2")]
    assert bonds == RING  # uncut input, callers remove the cut bonds themselves


def test_delete_bonds_symmetry_never_splits_the_molecule():
    """Cutting both C1-C2 and C3-C4 would split the ring in two: the group can not open mu=2."""
    valid = tp.valide_atoms_to_cut(RING, RING_MULT)
    removed, bonds = tp.delete_bonds_symmetry([("C1", "C2"), ("C3", "C4")], RING, 2, valid)
    assert removed == []
    assert bonds == RING


def test_delete_bonds_symmetry_skips_terminal_bonds():
    valid = tp.valide_atoms_to_cut(RING, RING_MULT)
    removed, _ = tp.delete_bonds_symmetry([("C1", "H1"), ("C2", "C3")], RING, 1, valid)
    assert removed == [("C2", "C3")]


def test_delete_bonds_does_not_mutate_input():
    valid = tp.valide_atoms_to_cut(RING, RING_MULT)
    before = list(RING)
    removed, remaining = tp.delete_bonds(RING, 1, valid)
    assert RING == before
    assert removed == [("C1", "C2")]
    assert remaining == [b for b in RING if b != ("C1", "C2")]
    assert tp.molecule_is_not_split(remaining)


def test_delete_bonds_warns_when_rings_can_not_be_opened(caplog):
    chain = [("C1", "C2"), ("C2", "C3")]
    valid = tp.valide_atoms_to_cut(chain, [("C1", 1), ("C2", 2), ("C3", 1)])
    with caplog.at_level(logging.WARNING):
        removed, remaining = tp.delete_bonds(chain, 1, valid)
    assert removed == []
    assert remaining == chain
    assert "could be cut" in caplog.text


def test_update_internal_coordinates_cyclic_removes_ics_through_the_cut_bond():
    angles = [("C2", "C1", "C4"), ("C1", "C2", "C3"), ("H1", "C1", "C2"), ("H1", "C1", "C4")]
    assert tp.update_internal_coordinates_cyclic([("C2", "C1")], angles) == [("H1", "C1", "C4")]


def test_update_internal_coordinates_cyclic_with_no_cut():
    angles = [("C2", "C1", "C4")]
    assert tp.update_internal_coordinates_cyclic([], angles) == angles


# --- planar submolecules --------------------------------------------------------------


def test_remove_angles_keeps_mult_minus_one_per_planar_center():
    angles = [("H1", "C1", "H2"), ("H1", "C1", "O"), ("H2", "C1", "O"), ("C1", "O", "H3")]
    assert tp.remove_angles(("C1", 3), angles) == [("H1", "C1", "O"), ("H2", "C1", "O"),
                                                    ("C1", "O", "H3")]


def test_param_planar_submolecule():
    """C1 planar (3 neighbours): 2 angles + 1 oop; C2..C4 non-planar, 2 neighbours: 1 angle each."""
    n_phi, n_gamma, _ = tp.get_param_planar_submolecule([("C1", 3)], RING_MULT, [], RING)
    assert (n_phi, n_gamma) == (2 + 3 * 1, 1)


def test_param_planar_submolecule_counts_neighbours_from_the_cut_bonds():
    """After cutting C1-C2, C1 has 2 neighbours although the intact multiplicity list says 3."""
    cut = [b for b in RING if b != ("C1", "C2")]
    n_phi, n_gamma, _ = tp.get_param_planar_submolecule([("C1", 3)], RING_MULT, [], cut)
    # C1: 2 neighbours -> 1 angle, 0 oop; C2 now terminal -> 0; C3, C4: 1 each
    assert (n_phi, n_gamma) == (1 + 0 + 1 + 1, 0)


# --- distributing intermolecular ICs over the submolecule sets ------------------------


def shared_list_sets():
    """Two sets sharing their angle list, as left by the shallow copies of the bond step."""
    angles = ["a1", "a2"]
    return [{0: {"angles": angles, "bonds": ["b1"]}, 1: {"angles": angles, "bonds": ["b2"]}}]


def test_distribute_elements_gives_every_set_exactly_the_requested_number():
    """
    Two of three intermolecular angles per set: each result has its own 2 + 2 angles. The
    combination lists used to be extended in place, so later sets inherited the ICs of earlier
    ones (methanol dimer: 17 instead of 15 angles, 2 redundancies per set).
    """
    # two submolecule sets with different angles, like the three angle subsets of one methanol
    sets = [{0: {"angles": ["a1", "a2"]}, 1: {"angles": ["a1", "a3"]}}]
    result = tp.distribute_elements(sets, ["X1", "X2", "X3"], "angles", 2)[0]
    assert len(result) == 2 * 3
    combos = {frozenset(c) for c in (("X1", "X2"), ("X1", "X3"), ("X2", "X3"))}
    for own in (["a1", "a2"], ["a1", "a3"]):
        mine = [set(s["angles"]) for s in result.values() if set(own) <= set(s["angles"])]
        assert [len(s) for s in mine] == [4, 4, 4]
        assert {frozenset(s - set(own)) for s in mine} == combos


def test_distribute_elements_does_not_change_its_input():
    for take in (2, 3):
        sets = shared_list_sets()
        tp.distribute_elements(sets, ["X1", "X2", "X3"], "angles", take)
        assert sets[0][0]["angles"] == ["a1", "a2"]


def test_add_element_to_all_entries_with_shared_lists():
    sets = shared_list_sets()
    tp.add_element_to_all_entries(sets, ["X1", "X2"], "angles", 1)
    assert [s["angles"] for s in sets[0].values()] == [["a1", "a2", "X1"]] * 2


def test_get_multiplicity():
    assert tp.get_multiplicity("C1", RING_MULT) == 3
    assert tp.get_multiplicity("X", RING_MULT) is None
