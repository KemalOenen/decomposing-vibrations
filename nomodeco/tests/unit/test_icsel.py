"""
icsel helpers: subset enumeration, symmetry grouping and the metric.
"""

import logging

import numpy as np
import pytest

from nomodeco.libraries import icsel

H2O_SPEC = {"equivalent_atoms": [["H1", "H2"]]}


# --- _flat_subsets_of_size ------------------------------------------------------------


def test_flat_subsets_target_zero_is_the_empty_selection():
    assert icsel._flat_subsets_of_size([], 0) == [[]]
    assert icsel._flat_subsets_of_size([["a"], ["b", "c"]], 0) == [[]]


def test_flat_subsets_no_groups():
    assert icsel._flat_subsets_of_size([], 2) == []


def test_flat_subsets_combines_whole_groups_only():
    groups = [["a"], ["b", "c"], ["d", "e"]]
    result = icsel._flat_subsets_of_size(groups, 3)
    assert sorted(map(sorted, result)) == [["a", "b", "c"], ["a", "d", "e"]]


def test_flat_subsets_unreachable_target():
    assert icsel._flat_subsets_of_size([["a", "b"], ["c", "d"]], 3) == []


def test_flat_subsets_is_capped(monkeypatch):
    monkeypatch.setattr(icsel, "_MAX_SUBSETS", 5)
    groups = [[i] for i in range(10)]
    assert len(icsel._flat_subsets_of_size(groups, 3)) == 5


# --- angle / dihedral subsets ---------------------------------------------------------


def test_angle_subsets_exact_size():
    symmetric = {("H1", "O", "H2"): [("H1", "O", "H2")]}
    assert icsel.get_angle_subsets(symmetric, 2, 1, 3, 1) == [[("H1", "O", "H2")]]


def test_angle_subsets_fall_back_to_one_redundant_angle():
    """Two groups of two can not give 3 angles; the symmetric 4-angle selection is returned."""
    symmetric = {"a1": ["a1", "a2"], "a2": ["a1", "a2"], "b1": ["b1", "b2"], "b2": ["b1", "b2"]}
    result = icsel.get_angle_subsets(symmetric, 0, 4, 0, 3)
    assert sorted(map(sorted, result)) == [["a1", "a2", "b1", "b2"]]


def test_dihedral_subsets_no_dihedrals_needed():
    assert icsel.get_dihedral_subsets({}, 0, 0, 0, 0) == [[]]


def test_dihedral_subsets_give_up_after_one_redundant():
    symmetric = {"d1": ["d1", "d2", "d3"]}
    assert icsel.get_dihedral_subsets(symmetric, 0, 0, 0, 1) == []


# --- out-of-plane subsets -------------------------------------------------------------


OOPS = [("C1", "H1", "H2", "C2"), ("C1", "H2", "H1", "C2"), ("C1", "C2", "H1", "H2"),
        ("C2", "H3", "H4", "C1"), ("C2", "H4", "H3", "C1"), ("C2", "C1", "H3", "H4")]


def test_oop_subsets_none_needed():
    assert icsel.get_oop_subsets(OOPS, 0) == [[]]


def test_oop_subsets_one_per_center():
    result = icsel.get_oop_subsets(OOPS, 2)
    assert len(result) == 3 * 3
    for subset in result:
        assert sorted(oop[0] for oop in subset) == ["C1", "C2"]


def test_oop_subsets_single_center_choices():
    result = icsel.get_oop_subsets(OOPS, 1)
    assert sorted(map(tuple, result)) == sorted((oop,) for oop in OOPS)


def test_oop_subsets_more_needed_than_centers():
    assert icsel.get_oop_subsets(OOPS, 3) == []


# --- symmetry grouping ----------------------------------------------------------------


def test_symm_bonds_groups_equivalent_bonds_in_either_orientation():
    bonds = [("O", "H1"), ("H2", "O")]
    groups = icsel.get_symm_bonds(bonds, H2O_SPEC)
    assert groups[("O", "H1")] == bonds
    assert groups[("H2", "O")] == bonds


def test_symm_angles_matches_reversed_angle():
    angles = [("H1", "O", "H2"), ("H2", "O", "H1")]
    groups = icsel.get_symm_angles(angles, H2O_SPEC)
    assert groups[("H1", "O", "H2")] == angles


def test_symm_angles_without_equivalent_atoms():
    angles = [("H1", "C", "H2"), ("H1", "C", "O")]
    groups = icsel.get_symm_angles(angles, {"equivalent_atoms": []})
    assert groups == {angles[0]: [angles[0]], angles[1]: [angles[1]]}


def test_symm_dihedrals_matches_reversed_dihedral():
    spec = {"equivalent_atoms": [["H1", "H2", "H3", "H4"], ["C1", "C2"]]}
    dihedrals = [("H1", "C1", "C2", "H3"), ("H4", "C2", "C1", "H2")]
    groups = icsel.get_symm_dihedrals(dihedrals, spec)
    assert groups[dihedrals[0]] == dihedrals


def test_number_terminal_bonds():
    assert icsel.number_terminal_bonds([("O", 2), ("H1", 1), ("H2", 1)]) == 2


# --- metric ---------------------------------------------------------------------------


CONTRIBUTIONS = np.array([[90.0, 10.0], [20.0, 80.0], [50.0, 50.0]])


def test_metric_is_mean_of_row_maxima():
    diag = np.full((3, 2), 0.5)
    assert icsel.Kemalian_metric(CONTRIBUTIONS, diag, 0, 0.0, 0.0, None) == pytest.approx(220 / 3)


def test_metric_penalties():
    diag = np.array([[1.2, 0.0], [0.5, 0.5], [0.0, 1.0]])  # one diagonal PED element > 1
    value = icsel.Kemalian_metric(CONTRIBUTIONS, diag, 2, 1.5, 3.0, None)
    assert value == pytest.approx(220 / 3 - 1.5 * 2 - 0.2 * 10 * 3.0)


def test_metric_rejects_strongly_negative_entries():
    matrix = CONTRIBUTIONS.copy()
    matrix[1, 1] = -1.5
    assert icsel.Kemalian_metric(matrix, np.zeros((3, 2)), 0, 0.0, 0.0, None) == 0


def test_metric_and_logging_metric_agree():
    diag = np.array([[1.2, 0.0], [0.5, 0.5], [0.0, 1.0]])
    log = logging.getLogger("test_icsel")
    for matrix in (CONTRIBUTIONS, -2 * CONTRIBUTIONS):
        assert icsel.Kemalian_metric_log(matrix, diag, 2, 1.5, 3.0, log) == pytest.approx(
            icsel.Kemalian_metric(matrix, diag, 2, 1.5, 3.0, None)
        )


# --- completeness ---------------------------------------------------------------------


def test_completeness():
    rng = np.random.default_rng(0)
    B = rng.normal(size=(3, 3))
    F_int = np.diag([1.0, 2.0, 3.0])
    F_cart = B.T @ F_int @ B
    assert icsel.test_completeness(F_cart, B, None, F_int)
    assert not icsel.test_completeness(F_cart + 1e-3, B, None, F_int)
