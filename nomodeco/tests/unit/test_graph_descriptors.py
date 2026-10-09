"""
graph_descriptors: global, per-submolecule and per-atom descriptors of the molecular graph, and
their section in the .out file.
"""

import io
import logging

import numpy as np
import pytest

from nomodeco.libraries import logfile
from nomodeco.libraries.graph_descriptors import graph_descriptors
from nomodeco.tests import molecules


def test_benzene_ring_matches_closed_forms():
    """The heavy-atom graph of benzene is the cycle C6."""
    d = graph_descriptors(molecules.benzene())
    assert d["cyclomatic_number"] == 1
    assert d["ring_sizes"] == [6]
    assert d["rings"] == [["C1", "C2", "C3", "C4", "C5", "C6"]]
    ring = d["submolecules"][0]["heavy"]
    n = 6
    assert ring["wiener"] == n ** 3 / 8                     # even cycle: n^3 / 8
    assert ring["randic"] == pytest.approx(n / 2)
    assert ring["zagreb_m1"] == ring["zagreb_m2"] == 4 * n
    assert ring["balaban_j"] == pytest.approx(2.0)
    assert ring["algebraic_connectivity"] == pytest.approx(2 - 2 * np.cos(2 * np.pi / n))
    assert ring["kirchhoff"] == pytest.approx((n ** 3 - n) / 12)
    assert ring["graph_energy"] == pytest.approx(8.0)
    assert (ring["diameter"], ring["radius"]) == (3, 3)


def test_water_dimer_counts_every_bond_kind():
    d = graph_descriptors(molecules.water_dimer())
    assert d["n_bonds"] == {"cov": 4, "h_bond": 1, "acc_don": 1}
    assert d["n_components"] == 2
    assert d["n_components_with_h_bonds"] == 1    # the h-bond holds the dimer together
    assert d["n_zero_laplacian"] == 2             # one zero eigenvalue per covalent component
    assert [s["formula"] for s in d["submolecules"]] == ["H2O", "H2O"]


def test_no_h_bond_entry_without_h_bonds():
    assert "n_components_with_h_bonds" not in graph_descriptors(molecules.h2o())


@pytest.mark.parametrize("name, formula", [
    ("h2o", "H2O"), ("nh3", "H3N"), ("co2", "CO2"), ("hcocn", "C2HNO"),
    ("acetyl_cyanide", "C3H3NO"), ("benzene", "C6H6"),
])
def test_formula_in_hill_order(name, formula):
    assert graph_descriptors(molecules.ALL[name]())["submolecules"][0]["formula"] == formula


def test_degrees_terminal_atoms_and_bridges():
    d = graph_descriptors(molecules.propyne())
    assert d["degree_histogram"] == {1: 4, 2: 2, 4: 1}
    assert d["n_terminal"] == 4 and d["n_branching"] == 1
    assert d["bridges"] == 6                      # acyclic: every bond is a bridge
    assert d["articulation_atoms"] == ["C1", "C2", "C3"]


def test_heavy_atom_graph_needs_two_connected_heavy_atoms():
    assert graph_descriptors(molecules.h2o())["submolecules"][0]["heavy"] is None


def test_atom_table():
    d = graph_descriptors(molecules.benzene())
    rows = {r[0]: r for r in d["atoms"]}
    assert [r[0] for r in d["atoms"]] == [a.symbol for a in molecules.benzene()]
    assert rows["C1"][1:3] == (3, True) and rows["H1"][1:3] == (1, False)


def test_out_section_is_written(mol_name):
    stream = io.StringIO()
    log = logging.getLogger(f"graph-descriptors-{mol_name}")
    log.handlers = [logging.StreamHandler(stream)]
    log.setLevel(logging.INFO)
    log.propagate = False
    logfile.write_graph_descriptors(log, graph_descriptors(molecules.ALL[mol_name]()))
    text = stream.getvalue()
    assert "Molecular Graph" in text and "Per-atom descriptors" in text
    assert "Submolecule 1" in text
