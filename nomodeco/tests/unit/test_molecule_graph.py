"""
Molecule.graph(): the molecular graph (atoms as nodes; covalent, hydrogen and acceptor-donor
bonds as edges) and the quantities derived from it.
"""

import networkx as nx
import numpy as np
import pytest

from nomodeco.tests import molecules


def test_nodes_carry_index_and_element(mol_name):
    mol = molecules.ALL[mol_name]()
    G = mol.graph()
    assert list(G.nodes) == mol.list_of_atom_symbols()
    for i, atom in enumerate(mol):
        assert G.nodes[atom.symbol]["index"] == i
        assert G.nodes[atom.symbol]["element"] == atom.symbol.rstrip("0123456789")


def test_covalent_bonds_match_the_degree_of_covalence(mol_name):
    mol = molecules.ALL[mol_name]()
    G = mol.graph()
    assert G.graph["cov_bonds"] == mol.covalent_bonds(mol.degree_of_covalance())
    assert sorted(map(sorted, mol.bond_graph("cov").edges)) == sorted(map(sorted, G.graph["cov_bonds"]))


def test_edge_kinds_and_degree_of_covalence():
    mol = molecules.water_dimer()
    G = mol.graph()
    assert len(G.graph["cov_bonds"]) == 4
    assert G.graph["h_bonds"] == [("H2", "O2")]
    for a, b in G.graph["cov_bonds"]:
        assert G.edges[a, b]["kind"] == "cov" and G.edges[a, b]["degofc"] > 0.75
    assert G.edges["H2", "O2"]["kind"] == "h_bond"
    assert 0.27 < G.edges["H2", "O2"]["degofc"] < 0.7
    for a, b in G.graph["acc_don_bonds"]:
        assert G.edges[a, b]["kind"] == "acc_don"


@pytest.mark.parametrize("name, components, mu, beta", [
    ("h2o", 1, 0, 0),
    ("benzene", 1, 1, 1),
    ("water_dimer", 2, 0, 0),
    ("hcn_h2o", 2, 0, 0),
])
def test_components_and_cycle_numbers(name, components, mu, beta):
    mol = molecules.ALL[name]()
    assert len(mol.graph().graph["components"]) == components
    assert mol.mu() == mu
    assert mol.beta() == beta


def test_detect_submolecules():
    mol = molecules.water_dimer()
    components, submolecule_bonds, symbols = mol.detect_submolecules()
    assert sorted(map(sorted, components)) == [["H1", "H2", "O1"], ["H3", "H4", "O2"]]
    assert sum(map(len, submolecule_bonds)) == 4
    assert [symbols[i] for i in range(2)] == components


def test_detect_submolecules_returns_copies():
    mol = molecules.water_dimer()
    components, _, _ = mol.detect_submolecules()
    components[0].add("X")
    assert all("X" not in c for c in mol.graph().graph["components"])


def test_graph_is_cached_and_frozen():
    mol = molecules.h2o()
    G = mol.graph()
    assert mol.graph() is G
    assert nx.is_frozen(G)
    with pytest.raises(nx.NetworkXError):
        G.add_edge("O", "X")


def test_graph_is_rebuilt_when_the_geometry_changes():
    mol = molecules.h2o()
    G = mol.graph()
    mol[1].coordinates = (0.0, 3.0, -0.4692)  # pull H1 away: O-H1 is no longer a bond
    G_new = mol.graph()
    assert G_new is not G
    assert ("O", "H1") not in G_new.graph["cov_bonds"]


def test_explicit_degofc_table_is_not_cached():
    mol = molecules.h2o()
    G = mol.graph(mol.degree_of_covalance())
    assert mol.graph() is not G


def test_adjacency_matrices():
    mol = molecules.water_dimer()
    labels = mol.list_of_atom_symbols()
    A = mol.covalent_adjacency_matrix()
    assert np.array_equal(A, A.T) and A.sum() == 2 * 4
    H = mol.hydrogen_adjacency_matrix()
    i, j = labels.index("H2"), labels.index("O2")
    assert H.sum() == 2 and H[i, j] == H[j, i] == 1
