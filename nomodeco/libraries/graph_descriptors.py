from __future__ import annotations
from collections import Counter

import networkx as nx
import numpy as np

HYDROGENS = {"H", "D"}
_TOL = 1e-10

def graph_descriptors(atoms) -> dict:
    """ 
    Computes several graph descriptors from a given atoms list
    """
    G = atoms.graph()
    C = atoms.bond_graph("cov")
    # Compute order of the graph
    order = [atom.symbol for atom in atoms]

    # all edge kinds: C only has the covalent ones
    kinds = Counter(d["kind"] for _, _, d in G.edges(data=True))
    degrees = dict(C.degree())
    n_comp = nx.number_connected_components(C)
    rings = nx.minimum_cycle_basis(C)

    lap = np.sort(nx.laplacian_spectrum(C))
    lap[np.abs(lap) < _TOL] = 0.0

    d = {
        "n_atoms": C.number_of_nodes(),
        "n_bonds": {"cov": kinds.get("cov", 0), "h_bond": kinds.get("h_bond", 0),
                    "acc_don": kinds.get("acc_don", 0)},
        "n_components": n_comp,
        # cyclomatic number b - a + c = number of independent rings
        "cyclomatic_number": C.number_of_edges() - C.number_of_nodes() + n_comp,
        "ring_sizes": sorted(len(r) for r in rings),
        "rings": [sorted(r, key=order.index) for r in rings],
        "degree_histogram": dict(sorted(Counter(degrees.values()).items())),
        "n_terminal": sum(1 for v in degrees.values() if v == 1),
        "n_branching": sum(1 for v in degrees.values() if v >= 3),
        "bridges": len(list(nx.bridges(C))),
        "articulation_atoms": sorted(nx.articulation_points(C), key=order.index),
        "laplacian_spectrum": lap,
        "n_zero_laplacian": int(np.sum(lap == 0.0)),
        "submolecules": [],
        "atoms": _atom_table(C, order),
    }
    # hydrogen-bond network: is the cluster held together by its h-bonds?
    if kinds.get("h_bond", 0):
        H = atoms.bond_graph("cov", "h_bond")
        d["n_components_with_h_bonds"] = nx.number_connected_components(H)

    for comp in sorted(nx.connected_components(C), key=lambda c: min(map(order.index, c))):
        d["submolecules"].append(_submolecule_descriptors(C.subgraph(comp), order))
    return d

def _submolecule_descriptors(S: nx.Graph, order) -> dict:
    atoms = sorted(S.nodes, key=order.index)
    heavy = S.subgraph(n for n in S if S.nodes[n]["element"] not in HYDROGENS)
    out = {
        "atoms": atoms,
        "formula": _formula(S),
        "full": _indices(S),
        "heavy": _indices(heavy) if heavy.number_of_nodes() > 1 and nx.is_connected(heavy) else None,
    }
    return out

def _indices(S: nx.Graph) -> dict:
    """Topological and spectral indices of a connected graph."""
    n, m = S.number_of_nodes(), S.number_of_edges()
    deg = dict(S.degree())
    res = {"n": n, "m": m}
    if n < 2:
        return res

    dist = dict(nx.all_pairs_shortest_path_length(S))
    wiener = sum(dist[u][v] for u in S for v in S) / 2
    lap = np.sort(nx.laplacian_spectrum(S))
    adj = np.linalg.eigvalsh(nx.to_numpy_array(S))
    mu_cyc = m - n + 1

    res.update({
        "diameter": nx.diameter(S),
        "radius": nx.radius(S),
        "wiener": wiener,
        "avg_path_length": wiener / (n * (n - 1) / 2),
        # Randic connectivity index chi = sum 1/sqrt(d_u d_v)
        "randic": sum(1 / np.sqrt(deg[u] * deg[v]) for u, v in S.edges),
        # Zagreb indices M1 = sum d^2, M2 = sum d_u d_v
        "zagreb_m1": sum(v ** 2 for v in deg.values()),
        "zagreb_m2": sum(deg[u] * deg[v] for u, v in S.edges),
        # Balaban J = m/(mu+1) * sum 1/sqrt(D_u D_v), D = distance sums
        "balaban_j": _balaban(S, dist, m, mu_cyc),
        "algebraic_connectivity": lap[1],
        # Kirchhoff index Kf = n * sum 1/lambda_i (resistance distance)
        "kirchhoff": n * float(np.sum(1 / lap[1:])),
        "spectral_radius": float(adj[-1]),
        # graph energy E = sum |lambda_i(A)|
        "graph_energy": float(np.sum(np.abs(adj))),
        "estrada": float(np.sum(np.exp(adj))),
    })
    return res


def _balaban(S, dist, m, mu_cyc):
    dsum = {u: sum(dist[u].values()) for u in S}
    return m / (mu_cyc + 1) * sum(1 / np.sqrt(dsum[u] * dsum[v]) for u, v in S.edges)


def _atom_table(C: nx.Graph, order) -> list:
    """Per-atom: degree, ring membership, eccentricity and centralities (within its submolecule)."""
    in_ring = set().union(*nx.cycle_basis(C)) if C.number_of_edges() else set()
    rows = []
    for comp in nx.connected_components(C):
        S = C.subgraph(comp)
        ecc = nx.eccentricity(S) if S.number_of_nodes() > 1 else {next(iter(S)): 0}
        btw = nx.betweenness_centrality(S, normalized=True)
        clo = nx.closeness_centrality(S)
        for a in S:
            rows.append((a, C.degree[a], a in in_ring, ecc[a], btw[a], clo[a]))
    return sorted(rows, key=lambda r: order.index(r[0]))


def _formula(S):
    # Hill order: C, then H, then alphabetical; without carbon everything alphabetical
    counts = Counter(S.nodes[n]["element"] for n in S)
    first = ["C", "H"] if "C" in counts else []
    elements = first + sorted(el for el in counts if el not in first)
    return "".join(f"{el}{counts[el] if counts[el] > 1 else ''}" for el in elements if el in counts)
