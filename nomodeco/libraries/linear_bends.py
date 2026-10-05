"""
Reference vectors for linear bends, aligned with the normal modes.

Runs once per molecule, before any B matrix is built; every b_matrix call reuses its result.
Bend vector of mode k for the unit M-O-N (center O):
    beta_k = P [(d_M - d_O)/lambda_M + (d_N - d_O)/lambda_N],  P = I - u u^T
The D element of a bend with reference w is (u x w) . beta_k.
"""

import numpy as np

HARTREE_BOHR_AMU_TO_CM1 = 5140.4981


def _bend_vector(xyz, triple, d):
    """beta for one Cartesian displacement d (N x 3) of the unit M-O-N."""
    M, O, N = triple
    lam_M = np.linalg.norm(xyz[M] - xyz[O])
    lam_N = np.linalg.norm(xyz[N] - xyz[O])
    u = (xyz[M] - xyz[O]) / lam_M
    beta = (d[M] - d[O]) / lam_M + (d[N] - d[O]) / lam_N
    return beta - np.dot(beta, u) * u


def _same_axis(xyz, t1, t2):
    u1 = xyz[t1[0]] - xyz[t1[1]]
    u2 = xyz[t2[0]] - xyz[t2[1]]
    offset = xyz[t2[1]] - xyz[t1[1]]
    return (np.linalg.norm(np.cross(u1, u2)) < 1e-6 * np.linalg.norm(u1) * np.linalg.norm(u2)
            and np.linalg.norm(np.cross(offset, u1)) < 1e-4 * np.linalg.norm(u1))


def linear_bend_references(atoms, linear_angles, L, eigenvalues, diag_reciprocal_square, idof,
                           tol_cm1=0.5):
    """
    Returns (L, refs):
        L     eigenvectors with every degenerate vibrational pair rotated onto t
        refs  {(M, O, N) label triple: (w1, w2)} for every linear triple
    """
    L = L.copy()
    xyz = np.array([a.coordinates for a in atoms], dtype=float)
    index = {a.symbol: i for i, a in enumerate(atoms)}
    triples = [tuple(index[a] for a in t) for t in dict.fromkeys(map(tuple, linear_angles))]
    if not triples:
        return L, {}

    # Step 1: vibrational modes as Cartesian displacements
    vib = range(L.shape[0] - idof, L.shape[0])

    def displacement(k):
        return (diag_reciprocal_square * L[:, k]).reshape(-1, 3)

    # Step 4: one direction t per axis, from the mode and triple with the largest |beta|
    axes = []
    for t in triples:
        for axis in axes:
            if _same_axis(xyz, axis[0], t):
                axis.append(t)
                break
        else:
            axes.append([t])
    t_of = {}
    for axis in axes:
        betas = [_bend_vector(xyz, t, displacement(k)) for t in axis for k in vib]
        best = max(betas, key=np.linalg.norm)
        for t in axis:
            t_of[t] = best / np.linalg.norm(best)

    # Step 5: rotate degenerate pairs so that beta_a || t and beta_b _|_ t
    freqs = np.sqrt(np.abs(eigenvalues)) * HARTREE_BOHR_AMU_TO_CM1
    for a in vib[:-1]:
        b = a + 1
        if abs(freqs[a] - freqs[b]) >= tol_cm1:
            continue
        def pair_weight(t):
            return (np.linalg.norm(_bend_vector(xyz, t, displacement(a))) ** 2
                    + np.linalg.norm(_bend_vector(xyz, t, displacement(b))) ** 2)
        triple = max(triples, key=pair_weight)
        t = t_of[triple]
        alpha = np.arctan2(t @ _bend_vector(xyz, triple, displacement(b)),
                           t @ _bend_vector(xyz, triple, displacement(a)))
        L_a, L_b = L[:, a].copy(), L[:, b].copy()
        L[:, a] = np.cos(alpha) * L_a + np.sin(alpha) * L_b
        L[:, b] = -np.sin(alpha) * L_a + np.cos(alpha) * L_b

    # Step 6: reference vectors w1 = t x u, w2 = t
    labels = [a.symbol for a in atoms]
    refs = {}
    for triple, t in t_of.items():
        M, O, _ = triple
        u = (xyz[M] - xyz[O]) / np.linalg.norm(xyz[M] - xyz[O])
        w1 = np.cross(t, u)
        refs[tuple(labels[i] for i in triple)] = (w1 / np.linalg.norm(w1), t)
    return L, refs
