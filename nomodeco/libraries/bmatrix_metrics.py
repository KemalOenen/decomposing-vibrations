import numpy as np


def condition_metrics(B, idof, rtol=1e-8):
    """Rank, completeness and condition number of the B matrix."""
    s = np.linalg.svd(B, compute_uv=False)
    rank = int(np.sum(s > rtol * s[0]))
    sigma_min_vib = s[idof - 1] if len(s) >= idof else 0.0
    s_next = s[idof] if len(s) > idof else 0.0
    return {
        "shape": B.shape,
        "idof": idof,
        "numerical_rank": rank,
        "n_redundant": B.shape[0] - rank,
        "complete": rank >= idof,
        "sigma_max": s[0],
        "sigma_min_vib": sigma_min_vib,
        "cond_vib": s[0] / sigma_min_vib if sigma_min_vib > 0 else np.inf,
        # large gap -> ICs cleanly span the vibrational space
        "spectral_gap": sigma_min_vib / s_next if s_next > 0 else np.inf,
        "singular_values": s,
    }


def rigid_body_residuals(B, coords):
    """|B @ t| for translations and rotations; should vanish for a correct B matrix."""
    coords = np.asarray(coords, dtype=float)
    N = coords.shape[0]
    r = coords - coords.mean(axis=0)
    vecs = []
    for k in range(3):  # translations
        t = np.zeros((N, 3))
        t[:, k] = 1.0
        vecs.append(t.ravel())
    for k in range(3):  # infinitesimal rotations
        axis = np.zeros(3)
        axis[k] = 1.0
        vecs.append(np.cross(axis, r).ravel())
    T = np.array(vecs).T  # (3N, 6)
    T /= np.linalg.norm(T, axis=0, keepdims=True) + 1e-300
    R = np.abs(B @ T)  # (n_int, 6)
    return {
        "max_trans_residual": R[:, :3].max(),
        "max_rot_residual": R[:, 3:].max(),
        "worst_ic_trans": int(R[:, :3].max(axis=1).argmax()),
        "worst_ic_rot": int(R[:, 3:].max(axis=1).argmax()),
    }


def row_collinearity(B, labels, thresh=0.95):
    """Cosine similarity between B rows; reports IC pairs with |cos| > thresh."""
    norms = np.linalg.norm(B, axis=1)
    Bn = B / np.where(norms > 0, norms, 1.0)[:, None]
    C = Bn @ Bn.T
    iu = np.triu_indices_from(C, k=1)
    mask = np.abs(C[iu]) > thresh
    pairs = sorted(
        ((labels[i], labels[j], C[i, j]) for i, j in zip(iu[0][mask], iu[1][mask])),
        key=lambda x: -abs(x[2]),
    )
    return {
        "collinear_pairs": pairs,
        "max_offdiag_cos": np.abs(C[iu]).max() if len(iu[0]) else 0.0,
    }


def bmatrix_metrics(B, idof, coords=None, ic_labels=None, rtol=1e-8, collinear_thresh=0.95):
    """Run all B matrix checks. Rigid body residuals are skipped if coords is None."""
    B = np.asarray(B, dtype=float)
    labels = ic_labels if ic_labels is not None else [str(i) for i in range(B.shape[0])]

    m = condition_metrics(B, idof, rtol)
    m.update(row_collinearity(B, labels, collinear_thresh))
    if coords is not None:
        m.update(rigid_body_residuals(B, coords))
    return m

if __name__ == "__main__":
    # Water: r1, r2, angle + a duplicated angle (should show up as redundant / collinear)
    coords = np.array([[0.0, 0.0, 0.117], [0.0, 0.757, -0.468], [0.0, -0.757, -0.468]])

    def internals(x):
        x = x.reshape(-1, 3)
        a, b = x[1] - x[0], x[2] - x[0]
        ang = np.arccos(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
        return np.array([np.linalg.norm(a), np.linalg.norm(b), ang, ang])

    # B matrix by central finite differences
    h = 1e-6
    x0 = coords.ravel()
    B = np.array([(internals(x0 + h * e) - internals(x0 - h * e)) / (2 * h)
                  for e in np.eye(x0.size)]).T

    metrics = bmatrix_metrics(B, idof=3, coords=coords, ic_labels=["r1", "r2", "a", "a_dup"])
    for key, val in metrics.items():
        print(f"{key:<20}: {val}")