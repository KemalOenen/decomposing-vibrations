"""
Completeness check for internal coordinate sets.

A set is complete if rank(B) = idof. Methods:
    1. svd       rank of B via singular values                (reference, slow)
    2. gram      smallest eigenvalue of (BQ)^T (BQ)            (per set)
    3. product   all sets of a LazyIcDict at once              (fast)
"""

import numpy as np

from nomodeco.libraries import bmatrix, icset_opt

KINDS = ["bonds", "angles", "linear valence angles", "out of plane angles", "dihedrals"]

RANK_TOL = 1e-8      # singular values below RANK_TOL * largest count as zero
GRAM_TOL = 1e-12     # eigenvalues of a Gram matrix are squared singular values


def numerical_rank(singular_values, tol=RANK_TOL):
    """  
    Number of singular vlaues above tol * largest singular value.
    """
    if len(singular_values) == 0 or singular_values[0] == 0:
        return 0
    return int(np.sum(singular_values > tol * singular_values[0]))

def shape_space_basis(atoms):
    """
    Q (3N x idof): orthonormal basis of all displacements that are not
    translations or rotations. Every B row lies in this space.
    """
    coords = np.array([atom.coordinates for atom in atoms], dtype=float)
    coords = coords - coords.mean(axis=0)
    n_atoms = len(atoms)

    # the 6 rigid-body motions as columns: 3 translations, 3 rotations
    rigid = np.zeros((3 * n_atoms, 6))
    for i in range(n_atoms):
        x, y, z = coords[i]
        rigid[3 * i:3 * i + 3, 0] = [1, 0, 0]
        rigid[3 * i:3 * i + 3, 1] = [0, 1, 0]
        rigid[3 * i:3 * i + 3, 2] = [0, 0, 1]
        rigid[3 * i:3 * i + 3, 3] = [0, -z, y]
        rigid[3 * i:3 * i + 3, 4] = [z, 0, -x]
        rigid[3 * i:3 * i + 3, 5] = [-y, x, 0]

    U, s, _ = np.linalg.svd(rigid)
    n_rigid = numerical_rank(s)          # 6, or 5 for a linear molecule
    return U[:, n_rigid:]                # the remaining columns span the shape space


def set_b_matrix(atoms, ic_set, bend_refs=None):
    """B matrix of one IC set."""
    return bmatrix.b_matrix(atoms, *[ic_set[k] for k in KINDS], 0, bend_refs)

# Method 1

def is_complete_svd(B, idof):
    """ 
    Check the completeness of a set of internal coordinates
    by computing the singular values
    """ 
    singular_values = np.linalg.svd(B, compute_uv=False)
    return numerical_rank(singular_values) == idof

# Method 2

def is_complete_gram(B, idof, Q):
    """  
    M = (BQ)^T (BQ) is a matrix of idof x idof. The set is complete
    if M is not singular
    """
    Bv = B @ Q
    eigenvalues = np.linalg.eigvalsh(Bv.T @ Bv)
    return eigenvalues[0] > GRAM_TOL * eigenvalues[-1]


# Method 3


def rows_for(atoms, Q, ics, kind):
    """B rows (projected onto the shape space) for a list of ICs of one kind, as a dict ic -> row."""
    lists = {k: [] for k in KINDS}
    lists[kind] = ics
    B = bmatrix.b_matrix(atoms, *[lists[k] for k in KINDS], 0)
    return {tuple(ic): row for ic, row in zip(ics, B @ Q)}


def unique(subsets):
    """All different ICs that occur in a list of subsets."""
    seen = {}
    for subset in subsets:
        for ic in subset:
            seen[tuple(ic)] = True
    return list(seen)


def subsets_by_size(subsets, index):
    """Subsets grouped by size: size -> (positions in subsets, row indices (n_subsets x size))."""
    groups = {}
    for i, subset in enumerate(subsets):
        groups.setdefault(len(subset), []).append((i, [index[tuple(ic)] for ic in subset]))
    return {
        size: (np.array([i for i, _ in g]), np.array([rows for _, rows in g], dtype=np.intp).reshape(len(g), size))
        for size, g in groups.items()
    }


def dihedrals_complete(R, groups, n_subsets, missing, scale=1.0):
    """
    For every dihedral subset: does it add exactly `missing` new dimensions?
    R = dihedral rows @ N, N (idof x missing) an orthonormal basis of the directions the prefix
    does not cover, so the subset matrix is only (size x missing). A direction counts if its
    singular value is above RANK_TOL * scale (scale = largest singular value of the prefix),
    the criterion of numerical_rank. One batched SVD per subset size.
    """
    result = np.zeros(n_subsets, dtype=bool)
    for size, (positions, rows) in groups.items():
        if size < missing:
            continue                                     # can never be enough
        s = np.linalg.svd(R[rows], compute_uv=False)     # (n_subsets, missing), descending
        result[positions] = s[:, missing - 1] > RANK_TOL * scale
    return result


def complete_mask_product(ic_dict, atoms, idof, bend_refs=None):
    """
    True/False for every set of a LazyIcDict, in the order of ic_dict[k].
    A set = prefix (bonds, linear angles, angle subset, oop subset) + dihedral subset.
    The prefix rank is computed once, the dihedral subsets only fill the rest.
    """
    Q = shape_space_basis(atoms)
    masks = []

    for bonds, linear_angles, angle_subsets, oop_subsets, dihedral_subsets in ic_dict._parts:
        # B rows of all ICs in this block, once
        base = bmatrix.b_matrix(atoms, bonds, [], linear_angles, [], [], 0, bend_refs) @ Q
        angle_rows = rows_for(atoms, Q, unique(angle_subsets), "angles")
        oop_rows = rows_for(atoms, Q, unique(oop_subsets), "out of plane angles")
        dihedral_rows = rows_for(atoms, Q, unique(dihedral_subsets), "dihedrals")
        D = np.array(list(dihedral_rows.values())).reshape(-1, idof)
        groups = subsets_by_size(dihedral_subsets, {ic: i for i, ic in enumerate(dihedral_rows)})

        mask = np.zeros((len(angle_subsets), len(oop_subsets), len(dihedral_subsets)), dtype=bool)

        for ia, angle_subset in enumerate(angle_subsets):
            for io, oop_subset in enumerate(oop_subsets):
                prefix = np.vstack(
                    [base]
                    + [angle_rows[tuple(ic)][None, :] for ic in angle_subset]
                    + [oop_rows[tuple(ic)][None, :] for ic in oop_subset]
                )

                # rank of the prefix; Vt (idof x idof) holds the covered directions first
                _, s, Vt = np.linalg.svd(prefix)
                rank = numerical_rank(s)
                missing = idof - rank

                if missing == 0:
                    mask[ia, io, :] = True
                    continue

                # orthonormal basis of the directions the prefix does NOT cover yet
                N = Vt[rank:].T

                mask[ia, io, :] = dihedrals_complete(
                    D @ N, groups, len(dihedral_subsets), missing, s[0]
                )

        masks.append(mask.ravel())       # same order as LazyIcDict: angle, oop, dihedral

    return np.concatenate(masks)

def complete_keys(ic_dict, atoms, idof, bend_refs=None, method="product"):
    """Keys of all complete sets."""
    if method == "off":
        return np.arange(len(ic_dict))

    if method == "product":
        if hasattr(ic_dict, "_parts"):
            return np.flatnonzero(complete_mask_product(ic_dict, atoms, idof, bend_refs))
        # plain dicts have no product structure: the threaded batched SVD is the fastest there
        return np.array(icset_opt.complete_set_keys(ic_dict, atoms, idof, bend_refs), dtype=int)

    # per set: svd or gram
    Q = shape_space_basis(atoms)
    keys = []
    for k in range(len(ic_dict)):
        B = set_b_matrix(atoms, ic_dict[k], bend_refs)
        if method == "svd":
            ok = is_complete_svd(B, idof)
        else:
            ok = is_complete_gram(B, idof, Q)
        if ok:
            keys.append(k)
    return np.array(keys, dtype=int)