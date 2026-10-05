import numpy as np
from numba import njit


@njit(cache=True)
def bond_length(a, b):
    d = b - a
    return np.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])


@njit(cache=True)
def normalized_bond_vector(a, b):
    d = b - a
    return d / np.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])


@njit(cache=True)
def bond_angle(a, b, c):
    ba = a - b
    bc = c - b
    cos_a = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    if cos_a > 1.0:
        cos_a = 1.0
    elif cos_a < -1.0:
        cos_a = -1.0
    return np.arccos(cos_a)


@njit(cache=True)
def torsion_angle(a, b, c, d):
    a1 = bond_angle(a, b, c)
    a2 = bond_angle(b, c, d)
    eba = normalized_bond_vector(b, a)
    edc = normalized_bond_vector(d, c)
    cos_t = (np.cos(a1) * np.cos(a2) - np.dot(eba, edc)) / (np.sin(a1) * np.sin(a2))
    if cos_t > 1.0:
        cos_t = 1.0
    elif cos_t < -1.0:
        cos_t = -1.0
    return np.arccos(cos_t)


@njit(cache=True)
def B_Matrix_Entry_BondLength(a, b):
    return normalized_bond_vector(a, b)


@njit(cache=True)
def B_Matrix_Entry_Angle_AtomB(A, B, C):
    ang = bond_angle(C, A, B)
    eAB = normalized_bond_vector(A, B)
    eAC = normalized_bond_vector(A, C)
    return (eAB * np.cos(ang) - eAC) / (bond_length(A, B) * np.sin(ang))


@njit(cache=True)
def B_Matrix_Entry_Angle_AtomC(A, B, C):
    ang = bond_angle(C, A, B)
    eAB = normalized_bond_vector(A, B)
    eAC = normalized_bond_vector(A, C)
    return (eAC * np.cos(ang) - eAB) / (bond_length(A, C) * np.sin(ang))


@njit(cache=True)
def B_Matrix_Entry_Angle_AtomA(A, B, C):
    return -(B_Matrix_Entry_Angle_AtomB(A, B, C) + B_Matrix_Entry_Angle_AtomC(A, B, C))


@njit(cache=True)
def _linear_bend_w(u, v):
    # Bend direction w (Bakken & Helgaker, JCP 117, 9160 (2002), Eq. 24).
    # w = u x v when the angle is not exactly linear: only a w perpendicular to both
    # u and v gives a rotation-invariant row. At 180 deg, cross u with a reference
    # vector; take whichever of [1,-1,1], [-1,1,1] is less parallel to u.
    w = np.cross(u, v)
    if np.dot(w, w) > 1e-12:
        return w / np.sqrt(np.dot(w, w))
    w1 = np.cross(u, np.array([1.0, -1.0, 1.0]))
    w2 = np.cross(u, np.array([-1.0, 1.0, 1.0]))
    if np.dot(w1, w1) < np.dot(w2, w2):
        w1 = w2
    return w1 / np.sqrt(np.dot(w1, w1))


@njit(cache=True)
def B_Matrix_Entry_LinearBend(M, O, N, second):
    # Linear bend M-O-N with center O (Bakken & Helgaker Eq. 25); returns (dq/dM, dq/dO, dq/dN).
    # The second bend uses w2 = u x w1, orthogonal to the first (already unit length).
    u = normalized_bond_vector(O, M)
    v = normalized_bond_vector(O, N)
    w = _linear_bend_w(u, v)
    if second:
        w = np.cross(u, w)
    dm = np.cross(u, w) / bond_length(O, M)
    dn = np.cross(w, v) / bond_length(O, N)
    return dm, -(dm + dn), dn


@njit(cache=True)
def B_Matrix_Entry_LinearBend_Ref(M, O, N, w):
    # Linear bend M-O-N with a given reference w (from linear_bends.linear_bend_references), Eq. 25
    u = normalized_bond_vector(O, M)
    v = normalized_bond_vector(O, N)
    dm = np.cross(u, w) / bond_length(O, M)
    dn = np.cross(w, v) / bond_length(O, N)
    return dm, -(dm + dn), dn


@njit(cache=True)
def B_Matrix_Entry_Torsion_AtomB(A, B, C, D):
    eAC = normalized_bond_vector(A, C)
    eBA = normalized_bond_vector(B, A)
    sin_BAC = np.sin(bond_angle(B, A, C))
    return np.cross(eAC, eBA) / (bond_length(A, B) * sin_BAC * sin_BAC)


@njit(cache=True)
def B_Matrix_Entry_Torsion_AtomD(A, B, C, D):
    eCA = normalized_bond_vector(C, A)
    eDC = normalized_bond_vector(D, C)
    sin_ACD = np.sin(bond_angle(A, C, D))
    return np.cross(eCA, eDC) / (bond_length(C, D) * sin_ACD * sin_ACD)


@njit(cache=True)
def B_Matrix_Entry_Torsion_AtomA(A, B, C, D):
    # Deduplicated: bond_angle(B,A,C) x3 -> once, bond_angle(A,C,D) x2 -> once,
    # cross(eBA,eAC) x2 -> once, cross(eDC,eCA) x2 -> once
    eBA = normalized_bond_vector(B, A)
    eAC = normalized_bond_vector(A, C)
    eDC = normalized_bond_vector(D, C)
    eCA = normalized_bond_vector(C, A)
    r_AB = bond_length(A, B)
    r_AC = bond_length(A, C)
    a_BAC = bond_angle(B, A, C)
    a_ACD = bond_angle(A, C, D)
    sin2_BAC = np.sin(a_BAC) ** 2
    cos_BAC = np.cos(a_BAC)
    sin2_ACD = np.sin(a_ACD) ** 2
    cos_ACD = np.cos(a_ACD)
    cross_BAC = np.cross(eBA, eAC)
    cross_DCA = np.cross(eDC, eCA)
    return (
        cross_BAC / (r_AB * sin2_BAC)
        - (cos_BAC / (r_AC * sin2_BAC)) * cross_BAC
        + (cos_ACD / (r_AC * sin2_ACD)) * cross_DCA  # Wilson s_2: + cos(phi3)/r23 * Y
    )


@njit(cache=True)
def B_Matrix_Entry_Torsion_AtomC(A, B, C, D):
    # Deduplicated: bond_angle(A,C,D) x3 -> once, bond_angle(B,A,C) x2 -> once,
    # cross(eDC,eCA) x2 -> once, cross(eBA,eAC) x2 -> once
    eDC = normalized_bond_vector(D, C)
    eCA = normalized_bond_vector(C, A)
    eBA = normalized_bond_vector(B, A)
    eAC = normalized_bond_vector(A, C)
    r_CD = bond_length(C, D)
    r_CA = bond_length(C, A)
    a_ACD = bond_angle(A, C, D)
    a_BAC = bond_angle(B, A, C)
    sin2_ACD = np.sin(a_ACD) ** 2
    cos_ACD = np.cos(a_ACD)
    sin2_BAC = np.sin(a_BAC) ** 2
    cos_BAC = np.cos(a_BAC)
    cross_DCA = np.cross(eDC, eCA)
    cross_BAC = np.cross(eBA, eAC)
    return (
        cross_DCA / (r_CD * sin2_ACD)
        - (cos_ACD / (r_CA * sin2_ACD)) * cross_DCA
        + (cos_BAC / (r_CA * sin2_BAC)) * cross_BAC  # Wilson s_3: + cos(phi2)/r23 * X
    )


@njit(cache=True)
def B_Matrix_Entry_OutOfPlane_AtomB(A, B, C, D):
    r_ab = bond_length(A, B)
    e_ab = normalized_bond_vector(A, B)
    e_ac = normalized_bond_vector(A, C)
    e_ad = normalized_bond_vector(A, D)
    phi_b = bond_angle(C, A, D)
    sin_phi_b = np.sin(phi_b)
    cross_cd = np.cross(e_ac, e_ad)
    sin_theta = np.dot(e_ab, cross_cd / sin_phi_b)
    if sin_theta > 1.0:
        sin_theta = 1.0
    elif sin_theta < -1.0:  # theta is signed; clamping to 0 broke all theta < 0 permutations
        sin_theta = -1.0
    theta = np.arcsin(sin_theta)
    if abs(theta) < 1e-10:
        return cross_cd / (sin_phi_b * r_ab)
    return (1.0 / r_ab) * (cross_cd / (np.cos(theta) * sin_phi_b) - np.tan(theta) * e_ab)


@njit(cache=True)
def B_Matrix_Entry_OutOfPlane_AtomC(A, B, C, D):
    r_ac = bond_length(A, C)
    e_ab = normalized_bond_vector(A, B)
    e_ac = normalized_bond_vector(A, C)
    e_ad = normalized_bond_vector(A, D)
    phi_b = bond_angle(C, A, D)
    phi_c = bond_angle(B, A, D)
    phi_d = bond_angle(B, A, C)
    sin_phi_b = np.sin(phi_b)
    cross_cd = np.cross(e_ac, e_ad)
    sin_theta = np.dot(e_ab, cross_cd / sin_phi_b)
    if sin_theta > 1.0:
        sin_theta = 1.0
    elif sin_theta < -1.0:  # theta is signed; clamping to 0 broke all theta < 0 permutations
        sin_theta = -1.0
    theta = np.arcsin(sin_theta)
    if abs(theta) < 1e-10:
        return (1.0 / r_ac) * (cross_cd / sin_phi_b) * (np.sin(phi_c) / sin_phi_b)
    return (1.0 / r_ac) * (cross_cd / sin_phi_b) * (
        (np.cos(phi_b) * np.cos(phi_c) - np.cos(phi_d)) / (np.cos(theta) * sin_phi_b ** 2)
    )


@njit(cache=True)
def B_Matrix_Entry_OutOfPlane_AtomD(A, B, C, D):
    r_ad = bond_length(A, D)
    e_ab = normalized_bond_vector(A, B)
    e_ac = normalized_bond_vector(A, C)
    e_ad = normalized_bond_vector(A, D)
    phi_b = bond_angle(C, A, D)
    phi_c = bond_angle(B, A, D)
    phi_d = bond_angle(B, A, C)
    sin_phi_b = np.sin(phi_b)
    cross_cd = np.cross(e_ac, e_ad)
    sin_theta = np.dot(e_ab, cross_cd / sin_phi_b)
    if sin_theta > 1.0:
        sin_theta = 1.0
    elif sin_theta < -1.0:  # theta is signed; clamping to 0 broke all theta < 0 permutations
        sin_theta = -1.0
    theta = np.arcsin(sin_theta)
    if abs(theta) < 1e-10:
        return (1.0 / r_ad) * (cross_cd / sin_phi_b) * (np.sin(phi_d) / sin_phi_b)
    return (1.0 / r_ad) * (cross_cd / sin_phi_b) * (
        (np.cos(phi_b) * np.cos(phi_d) - np.cos(phi_c)) / (np.cos(theta) * sin_phi_b ** 2)
    )


@njit(cache=True)
def B_Matrix_Entry_OutOfPlane_AtomA(A, B, C, D):
    return -(
        B_Matrix_Entry_OutOfPlane_AtomB(A, B, C, D)
        + B_Matrix_Entry_OutOfPlane_AtomC(A, B, C, D)
        + B_Matrix_Entry_OutOfPlane_AtomD(A, B, C, D)
    )


def b_matrix(atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof, bend_refs=None) -> np.ndarray:
    # bend_refs: {(M, O, N): (w1, w2)} from linear_bends.linear_bend_references; triples without
    # an entry use the geometric rule in B_Matrix_Entry_LinearBend
    n_atoms = len(atoms)
    coordinates = np.array([a.coordinates for a in atoms])
    atom_index = {a.symbol: i for i, a in enumerate(atoms)}
    n_internal = (
        len(bonds) + len(angles) + len(linear_angles) + len(out_of_plane) + len(dihedrals)
    )
    assert n_internal >= idof, (
        f"Wrong number of internal coordinates, n_internal ({n_internal}) should be >= {idof}."
    )
    matrix = np.zeros((n_internal, 3 * n_atoms))
    i_internal = 0
    # each linear triple appears twice (two orthogonal bends): first occurrence -> first bend,
    # later occurrence -> second bend; M-O-N and N-O-M are the same triple
    seen_linear_angles = set()

    for bond in bonds:
        index = [atom_index[a] * 3 for a in bond]
        coord = [coordinates[atom_index[a]] for a in bond]
        matrix[i_internal, index[0]:index[0]+3] = B_Matrix_Entry_BondLength(coord[1], coord[0])
        matrix[i_internal, index[1]:index[1]+3] = B_Matrix_Entry_BondLength(coord[0], coord[1])
        i_internal += 1

    for angle in angles:
        index = [atom_index[a] * 3 for a in angle]
        coord = [coordinates[atom_index[a]] for a in angle]
        matrix[i_internal, index[0]:index[0]+3] = B_Matrix_Entry_Angle_AtomB(coord[1], coord[0], coord[2])
        matrix[i_internal, index[1]:index[1]+3] = B_Matrix_Entry_Angle_AtomA(coord[1], coord[0], coord[2])
        matrix[i_internal, index[2]:index[2]+3] = B_Matrix_Entry_Angle_AtomC(coord[1], coord[0], coord[2])
        i_internal += 1

    for linear_angle in linear_angles:
        index = [atom_index[a] * 3 for a in linear_angle]
        coord = [coordinates[atom_index[a]] for a in linear_angle]
        key = min(tuple(linear_angle), tuple(linear_angle)[::-1])
        second = key in seen_linear_angles
        seen_linear_angles.add(key)
        # normal-mode aligned (w1, w2) if given; a reversed triple only flips the row's sign
        refs = None
        if bend_refs:
            refs = bend_refs.get(tuple(linear_angle)) or bend_refs.get(tuple(linear_angle)[::-1])
        if refs is not None:
            dm, do, dn = B_Matrix_Entry_LinearBend_Ref(coord[0], coord[1], coord[2], refs[1 if second else 0])
        else:
            dm, do, dn = B_Matrix_Entry_LinearBend(coord[0], coord[1], coord[2], second)
        matrix[i_internal, index[0]:index[0]+3] = dm
        matrix[i_internal, index[1]:index[1]+3] = do
        matrix[i_internal, index[2]:index[2]+3] = dn
        i_internal += 1

    for outofplane in out_of_plane:
        index = [atom_index[a] * 3 for a in outofplane]
        coord = [coordinates[atom_index[a]] for a in outofplane]
        # oop tuples are (center, wing, c, d) (Molecule.generate_out_of_plane); the
        # B_Matrix_Entry_OutOfPlane_* functions take A = center, B = wing atom
        matrix[i_internal, index[0]:index[0]+3] = B_Matrix_Entry_OutOfPlane_AtomA(coord[0], coord[1], coord[2], coord[3])
        matrix[i_internal, index[1]:index[1]+3] = B_Matrix_Entry_OutOfPlane_AtomB(coord[0], coord[1], coord[2], coord[3])
        matrix[i_internal, index[2]:index[2]+3] = B_Matrix_Entry_OutOfPlane_AtomC(coord[0], coord[1], coord[2], coord[3])
        matrix[i_internal, index[3]:index[3]+3] = B_Matrix_Entry_OutOfPlane_AtomD(coord[0], coord[1], coord[2], coord[3])
        i_internal += 1

    for dihedral in dihedrals:
        index = [atom_index[a] * 3 for a in dihedral]
        coord = [coordinates[atom_index[a]] for a in dihedral]
        matrix[i_internal, index[0]:index[0]+3] = B_Matrix_Entry_Torsion_AtomB(coord[1], coord[0], coord[2], coord[3])
        matrix[i_internal, index[1]:index[1]+3] = B_Matrix_Entry_Torsion_AtomA(coord[1], coord[0], coord[2], coord[3])
        matrix[i_internal, index[2]:index[2]+3] = B_Matrix_Entry_Torsion_AtomC(coord[1], coord[0], coord[2], coord[3])
        matrix[i_internal, index[3]:index[3]+3] = B_Matrix_Entry_Torsion_AtomD(coord[1], coord[0], coord[2], coord[3])
        i_internal += 1

    return matrix


if __name__ == "__main__":
    def test_b_matrix_speed():
        from time import time
        from molecule_class import Molecule

        mol = Molecule.from_xyz_file('/home/lme/decomposing-vibrations/test_calculations/xyz_tests/1citric_001.xyz')
        degofc = mol.degree_of_covalance()
        bonds = mol.covalent_bonds(degofc)
        angles, linear_angles = mol.generate_angles(bonds)
        dihedrals = mol.generate_dihedrals(bonds)
        dof = mol.idof_general()

        start_time = time()
        b_mat = b_matrix(mol, bonds, angles, linear_angles, [], dihedrals, dof)
        end_time = time()
        # First time with numba is the compilation time, subsequent calls should be much faster
        print(f"B-Matrix calculation time: {end_time - start_time:.4f} seconds")

    test_b_matrix_speed()
    test_b_matrix_speed()