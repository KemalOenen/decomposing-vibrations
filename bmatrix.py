from cmath import isclose
import numpy as np
import scipy

# General Chemistry-related functions

# Note: Vector points from Atom A to Atom B
def bond_length(Coordinates_AtomA,Coordinates_AtomB): #ANG
    return np.linalg.norm((Coordinates_AtomB - Coordinates_AtomA))

def normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB):
    return (Coordinates_AtomB - Coordinates_AtomA) / (bond_length(Coordinates_AtomA,Coordinates_AtomB))

# Angle for the Atoms in the conformation A-B-C. 
# Note: B is the central Atom, the bond vectors start from Atom B therefore
def bond_angle(Coordinates_AtomA, Coordinates_AtomB, Coordinates_AtomC): #RAD
    cosine_angle = np.clip((np.inner((Coordinates_AtomA - Coordinates_AtomB), 
                                     (Coordinates_AtomC - Coordinates_AtomB)))/
                           (bond_length(Coordinates_AtomA,Coordinates_AtomB) * 
                            (bond_length(Coordinates_AtomC,Coordinates_AtomB))), -1.0,1.0)
    return np.arccos(cosine_angle)

# Dihedral for the Atoms in the conformation 
def torsion_angle(Coordinates_AtomA, Coordinates_AtomB, Coordinates_AtomC, Coordinates_AtomD):
    cosine_torsion_angle = np.clip((np.cos(bond_angle(Coordinates_AtomA, 
                                                      Coordinates_AtomB, Coordinates_AtomC)) * 
                                    np.cos(bond_angle(Coordinates_AtomB, Coordinates_AtomC, Coordinates_AtomD)) 
                                    - np.dot(normalized_bond_vector(Coordinates_AtomB, Coordinates_AtomA),
                                             normalized_bond_vector(Coordinates_AtomD, Coordinates_AtomC)))/
                                   (np.sin(bond_angle(Coordinates_AtomA, Coordinates_AtomB, 
                                                      Coordinates_AtomC)) * 
                                    np.sin(bond_angle(Coordinates_AtomB, Coordinates_AtomC, 
                                                      Coordinates_AtomD))),-1.0,1.0)
    return np.arccos(cosine_torsion_angle)



# Functions for the construction of the B-Matrix
''''' 
Note: In this order, the vector points from Atom A to Atom B - that means: if you want to calculate the B-Matrix Entry for
Atom A then you need the NEGATIVE value of this function-call, for Atom B you just take the normal value
'''''
def B_Matrix_Entry_BondLength(Coordinates_AtomA,Coordinates_AtomB):
    return normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB)

''''' 
Note: For the Entries of the bending-part of the B-Matrix entry, we define the following geometry:
C-A-B or with vectors: C <- A -> BCoordinates_AtomC

Important the entries are different for the Central Atom and the Side Atoms

''''' 

def B_Matrix_Entry_Angle_AtomB(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC):
    return (normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB) * 
            np.cos(bond_angle(Coordinates_AtomC, Coordinates_AtomA, Coordinates_AtomB)) - 
            normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomC))/(
                bond_length(Coordinates_AtomA,Coordinates_AtomB) * 
                np.sin(bond_angle(Coordinates_AtomC, Coordinates_AtomA, Coordinates_AtomB)))

def B_Matrix_Entry_Angle_AtomC(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC):
    return (normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomC) * 
            np.cos(bond_angle(Coordinates_AtomC, Coordinates_AtomA, Coordinates_AtomB)) - 
            normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB))/(
                bond_length(Coordinates_AtomA,Coordinates_AtomC) * 
                np.sin(bond_angle(Coordinates_AtomC, Coordinates_AtomA, Coordinates_AtomB)))

def B_Matrix_Entry_Angle_AtomA(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC):
    return -(B_Matrix_Entry_Angle_AtomB(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC) + 
             B_Matrix_Entry_Angle_AtomC(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC))

''''' 
Note: For the Entries of the linear bending-part of the B-Matrix entries, we define the following geometry:
C-A-B or with vectors: C <- A -> B

There are two different linear angle bendings, i.e., one on the xz plane and one on the yz plane
The unitary vector u is obtained via rotating the vector A -> B counterclockweise by 90 degrees 
to the right and then dividing through
''''' 
#TODO: current linear valence angles only useful for degenerate linear valence angle modes - make more generic

def B_Matrix_Entry_LinearAngleFirstPlane_AtomB(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC):
    rotation_radians = np.pi/2
    rotation_axis = np.array([0,1,0])
    rotation_vector = rotation_radians * rotation_axis
    rotation = scipy.spatial.transform.Rotation.from_rotvec(rotation_vector)
    u = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB)
    u = rotation.apply(u)
    return -(u/bond_length(Coordinates_AtomA, Coordinates_AtomB))

def B_Matrix_Entry_LinearAngleFirstPlane_AtomC(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC):
    rotation_radians = np.pi/2
    rotation_axis = np.array([0,1,0])
    rotation_vector = rotation_radians * rotation_axis
    rotation = scipy.spatial.transform.Rotation.from_rotvec(rotation_vector)
    u = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB)
    u = rotation.apply(u)
    return -(u/bond_length(Coordinates_AtomA, Coordinates_AtomC))

def B_Matrix_Entry_LinearAngleFirstPlane_AtomA(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC):
    return -(B_Matrix_Entry_LinearAngleFirstPlane_AtomB(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC) +
    B_Matrix_Entry_LinearAngleFirstPlane_AtomC(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC))

def B_Matrix_Entry_LinearAngleSecondPlane_AtomB(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC):
    rotation_radians = np.pi/2
    rotation_axis = np.array([0,1,0])
    rotation_vector = rotation_radians * rotation_axis
    rotation = scipy.spatial.transform.Rotation.from_rotvec(rotation_vector)
    u = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB)
    u = rotation.apply(u)
    up = np.cross(normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB), u)
    return -(up/bond_length(Coordinates_AtomA, Coordinates_AtomB))

def B_Matrix_Entry_LinearAngleSecondPlane_AtomC(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC):
    rotation_radians = np.pi/2
    rotation_axis = np.array([0,1,0])
    rotation_vector = rotation_radians * rotation_axis
    rotation = scipy.spatial.transform.Rotation.from_rotvec(rotation_vector)
    u = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB)
    u = rotation.apply(u)
    up = np.cross(normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB), u)
    return -(up/bond_length(Coordinates_AtomA, Coordinates_AtomC))

def B_Matrix_Entry_LinearAngleSecondPlane_AtomA(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC):
    return -(B_Matrix_Entry_LinearAngleSecondPlane_AtomB(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC) +
    B_Matrix_Entry_LinearAngleSecondPlane_AtomC(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC))

''''' 
Note: For the Entries of the torsion-part of the B-Matrix entry, we define the following geometry:
B-A-C-D or with vectors: B <- A <-> C -> D

Important the entries are different for the 'Central Atoms' (A,C) and the 'Side Atoms' (B,D)

''''' 

#note that the Entries of the normalized bond vectors are simply the cross product here of the following kind:
# Vector from A to C (!) x Vector from B -> A(!)

def B_Matrix_Entry_Torsion_AtomB(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC,Coordinates_AtomD):
    return np.cross(normalized_bond_vector(Coordinates_AtomA,Coordinates_AtomC), 
                    normalized_bond_vector(Coordinates_AtomB,Coordinates_AtomA)) / (
                        bond_length(Coordinates_AtomA,Coordinates_AtomB) * 
                        np.square(np.sin(bond_angle(Coordinates_AtomB,Coordinates_AtomA,Coordinates_AtomC))))

#  Vector from C to A (!) x Vector from D -> C(!)
def B_Matrix_Entry_Torsion_AtomD(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC,Coordinates_AtomD):
    return np.cross(normalized_bond_vector(Coordinates_AtomC,Coordinates_AtomA), 
                    normalized_bond_vector(Coordinates_AtomD,Coordinates_AtomC)) / (
                        bond_length(Coordinates_AtomC,Coordinates_AtomD) * 
                        np.square(np.sin(bond_angle(Coordinates_AtomA,Coordinates_AtomC,Coordinates_AtomD))))


def B_Matrix_Entry_Torsion_AtomA(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC,Coordinates_AtomD):
    return (np.cross(normalized_bond_vector(Coordinates_AtomB, Coordinates_AtomA), 
                           normalized_bond_vector(Coordinates_AtomA,Coordinates_AtomC))/(bond_length(Coordinates_AtomA,Coordinates_AtomB) * 
              np.square(np.sin(bond_angle(Coordinates_AtomB,Coordinates_AtomA,Coordinates_AtomC))))) - (
                  (np.cos(bond_angle(Coordinates_AtomB,Coordinates_AtomA,Coordinates_AtomC)) / (
                      bond_length(Coordinates_AtomA,Coordinates_AtomC) * 
                      np.square(np.sin(bond_angle(Coordinates_AtomB,Coordinates_AtomA,Coordinates_AtomC))))) *
                  np.cross(normalized_bond_vector(Coordinates_AtomB, Coordinates_AtomA), 
                           normalized_bond_vector(Coordinates_AtomA,Coordinates_AtomC))) - (
                               (np.cos(bond_angle(Coordinates_AtomA,Coordinates_AtomC,Coordinates_AtomD)) / (
                                   bond_length(Coordinates_AtomA,Coordinates_AtomC) * 
                                   np.square(np.sin(bond_angle(Coordinates_AtomA,Coordinates_AtomC,
                                                               Coordinates_AtomD))))) * 
                               np.cross(normalized_bond_vector(Coordinates_AtomD, Coordinates_AtomC), 
                                        normalized_bond_vector(Coordinates_AtomC,Coordinates_AtomA))) 


def B_Matrix_Entry_Torsion_AtomC(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC,Coordinates_AtomD):
    return (np.cross(normalized_bond_vector(Coordinates_AtomD, Coordinates_AtomC),
                           normalized_bond_vector(Coordinates_AtomC,Coordinates_AtomA))/(bond_length(Coordinates_AtomC,Coordinates_AtomD) * 
              np.square(np.sin(bond_angle(Coordinates_AtomA,Coordinates_AtomC,Coordinates_AtomD))))) - (
                  (np.cos(bond_angle(Coordinates_AtomA,Coordinates_AtomC,Coordinates_AtomD)) / (
                      bond_length(Coordinates_AtomC,Coordinates_AtomA) * 
                    np.square(np.sin(bond_angle(Coordinates_AtomA,Coordinates_AtomC,Coordinates_AtomD))))) * 
                  np.cross(normalized_bond_vector(Coordinates_AtomD, Coordinates_AtomC),
                           normalized_bond_vector(Coordinates_AtomC,Coordinates_AtomA))) - (
                               (np.cos(bond_angle(Coordinates_AtomB,Coordinates_AtomA,Coordinates_AtomC)) / 
                    (bond_length(Coordinates_AtomC,Coordinates_AtomA) * 
                     np.square(np.sin(bond_angle(Coordinates_AtomB,Coordinates_AtomA,Coordinates_AtomC))))) *
                   np.cross(normalized_bond_vector(Coordinates_AtomB, Coordinates_AtomA), 
                            normalized_bond_vector(Coordinates_AtomA,Coordinates_AtomC))) 


''''' 
Note: For the Entries of the out-of-plane-part of the B-Matrix entry, we define the following geometry:
The atoms B, C, D are all bound to atom A and not to each other. 

The angle phi is defined as the angle of C-A-D. 

The angle theta is defined as the angle of A-B with the plane defined by A-C 
and A-D. 
It can be calculated by calculating the angle between:
    A-B and the normal vector (i.e. the cross product) of A-C and A-D

When handing in the out-of-plane: then do it in the Form (A,B,C,D)
The expressions are simplified for planar wages (i.e. theta = 0 degrees)
''''' 

def B_Matrix_Entry_OutOfPlane_AtomB(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC,Coordinates_AtomD):
    r_ab = bond_length(Coordinates_AtomA, Coordinates_AtomB)
    e_ab = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB)
    e_ac = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomC)
    e_ad = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomD)
    phi_b = bond_angle(Coordinates_AtomC, Coordinates_AtomA, Coordinates_AtomD)
    sin_theta = np.inner(e_ab, (np.cross(e_ac, e_ad)/np.sin(phi_b)))
    theta = np.arcsin(np.clip(sin_theta, 0, 1.0))
    
    if np.isclose(theta, 0):
        return np.cross(e_ac, e_ad)/(np.sin(phi_b)*r_ab)
    else:
        return (1/r_ab)*((np.cross(e_ac, e_ad)/(np.cos(theta)*np.sin(phi_b))) - np.tan(theta)*e_ab)

def B_Matrix_Entry_OutOfPlane_AtomC(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC,Coordinates_AtomD):
    r_ac = bond_length(Coordinates_AtomA, Coordinates_AtomC)
    e_ab = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB)
    e_ac = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomC)
    e_ad = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomD)
    phi_b = bond_angle(Coordinates_AtomC, Coordinates_AtomA, Coordinates_AtomD)
    phi_c = bond_angle(Coordinates_AtomB, Coordinates_AtomA, Coordinates_AtomD)
    phi_d = bond_angle(Coordinates_AtomB, Coordinates_AtomA, Coordinates_AtomC)
    sin_theta = np.inner(e_ab, (np.cross(e_ac, e_ad)/np.sin(phi_b)))
    theta = np.arcsin(np.clip(sin_theta, 0, 1.0))
    
    if np.isclose(theta, 0):
        return (1/r_ac)*((np.cross(e_ac, e_ad)/(np.sin(phi_b))) * (np.sin(phi_c)/np.sin(phi_b)))
    else:
        return (1/r_ac)*((np.cross(e_ac, e_ad)/(np.sin(phi_b))) * ((np.cos(phi_b)*np.cos(phi_c) - np.cos(phi_d))/(np.cos(theta)*np.square(np.sin(phi_b)))))
   
def B_Matrix_Entry_OutOfPlane_AtomD(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC,Coordinates_AtomD):
    r_ad = bond_length(Coordinates_AtomA, Coordinates_AtomD)
    e_ab = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomB)
    e_ac = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomC)
    e_ad = normalized_bond_vector(Coordinates_AtomA, Coordinates_AtomD)
    phi_b = bond_angle(Coordinates_AtomC, Coordinates_AtomA, Coordinates_AtomD)
    phi_c = bond_angle(Coordinates_AtomB, Coordinates_AtomA, Coordinates_AtomD)
    phi_d = bond_angle(Coordinates_AtomB, Coordinates_AtomA, Coordinates_AtomC)
    sin_theta = np.inner(e_ab, (np.cross(e_ac, e_ad)/np.sin(phi_b)))
    theta = np.arcsin(np.clip(sin_theta, 0, 1.0))
    
    if np.isclose(theta, 0):
        return (1/r_ad)*((np.cross(e_ac, e_ad)/(np.sin(phi_b))) * (np.sin(phi_d)/np.sin(phi_b)))
    else:
        return (1/r_ad)*((np.cross(e_ac, e_ad)/(np.sin(phi_b))) * ((np.cos(phi_b)*np.cos(phi_d) - np.cos(phi_c))/(np.cos(theta)*np.square(np.sin(phi_b)))))

def B_Matrix_Entry_OutOfPlane_AtomA(Coordinates_AtomA,Coordinates_AtomB,Coordinates_AtomC,Coordinates_AtomD):
    return -(B_Matrix_Entry_OutOfPlane_AtomB(Coordinates_AtomA,Coordinates_AtomB,
                                             Coordinates_AtomC,Coordinates_AtomD) 
             + B_Matrix_Entry_OutOfPlane_AtomC(Coordinates_AtomA,Coordinates_AtomB,
                                               Coordinates_AtomC,Coordinates_AtomD) 
             + B_Matrix_Entry_OutOfPlane_AtomD(Coordinates_AtomA,Coordinates_AtomB,
                                               Coordinates_AtomC,Coordinates_AtomD))

#TODO: current linear valence angles only useful for degenerate linear valence angle modes - make more generic
def b_matrix(atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof):
    n_atoms = len(atoms)
    coordinates = np.array([a.coordinates for a in atoms])
    atom_index = {a.symbol: i for i, a in enumerate(atoms)}
    n_internal = len(bonds) + len(angles) + len(linear_angles) + len(out_of_plane) + len(dihedrals)
    assert n_internal >= idof, \
        f"Wrong number of internal coordinates, n_internal ({n_internal}) should be >= {idof}."
    matrix = np.zeros((n_internal, 3*n_atoms))
    i_internal = 0
    n_used_linear_angles = 0
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
        if len(angles) == 0:
            if (n_used_linear_angles %2) == 0:
                matrix[i_internal, index[0]:index[0]+3] = B_Matrix_Entry_LinearAngleFirstPlane_AtomB(coord[1], coord[0], coord[2])
                matrix[i_internal, index[1]:index[1]+3] = B_Matrix_Entry_LinearAngleFirstPlane_AtomA(coord[1], coord[0], coord[2])
                matrix[i_internal, index[2]:index[2]+3] = B_Matrix_Entry_LinearAngleFirstPlane_AtomC(coord[1], coord[0], coord[2])
                i_internal += 1
                n_used_linear_angles += 1
            else:
                matrix[i_internal, index[0]:index[0]+3] = B_Matrix_Entry_LinearAngleSecondPlane_AtomB(coord[1], coord[0], coord[2])
                matrix[i_internal, index[1]:index[1]+3] = B_Matrix_Entry_LinearAngleSecondPlane_AtomA(coord[1], coord[0], coord[2])
                matrix[i_internal, index[2]:index[2]+3] = B_Matrix_Entry_LinearAngleSecondPlane_AtomC(coord[1], coord[0], coord[2])
                i_internal += 1
                n_used_linear_angles += 1
        else:
            matrix[i_internal, index[0]:index[0]+3] = B_Matrix_Entry_LinearAngleFirstPlane_AtomB(coord[1], coord[0], coord[2])
            matrix[i_internal, index[1]:index[1]+3] = B_Matrix_Entry_LinearAngleFirstPlane_AtomA(coord[1], coord[0], coord[2])
            matrix[i_internal, index[2]:index[2]+3] = B_Matrix_Entry_LinearAngleFirstPlane_AtomC(coord[1], coord[0], coord[2])
            i_internal += 1
    for outofplane in out_of_plane:
        index = [atom_index[a] * 3 for a in outofplane]
        coord = [coordinates[atom_index[a]] for a in outofplane]
        matrix[i_internal, index[0]:index[0]+3] = B_Matrix_Entry_OutOfPlane_AtomA(coord[1], coord[0], coord[2], coord[3])
        matrix[i_internal, index[1]:index[1]+3] = B_Matrix_Entry_OutOfPlane_AtomB(coord[1], coord[0], coord[2], coord[3])
        matrix[i_internal, index[2]:index[2]+3] = B_Matrix_Entry_OutOfPlane_AtomC(coord[1], coord[0], coord[2], coord[3])
        matrix[i_internal, index[3]:index[3]+3] = B_Matrix_Entry_OutOfPlane_AtomD(coord[1], coord[0], coord[2], coord[3])
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
