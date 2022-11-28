import itertools
from collections import Counter
from os import remove
from pickle import TRUE
import numpy as np
import pprint
import nomodeco

"""""
step 1: generate all IC sets that have 3N-6 ICs and have all bonds included - DONE
step 2: do a calculation run for every possibiliy - DONE
step 3: compute the B-Matrix for all possibilites and check for completeness, if not complete remove - DONE

step 4: order the decomposition according to new metric - NOT DONE
step 5: think of more methods to reduce the high dimensionaliy of new IC sets - NOT DONE
"""""

def Kemalian_metric():
    return 0

def avoid_double_oop(test_oop, used_out_of_plane) -> bool:
    if len(used_out_of_plane) == 0:
        return True
    for i in range(0, len(used_out_of_plane)):
        if set(test_oop).issubset(used_out_of_plane[i]):
            return False
    return True

def remove_enumeration(atom_list) -> list:
    atom_list = [list(tup) for tup in atom_list]
    for i in range(0, len(atom_list)):
        for j in range(0, len(atom_list[i])):
            atom_list[i][j] = ''.join(c for c in atom_list[i][j] if not c.isnumeric())
        atom_list[i] = tuple(atom_list[i])
    return atom_list

def remove_enumeration_tuple(atom_tuple) -> tuple:
    atom_list = list(atom_tuple)
    for i in range(0, len(atom_tuple)):
        atom_list[i] = ''.join(c for c in atom_list[i] if not c.isnumeric())
    return tuple(atom_list)

def get_symm_angles(angles):
    symmetric_angles = dict()
    symmetric_angles = {key:[] for (key, val) in Counter(remove_enumeration(angles)).items()}
    for i,key in itertools.product(range(len(angles)), symmetric_angles):
        if key == remove_enumeration_tuple(angles[i]):
            symmetric_angles[key].append(angles[i])
    return symmetric_angles

def get_angle_subsets(symmetric_angles,num_bonds, num_angles,idof) -> list:
    symmetric_angles_list, angles = [], []
    if (num_bonds + num_angles) <= idof:
        k = 1
    else:
        k = 0
    for symm_angles in symmetric_angles.keys():
        symmetric_angles_list.append(symmetric_angles[symm_angles])
    for i in range(1, len(symmetric_angles_list)+k):
        for angle_subset in itertools.combinations(symmetric_angles_list,i):
            flat_angle_subset = [item for sublist in angle_subset for item in sublist]
            angles.append(list(flat_angle_subset))
    return angles

def get_symm_dihedrals(dihedrals):
    symmetric_dihedrals = dict()
    symmetric_dihedrals = {key:[] for (key, val) in Counter(remove_enumeration(dihedrals)).items()}
    for i,key in itertools.product(range(len(dihedrals)), symmetric_dihedrals):
        if key == remove_enumeration_tuple(dihedrals[i]):
            symmetric_dihedrals[key].append(dihedrals[i])
    return symmetric_dihedrals


def reduce_dihedral_sets(symmetric_dihedrals):
    dihedrals = []
    for symm_dihedrals in symmetric_dihedrals.keys():
        dihedrals.append(symmetric_dihedrals[symm_dihedrals][0])
    return dihedrals

def test_completeness(CartesianF_Matrix, B, B_inv, InternalF_Matrix) -> bool:
    CartesianF_Matrix_check = np.transpose(B) @ InternalF_Matrix @ B
    if (np.allclose(CartesianF_Matrix_check, CartesianF_Matrix)) == True:
        return True
    else:
        return False

def matrix_norm(matrix, matrix_inv, p):
    return np.linalg.norm(matrix, p) * np.linalg.norm(matrix_inv, p) 

    # TODO: make approach more elegant? -> include filter(), to have less but useful combinations
    # more concrete: create IC subsets that are cleverly chosen; these criteria are:
    # TODO: the used oop per IC set should not be just changed in order; 
    # TODO: if there are angles that are the same, then they need to be included in the analysis       
    # TODO: think of reasonable cut off for maximum relative value of used dihedrals

def get_sets(n_atoms, idof, bonds, angles, linear_angles, out_of_plane, dihedrals):
    num_bonds = len(bonds)
    num_angles = len(angles)
    ic_dict = dict()

    k = 0
    symmetric_angles = get_symm_angles(angles)
    angle_subsets = get_angle_subsets(symmetric_angles, num_bonds, num_angles,idof)

    #TODO: think of this empirical value!!!!!!! make it relatve to the IC 
    if (len(dihedrals)/(num_bonds+num_angles+len(linear_angles))) > 0.40:
        symmetric_dihedrals = get_symm_dihedrals(dihedrals)
        dihedrals = reduce_dihedral_sets(symmetric_dihedrals)

    for i in range(0, (3*n_atoms) - idof):
        for j in range(0, len(angle_subsets)):
            for ic_subset in itertools.combinations(linear_angles + out_of_plane + dihedrals, idof - num_bonds - len(angle_subsets[j]) + i):
                used_linear_angles, used_out_of_plane, used_dihedrals = [], [], []
                for i in range(0, len(ic_subset)):
                    if ic_subset[i] in linear_angles:
                        used_linear_angles.append(ic_subset[i])
                    if ic_subset[i] in out_of_plane and avoid_double_oop(ic_subset[i], used_out_of_plane):
                        used_out_of_plane.append(ic_subset[i])
                    if ic_subset[i] in dihedrals:
                        used_dihedrals.append(ic_subset[i])
                if (num_bonds + len(angle_subsets[j]) + len(used_linear_angles) + len(used_out_of_plane) + len(used_dihedrals)) >= idof:
                        ic_dict[k] = {
                            "bonds" : bonds,
                            "angles" : angle_subsets[j],
                            "linear valence angles" : used_linear_angles,
                            "out of plane angles" : used_out_of_plane,
                            "dihedrals" : used_dihedrals
                        }
                        k +=1
            used_linear_angles, used_out_of_plane, used_dihedrals = [], [], []
    print(len(ic_dict), "internal coordinate sets (that should be tested) have been generated.")
    return ic_dict

"""
for i in range(0, (3*n_atoms) - idof):
            for ic_subset in itertools.combinations(angles + linear_angles + out_of_plane + dihedrals, idof - num_bonds + i):
                used_angles, used_linear_angles, used_out_of_plane, used_dihedrals = [], [], [], []
                for i in range(0, len(ic_subset)):
                    if ic_subset[i] in angles:
                        used_angles.append(ic_subset[i])
                    if ic_subset[i] in linear_angles:
                        used_linear_angles.append(ic_subset[i])
                    if ic_subset[i] in out_of_plane and avoid_double_oop(ic_subset[i], used_out_of_plane):
                        used_out_of_plane.append(ic_subset[i])
                    if ic_subset[i] in dihedrals:
                        used_dihedrals.append(ic_subset[i])
                if (num_bonds + len(used_angles) + len(used_linear_angles) + len(used_out_of_plane) + len(used_dihedrals)) >= idof:
                        ic_dict[k] = {
                            "bonds" : bonds,
                            "angles" : used_angles,
                            "linear valence angles" : used_linear_angles,
                            "out of plane angles" : used_out_of_plane,
                            "dihedrals" : used_dihedrals
                        }
                        k +=1
                used_angles, used_linear_angles, used_out_of_plane, used_dihedrals = [], [], [], []
    print(len(ic_dict), "internal coordinate sets (that should be tested) have been generated.")
"""