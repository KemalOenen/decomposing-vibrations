import itertools
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
    return atom_list

def ic_rules_angles(used_angles, angles) -> bool:
    used_angles = remove_enumeration(used_angles)
    angles = remove_enumeration(angles)

    angle_dict = {}
    for ang in angles:
        distinct_angle = tuple(ang)
        if distinct_angle in angle_dict:
            angle_dict[distinct_angle] += 1
        else:
            angle_dict[distinct_angle] = 1
    
    i = 0
    for test_ang in angle_dict.keys():
        if used_angles.count(list(test_ang)) == angle_dict.get(test_ang):
            i += 1
    
    if i == len(angle_dict):
        return True
    else:
        return False

def ic_rules(ic_subset):
    ic_subset = list(ic_subset)

    return True

def test_completeness(CartesianF_Matrix, B, B_inv, InternalF_Matrix) -> bool:
    CartesianF_Matrix_check = np.transpose(B) @ InternalF_Matrix @ B
    if (np.allclose(CartesianF_Matrix_check, CartesianF_Matrix)) == True:
        return True
    else:
        return False

def matrix_norm(matrix, matrix_inv, p):
    return np.linalg.norm(matrix, p) * np.linalg.norm(matrix_inv, p) 

def generate_all_possible_sets(n_atoms, idof, bonds, angles, linear_angles, out_of_plane, dihedrals):
    num_bonds = len(bonds)
    ic_dict = dict()

    # TODO: make approach more elegant? -> include filter(), to have less but useful combinations
    # more concrete: create IC subsets that are cleverly chosen; these criteria are:
    # TODO: the used oop per IC set should not be just changed in order; 
    # TODO: if there are angles that are the same, then they need to be included in the analysis       
    # TODO: think of reasonable cut off for maximum relative value of used dihedrals
    k = 0
    for i in range(0, (3*n_atoms) - idof):
        for ic_subset in filter(ic_rules,itertools.combinations(angles + linear_angles + out_of_plane + dihedrals, idof - num_bonds + i)):
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
                if ic_rules_angles(used_angles, angles):
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
    return ic_dict