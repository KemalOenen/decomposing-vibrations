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

def are_two_elements_same(tup1, tup2):
    return ((tup1[0] == tup2[0] and tup1[1] == tup2[1]) or
            (tup1[1] == tup2[1] and tup1[2] == tup2[2]) or
            (tup1[0] == tup2[0] and tup1[2] == tup2[2]))

def get_different_elements(tup1, tup2):
    differences = []
    for element1 in tup1:
        if element1 not in tup2:
            differences.append(element1)
    for element2 in tup2:
        if element2 not in tup1:
            differences.append(element2)
    return differences

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

def check_in_nested_list(check_list, nested_list):
    check = False
    for single_list in nested_list:
        if set(check_list).issubset(set(single_list)):
            check = True
    return check

def get_symm_angles(angles,specification):
    symmetric_angles = dict()
    symmetric_angles = {key:[] for (key, val) in Counter(angles).items()}
    
    #angles are the same if they share two atoms and
    #the other one is symmetry equivalent

    for i,key in itertools.product(range(len(angles)), symmetric_angles):
        if are_two_elements_same(key, angles[i]):
            symmetric_angles[key].append(angles[i])

    for key,val in symmetric_angles.items():
        i=0
        while i<len(val):
            ang = val[i]
            if key == ang:
                i += 1
                continue
            elif not check_in_nested_list(get_different_elements(ang,key), specification["equivalent_atoms"]):
                del val[i]
            elif check_in_nested_list(get_different_elements(ang,key), specification["equivalent_atoms"]):
                i += 1

    return symmetric_angles

def all_atoms_can_be_superimposed(test_angle, key_angle, nested_equivalent_atoms):
    return (test_angle[0] == key_angle[0] or check_in_nested_list([test_angle[0],key_angle[0]], nested_equivalent_atoms)) and (
            test_angle[1] == key_angle[1] or check_in_nested_list([test_angle[1],key_angle[1]], nested_equivalent_atoms)) and (
                    test_angle[2] == key_angle[2] or check_in_nested_list([test_angle[2],key_angle[2]], nested_equivalent_atoms))

def get_symm_angles2(angles,specification):
    symmetric_angles = dict()
    symmetric_angles = {key:[] for (key, val) in Counter(angles).items()}
    
    #angles are the same if the atoms can all be superimposed
    #on each other with symmetry operations

    for i,key in itertools.product(range(len(angles)), symmetric_angles):
        symmetric_angles[key].append(angles[i])

    for key,val in symmetric_angles.items():
        i=0
        while i<len(val):
            ang = val[i]
            if not all_atoms_can_be_superimposed(ang,key,specification["equivalent_atoms"]):
                del val[i]
            elif all_atoms_can_be_superimposed(ang,key,specification["equivalent_atoms"]):
                i +=1

    return symmetric_angles

#TODO: think about special case where num_bonds + num_angles > idof
#e.g.: make test if set has to many redundancies and if so -> break symmetry
def get_angle_subsets(symmetric_angles,num_bonds,num_angles,idof) -> list:
    symmetric_angles_list, angles = [], []
    #TODO: check for special case if numbond+numangle > idof, then break symmetry

    for ind_angle in symmetric_angles.keys():
        if symmetric_angles[ind_angle] not in symmetric_angles_list:
            symmetric_angles_list.append(symmetric_angles[ind_angle])
    print(symmetric_angles_list)
    for i in range(1,len(symmetric_angles_list)+1):
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


def reduce_dihedral_sets(symmetric_dihedrals, n_repr):
    dihedrals = []
    for symm_dihedrals in symmetric_dihedrals.keys():
        for i in range(0,n_repr):
            dihedrals.append(symmetric_dihedrals[symm_dihedrals][i])
    return dihedrals

def test_completeness(CartesianF_Matrix, B, B_inv, InternalF_Matrix) -> bool:
    CartesianF_Matrix_check = np.transpose(B) @ InternalF_Matrix @ B
    if (np.allclose(CartesianF_Matrix_check, CartesianF_Matrix)) == True:
        return True
    else:
        return False

def matrix_norm(matrix, matrix_inv, p):
    return np.linalg.norm(matrix, p) * np.linalg.norm(matrix_inv, p) 

#TODO:remove stupid angle subsets maybe even before
def get_sets(n_atoms, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, specification):
    num_bonds = len(bonds)
    num_angles = len(angles)
    ic_dict = dict()

    #symmetric_angles = get_symm_angles(angles,specification)
    symmetric_angles = get_symm_angles2(angles,specification)
    print(specification["equivalent_atoms"])
    pprint.pprint(symmetric_angles)
    angle_subsets = get_angle_subsets(symmetric_angles, num_bonds, num_angles,idof)
   # print(angle_subsets)

    k = 0
 
     #if specification["dihedral_reduction"][0] == "dihedral_reduction":
     #    symmetric_dihedrals = get_symm_dihedrals(dihedrals)
     #    dihedrals = reduce_dihedral_sets(symmetric_dihedrals,specification["dihedral_reduction"][1])
    
    for i in range(0, (3*n_atoms) - (idof+2)): # up to three reduncancies
        for j in range(0, len(angle_subsets)):
            for ic_subset in itertools.combinations(linear_angles + out_of_plane + dihedrals, idof - num_bonds - len(angle_subsets[j]) + i):
                used_linear_angles, used_out_of_plane, used_dihedrals = [], [], []
                for l in range(0, len(ic_subset)):
                    if ic_subset[l] in linear_angles:
                        used_linear_angles.append(ic_subset[l])
                    if ic_subset[l] in out_of_plane and avoid_double_oop(ic_subset[l], used_out_of_plane):
                        used_out_of_plane.append(ic_subset[l])
                    if ic_subset[l] in dihedrals:
                        used_dihedrals.append(ic_subset[l])
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

 #def get_sets(n_atoms, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, specification):
 #    num_bonds = len(bonds)
 #    num_angles = len(angles)
 #    ic_dict = dict()
 #
 #    k = 0
 # 
 #    if specification["dihedral_reduction"][0] == "dihedral_reduction":
 #        symmetric_dihedrals = get_symm_dihedrals(dihedrals)
 #        dihedrals = reduce_dihedral_sets(symmetric_dihedrals,specification["dihedral_reduction"][1])
 #    
 #    if specification["angle_symmetry"] == "angle_symmetry":
 #        symmetric_angles = get_symm_angles(angles)
 #        angle_subsets = get_angle_subsets(symmetric_angles, num_bonds, num_angles,idof)
 #        for i in range(0, (3*n_atoms) - idof):
 #            for j in range(0, len(angle_subsets)):
 #                for ic_subset in itertools.combinations(linear_angles + out_of_plane + dihedrals, idof - num_bonds - len(angle_subsets[j]) + i):
 #                    used_linear_angles, used_out_of_plane, used_dihedrals = [], [], []
 #                    for i in range(0, len(ic_subset)):
 #                        if ic_subset[i] in linear_angles:
 #                            used_linear_angles.append(ic_subset[i])
 #                        if ic_subset[i] in out_of_plane and avoid_double_oop(ic_subset[i], used_out_of_plane):
 #                            used_out_of_plane.append(ic_subset[i])
 #                        if ic_subset[i] in dihedrals:
 #                            used_dihedrals.append(ic_subset[i])
 #                    if (num_bonds + len(angle_subsets[j]) + len(used_linear_angles) + len(used_out_of_plane) + len(used_dihedrals)) >= idof:
 #                            ic_dict[k] = {
 #                                "bonds" : bonds,
 #                                "angles" : angle_subsets[j],
 #                                "linear valence angles" : used_linear_angles,
 #                                "out of plane angles" : used_out_of_plane,
 #                                "dihedrals" : used_dihedrals
 #                            }
 #                            k +=1
 #                used_linear_angles, used_out_of_plane, used_dihedrals = [], [], []
 #        print(len(ic_dict), "internal coordinate sets (that should be tested) have been generated.")
 #        return ic_dict
 #
 #    elif specification["angle_symmetry"] == "no-angle_symmetry":
 #        for i in range(0, (3*n_atoms) - idof):
 #                for ic_subset in itertools.combinations(angles + linear_angles + out_of_plane + dihedrals, idof - num_bonds + i):
 #                    used_angles, used_linear_angles, used_out_of_plane, used_dihedrals = [], [], [], []
 #                    for i in range(0, len(ic_subset)):
 #                        if ic_subset[i] in angles:
 #                            used_angles.append(ic_subset[i])
 #                        if ic_subset[i] in linear_angles:
 #                            used_linear_angles.append(ic_subset[i])
 #                        if ic_subset[i] in out_of_plane and avoid_double_oop(ic_subset[i], used_out_of_plane):
 #                            used_out_of_plane.append(ic_subset[i])
 #                        if ic_subset[i] in dihedrals:
 #                            used_dihedrals.append(ic_subset[i])
 #                    if (num_bonds + len(used_angles) + len(used_linear_angles) + len(used_out_of_plane) + len(used_dihedrals)) >= idof:
 #                            ic_dict[k] = {
 #                                "bonds" : bonds,
 #                                "angles" : used_angles,
 #                                "linear valence angles" : used_linear_angles,
 #                                "out of plane angles" : used_out_of_plane,
 #                                "dihedrals" : used_dihedrals
 #                            }
 #                            k +=1
 #                    used_angles, used_linear_angles, used_out_of_plane, used_dihedrals = [], [], [], []
 #        print(len(ic_dict), "internal coordinate sets (that should be tested) have been generated.")
 #        return ic_dict

 #def get_symm_angles(angles):
 #    symmetric_angles = dict()
 #    symmetric_angles = {key:[] for (key, val) in Counter(remove_enumeration(angles)).items()}
 #    for i,key in itertools.product(range(len(angles)), symmetric_angles):
 #        if key == remove_enumeration_tuple(angles[i]):
 #            symmetric_angles[key].append(angles[i])
 #    return symmetric_angles
 #
 #def get_angle_subsets(symmetric_angles,num_bonds, num_angles,idof) -> list:
 #    symmetric_angles_list, angles = [], []
 #    if (num_bonds + num_angles) <= idof:
 #        k = 1
 #    else:
 #        k = 0
 #    for symm_angles in symmetric_angles.keys():
 #        symmetric_angles_list.append(symmetric_angles[symm_angles])
 #    for i in range(1, len(symmetric_angles_list)+k):
 #        for angle_subset in itertools.combinations(symmetric_angles_list,i):
 #            flat_angle_subset = [item for sublist in angle_subset for item in sublist]
 #            angles.append(list(flat_angle_subset))
 #    return angles
 #
 #def get_symm_dihedrals(dihedrals):
 #    symmetric_dihedrals = dict()
 #    symmetric_dihedrals = {key:[] for (key, val) in Counter(remove_enumeration(dihedrals)).items()}
 #    for i,key in itertools.product(range(len(dihedrals)), symmetric_dihedrals):
 #        if key == remove_enumeration_tuple(dihedrals[i]):
 #            symmetric_dihedrals[key].append(dihedrals[i])
 #    return symmetric_dihedrals
 #
 #
 #def reduce_dihedral_sets(symmetric_dihedrals, n_repr):
 #    dihedrals = []
 #    for symm_dihedrals in symmetric_dihedrals.keys():
 #        for i in range(0,n_repr):
 #            dihedrals.append(symmetric_dihedrals[symm_dihedrals][i])
 #    return dihedrals
