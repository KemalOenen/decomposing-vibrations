import itertools
from collections import Counter
from os import remove
from pickle import TRUE
import numpy as np
import pandas as pd
import pprint
import nomodeco

"""""
step 1: generate all IC sets that have 3N-6 ICs and have all bonds included - DONE
step 2: do a calculation run for every possibiliy - DONE
step 3: compute the B-Matrix for all possibilites and check for completeness, if not complete remove - DONE

step 4: order the decomposition according to new metric - NOT DONE
step 5: think of more methods to reduce the high dimensionaliy of new IC sets - NOT DONE
"""""
#basically computes the MAD but I wanted a fancy name
def Kemalian_metric(Contribution_Matrix, intrinsic_frequencies, harmonic_frequencies):
    sum_abs = 0
    for i in range(len(harmonic_frequencies)):
        # print(np.dot(Contribution_Matrix[:,i],intrinsic_frequencies)-harmonic_frequencies[i])
        sum_abs += np.abs(np.dot(Contribution_Matrix[:,i],intrinsic_frequencies)-harmonic_frequencies[i])
    return sum_abs/len(harmonic_frequencies)

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

def all_atoms_can_be_superimposed_dihedral(test_dihedral, key_dihedral, nested_equivalent_atoms):
    return (test_dihedral[0] == key_dihedral[0] or check_in_nested_list([test_dihedral[0],key_dihedral[0]], nested_equivalent_atoms)) and (
            test_dihedral[1] == key_dihedral[1] or check_in_nested_list([test_dihedral[1],key_dihedral[1]], nested_equivalent_atoms)) and (
                    test_dihedral[2] == key_dihedral[2] or check_in_nested_list([test_dihedral[2],key_dihedral[2]], nested_equivalent_atoms)) and (
                            test_dihedral[3] == key_dihedral[3] or check_in_nested_list([test_dihedral[3],key_dihedral[3]], nested_equivalent_atoms))


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
def get_angle_subsets(symmetric_angles,num_bonds,num_angles,idof,n_phi) -> list:
    symmetric_angles_list, angles = [], []
    #TODO: check for special case if numbond+numangle > idof, then break symmetry

    for ind_angle in symmetric_angles.keys():
        if symmetric_angles[ind_angle] not in symmetric_angles_list:
            symmetric_angles_list.append(symmetric_angles[ind_angle])
    
    for i in range(1,len(symmetric_angles_list)+1):
        for angle_subset in itertools.combinations(symmetric_angles_list,i):
            flat_angle_subset = [item for sublist in angle_subset for item in sublist]
            if n_phi <= len(list(flat_angle_subset)) <= n_phi+2:
                angles.append(list(flat_angle_subset))
    return angles

def get_symm_dihedrals(dihedrals,specification):
    symmetric_dihedrals = dict()
    symmetric_dihedrals = {key:[] for (key, val) in Counter(dihedrals).items()}
   
    # symmetric dihedrals equally defined as in get_symm_angles2 --> make same function?
    for i, key in itertools.product(range(len(dihedrals)), symmetric_dihedrals):
        symmetric_dihedrals[key].append(dihedrals[i])

    for key,val in symmetric_dihedrals.items():
        i=0
        while i<len(val):
            dihedral = val[i]
            if not all_atoms_can_be_superimposed_dihedral(dihedral,key,specification["equivalent_atoms"]):
                del val[i]
            elif all_atoms_can_be_superimposed_dihedral(dihedral,key,specification["equivalent_atoms"]):
                i +=1
    return symmetric_dihedrals

def get_dihedral_subsets(symmetric_dihedrals,num_bonds,num_angles,idof,n_tau) -> list:
    symmetric_dihedrals_list, dihedrals = [], []
    #TODO: check for special cases

    for ind_dihedral in symmetric_dihedrals.keys():
        if symmetric_dihedrals[ind_dihedral] not in symmetric_dihedrals_list:
            symmetric_dihedrals_list.append(symmetric_dihedrals[ind_dihedral])
    for i in range(0,len(symmetric_dihedrals_list)+1):
        for dihedral_subset in itertools.combinations(symmetric_dihedrals_list,i):
            flat_dihedral_subset = [item for sublist in dihedral_subset for item in sublist]
            if n_tau <= len(list(flat_dihedral_subset)) <= n_tau+1:
                dihedrals.append(list(flat_dihedral_subset))
    return dihedrals

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

def check_evalue_f_matrix(reciprocal_square_massmatrix, B, B_inv, InternalF_Matrix):
    CartesianF_Matrix_check = np.transpose(B) @ InternalF_Matrix @ B
    evalue, evect = np.linalg.eigh(np.transpose(reciprocal_square_massmatrix) @ CartesianF_Matrix_check @ reciprocal_square_massmatrix)
    return evalue


def number_terminal_bonds(mult_list):
    number_of_terminal_bonds = 0
    for atom_and_mult in mult_list:
        if atom_and_mult[1] == 1:
            number_of_terminal_bonds += 1
    return number_of_terminal_bonds


def not_same_central_atom(list_oop_angles) -> bool:
    central_atoms = set()
    not_same_central_atom = True
    for oop_angle in list_oop_angles:
        if oop_angle[0] in central_atoms:
            not_same_central_atom = False
            break
        else:
            central_atoms.add(oop_angle[0])
    return not_same_central_atom

def matrix_norm(matrix, matrix_inv, p):
    return np.linalg.norm(matrix, p) * np.linalg.norm(matrix_inv, p) 


# This is the new get_sets method including Decius paper (https://doi.org/10.1063/1.1747158) on the selection of IC sets
#TODO: include symmetry considerations - currently we work highly non-redundant

def get_sets(idof,atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, specification):
    ic_dict = dict() 
    num_bonds = len(bonds)
    num_atoms = len(atoms)

    # @decision tree: linear
    if specification["linearity"] == "fully linear":
        ic_dict[0] = {
                "bonds" : bonds,
                "angles" : angles,
                "linear valence angles" : linear_angles,
                "out of plane angles" : out_of_plane,
                "dihedrals" : dihedrals} 
    
    # @decision tree: planar (acyclic, cyclic is smth missing and not linear submolecules) 
    if specification["planar"] == "yes" and not specification["linearity"] == "fully linear":
        num_of_red = 6*specification["mu"]
        a_1 = number_terminal_bonds(specification["multiplicity"])

        # set length of subsets
        n_r = num_bonds
        n_phi = 2*num_bonds - num_atoms
        n_gamma = 2*(num_bonds-num_atoms) + a_1
        n_tau = num_bonds - a_1
 
        symmetric_angles = get_symm_angles2(angles,specification)
        angle_subsets = get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)
        
        # the if statement ensures, that oop angles to the same central atom can not be in the same set
        oop_subsets = []
        for subset in itertools.combinations(out_of_plane, n_gamma):
            if not_same_central_atom(subset): 
                oop_subsets.append(list(subset))

        symmetric_dihedrals = get_symm_dihedrals(dihedrals,specification)
        dihedral_subsets = get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau)

        # TODO here: cyclic molecules have these redundancies that we need to stash

        k = 0
        for len_angles in range(0, len(angle_subsets)):
            for len_oop in range(0, len(oop_subsets)):
                for len_dihedrals in range(0, len(dihedral_subsets)):
                    ic_dict[k] = {
                            "bonds" : bonds,
                            "angles" : angle_subsets[len_angles],
                            "linear valence angles" : [] ,
                            "out of plane angles" : oop_subsets[len_oop],
                            "dihedrals" : dihedral_subsets[len_dihedrals] }
                    k +=1
 
    # @decision tree: general molecule (acyclic, cyclic is smth missing and not linear submolecules)
    if specification["planar"] == "no" and not specification["linearity"] == "fully linear":
        num_of_red = 6*specification["mu"]
        a_1 = number_terminal_bonds(specification["multiplicity"])

        # set length of subsets
        n_r = num_bonds
        n_phi = 4*num_bonds - 3*num_atoms + a_1
        n_tau = num_bonds - a_1

        symmetric_angles = get_symm_angles2(angles,specification)
        angle_subsets = get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)

        symmetric_dihedrals = get_symm_dihedrals(dihedrals,specification)
        dihedral_subsets = get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau)

        # TODO here: cyclic molecules have these redundancies that we need to stash

        k = 0
        for len_angles in range(0, len(angle_subsets)):
            for len_dihedrals in range(0, len(dihedral_subsets)):
                ic_dict[k] = {
                        "bonds" : bonds,
                        "angles" : angle_subsets[len_angles],
                        "linear valence angles" : [] ,
                        "out of plane angles" : [],
                        "dihedrals" : dihedral_subsets[len_dihedrals] }
                k +=1

    # @decision tree: planar molecules with linear submolecules
    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found":
        num_of_red = 6*specification["mu"]
        a_1 = number_terminal_bonds(specification["multiplicity"])
        l = specification["length of linear submolecule(s) l"]

        # find all combinations of x and y
        x_and_y = []
        for x in range(l):
            y = l - 1 -x
            x_and_y.append([x,y])

        # set length of subsets
        n_r = num_bonds
        n_phi = 2*num_bonds - num_atoms
        n_phi_prime = l - 1
        n_gamma = []
        n_tau = []
        
        for pair in x_and_y:
            n_gamma.append(2*(num_bonds-num_atoms) + a_1 - pair[0])
            n_tau.append(num_bonds - a_1 - pair[1])
 
        symmetric_angles = get_symm_angles2(angles,specification)
        angle_subsets = get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)
        
        # the if statement ensures, that oop angles to the same central atom can not be in the same set
        # for linear submolecules, multiple sets can be sadly formulated!

        oop_subsets = []
        for i in range(len(n_gamma)):
            for subset in itertools.combinations(out_of_plane, n_gamma[i]):
                if not_same_central_atom(subset): 
                    oop_subsets.append(list(subset))

        symmetric_dihedrals = get_symm_dihedrals(dihedrals,specification)
        
        dihedral_subsets = []
        for i in range(len(n_tau)):
            dihedral_subsets.append(list(
                get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau[i])))

        # TODO here: cyclic molecules have these redundancies that we need to stash

        k = 0
        for x_y in range(0,n_tau): # there as many taus and phis, and as x and y were ordered, the values are also ordered
            for len_angles in range(0, len(angle_subsets)):
                for len_oop in range(0, len(oop_subsets[x_y])):
                    for len_dihedrals in range(0, len(dihedral_subsets[x_y])):
                        ic_dict[k] = {
                                "bonds" : bonds,
                                "angles" : angle_subsets[len_angles],
                                "linear valence angles" : [] ,
                                "out of plane angles" : oop_subsets[x_y][len_oop],
                                "dihedrals" : dihedral_subsets[x_y][len_dihedrals] }
                        k +=1

    # @decision tree: general molecule and linear submolecules)
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found":
        num_of_red = 6*specification["mu"]
        a_1 = number_terminal_bonds(specification["multiplicity"])
        l = specification["length of linear submolecule(s) l"]

        # set length of subsets
        n_r = num_bonds
        n_phi = 4*num_bonds - 3*num_atoms + a_1
        n_phi_prime = l-1
        n_tau = num_bonds - a_1 - (l-1)

        symmetric_angles = get_symm_angles2(angles,specification)
        angle_subsets = get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)
 
        symmetric_lin_angles = get_symm_angles2(linear_angles,specification)
        lin_angle_subsets = get_angle_subsets(symmetric_lin_angles, len(bonds), len(angles),idof,n_phi_prime)

        symmetric_dihedrals = get_symm_dihedrals(dihedrals,specification)
        dihedral_subsets = get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau)

        # TODO here: cyclic molecules have these redundancies that we need to stash

        k = 0
        for len_angles in range(0, len(angle_subsets)):
            for len_lin_angles in range(0, len(lin_angle_subsets)):
                for len_dihedrals in range(0, len(dihedral_subsets)):
                    ic_dict[k] = {
                            "bonds" : bonds,
                            "angles" : angle_subsets[len_angles],
                            "linear valence angles" : lin_angle_subsets[len_lin_angles] ,
                            "out of plane angles" : [],
                            "dihedrals" : dihedral_subsets[len_dihedrals] }
                    k +=1


    print(len(ic_dict), "internal coordinate sets (that should be tested) have been generated.")
    return ic_dict

# The routine below includes new symmetry considerations but is not based on topological arguments; will be removed soon
#def get_sets(n_atoms, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, specification):
#    num_bonds = len(bonds)
#    num_angles = len(angles)
#    ic_dict = dict()
#
#    #use get_symm_angles2 as get_symm_angles has a "restricted" definition!
#    #symmetric_angles = get_symm_angles(angles,specification)
#    symmetric_angles = get_symm_angles2(angles,specification)
#    angle_subsets = get_angle_subsets(symmetric_angles, num_bonds, num_angles,idof)
#
#    k = 0

#    # Dihedrals are treated like this: the symmetric dihedrals are identified, then 
#    # it is checked how big the subsets are; take only dihedrals in a new reduced set
#    # which are included with their symmetrical counterparts or are together with the other IC's
#    # under a certain threshold 
#
#    symmetric_dihedrals = get_symm_dihedrals(dihedrals,specification)
#    dihedral_subsets = get_dihedral_subsets(symmetric_dihedrals, num_bonds, num_angles,idof)
#
#    #possible strategy for breaking symmetry
#    #dihedrals = reduce_dihedral_sets(symmetric_dihedrals,specification["dihedral_reduction"][1])
#    
#    #TODO: do something about this horrible for loop!!!!!!!!!!
#    for i in range(0, (3*n_atoms) - (idof+2)): # idof+2: up to three redundancies
#        for j in range(0, len(angle_subsets)):
#            for d in range(0, len(dihedral_subsets)):
#                if (0.4*len(angle_subsets[j]) >= len(dihedral_subsets[d])) and (idof - num_bonds - len(angle_subsets[j]) - len(dihedral_subsets[d]) + i) >= 0:
#                    for ic_subset in itertools.combinations(linear_angles + out_of_plane, idof - num_bonds - len(angle_subsets[j]) - len(dihedral_subsets[d]) + i):
#                        used_linear_angles, used_out_of_plane = [], []
#                        for l in range(0, len(ic_subset)):
#                            if ic_subset[l] in linear_angles:
#                                used_linear_angles.append(ic_subset[l])
#                            if ic_subset[l] in out_of_plane and avoid_double_oop(ic_subset[l], used_out_of_plane):
#                                used_out_of_plane.append(ic_subset[l])
#                        if (num_bonds + len(angle_subsets[j]) + len(used_linear_angles) + len(used_out_of_plane) + len(dihedral_subsets[d]) >= idof):
#                                ic_dict[k] = {
#                                    "bonds" : bonds,
#                                    "angles" : angle_subsets[j],
#                                    "linear valence angles" : used_linear_angles,
#                                    "out of plane angles" : used_out_of_plane,
#                                    "dihedrals" : dihedral_subsets[d]
#                                }
#                                k +=1
#                    used_linear_angles, used_out_of_plane = [], []
#    print(len(ic_dict), "internal coordinate sets (that should be tested) have been generated.")
#    return ic_dict
