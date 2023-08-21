'''''
This module contains the get_sets method as well as symmetry functions for selecting internal coordinates
'''''

import itertools
import logging
from collections import Counter
from os import remove
from pickle import TRUE
import numpy as np
import pandas as pd
import pprint
import topology
import nomodeco


#TODO: symmetry breaking needs to be revised completely -> currently only 3N-6 coordinate sets

def Kemalian_metric(ved_matrix):
    # axis = 0, when maximum of each column
    max_values = np.max(ved_matrix, axis=1)
    return np.mean(max_values)

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

def all_atoms_can_be_superimposed(test_angle, key_angle, nested_equivalent_atoms):
    return (test_angle[0] == key_angle[0] or check_in_nested_list([test_angle[0],key_angle[0]], nested_equivalent_atoms)) and (
            test_angle[1] == key_angle[1] or check_in_nested_list([test_angle[1],key_angle[1]], nested_equivalent_atoms)) and (
                    test_angle[2] == key_angle[2] or check_in_nested_list([test_angle[2],key_angle[2]], nested_equivalent_atoms))

def all_atoms_can_be_superimposed_dihedral(test_dihedral, key_dihedral, nested_equivalent_atoms):
    return (test_dihedral[0] == key_dihedral[0] or check_in_nested_list([test_dihedral[0],key_dihedral[0]], nested_equivalent_atoms)) and (
            test_dihedral[1] == key_dihedral[1] or check_in_nested_list([test_dihedral[1],key_dihedral[1]], nested_equivalent_atoms)) and (
                    test_dihedral[2] == key_dihedral[2] or check_in_nested_list([test_dihedral[2],key_dihedral[2]], nested_equivalent_atoms)) and (
                            test_dihedral[3] == key_dihedral[3] or check_in_nested_list([test_dihedral[3],key_dihedral[3]], nested_equivalent_atoms))


def get_symm_angles(angles,specification):
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

def get_angle_subsets(symmetric_angles,num_bonds,num_angles,idof,n_phi) -> list:
    symmetric_angles_list, angles = [], []

    for ind_angle in symmetric_angles.keys():
        if symmetric_angles[ind_angle] not in symmetric_angles_list:
            symmetric_angles_list.append(symmetric_angles[ind_angle])
    
    for i in range(1,len(symmetric_angles_list)+1):
        for angle_subset in itertools.combinations(symmetric_angles_list,i):
            flat_angle_subset = [item for sublist in angle_subset for item in sublist]
            if len(list(flat_angle_subset)) == n_phi:
                angles.append(list(flat_angle_subset))
    return angles

def get_symm_dihedrals(dihedrals,specification):
    symmetric_dihedrals = dict()
    symmetric_dihedrals = {key:[] for (key, val) in Counter(dihedrals).items()}
   
    # symmetric dihedrals equally defined as in get_symm_angles --> make same function?
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

def get_oop_subsets(out_of_plane, n_gamma):
    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if not_same_central_atom(subset):
                oop_subsets.append(list(subset))
    return oop_subsets

def get_dihedral_subsets(symmetric_dihedrals,num_bonds,num_angles,idof,n_tau) -> list:
    symmetric_dihedrals_list, dihedrals = [], []
    for ind_dihedral in symmetric_dihedrals.keys():
        if symmetric_dihedrals[ind_dihedral] not in symmetric_dihedrals_list:
            symmetric_dihedrals_list.append(symmetric_dihedrals[ind_dihedral])
    for i in range(0,len(symmetric_dihedrals_list)+1):
        for dihedral_subset in itertools.combinations(symmetric_dihedrals_list,i):
            flat_dihedral_subset = [item for sublist in dihedral_subset for item in sublist]
            if len(list(flat_dihedral_subset)) == n_tau:
                dihedrals.append(list(flat_dihedral_subset))
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


def get_sets(idof,atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, specification):
    ic_dict = dict() 
    num_bonds = len(bonds)
    num_atoms = len(atoms)

    num_of_red = 6*specification["mu"]
    # @decision tree: linear
    if specification["linearity"] == "fully linear":
        ic_dict = topology.fully_linear_molecule(ic_dict, bonds, angles, linear_angles, out_of_plane, dihedrals)
    
    # @decision tree: planar, acyclic and no linear submolecules 
    if specification["planar"] == "yes" and not specification["linearity"] == "linear submolecules found" and (num_of_red == 0):
        ic_dict = topology.planar_acyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, 
                dihedrals, num_bonds, num_atoms, number_terminal_bonds(specification["multiplicity"]), specification)

    # @decision tree: planar, cyclic and no linear submolecules 
    if specification["planar"] == "yes" and not specification["linearity"] == "linear submolecules found" and (num_of_red != 0):
        ic_dict = topology.planar_cyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, 
                dihedrals, num_bonds, num_atoms, number_terminal_bonds(specification["multiplicity"]), specification)

    # @decision tree: general molecule, acyclic and no linear submolecules
    if specification["planar"] == "no" and not specification["linearity"] == "linear submolecules found" and (num_of_red == 0):
        ic_dict = topology.general_acyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, 
                dihedrals, num_bonds, num_atoms, number_terminal_bonds(specification["multiplicity"]), specification)

    # @decision tree: general molecule, cyclic and no linear submolecules
    if specification["planar"] == "no" and not specification["linearity"] == "linear submolecules found" and (num_of_red != 0):
        ic_dict = topology.general_cyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, 
                dihedrals, num_bonds, num_atoms, num_of_red, number_terminal_bonds(specification["multiplicity"]), specification)

    
    # @decision tree: planar, acyclic molecules with linear submolecules
    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (num_of_red == 0):
        
        ic_dict = topology.planar_acyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane,
                dihedrals, num_bonds, num_atoms, number_terminal_bonds(specification["multiplicity"]), 
                specification["length of linear submolecule(s) l"], specification)

    # @decision tree: planar, cyclic molecules with linear submolecules
    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (num_of_red != 0):
        ic_dict = topology.planar_cyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane,
                dihedrals, num_bonds, num_atoms, number_terminal_bonds(specification["multiplicity"]), 
                specification["length of linear submolecule(s) l"], specification) 

    # @decision tree: general, acyclic molecule with linear submolecules
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found" and (num_of_red == 0): 
        ic_dict = topology.general_acyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane,
                dihedrals, num_bonds, num_atoms, number_terminal_bonds(specification["multiplicity"]), 
                specification["length of linear submolecule(s) l"], specification)


    # @decision tree: general, cyclic molecule with linear submolecules
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found":
        ic_dict = topology.general_cyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane,
                dihedrals, num_bonds, num_atoms, num_of_red, number_terminal_bonds(specification["multiplicity"]), 
                specification["length of linear submolecule(s) l"], specification)

    print(len(ic_dict), "internal coordinate sets were generated.")
    print("The optimal coordinate set will be determined...")
    return ic_dict

