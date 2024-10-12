
"""
This module contains the get_sets method as well as symmetry functions for selecting internal coordinates
"""

import itertools
import logging
from collections import Counter
import numpy as np
from nomodeco.libraries import topology


def Kemalian_metric(matrix, Diag_elements, counter, intfreq_penalty, intfc_penalty, args) -> float:
    """
    Given a matrix calculates the metric implemented in Nomodeco
    """
    # axis = 0, when maximum of each column
    max_values = np.max(matrix, axis=1)

    # set 0 if negative values
    min_values = np.min(matrix, axis=1)
    if np.any(min_values < -1):
        return 0

    # penalty for high fc values
    penalty1 = intfreq_penalty * counter

    # penalty for high fc values
    penalty2 = 0
    max_values_int_fc = np.max(Diag_elements, axis=1)
    for max_fc_value in max_values_int_fc:
        if max_fc_value > 1:
            penalty2 += ((max_fc_value - 1) / 0.1) * intfc_penalty

    return np.mean(max_values) - penalty1 - penalty2


def Kemalian_metric_log(matrix, Diag_elements, counter, intfreq_penalty, intfc_penalty, log):
    # axis = 0, when maximum of each column
    max_values = np.max(matrix, axis=1)

    # set 0 if negative values
    min_values = np.min(matrix, axis=1)
    if np.any(min_values < -1):
        log.info("Negative values in energy distribution matrix found!")
        return 0

    # penalty for high fc values
    penalty1 = intfreq_penalty * counter

    # penalty for high fc values
    penalty2 = 0
    max_values_int_fc = np.max(Diag_elements, axis=1)
    for max_fc_value in max_values_int_fc:
        if max_fc_value > 1:
            penalty2 += ((max_fc_value - 1) / 0.1) * intfc_penalty

    log.info("The following diagonalization parameter, penalty for asymmetric instrinsic frequencies"
             " and penalty for unphysical contributions has been determined: %s, %s, %s ",
             np.around(np.mean(max_values),2), penalty1, penalty2)
    return np.mean(max_values) - penalty1 - penalty2


def are_two_elements_same(tup1, tup2) -> bool:
    """
    For two given tuples checks if they are the same and return True/False
    """
    return ((tup1[0] == tup2[0] and tup1[1] == tup2[1]) or
            (tup1[1] == tup2[1] and tup1[2] == tup2[2]) or
            (tup1[0] == tup2[0] and tup1[2] == tup2[2]))


def get_different_elements(tup1, tup2) -> list:
    """
    Given two tuples calculates the difference in elements and outputs as list
    """
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


def all_atoms_can_be_superimposed_bond(test_bond, key_bond, nested_equivalent_atoms):
    return (test_bond[0] == key_bond[0] or check_in_nested_list([test_bond[0], key_bond[0]],
                                                                nested_equivalent_atoms)) and (
            test_bond[1] == key_bond[1] or check_in_nested_list([test_bond[1], key_bond[1]], nested_equivalent_atoms))


def all_atoms_can_be_superimposed(test_angle, key_angle, nested_equivalent_atoms):
    return (test_angle[0] == key_angle[0] or check_in_nested_list([test_angle[0], key_angle[0]],
                                                                  nested_equivalent_atoms)) and (
            test_angle[1] == key_angle[1] or check_in_nested_list([test_angle[1], key_angle[1]],
                                                                  nested_equivalent_atoms)) and (
            test_angle[2] == key_angle[2] or check_in_nested_list([test_angle[2], key_angle[2]],
                                                                  nested_equivalent_atoms))


def all_atoms_can_be_superimposed_dihedral(test_dihedral, key_dihedral, nested_equivalent_atoms):
    return (test_dihedral[0] == key_dihedral[0] or check_in_nested_list([test_dihedral[0], key_dihedral[0]],
                                                                        nested_equivalent_atoms)) and (
            test_dihedral[1] == key_dihedral[1] or check_in_nested_list([test_dihedral[1], key_dihedral[1]],
                                                                        nested_equivalent_atoms)) and (
            test_dihedral[2] == key_dihedral[2] or check_in_nested_list([test_dihedral[2], key_dihedral[2]],
                                                                        nested_equivalent_atoms)) and (
            test_dihedral[3] == key_dihedral[3] or check_in_nested_list([test_dihedral[3], key_dihedral[3]],
                                                                        nested_equivalent_atoms))


def get_symm_bonds(bonds, specification):
    symmetric_bonds = dict()
    symmetric_bonds = {key: [] for (key, val) in Counter(bonds).items()}

    for i, key in itertools.product(range(len(bonds)), symmetric_bonds):
        symmetric_bonds[key].append(bonds[i])

    for key, val in symmetric_bonds.items():
        i = 0
        while i < len(val):
            bond = val[i]
            if not (all_atoms_can_be_superimposed_bond(bond, key, specification[
                "equivalent_atoms"]) or all_atoms_can_be_superimposed_bond((bond[1], bond[0]), key,
                                                                           specification["equivalent_atoms"])):
                del val[i]
            elif (all_atoms_can_be_superimposed_bond(bond, key, specification[
                "equivalent_atoms"]) or all_atoms_can_be_superimposed_bond((bond[1], bond[0]), key,
                                                                           specification["equivalent_atoms"])):
                i += 1

    return symmetric_bonds


def get_bond_subsets(symmetric_bonds) -> list:
    symmetric_bonds_list = []

    for ind_bond in symmetric_bonds.keys():
        if symmetric_bonds[ind_bond] not in symmetric_bonds_list:
            symmetric_bonds_list.append(symmetric_bonds[ind_bond])

    return symmetric_bonds_list


def get_symm_angles(angles, specification):
    symmetric_angles = dict()
    symmetric_angles = {key: [] for (key, val) in Counter(angles).items()}

    # angles are the same if the atoms can all be superimposed
    # on each other with symmetry operations

    for i, key in itertools.product(range(len(angles)), symmetric_angles):
        symmetric_angles[key].append(angles[i])

    for key, val in symmetric_angles.items():
        i = 0
        while i < len(val):
            ang = val[i]
            if not (all_atoms_can_be_superimposed(ang, key,
                                                  specification["equivalent_atoms"]) or all_atoms_can_be_superimposed(
                    (ang[2], ang[1], ang[0]), key, specification["equivalent_atoms"])):
                del val[i]
            elif (all_atoms_can_be_superimposed(ang, key,
                                                specification["equivalent_atoms"]) or all_atoms_can_be_superimposed(
                    (ang[2], ang[1], ang[0]), key, specification["equivalent_atoms"])):
                i += 1
    return symmetric_angles


def get_angle_subsets(symmetric_angles, num_bonds, num_angles, idof, n_phi) -> list:
    symmetric_angles_list, angles = [], []

    for ind_angle in symmetric_angles.keys():
        if symmetric_angles[ind_angle] not in symmetric_angles_list:
            symmetric_angles_list.append(symmetric_angles[ind_angle])

    for i in range(1, len(symmetric_angles_list) + 1):
        for angle_subset in itertools.combinations(symmetric_angles_list, i):
            flat_angle_subset = [item for sublist in angle_subset for item in sublist]
            if len(list(flat_angle_subset)) == n_phi:
                angles.append(list(flat_angle_subset))

    # allow the inclusion of red
    # if you don't want that ==> uncomment
    if not angles:
        logging.info(
            "In order to obtain symmetry in the angles and hence intrinsic frequencies, inclusion of 1 redundant angle coordinate will be attempted")
        for i in range(1, len(symmetric_angles_list) + 1):
            for angle_subset in itertools.combinations(symmetric_angles_list, i):
                flat_angle_subset = [item for sublist in angle_subset for item in sublist]
                if len(list(flat_angle_subset)) == n_phi + 1:
                    angles.append(list(flat_angle_subset))

    if not angles:
        logging.info(
            "In order to obtain symmetry in the angles and hence intrinsic frequencies, inclusion of 2 redundant angle coordinates will be attempted")
        for i in range(1, len(symmetric_angles_list) + 1):
            for angle_subset in itertools.combinations(symmetric_angles_list, i):
                flat_angle_subset = [item for sublist in angle_subset for item in sublist]
                if len(list(flat_angle_subset)) == n_phi + 2:
                    angles.append(list(flat_angle_subset))
    return angles


def get_symm_dihedrals(dihedrals, specification):
    symmetric_dihedrals = dict()
    symmetric_dihedrals = {key: [] for (key, val) in Counter(dihedrals).items()}

    # symmetric dihedrals equally defined as in get_symm_angles --> make same function?
    for i, key in itertools.product(range(len(dihedrals)), symmetric_dihedrals):
        symmetric_dihedrals[key].append(dihedrals[i])

    for key, val in symmetric_dihedrals.items():
        i = 0
        while i < len(val):
            dihedral = val[i]
            if not (all_atoms_can_be_superimposed_dihedral(dihedral, key, specification["equivalent_atoms"]) or
                    all_atoms_can_be_superimposed_dihedral((dihedral[3], dihedral[2], dihedral[1], dihedral[0]), key,
                                                           specification["equivalent_atoms"])):
                del val[i]
            elif (all_atoms_can_be_superimposed_dihedral(dihedral, key, specification["equivalent_atoms"]) or
                  all_atoms_can_be_superimposed_dihedral((dihedral[3], dihedral[2], dihedral[1], dihedral[0]), key,
                                                         specification["equivalent_atoms"])):
                i += 1
    return symmetric_dihedrals


def get_oop_subsets(out_of_plane, n_gamma):
    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if not_same_central_atom(subset):
            oop_subsets.append(list(subset))
    return oop_subsets


def get_dihedral_subsets(symmetric_dihedrals, num_bonds, num_angles, idof, n_tau) -> list:
    symmetric_dihedrals_list, dihedrals = [], []
    for ind_dihedral in symmetric_dihedrals.keys():
        if symmetric_dihedrals[ind_dihedral] not in symmetric_dihedrals_list:
            symmetric_dihedrals_list.append(symmetric_dihedrals[ind_dihedral])
    for i in range(0, len(symmetric_dihedrals_list) + 1):
        for dihedral_subset in itertools.combinations(symmetric_dihedrals_list, i):
            flat_dihedral_subset = [item for sublist in dihedral_subset for item in sublist]
            if len(list(flat_dihedral_subset)) == n_tau:
                dihedrals.append(list(flat_dihedral_subset))

    # allow the inclusion of red
    # if you don't want that ==> uncomment
    if not dihedrals:
        logging.info(
            "In order to obtain symmetry in the dihedrals and hence intrinsic frequencies, inclusion of 1 redundant dihedral coordinate will be attempted")
        for i in range(0, len(symmetric_dihedrals_list) + 1):
            for dihedral_subset in itertools.combinations(symmetric_dihedrals_list, i):
                flat_dihedral_subset = [item for sublist in dihedral_subset for item in sublist]
                if len(list(flat_dihedral_subset)) == n_tau + 1:
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
    evalue, evect = np.linalg.eigh(
        np.transpose(reciprocal_square_massmatrix) @ CartesianF_Matrix_check @ reciprocal_square_massmatrix)
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

def get_sets(idof, out, atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, specification) -> dict:
    """
    get_sets is the decicion tree of nomodeco where for a given molecular topology all possible sets are generated

    Attributes:
        idof:
            a integer with the vibrational degrees of freedom
        out:
            the output file of nomodeco
        atoms:
            a object of the molecule class
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing linear angles
        out_of_plane:
            a list of tuples containing oop's
        dihedrals:
            a list of tuples containing dihedrals
        specficiation:
            the specification dictionary of nomodeco
    
    Returns:
        a dictionary of all the possible IC sets this is generated using the Topology module
    """ 
    ic_dict = dict()
    num_bonds = len(bonds)
    num_atoms = len(atoms)

    num_of_red = 6 * specification["mu"]
 
    # @decision tree: linear
    if specification["linearity"] == "fully linear" and specification["intermolecular"] == "no":
        ic_dict = topology.fully_linear_molecule(ic_dict, bonds, angles, linear_angles, out_of_plane, dihedrals)

    # @decision tree: planar, acyclic and no linear submolecules 
    if specification["planar"] == "yes" and not specification["linearity"] == "linear submolecules found" and (
            num_of_red == 0) and specification["intermolecular"] == "no":
        ic_dict = topology.planar_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                             dihedrals, num_bonds, num_atoms,
                                                             number_terminal_bonds(specification["multiplicity"]),
                                                             specification)

    # @decision tree: planar, cyclic and no linear submolecules 
    if specification["planar"] == "yes" and not specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "no":
        ic_dict = topology.planar_cyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification)

    # @decision tree: general molecule, acyclic and no linear submolecules
    if specification["planar"] == "no" and not specification["linearity"] == "linear submolecules found" and (
            num_of_red == 0) and specification["intermolecular"] == "no":
        ic_dict = topology.general_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                              dihedrals, num_bonds, num_atoms,
                                                              number_terminal_bonds(specification["multiplicity"]),
                                                              specification)

    # @decision tree: general molecule, cyclic and no linear submolecules
    if specification["planar"] == "no" and not specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "no":
        ic_dict = topology.general_cyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                             dihedrals, num_bonds, num_atoms, num_of_red,
                                                             number_terminal_bonds(specification["multiplicity"]),
                                                             specification)

    # @decision tree: planar, acyclic molecules with linear submolecules
    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (
            num_of_red == 0) and specification["intermolecular"] == "no":
        ic_dict = topology.planar_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                           dihedrals, num_bonds, num_atoms,
                                                           number_terminal_bonds(specification["multiplicity"]),
                                                           specification["length of linear submolecule(s) l"],
                                                           specification)

    # @decision tree: planar, cyclic molecules with linear submolecules
    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "no":
        ic_dict = topology.planar_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                          dihedrals, num_bonds, num_atoms,
                                                          number_terminal_bonds(specification["multiplicity"]),
                                                          specification["length of linear submolecule(s) l"],
                                                          specification)

        # @decision tree: general, acyclic molecule with linear submolecules
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found" and (
            num_of_red == 0) and specification["intermolecular"] == "no":
        ic_dict = topology.general_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification["length of linear submolecule(s) l"],
                                                            specification)

    # @decision tree: general, cyclic molecule with linear submolecules
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "no":
        ic_dict = topology.general_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                           dihedrals, num_bonds, num_atoms, num_of_red,
                                                           number_terminal_bonds(specification["multiplicity"]),
                                                           specification["length of linear submolecule(s) l"],
                                                           specification)


# For intermolecular complexes through determination of connectivity c we can manipulate our specification
# Therefore all the cases on the top get duplicated for the intermolecular systems

# General Molecules:
   
    # This is already done and working
    if specification["planar"] == "no" and specification["linearity"] == "not linear" and (num_of_red != 0) and specification["intermolecular"] == "yes":
       ic_dict = topology.intermolecular_general_cyclic_nolinsub(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms, num_of_red,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification)

    # This is already done and working 
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found" and (num_of_red == 0) and specification["intermolecular"] == "yes":
       ic_dict = topology.intermolecular_general_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification["length of linear submolecule(s) l"],
                                                            specification)
    # This is already done and woring
    if specification["planar"] == "no" and not specification["linearity"] == "linear submolecules found" and (
            num_of_red == 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_general_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                              dihedrals, num_bonds, num_atoms,
                                                              number_terminal_bonds(specification["multiplicity"]),
                                                              specification)
    
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_general_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                           dihedrals, num_bonds, num_atoms, num_of_red,
                                                           number_terminal_bonds(specification["multiplicity"]),
                                                           specification["length of linear submolecule(s) l"],
                                                           specification)
# Linear Molecules:
    
    # basically done just needs a example
    if specification["linearity"] == "fully linear" and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_fully_linear_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                              dihedrals, num_bonds, num_atoms,
                                                              number_terminal_bonds(specification["multiplicity"]),
                                                              specification)

# Planar Molecules


    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_planar_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                          dihedrals, num_bonds, num_atoms,
                                                          number_terminal_bonds(specification["multiplicity"]),
                                                          specification["length of linear submolecule(s) l"],
                                                          specification)
 
    if specification["planar"] == "yes" and specification["linearity"] == "not linear" and (
            num_of_red != 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_planar_cyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification)
    # Here we need the specification of x+y = l-1 with kemal
    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (num_of_red == 0) and specification["intermolecular"] == "yes":
       ic_dict = topology.intermolecular_planar_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                           dihedrals, num_bonds, num_atoms,
                                                           number_terminal_bonds(specification["multiplicity"]),
                                                           specification["length of linear submolecule(s) l"],
                                                           specification)


    if specification["planar"] == "yes" and not specification["linearity"] == "linear submolecules found" and (num_of_red == 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_planar_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                             dihedrals, num_bonds, num_atoms,
                                                             number_terminal_bonds(specification["multiplicity"]),
                                                             specification)

 


    print(len(ic_dict), "internal coordinate sets were generated.")
    print("The optimal coordinate set will be determined...")
    return ic_dict
