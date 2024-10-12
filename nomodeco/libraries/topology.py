"""
This module contains the decision tree derived by Decius' work on complete sets and redundancies among
small vibrational coordinates
(https://doi.org/10.1063/1.1747158)
"""

import networkx as nx
import itertools
import logging
import string
import pandas as pd
import pymatgen.core as mg
import re
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
import os
import matplotlib.pyplot as plt
from nomodeco.libraries import logfile
from nomodeco.libraries import specifications
from nomodeco.libraries import icsel
from nomodeco.libraries.nomodeco_classes import Molecule
from nomodeco.libraries import arguments

"""
Shared Variables for the Intermolecular Functions
"""
# General information about the structure
atoms_list = os.getenv("atoms")

Total_IC_dict = os.getenv("Total_IC_dict")

"""
Arguments for IC Comb
"""


"""
GENERAL PURPOSE FUNCTIONS
"""
out = None


def molecule_is_not_split(bonds) -> bool:
    """
    Uses Networkx to check if for a given list of bonds the molecular structure is connected or split
    """
    G = nx.Graph()
    G.add_edges_from(bonds)
    if len(list(nx.connected_components(G))) == 1:
        return True
    else:
        return False


def strip_numbers(string) -> str:
    """
    Removes the digits from a string e.q H2 -> H
    """
    return "".join([char for char in string if not char.isdigit()])


def eliminate_symmetric_tuples(list_tuples) -> list:
    """
    Eliminates symmetric tuples in a list of tuples e.q (H2,O,H1) (H1,O,H2)
    """
    new_list = []
    seen_tuples = set()
    for tp1 in list_tuples:
        reversed_tp1 = tuple(reversed(tp1))
        if reversed_tp1 not in seen_tuples:
            new_list.append(tp1)
            seen_tuples.add(tp1)
    return new_list


def valide_atoms_to_cut(bonds, multiplicity_list) -> list:
    """
    Checks for a given list of bonds if the multiplicity of the atoms is greater then two. If this is the case return the atoms
    """
    # first of all we will identify where cutting bonds is even making sense to a first degree
    valide_atoms = []
    for tup in multiplicity_list:
        if tup[1] >= 2:
            valide_atoms.append(tup[0])
    return valide_atoms


def bonds_are_in_valide_atoms(symmetric_bond_group, valide_atoms):
    for bond in symmetric_bond_group:
        if not (bond[0] in valide_atoms and bond[1] in valide_atoms):
            return False
    return True


def delete_bonds_symmetry(symmetric_bond_group, bonds, mu, valide_atoms):
    removed_bonds = []
    while mu > 0:
        if not symmetric_bond_group:
            return [], bonds + removed_bonds
        # cut the bonds
        for bond in symmetric_bond_group:
            if bond[0] in valide_atoms and bond[1] in valide_atoms:
                removed_bonds.append(bond)
                bonds.remove(bond)
                symmetric_bond_group.remove(bond)
                break

        # we will check if the molecule is not split;
        if molecule_is_not_split(bonds):
            mu -= 1
        else:
            # if we can not cut bonds out in the symmetric group
            bonds.append(removed_bonds[-1])
            removed_bonds.pop()

    return removed_bonds, bonds + removed_bonds


def delete_bonds(bonds, mu, valide_atoms):
    removed_bonds = []
    while mu > 0:
        # cut the bonds
        for bond in bonds:
            if bond[0] in valide_atoms and bond[1] in valide_atoms:
                removed_bonds.append(bond)
                bonds.remove(bond)
                break

        # we will check if the molecule is not split;
        if molecule_is_not_split(bonds):
            mu -= 1
        else:
            bonds.append(removed_bonds[-1])
            removed_bonds.pop()

    return removed_bonds, bonds


# TODO: rename method
# TODO: possible bug here for cyclic systems, when using [:]!
def update_internal_coordinates_cyclic(removed_bonds, ic_list) -> list:
    """
    Given a set of ICs and a removed bond, eliminates this particular IC out of the list

    Attributes:
        removed_bonds: list
            a list of tuples containing the bonds to remove
        ic_list: list
            a list of internal coordinates of a specific type
    """
    ic_list_dup = ic_list[:]
    for bond in removed_bonds:
        for ic in ic_list_dup[:]:
            if bond[0] in ic and bond[1] in ic:
                ic_list_dup.remove(ic)
    return ic_list_dup


def find_common_index_with_most_subsets(list1, list2):
    combined_lengths = [
        len(sublist1) + len(sublist2) for sublist1, sublist2 in zip(list1, list2)
    ]
    max_combined_length = max(combined_lengths)
    common_index = combined_lengths.index(max_combined_length)
    return common_index


def find_common_index_with_least_subsets(list1, list2):
    combined_lengths = [
        len(sublist1) + len(sublist2) for sublist1, sublist2 in zip(list1, list2)
    ]
    max_combined_length = min(combined_lengths)
    common_index = combined_lengths.index(max_combined_length)
    return common_index


def remove_angles(atom_and_mult, angles):
    num_angles_to_be_removed = atom_and_mult[1] - 2
    for angle in angles:
        if angle[1] == atom_and_mult[0] and num_angles_to_be_removed >= 1:
            angles.remove(angle)
            num_angles_to_be_removed -= 1
    return angles


def get_param_planar_submolecule(planar_subunits_list, multiplicity_list, angles):
    n_phi = 0
    n_gamma = 0
    for atom_and_mult in multiplicity_list:
        if atom_and_mult[1] > 1:
            if atom_and_mult in planar_subunits_list:
                n_phi += atom_and_mult[1] - 1
                n_gamma += atom_and_mult[1] - 2
                angles = remove_angles(atom_and_mult, angles)
            else:
                n_phi += 2 * atom_and_mult[1] - 3
    return n_phi, n_gamma, angles


def get_multiplicity(atom_name, multiplicity_list):
    for atom_and_mult in multiplicity_list:
        if atom_name == atom_and_mult[0]:
            return atom_and_mult[1]


"""''
LINEAR SYSTEMS
""" ""


def fully_linear_molecule(
    ic_dict, bonds, angles, linear_angles, out_of_plane, dihedrals
) -> dict:
    """
    Generates the IC sets for the covalent fully linear conformation of a structure

    Attributes: 
        ic_dict:
            a dictionary containing IC_sets
        bonds:
            a list of tuples with bonds
        angles:
            a list of tuples with angles
        linear_angles:
            a list of tuples with linear angles
        out_of_plane:
            a list of tuples with oop's
        dihedrals:
            a list of tuples with dihedrals
    
    Returns:
        a IC dictionary where each entry is a valid IC set for the further analysis
    """
    # purely linear molecules do not have oop, one can define dihedrals, but they are not significant as
    # the intrinsic frequency equals 0 for them
    ic_dict[0] = {
        "bonds": bonds,
        "angles": angles,
        "linear valence angles": linear_angles,
        "out of plane angles": out_of_plane,
        "dihedrals": [],
    }
    return ic_dict


"""'
PLANAR SYSTEMS
""" ""


def planar_acyclic_nolinunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    specification,
) -> dict:
    """
    Generates IC sets for the planar acyclic case without any linear subunits

    Attributes:
        ic_dict:
            a dictionary with IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples with bonds
        angles:
            a list of tuples with angles
        linear_angles: 
            a list of tuples with linear angles
        out_of_plane:
            a list of tuples with oop's
        dihedrals:
            a list of tuples with dihedrals
        num_bonds:
            a integer containing the number of atoms
        a_1:
            a integer containing the number of terminal atoms
        specification:
            the specification for the given molecule (see specifcation documentation)

    Returns:
        all possible IC sets for this conformation in a dictionary, where each entry is a IC set

    """
    # set length of subsets
    n_r = num_bonds
    n_phi = 2 * num_bonds - num_atoms
    n_gamma = 2 * (num_bonds - num_atoms) + a_1
    n_tau = num_bonds - a_1

    # remove angles that are at specific oop spots
    oop_central_atoms = []
    for oop in out_of_plane:
        if oop[0] not in oop_central_atoms:
            oop_central_atoms.append(oop[0])
    for oop_central_atom in oop_central_atoms:
        angles = remove_angles(
            (
                oop_central_atom,
                get_multiplicity(oop_central_atom, specification["multiplicity"]),
            ),
            angles,
        )

    symmetric_angles = icsel.get_symm_angles(angles, specification)
    angle_subsets = icsel.get_angle_subsets(
        symmetric_angles, len(bonds), len(angles), idof, n_phi
    )

    # in cases of restrictive symmetry, we need to break angle symmetry
    if len(angle_subsets) == 0:
        logging.warning(
            "In order to obtain angle subsets, symmetry needs to be broken!"
        )
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset))

            # the if statement ensures, that oop angles to the same central atom can not be in the same set
    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset):
            oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals, specification)
    dihedral_subsets = icsel.get_dihedral_subsets(
        symmetric_dihedrals, len(bonds), len(angles), idof, n_tau
    )

    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning(
            "In order to obtain dihedral subsets, symmetry needs to be broken!"
        )
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))

    k = 0
    for len_angles in range(0, len(angle_subsets)):
        for len_oop in range(0, len(oop_subsets)):
            for len_dihedrals in range(0, len(dihedral_subsets)):
                ic_dict[k] = {
                    "bonds": bonds,
                    "angles": angle_subsets[len_angles],
                    "linear valence angles": [],
                    "out of plane angles": oop_subsets[len_oop],
                    "dihedrals": dihedral_subsets[len_dihedrals],
                }
                k += 1

    return ic_dict


# in the cyclic cases, the molecule is rendered acyclic (by removing random but possible mu bonds) and then the acylic function is called
def planar_cyclic_nolinunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the planar,cyclic case without linear submolecules
    
    Attributes:
        ic_dict:
            a dictionary with IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples with bonds
        angles:
            a list of tuples with angles
        linear_angles: 
            a list of tuples with linear angles
        out_of_plane:
            a list of tuples with oop's
        dihedrals:
            a list of tuples with dihedrals
        num_bonds:
            a integer containing the number of atoms
        a_1:
            a integer containing the number of terminal atoms
        specification:
            the specification for the given molecule (see specifcation documentation) 
    
    Returns:
        a dictionary where each entry is a valid IC set, note that in this case bonds are cut
    """
    # remove bonds without destroying the molecule
    # if there are several classes of symmetric bonds, we need to remove the corresponding
    # symmetric bonds

    symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
    symmetric_bonds_list = icsel.get_bond_subsets(symmetric_bonds)
    valide_atoms = valide_atoms_to_cut(bonds, specification["multiplicity"])
    ic_dict_list = []
    removed_bonds = []
    for symmetric_bond_group in symmetric_bonds_list:
        if len(symmetric_bond_group) >= specification[
            "mu"
        ] and bonds_are_in_valide_atoms(symmetric_bond_group, valide_atoms):
            removed_bonds, bonds = delete_bonds_symmetry(
                symmetric_bond_group, bonds, specification["mu"], valide_atoms
            )
        if not removed_bonds:
            continue

        # update bonds, angles, oop, and dihedrals to not include the coordinates that were removed
        bonds_updated = update_internal_coordinates_cyclic(removed_bonds, bonds)
        angles_updated = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane_updated = update_internal_coordinates_cyclic(
            removed_bonds, out_of_plane
        )
        dihedrals_updated = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(
            out,
            bonds_updated,
            angles_updated,
            linear_angles,
            out_of_plane_updated,
            dihedrals_updated,
            specification,
        )

        # we need to do some pre-calc of the symmetric angles and dihedrals etc. sadly
        # so that we do not sample a subspace, which is not feasible

        symmetric_angles = icsel.get_symm_angles(angles_updated, specification)
        angle_subsets = icsel.get_angle_subsets(
            symmetric_angles,
            len(bonds_updated),
            len(angles_updated),
            idof,
            2 * len(bonds_updated) - num_atoms,
        )
        if len(angle_subsets) == 0:
            logging.warning(
                "For this rendered molecule, angle symmetry can not be considered and hence this subspace of internal coordinates will be skipped"
            )
            continue

        symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals_updated, specification)
        dihedral_subsets = icsel.get_dihedral_subsets(
            symmetric_dihedrals,
            len(bonds_updated),
            len(angles_updated),
            idof,
            len(bonds_updated) - a_1,
        )
        if len(dihedral_subsets) == 0:
            logging.warning(
                "For this rendered molecule, dihedral symmetry can not be considered and hence this subspace of internal coordintes will be skipped"
            )
            continue

        # call the acyclic version

        ic_dict_list.append(
            planar_acyclic_nolinunit_molecule(
                dict(),
                out,
                idof,
                bonds_updated,
                angles_updated,
                linear_angles,
                out_of_plane_updated,
                dihedrals_updated,
                len(bonds_updated),
                num_atoms,
                a_1,
                specification,
            )
        )
        removed_bonds = []

    # if we can't cut according to symmetry, do random cutting
    # cut symmetry out if you want, by commenting everyting out
    if not ic_dict_list:
        removed_bonds, bonds = delete_bonds(bonds, specification["mu"], valide_atoms)
        angles = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(
            out, bonds, angles, linear_angles, out_of_plane, dihedrals, specification
        )

        ic_dict = planar_acyclic_nolinunit_molecule(
            ic_dict,
            out,
            idof,
            bonds,
            angles,
            linear_angles,
            out_of_plane,
            dihedrals,
            len(bonds),
            num_atoms,
            a_1,
            specification,
        )
        return ic_dict

    else:
        ic_dict = dict()
        new_key = 0
        for dictionary in ic_dict_list:
            for key, value in dictionary.copy().items():
                ic_dict[new_key] = value
                new_key += 1
        return ic_dict


def planar_acyclic_linunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    l,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the planar,acyclic case with linear submoleculesd

    Attributes:
        ic_dict:
            a dictionary with IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples with bonds
        angles:
            a list of tuples with angles
        linear_angles: 
            a list of tuples with linear angles
        out_of_plane:
            a list of tuples with oop's
        dihedrals:
            a list of tuples with dihedrals
        num_bonds:
            a integer containing the number of atoms
        a_1:
            a integer containing the number of terminal atoms
        l:
            a integer containing the number of linear bonds
        specification:
            the specification for the given molecule (see specifcation documentation)

    Returns:
        a dictionary where each entry is a valid IC set
    """
    # set length of subsets
    # IMPORTANT: there is a distinction to Decius work -> Decius counts one l.A. in n_phi and one in n_phi'!

    n_r = num_bonds
    n_phi = 2 * num_bonds - num_atoms - (l - 1)
    n_phi_prime = 2 * (l - 1)

    # remove angles that are at specific oop spots
    oop_central_atoms = []
    for oop in out_of_plane:
        if oop[0] not in oop_central_atoms:
            oop_central_atoms.append(oop[0])
    for oop_central_atom in oop_central_atoms:
        angles = remove_angles(
            (
                oop_central_atom,
                get_multiplicity(oop_central_atom, specification["multiplicity"]),
            ),
            angles,
        )

    n_gamma = 2 * (num_bonds - num_atoms) + a_1
    n_tau = num_bonds - a_1

    # before computing the number of ICs we will remove all oop that are associated with this linear angle
    # also remove dihedrals,if they are terminal ==> we will then also reset the number of internals

    linear_bonds = specifications.get_linear_bonds(linear_angles)

    # remove dihedrals if terminal
    for linear_bond in linear_bonds:
        if (
            get_multiplicity(linear_bond[0], specification["multiplicity"]) == 1
            or get_multiplicity(linear_bond[1], specification["multiplicity"]) == 1
        ):
            dihedrals = update_internal_coordinates_cyclic([linear_bond], dihedrals)
            n_tau -= l - 1

    # corret n_gamma if it defined for a linear submolecule
    # as we have 3 oop angles per central unit we need to divide by 3!
    out_of_plane_updated = update_internal_coordinates_cyclic(
        linear_bonds, out_of_plane
    )
    n_gamma = n_gamma - ((len(out_of_plane) - len(out_of_plane_updated)) // 3)
    out_of_plane = out_of_plane_updated

    logfile.write_logfile_updatedICs_linunit(out, out_of_plane, dihedrals)

    symmetric_angles = icsel.get_symm_angles(angles, specification)
    angle_subsets = icsel.get_angle_subsets(
        symmetric_angles, len(bonds), len(angles), idof, n_phi
    )

    # symmetry breaking if needed
    if len(angle_subsets) == 0:
        logging.warning(
            "In order to obtain angle subsets, symmetry needs to be broken!"
        )
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset))

    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset):
            oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals, specification)
    dihedral_subsets = icsel.get_dihedral_subsets(
        symmetric_dihedrals, len(bonds), len(angles), idof, n_tau
    )

    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning(
            "In order to obtain dihedral subsets, symmetry needs to be broken!"
        )
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))

    k = 0
    for len_angles in range(0, len(angle_subsets)):
        for len_oop in range(0, len(oop_subsets)):
            for len_dihedrals in range(0, len(dihedral_subsets)):
                ic_dict[k] = {
                    "bonds": bonds,
                    "angles": angle_subsets[len_angles],
                    "linear valence angles": linear_angles,
                    "out of plane angles": oop_subsets[len_oop],
                    "dihedrals": dihedral_subsets[len_dihedrals],
                }
                k += 1

    return ic_dict


def planar_cyclic_linunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    l,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the planar,cyclic case with linear submolecules

    Attributes:
        ic_dict:
            a dictionary with IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples with bonds
        angles:
            a list of tuples with angles
        linear_angles: 
            a list of tuples with linear angles
        out_of_plane:
            a list of tuples with oop's
        dihedrals:
            a list of tuples with dihedrals
        num_bonds:
            a integer containing the number of atoms
        a_1:
            a integer containing the number of terminal atoms
        l:
            a integer containing the number of linear bonds
        specification:
            the specification for the given molecule (see specifcation documentation)
    
    Returns:
        a dictionary where each entry is a valid IC set
    """
    # remove bonds without destroying the molecule
    # if there are several classes of symmetric bonds, we need to remove the corresponding
    # symmetric bonds

    symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
    symmetric_bonds_list = icsel.get_bond_subsets(symmetric_bonds)
    valide_atoms = valide_atoms_to_cut(bonds, specification["multiplicity"])
    ic_dict_list = []
    removed_bonds = []
    for symmetric_bond_group in symmetric_bonds_list:
        if len(symmetric_bond_group) >= specification[
            "mu"
        ] and bonds_are_in_valide_atoms(symmetric_bond_group, valide_atoms):
            removed_bonds, bonds = delete_bonds_symmetry(
                symmetric_bond_group, bonds, specification["mu"], valide_atoms
            )
        if not removed_bonds:
            continue

        # update bonds, angles, oop, and dihedrals to not include the coordinates that were removed
        bonds_updated = update_internal_coordinates_cyclic(removed_bonds, bonds)
        angles_updated = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane_updated = update_internal_coordinates_cyclic(
            removed_bonds, out_of_plane
        )
        dihedrals_updated = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(
            out,
            bonds_updated,
            angles_updated,
            linear_angles,
            out_of_plane_updated,
            dihedrals_updated,
            specification,
        )

        # we need to do some pre-calc of the symmetric angles and dihedrals etc. sadly
        # so that we do not sample a subspace, which is not feasible

        symmetric_angles = icsel.get_symm_angles(angles_updated, specification)
        angle_subsets = icsel.get_angle_subsets(
            symmetric_angles,
            len(bonds_updated),
            len(angles_updated),
            idof,
            2 * len(bonds_updated) - num_atoms,
        )
        if len(angle_subsets) == 0:
            logging.warning(
                "For this rendered molecule, angle symmetry can not be considered and hence this subspace of internal coordinates will be skipped"
            )
            continue

        symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals_updated, specification)
        dihedral_subsets = icsel.get_dihedral_subsets(
            symmetric_dihedrals,
            len(bonds_updated),
            len(angles_updated),
            idof,
            len(bonds_updated) - a_1,
        )
        if len(dihedral_subsets) == 0:
            logging.warning(
                "For this rendered molecule, dihedral symmetry can not be considered and hence this subspace of internal coordintes will be skipped"
            )
            continue

        # call the acyclic version #TODO Help there was a l missing
        ic_dict_list.append(
            planar_acyclic_linunit_molecule(
                dict(),
                out,
                idof,
                bonds_updated,
                angles_updated,
                linear_angles,
                out_of_plane_updated,
                dihedrals_updated,
                len(bonds_updated),
                num_atoms,
                a_1,
                l,
                specification,
            )
        )
        removed_bonds = []

    # if we can't cut according to symmetry, do random cutting
    # cut symmetry out if you want, by commenting everyting out
    if not ic_dict_list:
        removed_bonds, bonds = delete_bonds(bonds, specification["mu"], valide_atoms)
        angles = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(
            out, bonds, angles, linear_angles, out_of_plane, dihedrals, specification
        )

        ic_dict = planar_acyclic_linunit_molecule(
            ic_dict,
            out,
            idof,
            bonds,
            angles,
            linear_angles,
            out_of_plane,
            dihedrals,
            len(bonds),
            num_atoms,
            a_1,
            l,
            specification,
        )
        return ic_dict

    else:
        ic_dict = dict()
        new_key = 0
        for dictionary in ic_dict_list:
            for key, value in dictionary.copy().items():
                ic_dict[new_key] = value
                new_key += 1
        return ic_dict


"""''
GENERAL SYSTEMS
""" ""


def general_acyclic_nolinunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the general acyclic case without linear submolecules

    Attributes:
        ic_dict:
            a dictionary with IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples with bonds
        angles:
            a list of tuples with angles
        linear_angles: 
            a list of tuples with linear angles
        out_of_plane:
            a list of tuples with oop's
        dihedrals:
            a list of tuples with dihedrals
        num_bonds:
            a integer containing the number of atoms
        num_atoms:
            contains the number of bonds of the structure
        a_1:
            a integer containing the number of terminal atoms
        specification:
            the specification for the given molecule (see specifcation documentation)
    Returns:
        a dictionary where each entry is a IC set
    """
    # set length of subsets
    n_r = num_bonds
    n_phi = 4 * num_bonds - 3 * num_atoms + a_1
    n_gamma = 0
    planar_subunits_list = specification["planar submolecule(s)"]
    n_tau = num_bonds - a_1

    # if planar subunits exist, we need to do 2 things: change n_phi and n_gamma;
    # remove angles at the specified coordinate, as we else would have linear dependencies
    if len(planar_subunits_list) != 0:
        n_phi, n_gamma, angles = get_param_planar_submolecule(
            planar_subunits_list, specification["multiplicity"], angles
        )

    symmetric_angles = icsel.get_symm_angles(angles, specification)
    angle_subsets = icsel.get_angle_subsets(
        symmetric_angles, len(bonds), len(angles), idof, n_phi
    )
    if len(angle_subsets) == 0:
        logging.warning(
            "In order to obtain angle subsets, symmetry needs to be broken!"
        )
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset))

    oop_subsets = []

    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset):
            oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals, specification)
    dihedral_subsets = icsel.get_dihedral_subsets(
        symmetric_dihedrals, len(bonds), len(angles), idof, n_tau
    )

    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning(
            "In order to obtain dihedral subsets, symmetry needs to be broken!"
        )
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))

    k = 0
    for len_angles in range(0, len(angle_subsets)):
        for len_oop in range(0, len(oop_subsets)):
            for len_dihedrals in range(0, len(dihedral_subsets)):
                ic_dict[k] = {
                    "bonds": bonds,
                    "angles": angle_subsets[len_angles],
                    "linear valence angles": [],
                    "out of plane angles": oop_subsets[len_oop],
                    "dihedrals": dihedral_subsets[len_dihedrals],
                }
                k += 1

    return ic_dict


def general_cyclic_nolinunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    num_of_red,
    a_1,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the general cyclic case without linear submolecules

    Attributes:
        ic_dict:
            a dictionary with IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples with bonds
        angles:
            a list of tuples with angles
        linear_angles: 
            a list of tuples with linear angles
        out_of_plane:
            a list of tuples with oop's
        dihedrals:
            a list of tuples with dihedrals
        num_bonds:
            a integer containing the number of bonds
        a_1:
            a integer containing the number of terminal atoms
        num_atoms:
            a integer containing the number of atoms
        num_of_red:
            a integer containg the redundancies of a given molecular structure
        specification:
            the specification for the given molecule (see specifcation documentation)
    Returns:
        a dictionary where each entry is a valid IC set
    """
    # remove bonds without destroying the molecule
    # if there are several classes of symmetric bonds, we need to remove the corresponding
    # symmetric bonds
    symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
    symmetric_bonds_list = icsel.get_bond_subsets(symmetric_bonds)
    valide_atoms = valide_atoms_to_cut(bonds, specification["multiplicity"])
    ic_dict_list = []
    removed_bonds = []
    for symmetric_bond_group in symmetric_bonds_list:
        if len(symmetric_bond_group) >= specification[
            "mu"
        ] and bonds_are_in_valide_atoms(symmetric_bond_group, valide_atoms):
            removed_bonds, bonds = delete_bonds_symmetry(
                symmetric_bond_group, bonds, specification["mu"], valide_atoms
            )
        if not removed_bonds:
            continue

        # update bonds, angles, oop, and dihedrals to not include the coordinates that were removed
        bonds_updated = update_internal_coordinates_cyclic(removed_bonds, bonds)
        angles_updated = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane_updated = update_internal_coordinates_cyclic(
            removed_bonds, out_of_plane
        )
        dihedrals_updated = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(
            out,
            bonds_updated,
            angles_updated,
            linear_angles,
            out_of_plane_updated,
            dihedrals_updated,
            specification,
        )

        # we need to do some pre-calc of the symmetric angles and dihedrals etc. sadly
        # so that we do not sample a subspace, which is not feasible

        symmetric_angles = icsel.get_symm_angles(angles_updated, specification)
        angle_subsets = icsel.get_angle_subsets(
            symmetric_angles,
            len(bonds_updated),
            len(angles_updated),
            idof,
            2 * len(bonds_updated) - num_atoms,
        )
        if len(angle_subsets) == 0:
            logging.warning(
                "For this rendered molecule, angle symmetry can not be considered and hence this subspace of internal coordinates will be skipped"
            )
            continue

        symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals_updated, specification)
        dihedral_subsets = icsel.get_dihedral_subsets(
            symmetric_dihedrals,
            len(bonds_updated),
            len(angles_updated),
            idof,
            len(bonds_updated) - a_1,
        )
        if len(dihedral_subsets) == 0:
            logging.warning(
                "For this rendered molecule, dihedral symmetry can not be considered and hence this subspace of internal coordintes will be skipped"
            )
            continue

        # call the acyclic version

        ic_dict_list.append(
            general_acyclic_nolinunit_molecule(
                dict(),
                out,
                idof,
                bonds_updated,
                angles_updated,
                linear_angles,
                out_of_plane_updated,
                dihedrals_updated,
                len(bonds_updated),
                num_atoms,
                a_1,
                specification,
            )
        )
        removed_bonds = []

    # if we can't cut according to symmetry, do random cutting
    # cut symmetry out if you want, by commenting everything out
    if not ic_dict_list:
        removed_bonds, bonds = delete_bonds(bonds, specification["mu"], valide_atoms)
        angles = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(
            out, bonds, angles, linear_angles, out_of_plane, dihedrals, specification
        )

        ic_dict = general_acyclic_nolinunit_molecule(
            ic_dict,
            out,
            idof,
            bonds,
            angles,
            linear_angles,
            out_of_plane,
            dihedrals,
            len(bonds),
            num_atoms,
            a_1,
            specification,
        )
        return ic_dict

    else:
        ic_dict = dict()
        new_key = 0
        for dictionary in ic_dict_list:
            for key, value in dictionary.copy().items():
                ic_dict[new_key] = value
                new_key += 1
        return ic_dict


def general_acyclic_linunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    l,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the general,acyclic case where linear submolecules are in the structure

    Attributes:
        ic_dict:
            a dictionary with IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples with bonds
        angles:
            a list of tuples with angles
        linear_angles: 
            a list of tuples with linear angles
        out_of_plane:
            a list of tuples with oop's
        dihedrals:
            a list of tuples with dihedrals
        num_bonds:
            a integer containing the number of bonds
        num_atoms:
            a integer containing the number of atoms
        l:
            a integer with the number of linear bonds
        a_1:
            a integer containing the number of terminal atoms
        specification:
            the specification for the given molecule (see specifcation documentation)
    """
    # set length of subsets
    n_r = num_bonds
    n_phi = 4 * num_bonds - 3 * num_atoms + a_1 - (l - 1)
    n_gamma = 0
    planar_subunits_list = specification["planar submolecule(s)"]
    n_phi_prime = 2 * (l - 1)
    n_tau = num_bonds - a_1 - (l - 1)

    # occurs for SF6
    if n_tau < 0 or n_phi < 0:
        logging.warning(
            "Due to high number of linear angles, the topology conditions can not be considered. Linear angles will be removed from the analysis to do so."
        )
        ic_dict = general_acyclic_nolinunit_molecule(
            ic_dict,
            idof,
            bonds,
            angles,
            [],
            out_of_plane,
            dihedrals,
            num_bonds,
            num_atoms,
            a_1,
            specification,
        )
        return ic_dict

    # if planar subunits exist, we need to do 2 things: change n_phi and n_gamma;
    # remove angles at the specified coordinate, as we else would have linear dependencies
    if len(planar_subunits_list) != 0:
        n_phi, n_gamma, angles = get_param_planar_submolecule(
            planar_subunits_list, specification["multiplicity"], angles
        )
        # correct n_phi because we lose (l-1) DOF
        n_phi = n_phi - (l - 1)

    symmetric_angles = icsel.get_symm_angles(angles, specification)
    angle_subsets = icsel.get_angle_subsets(
        symmetric_angles, len(bonds), len(angles), idof, n_phi
    )
    if len(angle_subsets) == 0:
        logging.warning(
            "In order to obtain angle subsets, symmetry needs to be broken!"
        )
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset))

            # before computing the number of ICs we will remove all oop that are associated with this linear angle
    # also remove dihedrals,if they are terminal

    linear_bonds = specifications.get_linear_bonds(linear_angles)
    out_of_plane = update_internal_coordinates_cyclic(linear_bonds, out_of_plane)
    for linear_bond in linear_bonds:
        if (
            get_multiplicity(linear_bond[0], specification["multiplicity"]) == 1
            or get_multiplicity(linear_bond[1], specification["multiplicity"]) == 1
        ):
            dihedrals = update_internal_coordinates_cyclic([linear_bond], dihedrals)

    logfile.write_logfile_updatedICs_linunit(out, out_of_plane, dihedrals)

    # Uncomment, if you want to sample internal coordinates as well
    # Beware: This will lead to high combinatorics!
    # symmetric_lin_angles = icsel.get_symm_angles(linear_angles,specification)
    # lin_angle_subsets = icsel.get_angle_subsets(symmetric_lin_angles, len(bonds), len(angles),idof,n_phi_prime)

    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset):
            oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals, specification)
    dihedral_subsets = icsel.get_dihedral_subsets(
        symmetric_dihedrals, len(bonds), len(angles), idof, n_tau
    )

    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning(
            "In order to obtain dihedral subsets, symmetry needs to be broken!"
        )
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))

    k = 0
    for len_angles in range(0, len(angle_subsets)):
        for len_oop in range(0, len(oop_subsets)):
            for len_dihedrals in range(0, len(dihedral_subsets)):
                ic_dict[k] = {
                    "bonds": bonds,
                    "angles": angle_subsets[len_angles],
                    "linear valence angles": linear_angles,
                    "out of plane angles": oop_subsets[len_oop],
                    "dihedrals": dihedral_subsets[len_dihedrals],
                }
                k += 1

    return ic_dict


def general_cyclic_linunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    num_of_red,
    a_1,
    l,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the general, cyclic case with linear submolecules

    Attributes:
        ic_dict:
            a dictionary with IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples with bonds
        angles:
            a list of tuples with angles
        linear_angles: 
            a list of tuples with linear angles
        out_of_plane:
            a list of tuples with oop's
        dihedrals:
            a list of tuples with dihedrals
        num_bonds:
            a integer containing the number of bonds
        num_atoms:
            a integer containing the number of atoms
        a_1:
            a integer containing the number of terminal atoms
        l:
            a integer containing the number of linear bonds
        specification:
            the specification for the given molecule (see specifcation documentation)
    Returns:
        a dictionary where each entry is a valid IC set
    """
    # remove bonds without destroying the molecule
    # if there are several classes of symmetric bonds, we need to remove the corresponding
    # symmetric bonds

    symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
    symmetric_bonds_list = icsel.get_bond_subsets(symmetric_bonds)
    valide_atoms = valide_atoms_to_cut(bonds, specification["multiplicity"])
    ic_dict_list = []
    removed_bonds = []
    for symmetric_bond_group in symmetric_bonds_list:
        if len(symmetric_bond_group) >= specification[
            "mu"
        ] and bonds_are_in_valide_atoms(symmetric_bond_group, valide_atoms):
            removed_bonds, bonds = delete_bonds_symmetry(
                symmetric_bond_group, bonds, specification["mu"], valide_atoms
            )
        if not removed_bonds:
            continue

        # update bonds, angles, oop, and dihedrals to not include the coordinates that were removed
        bonds_updated = update_internal_coordinates_cyclic(removed_bonds, bonds)
        angles_updated = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane_updated = update_internal_coordinates_cyclic(
            removed_bonds, out_of_plane
        )
        dihedrals_updated = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(
            out,
            bonds_updated,
            angles_updated,
            linear_angles,
            out_of_plane_updated,
            dihedrals_updated,
            specification,
        )

        # we need to do some pre-calc of the symmetric angles and dihedrals etc. sadly
        # so that we do not sample a subspace, which is not feasible

        symmetric_angles = icsel.get_symm_angles(angles_updated, specification)
        angle_subsets = icsel.get_angle_subsets(
            symmetric_angles,
            len(bonds_updated),
            len(angles_updated),
            idof,
            2 * len(bonds_updated) - num_atoms,
        )
        if len(angle_subsets) == 0:
            logging.warning(
                "For this rendered molecule, angle symmetry can not be considered and hence this subspace of internal coordinates will be skipped"
            )
            continue

        symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals_updated, specification)
        dihedral_subsets = icsel.get_dihedral_subsets(
            symmetric_dihedrals,
            len(bonds_updated),
            len(angles_updated),
            idof,
            len(bonds_updated) - a_1,
        )
        if len(dihedral_subsets) == 0:
            logging.warning(
                "For this rendered molecule, dihedral symmetry can not be considered and hence this subspace of internal coordintes will be skipped"
            )
            continue

        # call the acyclic version

        ic_dict_list.append(
            general_acyclic_linunit_molecule(
                dict(),
                out,
                idof,
                bonds_updated,
                angles_updated,
                linear_angles,
                out_of_plane_updated,
                dihedrals_updated,
                len(bonds_updated),
                num_atoms,
                a_1,
                specification,
            )
        )
        removed_bonds = []

    # if we can't cut according to symmetry, do random cutting
    # cut symmetry out if you want, by commenting everyting out
    if not ic_dict_list:
        removed_bonds, bonds = delete_bonds(bonds, specification["mu"], valide_atoms)
        angles = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(
            out, bonds, angles, linear_angles, out_of_plane, dihedrals, specification
        )

        ic_dict = general_acyclic_linunit_molecule(
            ic_dict,
            out,
            idof,
            bonds,
            angles,
            linear_angles,
            out_of_plane,
            dihedrals,
            len(bonds),
            num_atoms,
            a_1,
            specification,
        )
        return ic_dict

    else:
        ic_dict = dict()
        new_key = 0
        for dictionary in ic_dict_list:
            for key, value in dictionary.copy().items():
                ic_dict[new_key] = value
                new_key += 1
        return ic_dict


"""''
Intermolecular Systems
""" ""
# Key goal for intermolecular systems is to cut the molecule into different submolecules
# For this we need to split up the molecule into the sperate parts and then use the topology descirption


def detect_submolecules(bonds):
    molecular_graph = nx.Graph()
    molecular_graph.add_edges_from(bonds)

    connected_components = list(nx.connected_components(molecular_graph))

    # First function extracts the Submolecules as a list of bonds
    submolecules = []
    for component in connected_components:
        subgraph = molecular_graph.subgraph(component)
        submolecule_bonds = list(subgraph.edges)
        submolecules.append(submolecule_bonds)

    return connected_components, submolecules


def extract_atoms_of_submolecules(connected_components, atoms_list):
    submolecule_atoms = []
    for component in connected_components:
        atom_objects = [
            atom_obj for atom_obj in atoms_list if atom_obj.symbol in component
        ]
        submolecule_atoms.append(atom_objects)
    return submolecule_atoms


# The combine dictionary function has the following task:
# Evaluate the smaller sub dictionary and make a element wise combination


def combine_dictionaries(dict_list) -> list:
    """
    Combines two submolecule dictionaries where in each dictionary all possible IC sets are contained. This combiniation is based on combinatorics

    The length of each dictionary determines how the combination takes place

    Attributes:
        dict_list:
            a list of dictionarys
    Returns:
        a list that contains the combined IC sets for both covalent submolecules
    """
    # First we determine the length of the list
    length_of_dict = []
    for d in dict_list:
        length_of_dict.append(len(d.items()))
    result = {}
    # FIRST-CASE:
    # Two Submolecules and one Dictionary of Length 1
    if (
        len(length_of_dict) == 2
        and min(length_of_dict) == 1
        and max(length_of_dict) != 1
    ):
        # In this case we evaluate the position of the first dictionary and then merge it entry wise
        length_of_small_dictionary = min(length_of_dict)
        length_of_big_dictionary = max(length_of_dict)
        small_dictionary = dict_list[length_of_dict.index(length_of_small_dictionary)]
        big_dictionary = dict_list[length_of_dict.index(length_of_big_dictionary)]
        idx = 0
        inner_dict_result = {}
        for outer_key in big_dictionary.keys():
            inner_dict_result["bonds"] = (
                big_dictionary[outer_key]["bonds"] + small_dictionary[0]["bonds"]
            )
            inner_dict_result["angles"] = (
                big_dictionary[outer_key]["angles"] + small_dictionary[0]["angles"]
            )
            inner_dict_result["linear valence angles"] = (
                big_dictionary[outer_key]["linear valence angles"]
                + small_dictionary[0]["linear valence angles"]
            )
            inner_dict_result["out of plane angles"] = (
                big_dictionary[outer_key]["out of plane angles"]
                + small_dictionary[0]["out of plane angles"]
            )
            inner_dict_result["dihedrals"] = (
                big_dictionary[outer_key]["dihedrals"]
                + small_dictionary[0]["dihedrals"]
            )
            result.setdefault(idx, inner_dict_result)
            inner_dict_result = (
                {}
            )  # After the values are unloaded we set the dict to {} i think this is because of some weird python stuff
            idx += 1
    # SECOND SPECIAL CASE
    # Both Submolecule_Dicts have length 1
    if (
        len(length_of_dict) == 2
        and min(length_of_dict) == 1
        and max(length_of_dict) == 1
    ):
        # we can arbitrary choose the small and the big one
        small_dictionary = dict_list[0]
        big_dictionary = dict_list[1]
        inner_dict_result = {}
        idx = 0
        for outer_key in big_dictionary.keys():
            for inner_key in big_dictionary[outer_key].keys():
                inner_dict_result[inner_key] = (
                    big_dictionary[outer_key][inner_key]
                    + small_dictionary[0][inner_key]
                )
                result.setdefault(0, inner_dict_result)
    # Third Case
    # Both Submolecule_Dicts have the same length but bigger then one
    # For example if both have 3 then we have 3*3 choises
    if (
        len(length_of_dict) == 2
        and min(length_of_dict) == max(length_of_dict)
        and min(length_of_dict) != 1
    ):
        first_dictionary = dict_list[0]
        second_dictionary = dict_list[1]
        combinations = len(first_dictionary.items()) * len(second_dictionary.items())
        idx = 0
        new_dict = {}  # new dict with the structure gets created
        subindex = 1
        while idx < combinations:
            copy_of_the_dictionary = second_dictionary[subindex - 1].copy()
            if idx < subindex * len(second_dictionary.items()):
                new_dict.setdefault(idx, copy_of_the_dictionary)
                idx += 1
            else:
                subindex += 1
        # This while loop gives us the structure now we need to combine the elements
        # now we have combine it with the first dictionary element wise
        idx = 0
        subindex = 1
        while idx < len(new_dict.items()):
            if idx < subindex * len(second_dictionary.items()):
                new_dict[idx]["bonds"] = (
                    new_dict[idx]["bonds"] + first_dictionary[subindex - 1]["bonds"]
                )
                new_dict[idx]["angles"] = (
                    new_dict[idx]["angles"] + first_dictionary[subindex - 1]["angles"]
                )
                new_dict[idx]["linear valence angles"] = (
                    new_dict[idx]["linear valence angles"]
                    + first_dictionary[subindex - 1]["linear valence angles"]
                )
                new_dict[idx]["out of plane angles"] = (
                    new_dict[idx]["out of plane angles"]
                    + first_dictionary[subindex - 1]["out of plane angles"]
                )
                new_dict[idx]["dihedrals"] = (
                    new_dict[idx]["dihedrals"]
                    + first_dictionary[subindex - 1]["dihedrals"]
                )
                idx += 1
            else:
                subindex += 1
        result = new_dict  # chefschmee
    # Fourth Case:
    # We get three submolecules all with length one
    if (
        len(length_of_dict) == 3
        and min(length_of_dict) == 1
        and max(length_of_dict) == 1
    ):
        first_dictionary = dict_list[0]
        second_dictionary = dict_list[1]
        third_dictionary = dict_list[2]
        idx = 0
        inner_dict_result = {}
        for outer_key in first_dictionary.keys():
            for inner_key in first_dictionary[outer_key].keys():
                inner_dict_result[inner_key] = (
                    first_dictionary[outer_key][inner_key]
                    + second_dictionary[0][inner_key]
                    + third_dictionary[0][inner_key]
                )
                result.setdefault(0, inner_dict_result)
    # TODO there are cases missing for every bigger complexes
    if (
        len(length_of_dict) == 4
        and min(length_of_dict) == 1
        and max(length_of_dict) == 1
    ):
        first_dictionary = dict_list[0]
        second_dictionary = dict_list[1]
        third_dictionary = dict_list[2]
        fourth_dictionary = dict_list[3]
        idx = 0
        inner_dict_result = {}
        for outer_key in first_dictionary.keys():
            for inner_key in first_dictionary[outer_key].keys():
                inner_dict_result[inner_key] = (
                    first_dictionary[outer_key][inner_key]
                    + second_dictionary[0][inner_key]
                    + third_dictionary[0][inner_key]
                    + fourth_dictionary[0][inner_key]
                )
                result.setdefault(0, inner_dict_result)

    # First generalisation: if the lenght of one Dictionary is bigger than the other one we can first create the template then fill them up
    if min(length_of_dict) != max(length_of_dict) and len(length_of_dict) == 2:
        min_idx = pd.Series(length_of_dict).idxmin()
        max_idx = pd.Series(length_of_dict).idxmax()
        small_dict = dict_list[min_idx]
        big_dict = dict_list[max_idx]
        subindex = 0  # Both a index and a subindex are used for the combination
        idx = 0
        result = {}
        length_new_dictionary = min(length_of_dict) * max(length_of_dict)
        while idx < length_new_dictionary:
            copy_of_big_dict = big_dict[subindex].copy()
            result.setdefault(idx, copy_of_big_dict)
            idx += 1
            if subindex < len(big_dict) - 1:
                subindex += 1
            elif subindex >= len(big_dict) - 1:
                subindex = 0
        # We now have our dictionary structure next we create n
        idx = 0
        subindex = 0
        while idx < len(result):
            result[idx]["bonds"] = result[idx]["bonds"] + small_dict[subindex]["bonds"]
            result[idx]["angles"] = (
                result[idx]["angles"] + small_dict[subindex]["angles"]
            )
            result[idx]["linear valence angles"] = (
                result[idx]["linear valence angles"]
                + small_dict[subindex]["linear valence angles"]
            )
            result[idx]["out of plane angles"] = (
                result[idx]["out of plane angles"]
                + small_dict[subindex]["out of plane angles"]
            )
            result[idx]["dihedrals"] = (
                result[idx]["dihedrals"] + small_dict[subindex]["dihedrals"]
            )
            idx += 1
            if subindex < len(small_dict) - 1:
                subindex += 1
            elif subindex >= len(small_dict) - 1:
                subindex = 0
    return [result]


def add_element_to_all_entries(dict_list, element, dict_key, amount) -> list:
    """
    Adds an element to all entries in a given dictionary.

    Attributes:
        dict_list:
            a list of dictionaries where the dictionaries are IC sets
        element:
            a element which gets added to a specific position
        dict_key:
            a string which determines where the element gets added
        amount:
            determines how often this element gets added to a specific position
    Returns:
        a list of dictionaries
    """
    # Specify how much elements will be added
    number_elements = element[:amount]
    # Ite rate through each dictionary in the list
    for d in dict_list:
        for key, value in d.items():
            # Append the new bond to the bonds list
            if dict_key in value:
                value[dict_key].extend(number_elements)
            else:
                value[dict_key] = number_elements.copy()
    return dict_list


# Alternative way would distributing the degrees of freedom over all the possible lists
# The Distribute Elements function will be the main tool for building up the dictionary

# the divide_chunks function therefore slices the list


def divide_chunks(l, n):

    # looping till length 1
    for i in range(0, len(l), n):
        yield l[i : i + n]


def distribute_elements(dict_list, elements, dict_key, number_of_elements) -> list:
    """
    Distributes elements to a given IC set also generating all possible combinations for the even distribution

    Here all combinatoric cases are considered for example if the there is one IC set and one needs to distribute two IC coordinates out of a list of two
    the total amount of IC sets in the end is two

    Attributes:
        dict_list:
            a list of dictionaries with possible IC sets
        elements: 
            a list of possible elements that get distributed
        dict_key:
            a str that determines the key where the elements get added
        number_of_elements: 
            a str with the number of elements that are taken out of the elements list
    
    Returns:
        a list of dictionaries where each entry in the dictionary is a possible IC set
    """


    # The dict list will always be a list of dict with exactly one dict in it
    dictionary_inner = dict_list[0]
    # we then define a copy which we we need when building up a new dictionary entry
    copy_of_the_dictionary = dictionary_inner[0].copy()

    # Next a index for iteration
    idx = 0

    # and the new_dict where the structure is formed, and the end result as a variable
    new_dict = {}
    result = {}

    # This print helps us to define the case and have a check for testing purposed
    #  print("Dictionary - Key currently used:", dict_key)
    #  print("Length of elements list:", len(elements))
    #  print("Number of elements to take out", number_of_elements)
    #  print("Length of the actual dictionary", len(dictionary_inner.items()))

    # Non-defining case when the number of elements is zero
    if number_of_elements == 0:
        return dict_list

    # here we evaluate the length of the elements and let the index run until, the length of currenct dictionary * length of elements is reaches
    if number_of_elements <= 1 and len(elements) <= (len(dictionary_inner.items())):
        # we also define a subindex so that the dictionary gets duplicated the right amount of times
        subindex = 0
        length_new_dictionary = len(dictionary_inner.items()) * len(elements)
        while idx < length_new_dictionary:
            copy_of_the_dictionary = dictionary_inner[subindex].copy()
            for key, value in copy_of_the_dictionary.items():
                if key == dict_key:
                    new_dict.setdefault(idx, copy_of_the_dictionary)
                    idx += 1
                    if subindex < len(dictionary_inner.items()) - 1:
                        subindex += 1
                    elif subindex > len(dictionary_inner.items()) - 1:
                        subindex = 0
        duplicated_elements = list(
            itertools.islice(itertools.cycle(elements), length_new_dictionary)
        )

        for outer_key, inner_dict in new_dict.items():
            result[outer_key] = inner_dict.copy()
            copy_elements = inner_dict[dict_key].copy()
            result[outer_key][dict_key] = [duplicated_elements.pop(0)]
            if len(copy_elements) > 0:
                result[outer_key][dict_key].extend(copy_elements)

    # FIRST-CASE: number_of_elements = 1 and len(elements) > dictionary entries
    # Only difference here is that new entries get created
    # With this variant alone for example the water dimer can be calculated
    if number_of_elements <= 1 and len(elements) > len(dictionary_inner.items()):

        subindex = 0
        length_new_dictionary = len(dictionary_inner.items()) * len(elements)
        while idx < length_new_dictionary:
            copy_of_the_dictionary = dictionary_inner[subindex].copy()
            for key, value in copy_of_the_dictionary.items():
                if key == dict_key:
                    new_dict.setdefault(idx, copy_of_the_dictionary)
                    idx += 1
                    if subindex < len(dictionary_inner.items()) - 1:
                        subindex += 1
                    elif subindex > len(dictionary_inner.items()) - 1:
                        subindex = 0
        duplicated_elements = list(
            itertools.islice(itertools.cycle(elements), length_new_dictionary)
        )

        for outer_key, inner_dict in new_dict.items():
            result[outer_key] = inner_dict.copy()
            copy_elements = inner_dict[dict_key].copy()
            result[outer_key][dict_key] = [duplicated_elements.pop(0)]
            if len(copy_elements) > 0:
                result[outer_key][dict_key].extend(copy_elements)

    # in this case we slice out the elements and create packages which then get added
    if number_of_elements > 1:
        # Second-CASE .1 the number of elements is equal to the len of the list
        # In this case we dont_need to take into account any combinatoris, because there is just one combination to choose out of
        if number_of_elements == len(elements) and len(elements) <= (
            len(dictionary_inner.items())
        ):
            while idx < len(dictionary_inner.items()):
                copy_of_the_dictionary = dictionary_inner[idx].copy()
                for key, value in copy_of_the_dictionary.items():
                    if key == dict_key:
                        new_dict.setdefault(idx, copy_of_the_dictionary)
                        idx += 1

            for outer_key, inner_dict in new_dict.items():
                result[outer_key] = inner_dict.copy()
                for element in elements:  # maybe do this with zip, but is fine for now
                    if element not in result[outer_key][dict_key]:
                        result[outer_key][dict_key].extend(elements)
        # Second-Case .2 if the length of the elements is bigger then the dictionary we create new elements
        # TODO i think this case is pointess and can be combined with the upper one
        if number_of_elements == len(elements) and len(elements) >= (
            len(dictionary_inner.items())
        ):
            while idx < len(dictionary_inner.items()):
                copy_of_the_dictionary = dictionary_inner[idx].copy()
                for key, value in copy_of_the_dictionary.items():
                    if key == dict_key:
                        new_dict.setdefault(idx, copy_of_the_dictionary)
                        idx += 1

            for outer_key, inner_dict in new_dict.items():
                result[outer_key] = inner_dict.copy()
                for element in elements:  # maybe do this with zip, but is fine for now
                    if element not in result[outer_key][dict_key]:
                        result[outer_key][dict_key].extend(elements)

        if number_of_elements != len(elements):
            # This is our fist combinatoric case, we have now n choose k and need to evaluate the possibilities
            # for example 8 choose 2 gives us 28 possible combinations which we have to consider
            # We then have to evaluate the number of possible combinations via multiplication with the length of the list
            element_combinations = [
                list(perm)
                for perm in itertools.combinations(elements, number_of_elements)
            ]
            possible_combinations = len(element_combinations)
            length_of_the_dictionary = len(dictionary_inner.items())
            combinations = possible_combinations * length_of_the_dictionary
            # And in order to do this pop method we have to duplicate the list n times --> n is the length of the dictionar<
            duplicated_element_combinations = [
                tup
                for _ in range(length_of_the_dictionary)
                for tup in element_combinations
            ]
            subindex = 1  # we let the subindex rise according to the number of dictionarys thus  far
            while idx < combinations:
                # Now we have to take pairs according to the length of the dictionary = 84 / 3 = 28 per dictionary
                copy_of_the_dictionary = dictionary_inner[subindex - 1].copy()
                for key, value in copy_of_the_dictionary.items():
                    if key == dict_key:
                        new_dict.setdefault(idx, copy_of_the_dictionary)
                        if idx < possible_combinations * subindex:
                            idx += 1
                        else:
                            subindex += 1
            # The dictionary therefore needs the length of the combinations
            for outer_key, inner_dict in new_dict.items():
                result[outer_key] = inner_dict.copy()
                # Fix elements to append
                copy_elements = inner_dict[dict_key].copy()
                if length_of_the_dictionary == 1:
                    result[outer_key][dict_key] = element_combinations.pop(0)
                elif length_of_the_dictionary > 1:
                    result[outer_key][dict_key] = duplicated_element_combinations.pop(0)
                # this is a safety function because for some strange reason the elements get copied 3 times ??
                if len(copy_elements) > 0:
                    unique_items = [
                        item
                        for item in copy_elements
                        if item not in result[outer_key][dict_key]
                    ]
                    result[outer_key][dict_key] = list(
                        itertools.chain(result[outer_key][dict_key], unique_items)
                    )  # at this point just ask Kemal
    return [result]


""" 
General Intermolecular Systems
"""
# In the general systems the oop angles are normally not considered


def intermolecular_general_acyclic_linunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    l,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the intermolecular, general acyclic case without linear submolecules.

    Attributes:
        ic_dict:
            a dictionary containing valid IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing angles
        out_of_plane:
            a list of tupels containing oop's
        dihedrals:
            a list of tupels containing dihedrals
        a_1:
            a integer with the number of terminal atoms
        l:
            a integer with the number of linear bonds
        num_bonds:
            a integer with the number of bonds 
        num_atoms:
            a integer with the number of atoms
        specification:
            the specification dictionary of the given molecule 

    Returns:
        a dictionary where each entry is a possible IC set
    """

    args = arguments.get_args()

    n_r = num_bonds
    n_phi = 4 * num_bonds - 3 * num_atoms + a_1 - (l - 1)
    n_gamma = 0
    planar_subunits_list = specification["planar submolecule(s)"]
    n_phi_prime = 2 * (l - 1)
    n_tau = num_bonds - a_1 - (l - 1)

    # If planar submolecules exist we have to change our definition of n_phi and n_gamma
    # remove angles at the specified coordinate, as we else would have linear dependencies
    if len(planar_subunits_list) != 0:
        n_phi, n_gamma, angles = get_param_planar_submolecule(
            planar_subunits_list, specification["multiplicity"], angles
        )

    # Open up the dictionary
    ic_dict_list = []
    ic_dic_list_1d_complex = []

    # first we determine the two disconnected submolecules
    connected_components, submolecules_bond_list, _ = atoms_list.detect_submolecules()
    submolecules_atoms_coordinates = extract_atoms_of_submolecules(
        connected_components, atoms_list
    )
    summarized_submolecule_ic_dict = []
    # every submolecule gets its own specification

    n = 1
    for submolecule_atom, submolecule_bonds in zip(
        submolecules_atoms_coordinates, submolecules_bond_list
    ):
        # Step 1: Point Group

        molecule = mg.Molecule(
            [strip_numbers(atom.symbol) for atom in submolecule_atom],
            [atom.coordinates for atom in submolecule_atom],
        )

        # in pymatgen we need to manualy adjust the mass when using deuterium
        indices_mg = []
        # we retrieve the indices where these atoms occure
        for i, atom in enumerate(molecule):
            atom_swapped = re.sub(r"[^a-zA-Z]", "", str(atom.species))
            if atom_swapped == "D":
                indices_mg.append(i)
        for indices in indices_mg:
            molecule[indices].species = {"H": 1.0}  # keep species as hydrogen
            molecule[indices].properties = {"atomic_mass": 2.014}
        # calculate spin reduction
        spin_reduction = len(
            indices_mg
        )  # determines how much we need to reduce the spin with this swapping trick
        molecule.set_charge_and_spin(0, molecule.spin_multiplicity - spin_reduction)

        molecule_pg = PointGroupAnalyzer(molecule)
        point_group_sch = molecule_pg.sch_symbol

        # Transform into a Molecule Class
        submolecule_atom = Molecule(submolecule_atom)

        # Step 1.1: Calculate the Connectivity C and pass it to specification (just a quick fix)

        molecular_graph = submolecule_atom.graph_rep()
        specifications.connectivity_c = submolecule_atom.count_connected_components(
            molecular_graph
        )

        # Step 2: Initialize Internal Coordinates:
        sub_angles, sub_linear_angles = submolecule_atom.generate_angles(
            submolecule_bonds
        )
        sub_dihedrals = submolecule_atom.generate_dihedrals(submolecule_bonds)

        # Step 3: Give the Specification for each molecule
        specification_submolecule = dict()
        specification_submolecule = specifications.calculation_specification(
            specification,
            submolecule_atom,
            molecule_pg,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
        )

        # Step 4: According to the specification create OOP ICs
        if specification["planar"] == "yes":
            sub_oop = submolecule_atom.generate_out_of_plane(submolecule_bonds)
        elif (
            specification["planar"] == "no"
            and specification["planar submolecule(s)"] == []
        ):
            sub_oop = []
        elif specification["planar"] == "no" and not (
            specification["planar submolecule(s)"] == []
        ):
            sub_oop = atoms.generate_oop_planar_subunits(
                specification["planar submolecule(s)"]
            )

        # Step 5: Calculate IDOF for the given Submolecules:
        idof_submolecule = 0
        if specification_submolecule["linearity"] == "fully linear":
            idof_submolecule = 3 * len(submolecule_atom) - 5
        else:
            idof_submolecule = 3 * len(submolecule_atom) - 6

        # Command Line Message about IC Set generation
        print(f"IC Set generation for submolecule {n}")

        # Additional Logging Information about Submolecule

        n += 1

        # Step 5: With the Given Specification we can calculate the submolecule ic_dict:
        ic_dict_submolecule = icsel.get_sets(
            idof_submolecule,
            out,
            submolecule_atom,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
            sub_oop,
            sub_dihedrals,
            specification_submolecule,
        )
        summarized_submolecule_ic_dict.append(ic_dict_submolecule)

    # logging for generation of submolecule IC sets

    logfile.write_logfile_submolecule_treatment(
        out, summarized_submolecule_ic_dict, connected_components
    )

    # Combine the Dictionary with all possible combinations:

    fixed_ic_dict = combine_dictionaries(summarized_submolecule_ic_dict)

    # Now we evaluate the number of internal coordinates we have so far
    bonds_length = []
    angles_length = []
    linear_angles_length = []
    dihedrals_length = []
    for d in fixed_ic_dict:
        for key, value in d.items():
            bonds_length.append(len(value.get("bonds", [])))
            angles_length.append(len(value.get("angles", [])))
            linear_angles_length.append(len(value.get("linear valence angles", [])))
            dihedrals_length.append(len(value.get("dihedrals", [])))

    # Now we can specify how much of the different intermolecular coordinates we need in our sets:

    ic_bonds_needed = n_r - max(bonds_length)
    ic_angles_needed = n_phi - max(angles_length)
    ic_linear_angles_needed = n_phi_prime - max(linear_angles_length)
    ic_dihedrals_needed = n_tau - max(dihedrals_length)

    # args comb controles which ics are used

    intermolecular_bonds = Total_IC_dict["h_bond"]
    intermolecular_angles = Total_IC_dict["h_bond_angles"]
    intermolecular_linear_angles = Total_IC_dict["h_bond_linear_angles"]
    intermolecular_dihedrals = Total_IC_dict["h_bond_dihedrals"]

    if args.comb == 2:
        intermolecular_bonds = intermolecular_bonds + Total_IC_dict["acc_don"]
        intermolecular_angles = intermolecular_angles + Total_IC_dict["acc_don_angles"]
        intermolecular_linear_angles = (
            intermolecular_linear_angles + Total_IC_dict["acc_don_linear_angles"]
        )

    if args.comb == 3:
        intermolecular_bonds = intermolecular_bonds + Total_IC_dict["acc_don"]
        intermolecular_angles = intermolecular_angles + Total_IC_dict["acc_don_angles"]
        intermolecular_linear_angles = (
            intermolecular_linear_angles + Total_IC_dict["acc_don_linear_angles"]
        )
        intermolecular_dihedrals = (
            intermolecular_dihedrals + Total_IC_dict["acc_don_dihedrals"]
        )

    ic_dict_list = distribute_elements(
        fixed_ic_dict, intermolecular_bonds, "bonds", ic_bonds_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list, intermolecular_angles, "angles", ic_angles_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        intermolecular_linear_angles,
        "linear valence angles",
        ic_linear_angles_needed,
    )
    ic_dict_list = distribute_elements(
        ic_dict_list, intermolecular_dihedrals, "dihedrals", ic_dihedrals_needed
    )

    ic_dict = dict()
    new_key = 0
    for dictionary in ic_dict_list:
        for key, value in dictionary.copy().items():
            ic_dict[new_key] = value
            new_key += 1
    return ic_dict


def intermolecular_general_acyclic_nolinunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the intermolecular, general, acyclic case without linear submolecules
    Attributes:
        ic_dict:
            a dictionary containing valid IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing angles
        out_of_plane:
            a list of tupels containing oop's
        dihedrals:
            a list of tupels containing dihedrals
        a_1:
            a integer with the number of terminal atoms
        num_bonds:
            a integer with the number of bonds 
        num_atoms:
            a integer with the number of atoms
        specification:
            the specification dictionary of the given molecule 
    
    Returns:
        a dictionary where each entry is a possible IC set
    """
    args = arguments.get_args()

    n_r = num_bonds
    n_phi = 4 * num_bonds - 3 * num_atoms + a_1
    n_gamma = 0
    planar_subunits_list = specification["planar submolecule(s)"]
    n_tau = num_bonds - a_1

    # If planar submolecules exist we have to change our definition of n_phi and n_gamma
    # remove angles at the specified coordinate, as we else would have linear dependencies
    if len(planar_subunits_list) != 0:
        n_phi, n_gamma, angles = get_param_planar_submolecule(
            planar_subunits_list, specification["multiplicity"], angles
        )

    # Open up the dictionary
    ic_dict_list = []
    ic_dic_list_1d_complex = []
    # first we determine the two disconnected submolecules
    connected_components, submolecules_bond_list, _ = atoms_list.detect_submolecules()

    submolecules_atoms_coordinates = extract_atoms_of_submolecules(
        connected_components, atoms_list
    )

    summarized_submolecule_ic_dict = []

    n = 1
    for submolecule_atom, submolecule_bonds in zip(
        submolecules_atoms_coordinates, submolecules_bond_list
    ):
        # Step 1: Point Group

        molecule = mg.Molecule(
            [strip_numbers(atom.symbol) for atom in submolecule_atom],
            [atom.coordinates for atom in submolecule_atom],
        )

        # in pymatgen we need to manualy adjust the mass when using deuterium
        indices_mg = []
        # we retrieve the indices where these atoms occure
        for i, atom in enumerate(molecule):
            atom_swapped = re.sub(r"[^a-zA-Z]", "", str(atom.species))
            if atom_swapped == "D":
                indices_mg.append(i)
        for indices in indices_mg:
            molecule[indices].species = {"H": 1.0}  # keep species as hydrogen
            molecule[indices].properties = {"atomic_mass": 2.014}
        # calculate spin reduction
        spin_reduction = len(
            indices_mg
        )  # determines how much we need to reduce the spin with this swapping trick
        molecule.set_charge_and_spin(0, molecule.spin_multiplicity - spin_reduction)

        molecule_pg = PointGroupAnalyzer(molecule)
        point_group_sch = molecule_pg.sch_symbol

        # Conversion into Molecule Class

        submolecule_atom = Molecule(submolecule_atom)

        # Step 1.1: Calculate the Connectivity C and pass it to specification

        molecular_graph = submolecule_atom.graph_rep()
        specifications.connectivity_c = submolecule_atom.count_connected_components(
            molecular_graph
        )

        # Step 2: Initialize Internal Coordinates:
        sub_angles, sub_linear_angles = submolecule_atom.generate_angles(
            submolecule_bonds
        )
        sub_dihedrals = submolecule_atom.generate_dihedrals(submolecule_bonds)

        # Step 3: Give the Specification for each molecule
        specification_submolecule = dict()
        specification_submolecule = specifications.calculation_specification(
            specification,
            submolecule_atom,
            molecule_pg,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
        )

        # Step 4: According to the specification create OOP ICs
        if specification["planar"] == "yes":
            sub_oop = submolecule_atom.generate_out_of_plane(submolecule_bonds)
        elif (
            specification["planar"] == "no"
            and specification["planar submolecule(s)"] == []
        ):
            sub_oop = []
        elif specification["planar"] == "no" and not (
            specification["planar submolecule(s)"] == []
        ):
            sub_oop = atoms.generate_oop_planar_subunits(
                specification["planar submolecule(s)"]
            )

        # Step 5: Calculate IDOF for the given Submolecules:
        idof_submolecule = 0
        if specification_submolecule["linearity"] == "fully linear":
            idof_submolecule = 3 * len(submolecule_atom) - 5
        else:
            idof_submolecule = 3 * len(submolecule_atom) - 6

        # Command Line Message about IC Set generation
        print(f"IC Set generation for submolecule {n}")
        n += 1

        # Step 5: With the Given Specification we can calculate the submolecule ic_dict:
        ic_dict_submolecule = icsel.get_sets(
            idof_submolecule,
            out,
            submolecule_atom,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
            sub_oop,
            sub_dihedrals,
            specification_submolecule,
        )
        summarized_submolecule_ic_dict.append(ic_dict_submolecule)

    # logging for generation of submolecule IC sets

    logfile.write_logfile_submolecule_treatment(
        out, summarized_submolecule_ic_dict, connected_components
    )

    # What we now do with this list of dictionary is to fix the internal coordinates of the two submolecules, now we need to evaluate the missing coordinates for a complete set
    fixed_ic_dict = combine_dictionaries(summarized_submolecule_ic_dict)

    # Now we evaluate the number of internal coordinates we have so far
    bonds_length = []
    angles_length = []
    linear_angles_length = []
    dihedrals_length = []
    for d in fixed_ic_dict:
        for key, value in d.items():
            bonds_length.append(len(value.get("bonds", [])))
            angles_length.append(len(value.get("angles", [])))
            linear_angles_length.append(len(value.get("linear valence angles", [])))
            dihedrals_length.append(len(value.get("dihedrals", [])))

    # Now we can specify how much of the different intermolecular coordinates we need in our sets:

    ic_bonds_needed = n_r - max(bonds_length)
    ic_angles_needed = n_phi - max(angles_length)
    ic_dihedrals_needed = n_tau - max(dihedrals_length)

    intermolecular_bonds = Total_IC_dict["h_bond"]
    intermolecular_angles = Total_IC_dict["h_bond_angles"]
    intermolecular_linear_angles = Total_IC_dict["h_bond_linear_angles"]
    intermolecular_dihedrals = Total_IC_dict["h_bond_dihedrals"]

    if args.comb == 2:
        intermolecular_bonds = intermolecular_bonds + Total_IC_dict["acc_don"]
        intermolecular_angles = intermolecular_angles + Total_IC_dict["acc_don_angles"]
        intermolecular_linear_angles = (
            intermolecular_linear_angles + Total_IC_dict["acc_don_linear_angles"]
        )

    if args.comb == 3:
        intermolecular_bonds = intermolecular_bonds + Total_IC_dict["acc_don"]
        intermolecular_angles = intermolecular_angles + Total_IC_dict["acc_don_angles"]
        intermolecular_linear_angles = (
            intermolecular_linear_angles + Total_IC_dict["acc_don_linear_angles"]
        )
        intermolecular_dihedrals = (
            intermolecular_dihedrals + Total_IC_dict["acc_don_dihedrals"]
        )

    ic_dict_list = distribute_elements(
        fixed_ic_dict, intermolecular_bonds, "bonds", ic_bonds_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list, intermolecular_angles, "angles", ic_angles_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list, intermolecular_dihedrals, "dihedrals", ic_dihedrals_needed
    )
    ic_dict = dict()

    new_key = 0
    for dictionary in ic_dict_list:
        for key, value in dictionary.copy().items():
            ic_dict[new_key] = value
            new_key += 1
    return ic_dict


def intermolecular_general_cyclic_nolinsub(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    num_of_red,
    a_1,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the intermolecular, general cyclic case where no linear submolecules are contained
    Attributes:
        ic_dict:
            a dictionary containing valid IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing angles
        out_of_plane:
            a list of tupels containing oop's
        dihedrals:
            a list of tupels containing dihedrals
        a_1:
            a integer with the number of terminal atoms
        num_of_red:
            a integer with the number of redundancies of the molecular structure
        num_bonds:
            a integer with the number of bonds 
        num_atoms:
            a integer with the number of atoms
        specification:
            the specification dictionary of the given molecule

    Returns:
        a dictionary where each entry is a valid IC set
    """

    args = arguments.get_args()

    # Open up the dictionary
    ic_dict_list = []
    ic_dic_list_1d_complex = []
    # first we determine the two disconnected submolecules
    connected_components, submolecules_bond_list, _ = atoms_list.detect_submolecules()
    submolecules_atoms_coordinates = extract_atoms_of_submolecules(
        connected_components, atoms_list
    )
    summarized_submolecule_ic_dict = []

    n = 1
    for submolecule_atom, submolecule_bonds in zip(
        submolecules_atoms_coordinates, submolecules_bond_list
    ):
        # Step 1: Point Group

        molecule = mg.Molecule(
            [strip_numbers(atom.symbol) for atom in submolecule_atom],
            [atom.coordinates for atom in submolecule_atom],
        )
        molecule_pg = PointGroupAnalyzer(molecule)
        point_group_sch = molecule_pg.sch_symbol

        # Convert into Molecule Class
        submolecule_atom = Molecule(submolecule_atom)

        # Step 1.1: Calculate the Connectivity C and pass it to specification (just a quick fix)

        molecular_graph = submolecule_atom.graph_rep()
        specifications.connectivity_c = submolecule_atom.count_connected_components(
            molecular_graph
        )

        # Step 2: Initialize Internal Coordinates:
        sub_angles, sub_linear_angles = submolecule_atom.generate_angles(
            submolecule_bonds
        )
        sub_dihedrals = submolecule_atom.generate_dihedrals(submolecule_bonds)

        # Step 3: Give the Specification for each molecule
        specification_submolecule = dict()
        specification_submolecule = specifications.calculation_specification(
            specification_submolecule,
            submolecule_atom,
            molecule_pg,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
        )

        # Step 4: With the specification generate OOP ICs

        if specification["planar"] == "yes":
            sub_oop = submolecule_atoms.generate_out_of_plane(submolecule_bonds)
        elif (
            specification["planar"] == "no"
            and specification["planar submolecule(s)"] == []
        ):
            sub_oop = []
        elif specification["planar"] == "no" and not (
            specification["planar submolecule(s)"] == []
        ):
            sub_oop = atoms.generate_oop_planar_subunits(
                specification["planar submolecule(s)"]
            )

        # Step 5: Calculate IDOF for the given Submolecules:
        idof_submolecule = 0
        if specification_submolecule["linearity"] == "fully linear":
            idof_submolecule = 3 * len(submolecule_atom) - 5
        else:
            idof_submolecule = 3 * len(submolecule_atom) - 6

        # Command Line Message about IC Set generation
        print(f"IC Set generation for submolecule {n}")
        n += 1

        # Step 5: With the Given Specification we can calculate the submolecule ic_dict:
        ic_dict_submolecule = icsel.get_sets(
            idof_submolecule,
            out,
            submolecule_atom,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
            sub_oop,
            sub_dihedrals,
            specification_submolecule,
        )
        summarized_submolecule_ic_dict.append(ic_dict_submolecule)

    # logging for generation of submolecule IC sets

    logfile.write_logfile_submolecule_treatment(
        out, summarized_submolecule_ic_dict, connected_components
    )

    # What we now do with this list of dictionary is to fix the internal coordinates of the two submolecules, now we need to evaluate the missing coordinates for a complete set

    fixed_ic_dict = combine_dictionaries(summarized_submolecule_ic_dict)

    # Redundancies mu get taken out of the intermolecular coordinates taking symmetry into account
    n_r = num_bonds - specification["mu"]
    n_phi = 4 * (num_bonds - specification["mu"]) - 3 * num_atoms + a_1
    n_tau = (num_bonds - specification["mu"]) - a_1

    # First we delete the h-bond coordinate then the corresponding acc_don_coordinate

    symmetric_intermolecular_bonds = icsel.get_symm_bonds(
        Total_IC_dict["h_bond"], specification
    )
    symmetric_intermolecular_bonds_list = icsel.get_bond_subsets(
        symmetric_intermolecular_bonds
    )
    valide_atoms = valide_atoms_to_cut(
        Total_IC_dict["h_bond"], specification["multiplicity"]
    )
    removed_bonds = []

    # initilize a bond dict for removal of acc_don bond
    bond_dict = atoms_list.bond_dict(
        Total_IC_dict["cov_bond"] + Total_IC_dict["h_bond"]
    )

    for symmetric_bond_group in symmetric_intermolecular_bonds_list:
        if len(symmetric_bond_group) >= specification[
            "mu"
        ] and bonds_are_in_valide_atoms(symmetric_bond_group, valide_atoms):
            removed_bonds, random_var = delete_bonds_symmetry(
                symmetric_bond_group,
                Total_IC_dict["h_bond"],
                specification["mu"],
                valide_atoms,
            )
        if not removed_bonds:
            continue

    # If the symmetric cut was possible we update the corresponding internal coordinates

    intermolecular_h_bonds_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond"]
    )
    intermolecular_h_angles_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond_angles"]
    )
    intermolecular_h_dihedrals_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond_dihedrals"]
    )

    # when no symmetric cutting is possible we cut a random bond
    # this is not ideal and leads to more possible solution of the decomposition process

    if len(removed_bonds) == 0:
        removed_bonds = Total_IC_dict["h_bond"][: specification["mu"]]

        intermolecular_h_bonds_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond"]
        )
        intermolecular_h_angles_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_angles"]
        )
        intermolecular_h_dihedrals_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_dihedrals"]
        )

    removed_acc_don_bonds = []
    # basically we use the hydrogen as a key to get our element out
    for removed_bond in removed_bonds:
        if removed_bond[0].strip(string.digits) == "H":
            removed_acc_don_bonds.append(bond_dict[removed_bond[0]])
        else:
            removed_acc_don_bonds.append(bond_dict[removed_bond[1]])

    intermolecular_acc_don_bonds_updated = update_internal_coordinates_cyclic(
        removed_acc_don_bonds, Total_IC_dict["acc_don"]
    )
    intermolecular_acc_don_angle_updated = update_internal_coordinates_cyclic(
        removed_acc_don_bonds, Total_IC_dict["acc_don_angles"]
    )
    intermolecular_acc_don_dihedrals_updated = update_internal_coordinates_cyclic(
        removed_acc_don_bonds, Total_IC_dict["acc_don_dihedrals"]
    )

    # evaluate the amount of internal coordinates fixed by the submolecules

    bonds_length = []
    angles_length = []
    linear_angles_length = []
    dihedrals_length = []
    oop_length = []
    for d in fixed_ic_dict:
        for key, value in d.items():
            bonds_length.append(len(value.get("bonds", [])))
            angles_length.append(len(value.get("angles", [])))
            linear_angles_length.append(len(value.get("linear valence angles", [])))
            dihedrals_length.append(len(value.get("dihedrals", [])))
            oop_length.append(len(value.get("out of plane angles", [])))

    intermolecular_bonds_needed = n_r - max(bonds_length)
    intermolecular_angles_needed = n_phi - max(angles_length)
    intermolecular_dihedrals_needed = n_tau - max(dihedrals_length)

    total_intermolecular_bonds = intermolecular_h_bonds_updated
    total_intermolecular_angles = intermolecular_h_angles_updated
    total_intermolecular_dihedrals = intermolecular_h_dihedrals_updated

    if args.comb == 2:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + intermolecular_acc_don_bonds_updated
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + intermolecular_acc_don_angle_updated
        )

    if args.comb == 3:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + intermolecular_acc_don_bonds_updated
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + intermolecular_acc_don_angle_updated
        )
        total_intermolecular_dihedrals = (
            total_intermolecular_dihedrals + intermolecular_acc_don_dihedrals_updated
        )

    ic_dict_list = distribute_elements(
        fixed_ic_dict, total_intermolecular_bonds, "bonds", intermolecular_bonds_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_angles,
        "angles",
        intermolecular_angles_needed,
    )

    # safety feature if not enough dihedrals are present include all acc_don_dihedrals
    if len(total_intermolecular_dihedrals) < intermolecular_dihedrals_needed:
        total_intermolecular_dihedrals = (
            intermolecular_h_dihedrals_updated + Total_IC_dict["acc_don_dihedrals"]
        )
        ic_dict_list = distribute_elements(
            ic_dict_list,
            total_intermolecular_dihedrals,
            "dihedrals",
            intermolecular_dihedrals_needed,
        )
    else:
        ic_dict_list = distribute_elements(
            ic_dict_list,
            total_intermolecular_dihedrals,
            "dihedrals",
            intermolecular_dihedrals_needed,
        )

    ic_dict = dict()

    new_key = 0
    for dictionary in ic_dict_list:
        for key, value in dictionary.copy().items():
            ic_dict[new_key] = value
            new_key += 1
    return ic_dict


def intermolecular_general_cyclic_linunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    num_of_red,
    a_1,
    l,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the intermolecular, general, cyclic case with linear submolecules

    Attributes:
        ic_dict:
            a dictionary containing valid IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing angles
        out_of_plane:
            a list of tupels containing oop's
        dihedrals:
            a list of tupels containing dihedrals
        a_1:
            a integer with the number of terminal atoms
        l:
            a integer with the number of linear bonds
        num_bonds:
            a integer with the number of bonds 
        num_of_red:
            a integer with the number of redundancies
        num_atoms:
            a integer with the number of atoms
        specification:
            the specification dictionary of the given molecule 

    Returns:
        a dictionary where each entry of the dictionary is a valid IC set
    """

    args = arguments.get_args()

    # Open up the dictionary
    ic_dict_list = []
    ic_dic_list_1d_complex = []

    # first we determine the two disconnected submolecules

    connected_components, submolecules_bond_list, _ = atoms_list.detect_submolecules()
    submolecules_atoms_coordinates = extract_atoms_of_submolecules(
        connected_components, atoms_list
    )

    summarized_submolecule_ic_dict = []
    # now the hard part every submolecule gets its own specification:
    n = 1
    for submolecule_atom, submolecule_bonds in zip(
        submolecules_atoms_coordinates, submolecules_bond_list
    ):

        # Step 1: Point Group

        molecule = mg.Molecule(
            [strip_numbers(atom.symbol) for atom in submolecule_atom],
            [atom.coordinates for atom in submolecule_atom],
        )
        molecule_pg = PointGroupAnalyzer(molecule)
        point_group_sch = molecule_pg.sch_symbol

        submolecule_atom = Molecule(submolecule_atom)
        # Step 1.1: Calculate the Connectivity C and pass it to specification (just a quick fix)

        molecular_graph = submolecule_atom.graph_rep()
        specifications.connectivity_c = submolecule_atom.count_connected_components(
            molecular_graph
        )

        # Step 2: Initialize Internal Coordinates:
        sub_angles, sub_linear_angles = submolecule_atom.generate_angles(
            submolecule_bonds
        )
        sub_dihedrals = submolecule_atom.generate_dihedrals(submolecule_bonds)

        # Step 3: Give the Specification for each molecule
        specification_submolecule = dict()
        specification_submolecule = specifications.calculation_specification(
            specification_submolecule,
            submolecule_atom,
            molecule_pg,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
        )
        # Step 4: With the specification generate OOP ICs
        if specification["planar"] == "yes":
            sub_oop = submolecule_atom.generate_out_of_plane(submolecule_bonds)
        elif (
            specification["planar"] == "no"
            and specification["planar submolecule(s)"] == []
        ):
            sub_oop = []
        elif specification["planar"] == "no" and not (
            specification["planar submolecule(s)"] == []
        ):
            sub_oop = submolecule_atom.generate_oop_planar_subunits(
                submolecule_bonds, specification["planar submolecule(s)"]
            )

        # Step 5: Calculate IDOF for the given Submolecules:
        idof_submolecule = 0
        if specification_submolecule["linearity"] == "fully linear":
            idof_submolecule = 3 * len(submolecule_atom) - 5
        else:
            idof_submolecule = 3 * len(submolecule_atom) - 6

        # Command Line Message about IC Set generation
        print(f"IC Set generation for submolecule {n}")
        n += 1

        # Step 5: With the Given Specification we can calculate the submolecule ic_dict:
        ic_dict_submolecule = icsel.get_sets(
            idof_submolecule,
            out,
            submolecule_atom,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
            sub_oop,
            sub_dihedrals,
            specification_submolecule,
        )

        summarized_submolecule_ic_dict.append(ic_dict_submolecule)

    # logging for generation of submolecule IC sets

    logfile.write_logfile_submolecule_treatment(
        out, summarized_submolecule_ic_dict, connected_components
    )

    # What we now do with this list of dictionary is to fix the internal coordinates of the two submolecules, now we need to evaluate the missing coordinates for a complete set
    fixed_ic_dict = combine_dictionaries(summarized_submolecule_ic_dict)

    symmetric_intermolecular_bonds = icsel.get_symm_bonds(
        Total_IC_dict["h_bond"], specification
    )
    symmetric_intermolecular_bonds_list = icsel.get_bond_subsets(
        symmetric_intermolecular_bonds
    )
    valide_atoms = valide_atoms_to_cut(
        Total_IC_dict["h_bond"], specification["multiplicity"]
    )
    removed_bonds = []

    # initialize a bond dictionary for removal of acc don bond
    bond_dict = atoms_list.bond_dict(
        Total_IC_dict["cov_bond"] + Total_IC_dict["h_bond"]
    )

    # Beta is used to determine if there are intramolecular redundancies present
    if specification["beta"] == 0:
        for symmetric_bond_group in symmetric_intermolecular_bonds_list:
            if len(symmetric_bond_group) >= specification[
                "mu"
            ] and bonds_are_in_valide_atoms(symmetric_bond_group, valide_atoms):
                removed_bonds, random_var = delete_bonds_symmetry(
                    symmetric_bond_group,
                    Total_IC_dict["h_bond"],
                    specification["mu"],
                    valide_atoms,
                )
            if not removed_bonds:
                continue

        intermolecular_h_bonds_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond"]
        )
        intermolecular_h_angles_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_angles"]
        )
        intermolecular_h_linear_angles_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_linear_angles"]
        )
        intermolecular_h_dihedrals_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_dihedrals"]
        )

        # If the symmetric cutting didnt work one can simply take the first random cut
        if len(removed_bonds) == 0:
            removed_bonds = Total_IC_dict["h_bond"][: specficiation["mu"]]

            intermolecular_h_bonds_updated = update_internal_coordinates_cyclic(
                removed_bonds, Total_IC_dict["h_bond"]
            )
            intermolecular_h_angles_updated = update_internal_coordinates_cyclic(
                removed_bonds, Total_IC_dict["h_bond_angles"]
            )
            intermolecular_h_linear_angles_updated = update_internal_coordinates_cyclic(
                removed_bonds, Total_IC_dict["h_bond_linear_angles"]
            )
            intermolecular_h_dihedrals_updated = update_internal_coordinates_cyclic(
                removed_bonds, Total_IC_dict["h_bond_dihedrals"]
            )

        removed_acc_don_bonds = []
        # basically we use the hydrogen as a key to get our element out
        for removed_bond in removed_bonds:
            if removed_bond[0].strip(string.digits) == "H":
                removed_acc_don_bonds.append(bond_dict[removed_bond[0]])
            else:
                removed_acc_don_bonds.append(bond_dict[removed_bond[1]])

        intermolecular_acc_don_bonds_updated = update_internal_coordinates_cyclic(
            removed_acc_don_bonds, Total_IC_dict["acc_don"]
        )
        intermolecular_acc_don_angle_updated = update_internal_coordinates_cyclic(
            removed_acc_don_bonds, Total_IC_dict["acc_don_angles"]
        )
        intermolecular_acc_don_linear_angle_updated = (
            update_internal_coordinates_cyclic(
                removed_acc_don_bonds, Total_IC_dict["acc_don_linear_angles"]
            )
        )
        intermolecular_acc_don_dihedrals_updated = update_internal_coordinates_cyclic(
            removed_acc_don_bonds, Total_IC_dict["acc_don_dihedrals"]
        )

        # with the redundancies taken out we generate the total ICs for the appending algorithm

        total_intermolecular_bonds = intermolecular_h_bonds_updated
        total_intermolecular_angles = intermolecular_h_angles_updated
        total_intermolecular_linear_angles = intermolecular_h_linear_angles_updated
        total_intermolecular_dihedrals = intermolecular_h_dihedrals_updated

        if args.comb == 2:
            total_intermolecular_bonds = (
                total_intermolecular_bonds + intermolecular_acc_don_bonds_updated
            )
            total_intermolecular_angles = (
                total_intermolecular_angles + intermolecular_acc_don_angle_updated
            )

        if args.comb == 3:
            total_intermolecular_bonds = (
                total_intermolecular_bonds + intermolecular_acc_don_bonds_updated
            )
            total_intermolecular_angles = (
                total_intermolecular_angles + intermolecular_acc_don_angle_updated
            )
            total_intermolecular_dihedrals = (
                total_intermolecular_dihedrals
                + intermolecular_acc_don_dihedrals_updated
            )

    # when beta is not equal to zero we know the redundancies have been taken out of the intramolecular coordinates
    else:

        total_intermolecular_bonds = Total_IC_dict["h_bond"]
        total_intermolecular_angles = Total_IC_dict["h_bond_angles"]
        total_intermolecular_linear_angles = Total_IC_dict["h_bond_linear_angles"]
        total_intermolecular_dihedrals = Total_IC_dict["h_bond_dihedrals"]

        if args.comb == 2:
            total_intermolecular_bonds = (
                total_intermolecular_bonds + Total_IC_dict["acc_don"]
            )

            total_intermolecular_angles = (
                total_intermolecular_angles + Total_IC_dict["acc_don_angles"]
            )
            total_intermolecular_linear_angles = (
                total_intermolecular_angles + Total_IC_dict["acc_don_linear_angles"]
            )

        if args.comb == 3:
            total_intermolecular_bonds = (
                total_intermolecular_bonds + Total_IC_dict["acc_don"]
            )
            total_intermolecular_angles = (
                total_intermolecular_angles + Total_IC_dict["acc_don_angles"]
            )
            total_intermolecular_linear_angles = (
                total_intermolecular_angles + Total_IC_dict["acc_don_linear_angles"]
            )
            total_intermolecular_dihedrals = (
                total_intermolecular_dihedrals + Total_IC_dict["acc_don_dihedrals"]
            )

    # evaluate the amount of internal coordinates fixed by the submolecules

    bonds_length = []
    angles_length = []
    linear_angles_length = []
    dihedrals_length = []
    oop_length = []
    for d in fixed_ic_dict:
        for key, value in d.items():
            bonds_length.append(len(value.get("bonds", [])))
            angles_length.append(len(value.get("angles", [])))
            linear_angles_length.append(len(value.get("linear valence angles", [])))
            dihedrals_length.append(len(value.get("dihedrals", [])))
            oop_length.append(len(value.get("out of plane angles", [])))

    if len(removed_bonds) != 0:
        linear_bonds = specifications.get_linear_bonds(linear_angles)
        swapped_removed_bonds = removed_bonds[0][::-1]
        # for n_r we just substract the mu

        if removed_bonds[0] in linear_bonds or swapped_removed_bonds in linear_bonds:
            print("The bond removed", removed_bonds, "was linear")
            # Decius for this topology : n_phi = 4 * b - 3a - a_1 - (l-1)
            # because the bond removed was linear we also reduce the number l by two because the whole linear system is destroyed
            n_phi = (
                4 * (num_bonds - specification["mu"])
                - 3 * num_atoms
                + a_1
                - ((l - 2 * specification["mu"]) - 1)
            )
            # for n_phi_prime l gets also reduced by two
            n_phi_prime = 2 * ((l - 2 * specification["mu"]) - 1)
            n_tau = (
                (num_bonds - specification["mu"])
                - a_1
                - (l - 2 * specification["mu"] - 1)
            )
    else:
        n_phi = 4 * (num_bonds - specification["mu"]) - 3 * num_atoms + a_1 - (l - 1)
        n_phi_prime = 2 * (l - 1)
        n_tau = (num_bonds - specification["mu"]) - a_1 - (l - 1)

    n_r = num_bonds - specification["mu"]

    intermolecular_bonds_needed = n_r - max(bonds_length)
    intermolecular_angles_needed = n_phi - max(angles_length)
    intermolecular_linear_angles_needed = n_phi_prime - max(linear_angles_length)
    intermolecular_dihedrals_needed = n_tau - max(dihedrals_length)

    # With the cleaned up coordinates we can define out total ICs for the appending algorithm:

    # Using the distribute elements function we append the missing ICs out of the updated intermolecular ones

    ic_dict_list = distribute_elements(
        fixed_ic_dict, total_intermolecular_bonds, "bonds", intermolecular_bonds_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_angles,
        "angles",
        intermolecular_angles_needed,
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_linear_angles,
        "linear valence angles",
        intermolecular_linear_angles_needed,
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_dihedrals,
        "dihedrals",
        intermolecular_dihedrals_needed,
    )

    ic_dict = dict()

    new_key = 0
    for dictionary in ic_dict_list:
        for key, value in dictionary.copy().items():
            ic_dict[new_key] = value
            new_key += 1
    return ic_dict


"""
Planar Molecules
"""


def intermolecular_planar_acyclic_linunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    l,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the intermolecular, planar, acyclic case where linear submolecules are contained

    Attributes:
        ic_dict:
            a dictionary containing valid IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing angles
        out_of_plane:
            a list of tupels containing oop's
        dihedrals:
            a list of tupels containing dihedrals
        a_1:
            a integer with the number of terminal atoms
        l:
            a integer with the number of linear bonds
        num_bonds:
            a integer with the number of bonds 
        num_atoms:
            a integer with the number of atoms
        specification:
            the specification dictionary of the given molecule 
    
    Returns:
        a dictionary where each entry is a possible IC set
    """

    args = arguments.get_args()

    # Open up the dictionary
    ic_dict_list = []
    ic_dic_list_1d_complex = []
    # first we determine the two disconnected submolecules
    connected_components, submolecules_bond_list, _ = atoms_list.detect_submolecules()
    submolecules_atoms_coordinates = extract_atoms_of_submolecules(
        connected_components, atoms_list
    )

    summarized_submolecule_ic_dict = []
    # now the hard part every submolecule gets its own specification:
    n = 1
    for submolecule_atom, submolecule_bonds in zip(
        submolecules_atoms_coordinates, submolecules_bond_list
    ):
        # Step 1: Point Group

        molecule = mg.Molecule(
            [strip_numbers(atom.symbol) for atom in submolecule_atom],
            [atom.coordinates for atom in submolecule_atom],
        )
        molecule_pg = PointGroupAnalyzer(molecule)
        point_group_sch = molecule_pg.sch_symbol

        submolecule_atom = Molecule(submolecule_atom)
        # Step 1.1: Calculate the Connectivity C and pass it to specification (just a quick fix)

        molecular_graph = submolecule_atom.graph_rep()
        specifications.connectivity_c = submolecule_atom.count_connected_components(
            molecular_graph
        )

        # Step 2: Initialize Internal Coordinates:
        sub_angles, sub_linear_angles = submolecule_atom.generate_angles(
            submolecule_bonds
        )
        sub_dihedrals = submolecule_atom.generate_dihedrals(submolecule_bonds)

        # Step 3: Give the Specification for each molecule
        specification_submolecule = dict()
        specification_submolecule = specifications.calculation_specification(
            specification,
            submolecule_atom,
            molecule_pg,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
        )

        # Step 4: With the specification generate OOP ICs

        if specification["planar"] == "yes":
            sub_oop = submolecule_atom.generate_out_of_plane(submolecule_bonds)
        elif (
            specification["planar"] == "no"
            and specification["planar submolecule(s)"] == []
        ):
            sub_oop = []
        elif specification["planar"] == "no" and not (
            specification["planar submolecule(s)"] == []
        ):
            sub_oop = atoms.generate_oop_planar_subunits(
                specification["planar submolecule(s)"]
            )

        # Step 5: Calculate IDOF for the given Submolecules:
        idof_submolecule = 0
        if specification_submolecule["linearity"] == "fully linear":
            idof_submolecule = 3 * len(submolecule_atom) - 5
        else:
            idof_submolecule = 3 * len(submolecule_atom) - 6

        # Command Line Message about IC Set generation
        print(f"IC Set generation for submolecule {n}")
        n += 1

        # Step 5: With the Given Specification we can calculate the submolecule ic_dict:
        ic_dict_submolecule = icsel.get_sets(
            idof_submolecule,
            out,
            submolecule_atom,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
            sub_oop,
            sub_dihedrals,
            specification_submolecule,
        )
        summarized_submolecule_ic_dict.append(ic_dict_submolecule)

    # logging for generation of submolecule IC sets

    logfile.write_logfile_submolecule_treatment(
        out, summarized_submolecule_ic_dict, connected_components
    )

    # What we now do with this list of dictionary is to fix the internal coordinates of the two submolecules, now we need to evaluate the missing coordinates for a complete set
    fixed_ic_dict = combine_dictionaries(summarized_submolecule_ic_dict)

    # Now we evaluate the number of internal coordinates we have so far
    bonds_length = []
    angles_length = []
    linear_angles_length = []
    dihedrals_length = []
    out_of_plane_length = []
    for d in fixed_ic_dict:
        for key, value in d.items():
            bonds_length.append(len(value.get("bonds", [])))
            angles_length.append(len(value.get("angles", [])))
            linear_angles_length.append(len(value.get("linear valence angles", [])))
            dihedrals_length.append(len(value.get("dihedrals", [])))
            out_of_plane_length.append(len(value.get("out of plane angles", [])))

    # Now we can specify how much of the different intermolecular coordinates we need in our sets:
    # According to Decius the out of plane anges are defined as n_gamma = 2*(b-a) + a_1 - x
    # And the Dihedral Angles are defined as n_tau = b- a_1 - y
    # in the end x + y = l - 1

    n_r = num_bonds
    n_phi = 2 * num_bonds - num_atoms - (l - 1)
    n_phi_prime = 2 * (l - 1)
    n_gamma = 2 * (num_bonds - num_atoms) + a_1
    n_tau = num_bonds - a_1

    # so therefore we have to update the whole definition
    linear_bonds = specifications.get_linear_bonds(linear_angles)

    # remove dihedrals if terminal

    for linear_bond in linear_bonds:
        if (
            get_multiplicity(linear_bond[0], specification["multiplicity"]) == 1
            or get_multiplicity(linear_bond[1], specification["multiplicity"]) == 1
        ):
            intermolecular_dihedrals_updated = update_internal_coordinates_cyclic(
                [linear_bond], Total_IC_dict["h_bond_dihedrals"]
            )
            n_tau -= l - 1
        else:
            # if there are no terminal dihedrals present
            intermolecular_dihedrals_updated = Total_IC_dict["h_bond_dihedrals"]

    # correct n_gamma if it defined for a linear submolecule
    # as we have 3 oop angles per central unit we need to divide by 3!

    intermolecular_out_of_plane_updated = update_internal_coordinates_cyclic(
        linear_bonds, Total_IC_dict["h_bond_oop"]
    )
    n_gamma = n_gamma - (
        (len(Total_IC_dict["h_bond_oop"]) - len(intermolecular_out_of_plane_updated))
        // 3
    )

    ic_bonds_needed = n_r - max(bonds_length)
    ic_angles_needed = n_phi - max(angles_length)
    ic_dihedrals_needed = n_tau - max(dihedrals_length)
    ic_linear_angles_needed = n_phi_prime - max(linear_angles_length)
    ic_out_of_plane_angles_needed = n_gamma - max(out_of_plane_length)

    total_intermolecular_bonds = Total_IC_dict["h_bond"]
    total_intermolecular_angles = Total_IC_dict["h_bond_angles"]
    total_intermolecular_linear_angles = Total_IC_dict["h_bond_linear_angles"]
    total_intermolecular_dihedrals = intermolecular_dihedrals_updated
    total_intermolecular_oop = intermolecular_out_of_plane_updated

    if args.comb == 2:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + Total_IC_dict["acc_don"]
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + Total_IC_dict["acc_don_angles"]
        )
        total_intermolecular_linear_angles = (
            total_intermolecular_linear_angles + Total_IC_dict["acc_don_linear_angles"]
        )

    if args.comb == 3:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + Total_IC_dict["acc_don"]
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + Total_IC_dict["acc_don_angles"]
        )
        total_intermolecular_linear_angles = (
            total_intermolecular_linear_angles + Total_IC_dict["acc_don_linear_angles"]
        )
        total_intermolecular_dihedrals = (
            total_intermolecular_dihedrals + Total_IC_dict["acc_don_dihedrals"]
        )

    ic_dict_list = distribute_elements(
        fixed_ic_dict, total_intermolecular_bonds, "bonds", ic_bonds_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list, total_intermolecular_angles, "angles", ic_angles_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_linear_angles,
        "linear valence angles",
        ic_linear_angles_needed,
    )
    ic_dict_list = distribute_elements(
        ic_dict_list, total_intermolecular_dihedrals, "dihedrals", ic_dihedrals_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_oop,
        "out of plane angles",
        ic_out_of_plane_angles_needed,
    )

    ic_dict = dict()
    new_key = 0
    for dictionary in ic_dict_list:
        for key, value in dictionary.copy().items():
            ic_dict[new_key] = value
            new_key += 1
    return ic_dict


def intermolecular_planar_cyclic_linunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    l,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the intermolecular, planar, cyclic case with linear submolecules

    Attributes:
        ic_dict:
            a dictionary containing valid IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing angles
        out_of_plane:
            a list of tupels containing oop's
        dihedrals:
            a list of tupels containing dihedrals
        a_1:
            a integer with the number of terminal atoms
        l:
            a integer with the number of linear bonds
        num_bonds:
            a integer with the number of bonds 
        num_atoms:
            a integer with the number of atoms
        specification:
            the specification dictionary of the given molecule 
    
    Returns:
        a dictionary where each entry is a possible IC set
    """

    args = arguments.get_args()

    # Open up the dictionary
    ic_dict_list = []
    ic_dic_list_1d_complex = []
    # first we determine the two disconnected submolecules
    connected_components, submolecules_bond_list, _ = atoms_list.detect_submolecules()
    submolecules_atoms_coordinates = extract_atoms_of_submolecules(
        connected_components, atoms_list
    )

    summarized_submolecule_ic_dict = []
    # now the hard part every submolecule gets its own specification:
    n = 1
    for submolecule_atom, submolecule_bonds in zip(
        submolecules_atoms_coordinates, submolecules_bond_list
    ):
        # Step 1: Point Group

        molecule = mg.Molecule(
            [strip_numbers(atom.symbol) for atom in submolecule_atom],
            [atom.coordinates for atom in submolecule_atom],
        )
        molecule_pg = PointGroupAnalyzer(molecule)
        point_group_sch = molecule_pg.sch_symbol

        # Convert Submolecule_Atom into Molecule Class
        submolecule_atom = Molecule(submolecule_atom)

        # Step 1.1: Calculate the Connectivity C and pass it to specification (just a quick fix)

        molecular_graph = submolecule_atom.graph_rep()
        specifications.connectivity_c = submolecule_atom.count_connected_components(
            molecular_graph
        )

        # Step 2: Initialize Internal Coordinates:
        sub_angles, sub_linear_angles = submolecule_atom.generate_angles(
            submolecule_bonds
        )
        sub_dihedrals = submolecule_atom.generate_dihedrals(submolecule_bonds)
        # Step 3: Give the Specification for each molecule
        specification_submolecule = dict()
        specification_submolecule = specifications.calculation_specification(
            specification_submolecule,
            submolecule_atom,
            molecule_pg,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
        )

        if specification["planar"] == "yes":
            sub_oop = submolecule_atom.generate_out_of_plane(submolecule_bonds)
        elif (
            specification["planar"] == "no"
            and specification["planar submolecule(s)"] == []
        ):
            sub_oop = []
        elif specification["planar"] == "no" and not (
            specification["planar submolecule(s)"] == []
        ):
            sub_oop = atoms.generate_oop_planar_subunits(
                specification["planar submolecule(s)"]
            )

        # Step 4: Calculate IDOF for the given Submolecules:
        idof_submolecule = 0
        if specification_submolecule["linearity"] == "fully linear":
            idof_submolecule = 3 * len(submolecule_atom) - 5
        else:
            idof_submolecule = 3 * len(submolecule_atom) - 6

        # Command Line Message about IC Set generation
        print(f"IC Set generation for submolecule {n}")
        n += 1

        # Step 5: With the Given Specification we can calculate the submolecule ic_dict:
        ic_dict_submolecule = icsel.get_sets(
            idof_submolecule,
            out,
            submolecule_atom,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
            sub_oop,
            sub_dihedrals,
            specification_submolecule,
        )
        summarized_submolecule_ic_dict.append(ic_dict_submolecule)

    # logging for generation of submolecule IC sets

    logfile.write_logfile_submolecule_treatment(
        out, summarized_submolecule_ic_dict, connected_components
    )

    # What we now do with this list of dictionary is to fix the internal coordinates of the two submolecules, now we need to evaluate the missing coordinates for a complete set
    fixed_ic_dict = combine_dictionaries(summarized_submolecule_ic_dict)

    bonds_length = []
    angles_length = []
    linear_angles_length = []
    dihedrals_length = []
    oop_length = []
    for d in fixed_ic_dict:
        for key, value in d.items():
            bonds_length.append(len(value.get("bonds", [])))
            angles_length.append(len(value.get("angles", [])))
            linear_angles_length.append(len(value.get("linear valence angles", [])))
            dihedrals_length.append(len(value.get("dihedrals", [])))
            oop_length.append(len(value.get("out of plane angles", [])))

    # n_r = b so therefore this can get appended out of the pool of intermolecular ics
    # n_phi = 2b-a-(l-1) here it generally makes a difference if the bond that gets deleted is of the linear type

    # before randomly cutting the redundancies in the intermolecular set we first check for symmetry in the intermolecular bonds

    symmetric_intermolecular_bonds = icsel.get_symm_bonds(
        Total_IC_dict["h_bond"], specification
    )
    symmetric_intermolecular_bonds_list = icsel.get_bond_subsets(
        symmetric_intermolecular_bonds
    )
    valide_atoms = valide_atoms_to_cut(
        Total_IC_dict["h_bond"], specification["multiplicity"]
    )
    removed_bonds = []

    # initialize a bond dictionary for removal of acc don bond
    bond_dict = atoms_list.bond_dict(
        Total_IC_dict["cov_bond"] + Total_IC_dict["h_bond"]
    )

    # if its possible we cut according to symmetry

    for symmetric_bond_group in symmetric_intermolecular_bonds_list:
        if len(symmetric_bond_group) >= specification[
            "mu"
        ] and bonds_are_in_valide_atoms(symmetric_bond_group, valide_atoms):
            removed_bonds, random_var = delete_bonds_symmetry(
                symmetric_bond_group,
                Total_IC_dict["h_bond"],
                specification["mu"],
                valide_atoms,
            )
        if not removed_bonds:
            continue

    intermolecular_h_bonds_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond"]
    )
    intermolecular_h_angles_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond_angles"]
    )
    intermolecular_h_dihedrals_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond_dihedrals"]
    )
    intermolecular_h_linear_angles_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond_linear_angles"]
    )

    # when no symmetric cutting is possible we cut a random bond
    # this is not ideal and leads to more possible solution of the decomposition process

    if len(removed_bonds) == 0:
        removed_bonds = Total_IC_dict["h_bond"][: specification["mu"]]

        intermolecular_h_bonds_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond"]
        )
        intermolecular_h_angles_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_angles"]
        )
        intermolecular_h_dihedrals_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_dihedrals"]
        )
        intermolecular_h_linear_angles_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_linear_angles"]
        )

    removed_acc_don_bonds = []
    # basically we use the hydrogen as a key to get our element out
    for removed_bond in removed_bonds:
        if removed_bond[0].strip(string.digits) == "H":
            removed_acc_don_bonds.append(bond_dict[removed_bond[0]])
        else:
            removed_acc_don_bonds.append(bond_dict[removed_bond[1]])

    intermolecular_acc_don_bonds_updated = update_internal_coordinates_cyclic(
        removed_acc_don_bonds, Total_IC_dict["acc_don"]
    )
    intermolecular_acc_don_angle_updated = update_internal_coordinates_cyclic(
        removed_acc_don_bonds, Total_IC_dict["acc_don_angles"]
    )
    intermolecular_acc_don_dihedrals_updated = update_internal_coordinates_cyclic(
        removed_acc_don_bonds, Total_IC_dict["acc_don_dihedrals"]
    )

    total_intermolecular_bonds = intermolecular_h_bonds_updated
    total_intermolecular_angles = intermolecular_h_angles_updated
    total_intermolecular_linear_angles = intermolecular_h_linear_angles_updated
    total_intermolecular_dihedrals = intermolecular_h_dihedrals_updated

    if args.comb == 2:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + intermolecular_acc_don_bonds_updated
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + intermolecular_acc_don_angle_updated
        )

    if args.comb == 3:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + intermolecular_acc_don_bonds_updated
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + intermolecular_acc_don_angle_updated
        )
        total_intermolecular_dihedrals = (
            total_intermolecular_dihedrals + intermolecular_acc_don_dihedrals_updated
        )

    # now we check our descriptions for the acyclic case according to decius first we also determine if the bond we removed was linear
    # this alters the decius description of the amount of angles needed because l can be reduced
    # in this case we also generate one terminal atom through the deletion of a bond

    linear_bonds = specifications.get_linear_bonds(linear_angles)
    swapped_removed_bonds = removed_bonds[0][::-1]

    n_r = num_bonds - specification["mu"]
    n_gamma = 2 * (num_bonds - num_atoms) + a_1
    n_tau = (num_bonds - specification["mu"]) - a_1

    if removed_bonds[0] in linear_bonds or swapped_removed_bonds in linear_bonds:
        n_phi = (
            2 * (num_bonds - specification["mu"])
            - num_atoms
            - (l - 2 * specification["mu"] - 1)
        )
        n_phi_prime = 2 * ((l - 2 * specification["mu"]) - 1)

    # what we can do immediatly is to find the number of bonds needed for a complete set also the normal angles provide no problem
    intermolecular_bonds_needed = n_r - max(bonds_length)
    intermolecular_angles_needed = n_phi - max(angles_length)
    intermolecular_linear_angles_needed = n_phi_prime - max(linear_angles_length)

    # remove dihedrals if terminal

    for linear_bond in linear_bonds:
        if (
            get_multiplicity(linear_bond[0], specification["multiplicity"]) == 1
            or get_multiplicity(linear_bond[1], specification["multiplicity"]) == 1
        ):
            intermolecular_dihedrals_updated = update_internal_coordinates_cyclic(
                [linear_bond], Total_IC_dict["h_bond_dihedrals"]
            )
            n_tau -= l - 1

    # correct n_gamma if it defined for a linear submolecule
    # as we have 3 oop angles per central unit we need to divide by 3!

    intermolecular_out_of_plane_updated = update_internal_coordinates_cyclic(
        linear_bonds, Total_IC_dict["h_bond_oop"]
    )
    n_gamma = n_gamma - (
        (len(Total_IC_dict["h_bond_oop"]) - len(intermolecular_out_of_plane_updated))
        // 3
    )

    # Calculate Dihedrals and OOP after Cleanup
    intermolecular_dihedrals_needed = n_tau - max(dihedrals_length)
    intermolecular_out_of_plane_needed = n_gamma - max(oop_length)

    # TODO For the Formic Acid Dimer no linear angles get appended, a set is created with additional dihedral and angle coordinates
    # AT THE MOMENT MANUAL GENERATION OF ALTERNATIVE SETS

    ic_dict_list = distribute_elements(
        fixed_ic_dict, total_intermolecular_bonds, "bonds", intermolecular_bonds_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_linear_angles,
        "linear valence angles",
        intermolecular_linear_angles_needed,
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_angles,
        "angles",
        intermolecular_angles_needed,
    )

    # TODO also a quickfix
    fixed_ICs = (
        intermolecular_bonds_needed
        + intermolecular_angles_needed
        + intermolecular_linear_angles_needed
    )

    if ((3 * num_atoms - 6) - fixed_ICs) < intermolecular_dihedrals_needed:
        intermolecular_dihedrals_needed = (3 * num_atoms - 6) - fixed_ICs

    ic_dict_list = distribute_elements(
        ic_dict_list, total_intermolecular_dihedrals, "dihedrals", 2
    )  # intermolecular_dihedrals_needed)

    ic_dict = dict()

    new_key = 0
    for dictionary in ic_dict_list:
        for key, value in dictionary.copy().items():
            ic_dict[new_key] = value
            new_key += 1
    return ic_dict


def intermolecular_planar_cyclic_nolinunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the intermolecular, planar, cyclic case without linear submolecules

    Attributes:
        ic_dict:
            a dictionary containing valid IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing angles
        out_of_plane:
            a list of tupels containing oop's
        dihedrals:
            a list of tupels containing dihedrals
        a_1:
            a integer with the number of terminal atoms
        l:
            a integer with the number of linear bonds
        num_bonds:
            a integer with the number of bonds 
        num_atoms:
            a integer with the number of atoms
        specification:
            the specification dictionary of the given molecule 
    
    Returns:    
        a dictionary containing all the possible IC sets for a given molecule
    """
    args = arguments.get_args()

    # Open up the dictionary
    ic_dict_list = []
    ic_dic_list_1d_complex = []
    # first we determine the two disconnected submolecules
    connected_components, submolecules_bond_list, _ = atoms_list.detect_submolecules()
    submolecules_atoms_coordinates = extract_atoms_of_submolecules(
        connected_components, atoms_list
    )

    summarized_submolecule_ic_dict = []

    n = 1
    for submolecule_atom, submolecule_bonds in zip(
        submolecules_atoms_coordinates, submolecules_bond_list
    ):
        # Step 1: Point Group
        molecule = mg.Molecule(
            [strip_numbers(atom.symbol) for atom in submolecule_atom],
            [atom.coordinates for atom in submolecule_atom],
        )
        molecule_pg = PointGroupAnalyzer(molecule)
        point_group_sch = molecule_pg.sch_symbol

        # Convert submolecule_atom into Molecule_Class
        submolecule_atom = Molecule(submolecule_atom)

        # Step 1.1: Calculate the Connectivity C and pass it to specification (just a quick fix)

        molecular_graph = submolecule_atom.graph_rep()
        specifications.connectivity_c = submolecule_atom.count_connected_components(
            molecular_graph
        )

        # Step 2: Initialize Internal Coordinates:
        sub_angles, sub_linear_angles = submolecule_atom.generate_angles(
            submolecule_bonds
        )
        sub_dihedrals = submolecule_atom.generate_dihedrals(submolecule_bonds)
        # Step 3: Give the Specification for each molecule
        specification_submolecule = dict()
        specification_submolecule = specifications.calculation_specification(
            specification_submolecule,
            submolecule_atom,
            molecule_pg,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
        )
        if specification["planar"] == "yes":
            sub_oop = submolecule_atom.generate_out_of_plane(submolecule_bonds)
        elif (
            specification["planar"] == "no"
            and specification["planar submolecule(s)"] == []
        ):
            sub_oop = []
        elif specification["planar"] == "no" and not (
            specification["planar submolecule(s)"] == []
        ):
            sub_oop = atoms.generate_oop_planar_subunits(
                specification["planar submolecule(s)"]
            )

        idof_submolecule = 0
        if specification_submolecule["linearity"] == "fully linear":
            idof_submolecule = 3 * len(submolecule_atom) - 5
        else:
            idof_submolecule = 3 * len(submolecule_atom) - 6

        # Command Line Message about IC Set generation
        print(f"IC Set generation for submolecule {n}")
        n += 1

        # Step 5: With the Given Specification we can calculate the submolecule ic_dict:
        ic_dict_submolecule = icsel.get_sets(
            idof_submolecule,
            out,
            submolecule_atom,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
            sub_oop,
            sub_dihedrals,
            specification_submolecule,
        )
        summarized_submolecule_ic_dict.append(ic_dict_submolecule)

    # logging for generation of submolecule IC sets

    logfile.write_logfile_submolecule_treatment(
        out, summarized_submolecule_ic_dict, connected_components
    )

    # What we now do with this list of dictionary is to fix the internal coordinates of the two submolecules, now we need to evaluate the missing coordinates for a complete set
    fixed_ic_dict = combine_dictionaries(summarized_submolecule_ic_dict)

    # next we evaluate the length of the internal coordinates fixed by this procedure
    bonds_length = []
    angles_length = []
    linear_angles_length = []
    dihedrals_length = []
    oop_length = []
    for d in fixed_ic_dict:
        for key, value in d.items():
            bonds_length.append(len(value.get("bonds", [])))
            angles_length.append(len(value.get("angles", [])))
            linear_angles_length.append(len(value.get("linear valence angles", [])))
            dihedrals_length.append(len(value.get("dihedrals", [])))
            oop_length.append(len(value.get("out of plane angles", [])))

    # before randomly cutting the redundancies in the intermolecular set we first check for symmetry in the intermolecular bonds
    symmetric_intermolecular_bonds = icsel.get_symm_bonds(
        Total_IC_dict["h_bond"], specification
    )
    symmetric_intermolecular_bonds_list = icsel.get_bond_subsets(
        symmetric_intermolecular_bonds
    )
    valide_atoms = valide_atoms_to_cut(
        Total_IC_dict["h_bond"], specification["multiplicity"]
    )
    removed_bonds = []

    # initialize a bond dictionary for removal of acc don bond
    bond_dict = atoms_list.bond_dict(
        Total_IC_dict["cov_bond"] + Total_IC_dict["h_bond"]
    )

    # if its possible we cut according to symmetry

    for symmetric_bond_group in symmetric_intermolecular_bonds_list:
        if len(symmetric_bond_group) >= specification[
            "mu"
        ] and bonds_are_in_valide_atoms(symmetric_bond_group, valide_atoms):
            removed_bonds, random_var = delete_bonds_symmetry(
                symmetric_bond_group,
                Total_IC_dict["h_bond"],
                specification["mu"],
                valide_atoms,
            )
        if not removed_bonds:
            continue

    intermolecular_h_bonds_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond"]
    )
    intermolecular_h_angles_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond_angles"]
    )
    intermolecular_h_dihedrals_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond_dihedrals"]
    )
    intermolecular_h_linear_angles_updated = update_internal_coordinates_cyclic(
        removed_bonds, Total_IC_dict["h_bond_linear_angles"]
    )

    # when no symmetric cutting is possible we cut a random bond
    # this is not ideal and leads to more possible solution of the decomposition process

    if len(removed_bonds) == 0:
        removed_bonds = Total_IC_dict["h_bond"][: specification["mu"]]

        intermolecular_h_bonds_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond"]
        )
        intermolecular_h_angles_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_angles"]
        )
        intermolecular_h_dihedrals_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_dihedrals"]
        )
        intermolecular_h_linear_angles_updated = update_internal_coordinates_cyclic(
            removed_bonds, Total_IC_dict["h_bond_linear_angles"]
        )

    removed_acc_don_bonds = []
    # basically we use the hydrogen as a key to get our element out
    for removed_bond in removed_bonds:
        if removed_bond[0].strip(string.digits) == "H":
            removed_acc_don_bonds.append(bond_dict[removed_bond[0]])
        else:
            removed_acc_don_bonds.append(bond_dict[removed_bond[1]])

    intermolecular_acc_don_bonds_updated = update_internal_coordinates_cyclic(
        removed_acc_don_bonds, Total_IC_dict["acc_don"]
    )
    intermolecular_acc_don_angle_updated = update_internal_coordinates_cyclic(
        removed_acc_don_bonds, Total_IC_dict["acc_don_angles"]
    )
    intermolecular_acc_don_dihedrals_updated = update_internal_coordinates_cyclic(
        removed_acc_don_bonds, Total_IC_dict["acc_don_dihedrals"]
    )

    # Here one can define which coordinates get added!
    total_intermolecular_bonds = (
        intermolecular_h_bonds_updated + intermolecular_acc_don_bonds_updated
    )
    total_intermolecular_angles = (
        intermolecular_h_angles_updated + intermolecular_acc_don_angle_updated
    )
    total_intermolecular_linear_angles = intermolecular_h_linear_angles_updated
    total_intermolecular_dihedrals = (
        intermolecular_h_dihedrals_updated + intermolecular_acc_don_dihedrals_updated
    )

    total_intermolecular_bonds = intermolecular_h_bonds_updated
    total_intermolecular_angles = intermolecular_h_angles_updated
    total_intermolecular_linear_angles = intermolecular_h_linear_angles_updated
    total_intermolecular_dihedrals = intermolecular_h_dihedrals_updated

    if args.comb == 2:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + intermolecular_acc_don_bonds_updated
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + intermolecular_acc_don_angle_updated
        )

    if args.comb == 3:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + intermolecular_acc_don_bonds_updated
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + intermolecular_acc_don_angle_updated
        )
        total_intermolecular_dihedrals = (
            total_intermolecular_dihedrals + intermolecular_acc_don_dihedrals_updated
        )

    # Evaluate the Internal Coordinates needed
    n_r = num_bonds - specification["mu"]
    n_phi = 2 * (num_bonds - specification["mu"]) - num_atoms
    n_gamma = 2 * ((num_bonds - specification["mu"]) - num_atoms) + a_1
    n_tau = (num_bonds - specification["mu"]) - a_1
    intermolecular_bonds_needed = n_r - max(bonds_length)
    intermolecular_angles_needed = n_phi - max(angles_length)
    intermolecular_oop_needed = n_gamma - max(oop_length)
    intermolecular_dihedrals_needed = n_tau - max(dihedrals_length)

    # We also need to make a safety feature if the amount of one intermolecular coordinate is not large enough
    # for this we take the heavy coordinate

    ic_dict_list = distribute_elements(
        fixed_ic_dict, total_intermolecular_bonds, "bonds", intermolecular_bonds_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_angles,
        "angles",
        intermolecular_angles_needed,
    )
    # If not enough dihedrals present just append acc don dihedrals, ignore deletion
    if len(total_intermolecular_dihedrals) < intermolecular_dihedrals_needed:
        total_intermolecular_dihedrals = (
            intermolecular_h_dihedrals_updated + Total_IC_dict["acc_don_dihedrals"]
        )
        ic_dict_list = distribute_elements(
            ic_dict_list,
            total_intermolecular_dihedrals,
            "dihedrals",
            intermolecular_dihedrals_needed,
        )
    else:
        ic_dict_list = distribute_elements(
            ic_dict_list,
            total_intermolecular_dihedrals,
            "dihedrals",
            intermolecular_dihedrals_needed,
        )

    new_key = 0
    for dictionary in ic_dict_list:
        for key, value in dictionary.copy().items():
            ic_dict[new_key] = value
            new_key += 1
    return ic_dict


def intermolecular_planar_acyclic_nolinunit_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    specification,
) -> dict:
    """
    Generates all possible IC sets for the intermolecular, planar, acyclic case without linear submolecules

    Attributes:
        ic_dict:
            a dictionary containing valid IC sets
        out:
            the output file of nomodeco
        idof:
            a integer with the vibrational degrees of freedom
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing angles
        out_of_plane:
            a list of tupels containing oop's
        dihedrals:
            a list of tupels containing dihedrals
        a_1:
            a integer with the number of terminal atoms
        l:
            a integer with the number of linear bonds
        num_bonds:
            a integer with the number of bonds 
        num_atoms:
            a integer with the number of atoms
        specification:
            the specification dictionary of the given molecule 
    
    Returns:
        a dictionary where each entry is a possible IC set
    """

    args = arguments.get_args()

    # Open up the dictionary

    ic_dict_list = []
    ic_dic_list_1d_complex = []
    # first we determine the two disconnected submolecules
    connected_components, submolecules_bond_list, _ = atoms_list.detect_submolecules()
    submolecules_atoms_coordinates = extract_atoms_of_submolecules(
        connected_components, atoms_list
    )

    summarized_submolecule_ic_dict = []
    # now the hard part every submolecule gets its own specification:
    n = 1
    for submolecule_atom, submolecule_bonds in zip(
        submolecules_atoms_coordinates, submolecules_bond_list
    ):
        # Step 1: Point Group

        molecule = mg.Molecule(
            [strip_numbers(atom.symbol) for atom in submolecule_atom],
            [atom.coordinates for atom in submolecule_atom],
        )
        molecule_pg = PointGroupAnalyzer(molecule)
        point_group_sch = molecule_pg.sch_symbol

        # Transform into Molecule Class
        submolecule_atom = Molecule(submolecule_atom)
        # Step 1.1: Calculate the Connectivity C and pass it to specification (just a quick fix)

        molecular_graph = submolecule_atom.graph_rep()
        specifications.connectivity_c = submolecule_atom.count_connected_components(
            molecular_graph
        )

        # Step 2: Initialize Internal Coordinates:
        sub_angles, sub_linear_angles = submolecule_atom.generate_angles(
            submolecule_bonds
        )
        sub_dihedrals = submolecule_atom.generate_dihedrals(submolecule_bonds)

        # Step 3: Give the Specification for each molecule
        specification_submolecule = dict()
        specification_submolecule = specifications.calculation_specification(
            specification,
            submolecule_atom,
            molecule_pg,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
        )

        if specification["planar"] == "yes":
            sub_oop = submolecule_atom.generate_out_of_plane(submolecule_bonds)
        elif (
            specification["planar"] == "no"
            and specification["planar submolecule(s)"] == []
        ):
            sub_oop = []
        elif specification["planar"] == "no" and not (
            specification["planar submolecule(s)"] == []
        ):
            sub_oop = atoms.generate_oop_planar_subunits(
                specification["planar submolecule(s)"]
            )

        # Step 4: Calculate IDOF for the given Submolecules:
        idof_submolecule = 0
        if specification_submolecule["linearity"] == "fully linear":
            idof_submolecule = 3 * len(submolecule_atom) - 5
        else:
            idof_submolecule = 3 * len(submolecule_atom) - 6

        # Command Line Message about IC Set generation
        print(f"IC Set generation for submolecule {n}")
        n += 1

        # Step 5: With the Given Specification we can calculate the submolecule ic_dict:
        ic_dict_submolecule = icsel.get_sets(
            idof_submolecule,
            out,
            submolecule_atom,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
            sub_oop,
            sub_dihedrals,
            specification_submolecule,
        )
        summarized_submolecule_ic_dict.append(ic_dict_submolecule)

    # logging for generation of submolecule IC sets

    logfile.write_logfile_submolecule_treatment(
        out, summarized_submolecule_ic_dict, connected_components
    )

    # What we now do with this list of dictionary is to fix the internal coordinates of the two submolecules, now we need to evaluate the missing coordinates for a complete set
    fixed_ic_dict = combine_dictionaries(summarized_submolecule_ic_dict)

    # Now we evaluate the number of internal coordinates we have so far
    bonds_length = []
    angles_length = []
    linear_angles_length = []
    dihedrals_length = []
    out_of_plane_length = []
    for d in fixed_ic_dict:
        for key, value in d.items():
            bonds_length.append(len(value.get("bonds", [])))
            angles_length.append(len(value.get("angles", [])))
            linear_angles_length.append(len(value.get("linear valence angles", [])))
            dihedrals_length.append(len(value.get("dihedrals", [])))
            out_of_plane_length.append(len(value.get("out of plane angles", [])))
    # Now we can specify how much of the different intermolecular coordinates we need in our sets:
    n_r = num_bonds
    n_phi = 2 * num_bonds - num_atoms
    n_gamma = 2 * (num_bonds - num_atoms) + a_1
    n_tau = num_bonds - a_1

    ic_bonds_needed = n_r - max(bonds_length)
    ic_angles_needed = n_phi - max(angles_length)
    ic_dihedrals_needed = n_tau - max(dihedrals_length)
    ic_out_of_plane_angles_needed = n_gamma - max(out_of_plane_length)

    # Here we define out total bonds using to Total_IC_dict

    # Again if the set size is not feasible the Acc_Don Coordinate can be commented out

    total_intermolecular_dihedrals = (
        Total_IC_dict["h_bond_dihedrals"] + Total_IC_dict["acc_don_dihedrals"]
    )
    total_intermolecular_oop = Total_IC_dict["h_bond_oop"]

    total_intermolecular_bonds = Total_IC_dict["h_bond"]
    total_intermolecular_angles = Total_IC_dict["h_bond_angles"]
    total_intermolecular_dihedrals = Total_IC_dict["h_bond_dihedrals"]
    total_intermolecular_oop = Total_IC_dict["h_bond_oop"]

    if args.comb == 2:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + Total_IC_dict["acc_don"]
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + Total_IC_dict["acc_don_angles"]
        )

    if args.comb == 3:
        total_intermolecular_bonds = (
            total_intermolecular_bonds + Total_IC_dict["acc_don"]
        )
        total_intermolecular_angles = (
            total_intermolecular_angles + Total_IC_dict["acc_don_angles"]
        )
        total_intermolecular_dihedrals = (
            total_intermolecular_dihedrals + Total_IC_dict["acc_don_dihedrals"]
        )

    ic_dict_list = distribute_elements(
        fixed_ic_dict, total_intermolecular_bonds, "bonds", ic_bonds_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list, total_intermolecular_angles, "angles", ic_angles_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list, total_intermolecular_dihedrals, "dihedrals", ic_dihedrals_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        total_intermolecular_oop,
        "out of plane angles",
        ic_out_of_plane_angles_needed,
    )

    ic_dict = dict()
    new_key = 0
    for dictionary in ic_dict_list:
        for key, value in dictionary.copy().items():
            ic_dict[new_key] = value
            new_key += 1
    return ic_dict


"""
Fully Linear
"""


def intermolecular_fully_linear_molecule(
    ic_dict,
    out,
    idof,
    bonds,
    angles,
    linear_angles,
    out_of_plane,
    dihedrals,
    num_bonds,
    num_atoms,
    a_1,
    specification,
):

    n_phi_prime = 2 * (num_atoms - 2)
    # for fully linear molecules Decius work gives us n_r = b and npi = 2(a-2)
    # yet again we first fix the two submolecules, evaluate the missing angles and append them out of the pool of linear angles
    covalent_bonds = list(set(bonds).difference(set(intermolecular_bonds)))
    # Open up the dictionary
    ic_dict_list = []
    ic_dic_list_1d_complex = []
    # first we determine the two disconnected submolecules
    connected_components, submolecules_bond_list = detect_submolecules(covalent_bonds)
    submolecules_atoms_coordinates = extract_atoms_of_submolecules(
        connected_components, atoms_list
    )

    summarized_submolecule_ic_dict = []
    # now the hard part every submolecule gets its own specification:
    for submolecule_atom, submolecule_bonds in zip(
        submolecules_atoms_coordinates, submolecules_bond_list
    ):
        # Step 1: Point Group

        molecule = mg.Molecule(
            [strip_numbers(atom.symbol) for atom in submolecule_atom],
            [atom.coordinates for atom in submolecule_atom],
        )
        molecule_pg = PointGroupAnalyzer(molecule)
        point_group_sch = molecule_pg.sch_symbol

        # Step 1.1: Calculate the Connectivity C and pass it to specification (just a quick fix)

        molecular_graph = dfs_connected.graph_rep(submolecule_bonds)
        specifications.connectivity_c = dfs_connected.count_connected_components(
            molecular_graph
        )

        # Step 2: Initialize Internal Coordinates:
        sub_angles, sub_linear_angles = alt_icgen.initialize_angles(
            atoms_list, submolecule_bonds
        )
        sub_dihedrals, sub_improper_dihedrals = alt_icgen.initialize_dihedrals(
            atoms_list, submolecule_bonds
        )
        sub_oop = alt_icgen.alt_oop(atoms_list, submolecule_bonds)

        # Step 3: Give the Specification for each molecule
        specification_submolecule = dict()
        specification_submolecule = specifications.calculation_specification(
            specification_submolecule,
            submolecule_atom,
            molecule_pg,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
        )

        # Step 4: Calculate IDOF for the given Submolecules:
        idof_submolecule = 0
        if specification_submolecule["linearity"] == "fully linear":
            idof_submolecule = 3 * len(submolecule_atom) - 5
        else:
            idof_submolecule = 3 * len(submolecule_atom) - 6

        # Step 5: With the Given Specification we can calculate the submolecule ic_dict:
        ic_dict_submolecule = icsel.get_sets(
            idof_submolecule,
            out,
            submolecule_atom,
            submolecule_bonds,
            sub_angles,
            sub_linear_angles,
            sub_oop,
            sub_dihedrals,
            specification_submolecule,
        )
        summarized_submolecule_ic_dict.append(ic_dict_submolecule)

    # logging for generation of submolecule IC sets

    logfile.write_logfile_submolecule_treatment(
        out, summarized_submolecule_ic_dict, connected_components
    )

    # What we now do with this list of dictionary is to fix the internal coordinates of the two submolecules, now we need to evaluate the missing coordinates for a complete set
    fixed_ic_dict = combine_dictionaries(summarized_submolecule_ic_dict)

    fixed_bonds = []
    fixed_linear_angles = []

    for d in fixed_ic_dict:
        for key, value in d.items():
            fixed_bonds.extend(value.get("bonds"))
            fixed_linear_angles.extend(value.get("linear valence angles"))

    linear_angles_needed = n_phi_prime - len(fixed_linear_angles)
    bonds_needed = num_bonds - len(fixed_bonds)

    ic_dict_list = distribute_elements(
        fixed_ic_dict, intermolecular_bonds, "bonds", bonds_needed
    )
    ic_dict_list = distribute_elements(
        ic_dict_list,
        intermolecular_linear_angles,
        "linear valence angles",
        linear_angles_needed,
    )

    ic_dict = dict()

    new_key = 0
    for dictionary in ic_dict_list:
        for key, value in dictionary.copy().items():
            ic_dict[new_key] = value
            new_key += 1
    return ic_dict
