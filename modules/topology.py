'''''
This module contains the decision tree derived by Decius' work on complete sets and redundancies among
small vibrational coordinates
(https://doi.org/10.1063/1.1747158)
'''''

import networkx as nx
import itertools
import logging
from . import logfile
from . import specifications
from . import icsel

'''''
GENERAL PURPOSE FUNCTIONS
'''''


def molecule_is_not_split(bonds):
    G = nx.Graph()
    G.add_edges_from(bonds)
    if len(list(nx.connected_components(G))) == 1:
        return True
    else:
        return False


def valide_atoms_to_cut(bonds, multiplicity_list):
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
def update_internal_coordinates_cyclic(removed_bonds, ic_list):
    ic_list_dup = ic_list[:]
    for bond in removed_bonds:
        for ic in ic_list_dup[:]:
            if bond[0] in ic and bond[1] in ic:
                ic_list_dup.remove(ic)
    return ic_list_dup


def find_common_index_with_most_subsets(list1, list2):
    combined_lengths = [len(sublist1) + len(sublist2) for sublist1, sublist2 in zip(list1, list2)]
    max_combined_length = max(combined_lengths)
    common_index = combined_lengths.index(max_combined_length)
    return common_index


def find_common_index_with_least_subsets(list1, list2):
    combined_lengths = [len(sublist1) + len(sublist2) for sublist1, sublist2 in zip(list1, list2)]
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
                n_phi += (2 * atom_and_mult[1] - 3)
    return n_phi, n_gamma, angles


def get_multiplicity(atom_name, multiplicity_list):
    for atom_and_mult in multiplicity_list:
        if atom_name == atom_and_mult[0]:
            return atom_and_mult[1]


'''''
LINEAR SYSTEMS
'''''


def fully_linear_molecule(ic_dict, bonds, angles, linear_angles, out_of_plane, dihedrals):
    # purely linear molecules do not have oop, one can define dihedrals, but they are not significant as
    # the intrinsic frequency equals 0 for them
    ic_dict[0] = {
        "bonds": bonds,
        "angles": angles,
        "linear valence angles": linear_angles,
        "out of plane angles": out_of_plane,
        "dihedrals": []}
    return ic_dict


# TODO: linunits need revision
'''''
PLANAR SYSTEMS
'''''


def planar_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds,
                                      num_atoms, a_1, specification):
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
        angles = remove_angles((oop_central_atom, get_multiplicity(oop_central_atom, specification["multiplicity"])),
                               angles)

    symmetric_angles = icsel.get_symm_angles(angles, specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles), idof, n_phi)

    # in cases of restrictive symmetry, we need to break angle symmetry
    if len(angle_subsets) == 0:
        logging.warning("In order to obtain angle subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset))

            # the if statement ensures, that oop angles to the same central atom can not be in the same set
    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset):
            oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals, specification)
    dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles), idof, n_tau)

    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning("In order to obtain dihedral subsets, symmetry needs to be broken!")
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
                    "dihedrals": dihedral_subsets[len_dihedrals]}
                k += 1

    return ic_dict


# in the cyclic cases, the molecule is rendered acyclic (by removing random but possible mu bonds) and then the acylic function is called
def planar_cyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds,
                                     num_atoms, a_1, specification):
    # remove bonds without destroying the molecule
    # if there are several classes of symmetric bonds, we need to remove the corresponding
    # symmetric bonds

    symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
    symmetric_bonds_list = icsel.get_bond_subsets(symmetric_bonds)
    valide_atoms = valide_atoms_to_cut(bonds, specification["multiplicity"])
    ic_dict_list = []
    removed_bonds = []
    for symmetric_bond_group in symmetric_bonds_list:
        if len(symmetric_bond_group) >= specification["mu"] and bonds_are_in_valide_atoms(
                symmetric_bond_group, valide_atoms):
            removed_bonds, bonds = delete_bonds_symmetry(symmetric_bond_group, bonds, specification["mu"], valide_atoms)
        if not removed_bonds:
            continue

        # update bonds, angles, oop, and dihedrals to not include the coordinates that were removed 
        bonds_updated = update_internal_coordinates_cyclic(removed_bonds, bonds)
        angles_updated = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane_updated = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals_updated = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(out, bonds_updated, angles_updated,
                                                linear_angles, out_of_plane_updated, dihedrals_updated, specification)

        # we need to do some pre-calc of the symmetric angles and dihedrals etc. sadly
        # so that we do not sample a subspace, which is not feasible

        symmetric_angles = icsel.get_symm_angles(angles_updated, specification)
        angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds_updated), len(angles_updated), idof,
                                                2 * len(bonds_updated) - num_atoms)
        if len(angle_subsets) == 0:
            logging.warning(
                "For this rendered molecule, angle symmetry can not be considered and hence this subspace of internal coordinates will be skipped")
            continue

        symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals_updated, specification)
        dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds_updated), len(angles_updated),
                                                      idof, len(bonds_updated) - a_1)
        if len(dihedral_subsets) == 0:
            logging.warning(
                "For this rendered molecule, dihedral symmetry can not be considered and hence this subspace of internal coordintes will be skipped")
            continue

        # call the acyclic version

        ic_dict_list.append(planar_acyclic_nolinunit_molecule(dict(), out, idof, bonds_updated, angles_updated,
                                                              linear_angles, out_of_plane_updated, dihedrals_updated,
                                                              len(bonds_updated), num_atoms, a_1, specification))
        removed_bonds = []

    # if we can't cut according to symmetry, do random cutting
    # cut symmetry out if you want, by commenting everyting out
    if not ic_dict_list:
        removed_bonds, bonds = delete_bonds(bonds, specification["mu"], valide_atoms)
        angles = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(out, bonds, angles,
                                                linear_angles, out_of_plane, dihedrals, specification)

        ic_dict = planar_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles,
                                                    linear_angles, out_of_plane, dihedrals, len(bonds), num_atoms, a_1,
                                                    specification)
        return ic_dict

    else:
        ic_dict = dict()
        new_key = 0
        for dictionary in ic_dict_list:
            for key, value in dictionary.copy().items():
                ic_dict[new_key] = value
                new_key += 1
        return ic_dict


def planar_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds,
                                    num_atoms, a_1, l, specification):
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
        angles = remove_angles((oop_central_atom, get_multiplicity(oop_central_atom, specification["multiplicity"])),
                               angles)

    n_gamma = 2 * (num_bonds - num_atoms) + a_1
    n_tau = num_bonds - a_1

    # before computing the number of ICs we will remove all oop that are associated with this linear angle
    # also remove dihedrals,if they are terminal ==> we will then also reset the number of internals

    linear_bonds = specifications.get_linear_bonds(linear_angles)

    # remove dihedrals if terminal
    for linear_bond in linear_bonds:
        if get_multiplicity(linear_bond[0], specification["multiplicity"]) == 1 or get_multiplicity(
                linear_bond[1], specification["multiplicity"]) == 1:
            dihedrals = update_internal_coordinates_cyclic([linear_bond], dihedrals)
            n_tau -= (l - 1)

    # corret n_gamma if it defined for a linear submolecule
    # as we have 3 oop angles per central unit we need to divide by 3!
    out_of_plane_updated = update_internal_coordinates_cyclic(linear_bonds, out_of_plane)
    n_gamma = n_gamma - ((len(out_of_plane) - len(out_of_plane_updated)) // 3)
    out_of_plane = out_of_plane_updated

    logfile.write_logfile_updatedICs_linunit(out, out_of_plane, dihedrals)

    symmetric_angles = icsel.get_symm_angles(angles, specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles), idof, n_phi)

    # symmetry breaking if needed
    if len(angle_subsets) == 0:
        logging.warning("In order to obtain angle subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset))

    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset):
            oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals, specification)
    dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles), idof, n_tau)

    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning("In order to obtain dihedral subsets, symmetry needs to be broken!")
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
                    "dihedrals": dihedral_subsets[len_dihedrals]}
                k += 1

    return ic_dict

def planar_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds,
                                   num_atoms, a_1, l, specification):
    # remove bonds without destroying the molecule
    # if there are several classes of symmetric bonds, we need to remove the corresponding
    # symmetric bonds

    symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
    symmetric_bonds_list = icsel.get_bond_subsets(symmetric_bonds)
    valide_atoms = valide_atoms_to_cut(bonds, specification["multiplicity"])
    ic_dict_list = []
    removed_bonds = []
    for symmetric_bond_group in symmetric_bonds_list:
        if len(symmetric_bond_group) >= specification["mu"] and bonds_are_in_valide_atoms(
                symmetric_bond_group, valide_atoms):
            removed_bonds, bonds = delete_bonds_symmetry(symmetric_bond_group, bonds, specification["mu"], valide_atoms)
        if not removed_bonds:
            continue

        # update bonds, angles, oop, and dihedrals to not include the coordinates that were removed 
        bonds_updated = update_internal_coordinates_cyclic(removed_bonds, bonds)
        angles_updated = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane_updated = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals_updated = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(out, bonds_updated, angles_updated,
                                                linear_angles, out_of_plane_updated, dihedrals_updated, specification)

        # we need to do some pre-calc of the symmetric angles and dihedrals etc. sadly
        # so that we do not sample a subspace, which is not feasible

        symmetric_angles = icsel.get_symm_angles(angles_updated, specification)
        angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds_updated), len(angles_updated), idof,
                                                2 * len(bonds_updated) - num_atoms)
        if len(angle_subsets) == 0:
            logging.warning(
                "For this rendered molecule, angle symmetry can not be considered and hence this subspace of internal coordinates will be skipped")
            continue

        symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals_updated, specification)
        dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds_updated), len(angles_updated),
                                                      idof, len(bonds_updated) - a_1)
        if len(dihedral_subsets) == 0:
            logging.warning(
                "For this rendered molecule, dihedral symmetry can not be considered and hence this subspace of internal coordintes will be skipped")
            continue

        # call the acyclic version

        ic_dict_list.append(planar_acyclic_linunit_molecule(dict(), out, idof, bonds_updated, angles_updated,
                                                            linear_angles, out_of_plane_updated, dihedrals_updated,
                                                            len(bonds_updated), num_atoms, a_1, specification))
        removed_bonds = []

    # if we can't cut according to symmetry, do random cutting
    # cut symmetry out if you want, by commenting everyting out
    if not ic_dict_list:
        removed_bonds, bonds = delete_bonds(bonds, specification["mu"], valide_atoms)
        angles = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(out, bonds, angles,
                                                linear_angles, out_of_plane, dihedrals, specification)

        ic_dict = planar_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles,
                                                  linear_angles, out_of_plane, dihedrals, len(bonds), num_atoms, a_1,
                                                  specification)
        return ic_dict

    else:
        ic_dict = dict()
        new_key = 0
        for dictionary in ic_dict_list:
            for key, value in dictionary.copy().items():
                ic_dict[new_key] = value
                new_key += 1
        return ic_dict


'''''
GENERAL SYSTEMS
'''''


def general_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds,
                                       num_atoms, a_1, specification):
    # set length of subsets
    n_r = num_bonds
    n_phi = 4 * num_bonds - 3 * num_atoms + a_1
    n_gamma = 0
    planar_subunits_list = specification["planar submolecule(s)"]
    n_tau = num_bonds - a_1

    # if planar subunits exist, we need to do 2 things: change n_phi and n_gamma;
    # remove angles at the specified coordinate, as we else would have linear dependencies
    if len(planar_subunits_list) != 0:
        n_phi, n_gamma, angles = get_param_planar_submolecule(planar_subunits_list, specification["multiplicity"],
                                                              angles)

    symmetric_angles = icsel.get_symm_angles(angles, specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles), idof, n_phi)
    if len(angle_subsets) == 0:
        logging.warning("In order to obtain angle subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset))

    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset):
            oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals, specification)
    dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles), idof, n_tau)

    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning("In order to obtain dihedral subsets, symmetry needs to be broken!")
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
                    "dihedrals": dihedral_subsets[len_dihedrals]}
                k += 1
    return ic_dict


def general_cyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds,
                                      num_atoms, num_of_red, a_1, specification):
    # remove bonds without destroying the molecule
    # if there are several classes of symmetric bonds, we need to remove the corresponding
    # symmetric bonds

    symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
    symmetric_bonds_list = icsel.get_bond_subsets(symmetric_bonds)
    valide_atoms = valide_atoms_to_cut(bonds, specification["multiplicity"])
    ic_dict_list = []
    removed_bonds = []
    for symmetric_bond_group in symmetric_bonds_list:
        if len(symmetric_bond_group) >= specification["mu"] and bonds_are_in_valide_atoms(
                symmetric_bond_group, valide_atoms):
            removed_bonds, bonds = delete_bonds_symmetry(symmetric_bond_group, bonds, specification["mu"], valide_atoms)
        if not removed_bonds:
            continue

        # update bonds, angles, oop, and dihedrals to not include the coordinates that were removed 
        bonds_updated = update_internal_coordinates_cyclic(removed_bonds, bonds)
        angles_updated = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane_updated = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals_updated = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(out, bonds_updated, angles_updated,
                                                linear_angles, out_of_plane_updated, dihedrals_updated, specification)

        # we need to do some pre-calc of the symmetric angles and dihedrals etc. sadly
        # so that we do not sample a subspace, which is not feasible

        symmetric_angles = icsel.get_symm_angles(angles_updated, specification)
        angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds_updated), len(angles_updated), idof,
                                                2 * len(bonds_updated) - num_atoms)
        if len(angle_subsets) == 0:
            logging.warning(
                "For this rendered molecule, angle symmetry can not be considered and hence this subspace of internal coordinates will be skipped")
            continue

        symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals_updated, specification)
        dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds_updated), len(angles_updated),
                                                      idof, len(bonds_updated) - a_1)
        if len(dihedral_subsets) == 0:
            logging.warning(
                "For this rendered molecule, dihedral symmetry can not be considered and hence this subspace of internal coordintes will be skipped")
            continue

        # call the acyclic version

        ic_dict_list.append(general_acyclic_nolinunit_molecule(dict(), out, idof, bonds_updated, angles_updated,
                                                               linear_angles, out_of_plane_updated, dihedrals_updated,
                                                               len(bonds_updated), num_atoms, a_1, specification))
        removed_bonds = []

    # if we can't cut according to symmetry, do random cutting
    # cut symmetry out if you want, by commenting everything out
    if not ic_dict_list:
        removed_bonds, bonds = delete_bonds(bonds, specification["mu"], valide_atoms)
        angles = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(out, bonds, angles,
                                                linear_angles, out_of_plane, dihedrals, specification)

        ic_dict = general_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles,
                                                     linear_angles, out_of_plane, dihedrals, len(bonds), num_atoms, a_1,
                                                     specification)
        return ic_dict

    else:
        ic_dict = dict()
        new_key = 0
        for dictionary in ic_dict_list:
            for key, value in dictionary.copy().items():
                ic_dict[new_key] = value
                new_key += 1
        return ic_dict


def general_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds,
                                     num_atoms, a_1, l, specification):
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
            "Due to high number of linear angles, the topology conditions can not be considered. Linear angles will be removed from the analysis to do so.")
        ic_dict = general_acyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, [], out_of_plane, dihedrals,
                                                     num_bonds, num_atoms, a_1, specification)
        return ic_dict

    # if planar subunits exist, we need to do 2 things: change n_phi and n_gamma; 
    # remove angles at the specified coordinate, as we else would have linear dependencies
    if len(planar_subunits_list) != 0:
        n_phi, n_gamma, angles = get_param_planar_submolecule(planar_subunits_list, specification["multiplicity"],
                                                              angles)
        # correct n_phi because we lose (l-1) DOF
        n_phi = n_phi - (l - 1)

    symmetric_angles = icsel.get_symm_angles(angles, specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles), idof, n_phi)
    if len(angle_subsets) == 0:
        logging.warning("In order to obtain angle subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset))

            # before computing the number of ICs we will remove all oop that are associated with this linear angle
    # also remove dihedrals,if they are terminal

    linear_bonds = specifications.get_linear_bonds(linear_angles)
    out_of_plane = update_internal_coordinates_cyclic(linear_bonds, out_of_plane)
    for linear_bond in linear_bonds:
        if get_multiplicity(linear_bond[0], specification["multiplicity"]) == 1 or get_multiplicity(
                linear_bond[1], specification["multiplicity"]) == 1:
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
    dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles), idof, n_tau)

    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning("In order to obtain dihedral subsets, symmetry needs to be broken!")
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
                    "dihedrals": dihedral_subsets[len_dihedrals]}
                k += 1

    return ic_dict


def general_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds,
                                    num_atoms, num_of_red, a_1, l, specification):
    # remove bonds without destroying the molecule
    # if there are several classes of symmetric bonds, we need to remove the corresponding
    # symmetric bonds

    symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
    symmetric_bonds_list = icsel.get_bond_subsets(symmetric_bonds)
    valide_atoms = valide_atoms_to_cut(bonds, specification["multiplicity"])
    ic_dict_list = []
    removed_bonds = []
    for symmetric_bond_group in symmetric_bonds_list:
        if len(symmetric_bond_group) >= specification["mu"] and bonds_are_in_valide_atoms(
                symmetric_bond_group, valide_atoms):
            removed_bonds, bonds = delete_bonds_symmetry(symmetric_bond_group, bonds, specification["mu"], valide_atoms)
        if not removed_bonds:
            continue

        # update bonds, angles, oop, and dihedrals to not include the coordinates that were removed 
        bonds_updated = update_internal_coordinates_cyclic(removed_bonds, bonds)
        angles_updated = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane_updated = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals_updated = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(out, bonds_updated, angles_updated,
                                                linear_angles, out_of_plane_updated, dihedrals_updated, specification)

        # we need to do some pre-calc of the symmetric angles and dihedrals etc. sadly
        # so that we do not sample a subspace, which is not feasible

        symmetric_angles = icsel.get_symm_angles(angles_updated, specification)
        angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds_updated), len(angles_updated), idof,
                                                2 * len(bonds_updated) - num_atoms)
        if len(angle_subsets) == 0:
            logging.warning(
                "For this rendered molecule, angle symmetry can not be considered and hence this subspace of internal coordinates will be skipped")
            continue

        symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals_updated, specification)
        dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds_updated), len(angles_updated),
                                                      idof, len(bonds_updated) - a_1)
        if len(dihedral_subsets) == 0:
            logging.warning(
                "For this rendered molecule, dihedral symmetry can not be considered and hence this subspace of internal coordintes will be skipped")
            continue

        # call the acyclic version

        ic_dict_list.append(general_acyclic_linunit_molecule(dict(), out, idof, bonds_updated, angles_updated,
                                                             linear_angles, out_of_plane_updated, dihedrals_updated,
                                                             len(bonds_updated), num_atoms, a_1, specification))
        removed_bonds = []

    # if we can't cut according to symmetry, do random cutting
    # cut symmetry out if you want, by commenting everyting out
    if not ic_dict_list:
        removed_bonds, bonds = delete_bonds(bonds, specification["mu"], valide_atoms)
        angles = update_internal_coordinates_cyclic(removed_bonds, angles)
        out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
        dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)

        logfile.write_logfile_updatedICs_cyclic(out, bonds, angles,
                                                linear_angles, out_of_plane, dihedrals, specification)

        ic_dict = general_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles,
                                                   linear_angles, out_of_plane, dihedrals, len(bonds), num_atoms, a_1,
                                                   specification)
        return ic_dict

    else:
        ic_dict = dict()
        new_key = 0
        for dictionary in ic_dict_list:
            for key, value in dictionary.copy().items():
                ic_dict[new_key] = value
                new_key += 1
        return ic_dict
