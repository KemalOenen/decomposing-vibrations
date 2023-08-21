'''''
This module contains the decision tree derived by Decius' work on complete sets and redundancies among
small vibrational coordinates
(https://doi.org/10.1063/1.1747158)
'''''

import pandas as pd
import numpy as np
import networkx as nx
import itertools
import logging
import logfile
import nomodeco
import specifications
import icsel

'''''
GENERAL PURPOSE FUNCTIONS
'''''

def molecule_is_not_split(bonds):
    G=nx.Graph()
    G.add_edges_from(bonds)
    if len(list(nx.connected_components(G))) == 1:
        return True
    else:
        return False

def delete_bonds(bonds, mu, multiplicity_list):
    removed_bonds = []
    while mu > 0:
        # first of all we will identify where cutting bonds is even making sense to a first degree
        valide_atoms = []
        for tup in multiplicity_list:
            if tup[1] >= 3:
                valide_atoms.append(tup[0])
        
        # then we will cut the bonds out
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

#TODO: rename method
def update_internal_coordinates_cyclic(removed_bonds, ic_list):
    for bond in removed_bonds:
        for ic in ic_list[:]:
            if bond[0] in ic and bond[1] in ic:
                ic_list.remove(ic)
    return ic_list

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
                n_gamma += atom_and_mult[1]-2
                angles = remove_angles(atom_and_mult, angles)
            else:
                n_phi += (2*atom_and_mult[1] - 3)
    return n_phi, n_gamma, angles

'''''
LINEAR SYSTEMS
'''''

def fully_linear_molecule(ic_dict, bonds, angles, linear_angles, out_of_plane, dihedrals):
# purely linear molecules do not have oop, one can define dihedrals, but they are not significant as 
# the intrinsic frequency equals 0 for them
    ic_dict[0] = {
            "bonds" : bonds,
            "angles" : angles,
            "linear valence angles" : linear_angles,
            "out of plane angles" : out_of_plane,
            "dihedrals" : []} 
    return ic_dict

'''''
PLANAR SYSTEMS
'''''

def planar_acyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, a_1, specification):

    # set length of subsets
    n_r = num_bonds
    n_phi = 2*num_bonds - num_atoms 
    n_gamma = 2*(num_bonds - num_atoms) + a_1
    n_tau = num_bonds - a_1

    symmetric_angles = icsel.get_symm_angles(angles, specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles), idof, n_phi)

    # in cases of restrictive symmetry, we need to break angle symmetry
    if len(angle_subsets) == 0:
        logging.warning("In order to gain angle subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset)) 

    # the if statement ensures, that oop angles to the same central atom can not be in the same set
    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset): 
            oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)
    dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau)

    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))

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
    
    return ic_dict

#in the cyclic cases, the molecule is rendered acyclic (by removing random but possible mu bonds) and then the acylic function is called
def planar_cyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, a_1, specification):

    # remove bonds without destroying the molecule
    removed_bonds, bonds = delete_bonds(bonds, specification["mu"], specification["multiplicity"])

    # update angles, oop, and dihedrals to not include the coordinates that were removed 
    angles = update_internal_coordinates_cyclic(removed_bonds, angles)
    out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
    dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)
    
    logfile.write_logfile_updatedICs_cyclic(bonds, angles, linear_angles, out_of_plane, dihedrals, specification)
    
    # call the acyclic version
    ic_dict = planar_acyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, 
            linear_angles, out_of_plane, dihedrals, len(bonds), num_atoms, a_1, specification)

    return ic_dict


def planar_acyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, a_1,l, specification):

    # before computing the number of ICs we will remove all oop that are associated with this linear angle
    # you could also remove dihedrals, but looking at Decius' work it is not advised
    
    linear_bonds = specifications.get_linear_bonds(linear_angles)
    out_of_plane = update_internal_coordinates_cyclic(linear_bonds, out_of_plane)

    logfile.write_logfile_updatedICs_linunit(out_of_plane, dihedrals)
    
    # find all combinations of x and y
    x_and_y = []
    for x in range(l+1):
        y = l - 1 - x
        if x >= 0 and y >= 0:
            x_and_y.append([x,y])
    
    # set length of subsets
    #IMPORTANT: there is a distinction to Decius work -> Decius counts one l.A. in n_phi and one in n_phi'!
   
    n_r = num_bonds
    n_phi = 2*num_bonds - num_atoms - (l-1)
    n_phi_prime = 2*(l - 1)
    n_gamma = []
    n_tau = []
    
    for pair in x_and_y:
        n_gamma.append(2*(num_bonds-num_atoms) + a_1 - pair[0])
        n_tau.append(num_bonds - a_1 - pair[1])

    symmetric_angles = icsel.get_symm_angles(angles,specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)

    # symmetry breaking if needed
    if len(angle_subsets) == 0:
        logging.warning("In order to gain angle subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset)) 

    #TODO: we actually include all linear angles anyways! So we not need this right?
    #symmetric_lin_angles = icsel.get_symm_angles(linear_angles,specification)
    #lin_angle_subsets = icsel.get_angle_subsets(symmetric_lin_angles, len(bonds), len(angles),idof,n_phi_prime)

    oop_subsets = []
    for i in range(len(n_gamma)):
        one_gamma_oop_subset = []
        for subset in itertools.combinations(out_of_plane, n_gamma[i]):
            if icsel.not_same_central_atom(subset): 
                one_gamma_oop_subset.append(list(subset))
        oop_subsets.append(list(one_gamma_oop_subset))

    '''''
    Here we need to choose one set for n_gamma and n_tau, otherwise it is computationally unfeasible
    '''''
    # if there are oop subsets which are already empty, ommit them, get the index, so that that n_tau is removed;
    # for the dihedrals that is not generally applicable, see below
    for index, one_gamma_oop_subset in enumerate(oop_subsets):
        if len(one_gamma_oop_subset) == 0:
            oop_subsets.remove(one_gamma_oop_subset)
            n_tau.pop(index)
            n_gamma.pop(index)

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)

    dihedral_subsets = []
    for i in range(len(n_tau)):
        dihedral_subsets.append(list(
            icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau[i])))

    # if all dihedral subsets are empty, then do symmetry breaking with the value that corresponds to n_gamma
    # with the most entries
    if all(not one_tau_dihedral_subset for one_tau_dihedral_subset in dihedral_subsets):
        index, longest_oop_subset = max(enumerate(oop_subsets), key=lambda x: len(x[1]))
        n_tau = n_tau[index]
        n_gamma = n_gamma[index]
        oop_subsets = oop_subsets[index]
        dihedral_subsets = []
        logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))

    # else we will choose the final n_tau and n_gamma to be the ones, where the highest number of individual subsets are!
    else:
        #index = find_common_index_with_most_subsets(oop_subsets, dihedral_subsets)
        index = find_common_index_with_least_subsets(oop_subsets, dihedral_subsets)
        n_tau = n_tau[index]
        n_gamma = n_gamma[index]
        oop_subsets = oop_subsets[index]
        dihedral_subsets = dihedral_subsets[index]

    k = 0
    for len_angles in range(0, len(angle_subsets)):
        for len_oop in range(0, len(oop_subsets)):
            for len_dihedrals in range(0, len(dihedral_subsets)):
                ic_dict[k] = {
                        "bonds" : bonds,
                        "angles" : angle_subsets[len_angles],
                        "linear valence angles" : linear_angles,
                        "out of plane angles" : oop_subsets[len_oop],
                        "dihedrals" : dihedral_subsets[len_dihedrals] }
                k +=1
    
    #print(ic_dict)
    return ic_dict


def planar_cyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, a_1,l, specification):

    # remove bonds without destroying the molecule
    removed_bonds, bonds = delete_bonds(bonds, specification["mu"], specification["multiplicity"])
 
    # update angles, oop, and dihedrals to not include the coordinates that were removed 
    angles = update_internal_coordinates_cyclic(removed_bonds, angles)
    out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
    dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)
    
    logfile.write_logfile_updatedICs_cyclic(bonds, angles, linear_angles, out_of_plane, dihedrals, specification)
    
    # call the acyclic version
    ic_dict = planar_acyclic_linunit_molecule(ic_dict, idof, bonds, angles, 
            linear_angles, out_of_plane, dihedrals, len(bonds), num_atoms, a_1, l, specification)

    return ic_dict


'''''
GENERAL SYSTEMS
'''''

def general_acyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, a_1, specification):

    # set length of subsets
    n_r = num_bonds
    n_phi = 4*num_bonds - 3*num_atoms + a_1
    n_gamma = 0
    planar_subunits_list = specification["planar submolecule(s)"]
    n_tau = num_bonds - a_1


    # if planar subunits exist, we need to do 2 things: change n_phi and n_gamma; 
    # remove angles at the specified coordinate, as we else would have linear dependencies
    if len(planar_subunits_list) != 0:
        n_phi, n_gamma, angles = get_param_planar_submolecule(planar_subunits_list, specification["multiplicity"], angles)

    symmetric_angles = icsel.get_symm_angles(angles,specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)
    if len(angle_subsets) == 0:
        logging.warning("In order to gain angle subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset)) 

    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset): 
            oop_subsets.append(list(subset))
    
    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)
    dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau)


    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))

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
    return ic_dict 

def general_cyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, num_of_red, a_1, specification):

    #remove bonds without destroying the molecule
    removed_bonds, bonds = delete_bonds(bonds, specification["mu"], specification["multiplicity"])

    # update angles, oop, and dihedrals to not include the coordinates that were removed 
    angles = update_internal_coordinates_cyclic(removed_bonds, angles)
    out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
    dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)
    
    logfile.write_logfile_updatedICs_cyclic(bonds, angles, linear_angles, out_of_plane, dihedrals, specification)
    
    # call the acyclic version
    ic_dict = general_acyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, 
            linear_angles, out_of_plane, dihedrals, len(bonds), num_atoms, a_1, specification)

    return ic_dict

def general_acyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, a_1,l, specification):

    # set length of subsets
    n_r = num_bonds
    n_phi = 4*num_bonds - 3*num_atoms + a_1 - (l-1)
    n_gamma = 0
    planar_subunits_list = specification["planar submolecule(s)"]
    n_phi_prime = 2*(l-1)
    n_tau = num_bonds - a_1 - (l-1)
   
   # TODO: for SF6, we will not use the linear angles in the analysis, as we would have faulty values for some parameters
   # so we need a better check here
    if n_tau < 0:
        ic_dict = general_acyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, [], out_of_plane, dihedrals, 
                num_bonds, num_atoms, a_1, specification)
        return ic_dict

    # if planar subunits exist, we need to do 2 things: change n_phi and n_gamma; 
    # remove angles at the specified coordinate, as we else would have linear dependencies
    if len(planar_subunits_list) != 0:
        n_phi, n_gamma, angles = get_param_planar_submolecule(planar_subunits_list, specification["multiplicity"], angles)
        # correct n_phi because we lose again (l-1) DOF
        n_phi = n_phi - (l-1)

    symmetric_angles = icsel.get_symm_angles(angles,specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)
    if len(angle_subsets) == 0:
        logging.warning("In order to gain angle subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset)) 
    
    #TODO: we actually include all linear angles anyways! So we not need this right?
    #symmetric_lin_angles = icsel.get_symm_angles(linear_angles,specification)
    #lin_angle_subsets = icsel.get_angle_subsets(symmetric_lin_angles, len(bonds), len(angles),idof,n_phi_prime)

    oop_subsets = []
    for subset in itertools.combinations(out_of_plane, n_gamma):
        if icsel.not_same_central_atom(subset): 
            oop_subsets.append(list(subset))    

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)
    dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau)
    
    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))

    k = 0
    for len_angles in range(0, len(angle_subsets)):
        for len_oop in range(0, len(oop_subsets)):
            for len_dihedrals in range(0, len(dihedral_subsets)):
                ic_dict[k] = {
                        "bonds" : bonds,
                        "angles" : angle_subsets[len_angles],
                        "linear valence angles" : linear_angles,
                        "out of plane angles" : oop_subsets[len_oop],
                        "dihedrals" : dihedral_subsets[len_dihedrals] }
                k +=1
    
    return ic_dict 


def general_cyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, num_of_red, a_1,l, specification):

    #remove bonds without destroying the molecule
    removed_bonds, bonds = delete_bonds(bonds, specification["mu"], specification["multiplicity"])

    # update angles, oop, and dihedrals to not include the coordinates that were removed 
    angles = update_internal_coordinates_cyclic(removed_bonds, angles)
    out_of_plane = update_internal_coordinates_cyclic(removed_bonds, out_of_plane)
    dihedrals = update_internal_coordinates_cyclic(removed_bonds, dihedrals)
    
    logfile.write_logfile_updatedICs_cyclic(bonds, angles, linear_angles, out_of_plane, dihedrals, specification)
    
    # call the acyclic version
    ic_dict = general_acyclic_linunit_molecule(ic_dict, idof, bonds, angles, 
            linear_angles, out_of_plane, dihedrals, len(bonds), num_atoms, a_1,l, specification)

    return ic_dict
