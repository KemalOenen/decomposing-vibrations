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

def update_internal_coordinates_cyclic(removed_bonds, ic_list):
    for bond in removed_bonds:
        for ic in ic_list[:]:
            if bond[0] in ic and bond[1] in ic:
                ic_list.remove(ic)
    return ic_list

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

    symmetric_angles = icsel.get_symm_angles(angles,specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)
    
    if len(angle_subsets) == 0:
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset)) 
        for subset in itertools.combinations(angles, n_phi+1):
            angle_subsets.append(list(subset))
        for subset in itertools.combinations(angles, n_phi+2):
            angle_subsets.append(list(subset))

    symmetric_lin_angles = icsel.get_symm_angles(linear_angles,specification)
    lin_angle_subsets = icsel.get_angle_subsets(symmetric_lin_angles, len(bonds), len(angles),idof,n_phi_prime)
    
    oop_subsets = []
    for i in range(len(n_gamma)):
        for subset in itertools.combinations(out_of_plane, n_gamma[i]):
            if icsel.not_same_central_atom(subset): 
                oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)
    
    dihedral_subsets = []
    for i in range(len(n_tau)):
        dihedral_subsets.append(list(
            icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau[i])))

    for i in range(len(n_tau)):
        if n_tau[i] != 0 and len(dihedral_subsets) == 0:
                logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
                for subset in itertools.combinations(dihedrals, n_tau[i]):
                    dihedral_subsets.append(list(subset))
                for subset in itertools.combinations(dihedrals, n_tau[i]+1):
                    dihedral_subsets.append(list(subset))
    
    k = 0
    for x_y in range(0,len(n_tau)): # there as many taus and phis, and as x and y were ordered, the values are also ordered
        for len_angles in range(0, len(angle_subsets)):
            for len_oop in range(0, len(oop_subsets[x_y])):
                for len_lin_angles in range(0, len(lin_angle_subsets)):
                    for len_dihedrals in range(0, len(dihedral_subsets[x_y])):
                        ic_dict[k] = {
                                "bonds" : bonds,
                                "angles" : angle_subsets[len_angles],
                                "linear valence angles" : linear_angles[len_lin_angles] ,
                                "out of plane angles" : oop_subsets[x_y][len_oop],
                                "dihedrals" : dihedral_subsets[x_y][len_dihedrals] }
                        k +=1

    return ic_dict


def planar_cyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, a_1,l, specification):

    # find all combinations of x, y AND a,b
    x_and_y = []
    for x in range(l):
        y = l - 1 - x
        x_and_y.append([x,y])
    
    a_and_b = []
    for a in range(3*specification["mu"]):
        b = 3*specification["mu"] - a
        if a >= 0 and b >= 0:
            a_and_b.append([a, b])
            
    # set length of subsets
    n_r = num_bonds
    n_phi = 2*num_bonds - num_atoms
    n_phi_prime = l - 1
    n_gamma = []
    n_tau = []

    for pair1 in a_and_b:
        for pair2 in x_and_y:
            n_gamma.append(2*(num_bonds-num_atoms) + a_1 - pair1[0] - pair2[0])
            n_tau.append(num_bonds - a_1 - pair1[1] - pair2[1])


    symmetric_angles = icsel.get_symm_angles(angles,specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)

    # case, when angle symmetry needs to be broken
    if len(angle_subsets) == 0:
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset)) 
        for subset in itertools.combinations(angles, n_phi+1):
            angle_subsets.append(list(subset))
        for subset in itertools.combinations(angles, n_phi+2):
            angle_subsets.append(list(subset))
   
    symmetric_lin_angles = icsel.get_symm_angles(linear_angles,specification)
    lin_angle_subsets = icsel.get_angle_subsets(symmetric_lin_angles, len(bonds), len(angles),idof,n_phi_prime)


    oop_subsets = []
    for i in range(len(n_gamma)):
        for subset in itertools.combinations(out_of_plane, n_gamma[i]):
            if icsel.not_same_central_atom(subset): 
                oop_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)
    
    dihedral_subsets = []
    for i in range(len(n_tau)):
        dihedral_subsets.append(list(
            icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau[i])))

    for i in range(len(n_tau)):
        if n_tau[i] != 0 and len(dihedral_subsets) == 0:
                logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
                for subset in itertools.combinations(dihedrals, n_tau[i]):
                    dihedral_subsets.append(list(subset))
                for subset in itertools.combinations(dihedrals, n_tau[i]+1):
                    dihedral_subsets.append(list(subset))

    k = 0
    for x_y in range(0,len(n_tau)): # there as many taus and phis, and as x and y were ordered, the values are also ordered
        for len_angles in range(0, len(angle_subsets)):
            for len_oop in range(0, len(oop_subsets[x_y])):
                for len_lin_angles in range(0, len(lin_angle_subsets)):
                    for len_dihedrals in range(0, len(dihedral_subsets[x_y])):
                        ic_dict[k] = {
                                "bonds" : bonds,
                                "angles" : angle_subsets[len_angles],
                                "linear valence angles" : linear_angles[len_lin_angles] ,
                                "out of plane angles" : oop_subsets[x_y][len_oop],
                                "dihedrals" : dihedral_subsets[x_y][len_dihedrals] }
                        k +=1

    return ic_dict

'''''
GENERAL SYSTEMS
'''''
#TODO: planar subunits

def general_acyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, a_1, specification):

    # set length of subsets
    n_r = num_bonds
    n_phi = 4*num_bonds - 3*num_atoms + a_1
    n_tau = num_bonds - a_1

    symmetric_angles = icsel.get_symm_angles(angles,specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)
    if len(angle_subsets) == 0:
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset)) 
        for subset in itertools.combinations(angles, n_phi+1):
            angle_subsets.append(list(subset))
        for subset in itertools.combinations(angles, n_phi+2):
            angle_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)
    dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau)


    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))
        for subset in itertools.combinations(dihedrals, n_tau+1):
            dihedral_subsets.append(list(subset))
    
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
    return ic_dict 

def general_cyclic_nolinunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, num_of_red, a_1, specification):

    # set length of subsets
    n_r = num_bonds
    
    a_and_b = []
    for a in range(num_of_red):
        b = num_of_red - a
        if a >= 0 and b >= 0: 
            a_and_b.append([a,b])
    
    n_phi = []
    n_tau = []

    for pair in a_and_b:
        n_phi.append(4*num_bonds - 3*num_atoms + a_1 - pair[0])
        n_tau.append(num_bonds - a_1 - pair[1])

    symmetric_angles = icsel.get_symm_angles(angles,specification) 
    angle_subsets = []
    for i in range(len(n_phi)):
        angle_subsets.append(list(
            icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi[i])))

    for i in range(len(n_phi)):
        if len(angle_subsets) == 0:
            for subset in itertools.combinations(angles, n_phi[i]):
                angle_subsets.append(list(subset)) 
            for subset in itertools.combinations(angles, n_phi[i]+1):
                angle_subsets.append(list(subset))
            for subset in itertools.combinations(angles, n_phi[i]+2):
                angle_subsets.append(list(subset))

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)
    dihedral_subsets = []
    for i in range(len(n_tau)):
        dihedral_subsets.append(list(
            icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau[i])))


    # special case where symmetry of dihedrals needs to be broken
    for i in range(len(n_tau)):
        if n_tau[i] != 0 and len(dihedral_subsets) == 0:
                logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
                for subset in itertools.combinations(dihedrals, n_tau[i]):
                    dihedral_subsets.append(list(subset))
                for subset in itertools.combinations(dihedrals, n_tau[i]+1):
                    dihedral_subsets.append(list(subset))
    
    k = 0
    for a_b in range(0, len(n_tau)):
        for len_angles in range(0, len(angle_subsets[a_b])):
            for len_dihedrals in range(0, len(dihedral_subsets[a_b])):
                ic_dict[k] = {
                        "bonds" : bonds,
                        "angles" : angle_subsets[a_b][len_angles],
                        "linear valence angles" : [] ,
                        "out of plane angles" : [],
                        "dihedrals" : dihedral_subsets[a_b][len_dihedrals] }
                k +=1
    return ic_dict



def general_acyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, a_1,l, specification):

    # set length of subsets
    n_r = num_bonds
    n_phi = 4*num_bonds - 3*num_atoms + a_1
    n_phi_prime = l-1
    n_tau = num_bonds - a_1 - (l-1)

    symmetric_angles = icsel.get_symm_angles(angles,specification)
    angle_subsets = icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi)

    if len(angle_subsets) == 0:
        for subset in itertools.combinations(angles, n_phi):
            angle_subsets.append(list(subset)) 
        for subset in itertools.combinations(angles, n_phi+1):
            angle_subsets.append(list(subset))
        for subset in itertools.combinations(angles, n_phi+2):
            angle_subsets.append(list(subset))
    
    symmetric_lin_angles = icsel.get_symm_angles(linear_angles,specification)
    lin_angle_subsets = icsel.get_angle_subsets(symmetric_lin_angles, len(bonds), len(angles),idof,n_phi_prime)

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)
    dihedral_subsets = icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau)
    
    # special case where symmetry of dihedrals needs to be broken
    if n_tau != 0 and len(dihedral_subsets) == 0:
        logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
        for subset in itertools.combinations(dihedrals, n_tau):
            dihedral_subsets.append(list(subset))
        for subset in itertools.combinations(dihedrals, n_tau+1):
            dihedral_subsets.append(list(subset))


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
    
    return ic_dict 


def general_cyclic_linunit_molecule(ic_dict, idof, bonds, angles, linear_angles, out_of_plane, dihedrals, num_bonds, num_atoms, num_of_red, a_1,l, specification):

    # set length of subsets
    n_r = num_bonds
    
    a_and_b = []
    for a in range(num_of_red):
        b = num_of_red - a
        if a >= 0 and b >= 0: 
            a_and_b.append([a,b])
    
    n_phi_prime = l-1
    n_phi = []
    n_tau = []

    for pair in a_and_b:
        n_phi.append(4*num_bonds - 3*num_atoms + a_1 - pair[0])
        n_tau.append(num_bonds - a_1 - pair[1] - (l-1))

    symmetric_angles = icsel.get_symm_angles(angles,specification) 
    angle_subsets = []
    for i in range(len(n_phi)):
        angle_subsets.append(list(
            icsel.get_angle_subsets(symmetric_angles, len(bonds), len(angles),idof,n_phi[i])))

    for i in range(len(n_phi)):
        if len(angle_subsets) == 0:
            for subset in itertools.combinations(angles, n_phi[i]):
                angle_subsets.append(list(subset)) 
            for subset in itertools.combinations(angles, n_phi[i]+1):
                angle_subsets.append(list(subset))
            for subset in itertools.combinations(angles, n_phi[i]+2):
                angle_subsets.append(list(subset))

    symmetric_lin_angles = icsel.get_symm_angles(linear_angles,specification)
    lin_angle_subsets = icsel.get_angle_subsets(symmetric_lin_angles, len(bonds), len(angles),idof,n_phi_prime)

    symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals,specification)
    dihedral_subsets = []
    for i in range(len(n_tau)):
        dihedral_subsets.append(list(
            icsel.get_dihedral_subsets(symmetric_dihedrals, len(bonds), len(angles),idof,n_tau[i])))


    # special case where symmetry of dihedrals needs to be broken
    for i in range(len(n_tau)):
        if n_tau[i] != 0 and len(dihedral_subsets) == 0:
                logging.warning("In order to gain dihedral subsets, symmetry needs to be broken!")
                for subset in itertools.combinations(dihedrals, n_tau[i]):
                    dihedral_subsets.append(list(subset))
                for subset in itertools.combinations(dihedrals, n_tau[i]+1):
                    dihedral_subsets.append(list(subset))
    
    k = 0
    for a_b in range(0, len(n_tau)):
        for len_angles in range(0, len(angle_subsets[a_b])):
            for len_dihedrals in range(0, len(dihedral_subsets[a_b])):
                for len_lin_angles in range(0, len(lin_angle_subsets)):
                    ic_dict[k] = {
                            "bonds" : bonds,
                            "angles" : angle_subsets[a_b][len_angles],
                            "linear valence angles" : lin_angle_subsets[len_lin_angles] ,
                            "out of plane angles" : [],
                            "dihedrals" : dihedral_subsets[a_b][len_dihedrals] }
                    k +=1
    return ic_dict
