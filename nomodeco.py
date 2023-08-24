from __future__ import annotations

import itertools
from typing import NamedTuple
from typing import Iterable
from sklearn.preprocessing import normalize
import pprint
import string
import os
import numpy as np
import pandas as pd
import argparse
import logging
import time
import pymatgen.core as mg
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from mendeleev.fetch import fetch_table


import icgen
import icsel
import bmatrix
import logfile
import molpro_parser
import specifications
import icset_opt

class Atom(NamedTuple):
    symbol: str
    coordinates: tuple

#TODO: isotopes for all elements with command line input
def get_mass_information() -> pd.DataFrame:
    df = fetch_table('elements')
    mass_info = df.loc[:, ['symbol', 'atomic_weight']]
    deuterium_info = pd.DataFrame({'symbol': ['D'], 'atomic_weight': [2.014102]})
    mass_info = pd.concat([mass_info, deuterium_info])
    mass_info.set_index("symbol", inplace=True)
    return mass_info

def reciprocal_square_massvector(atoms):
    n_atoms = len(atoms)
    diag_reciprocal_square = np.zeros(3*n_atoms)
    MASS_INFO = get_mass_information()
    for i in range(0, n_atoms):
        diag_reciprocal_square[3*i:3*i+3] = 1/np.sqrt(MASS_INFO.loc[atoms[i].symbol.strip(string.digits)])
    return diag_reciprocal_square

def reciprocal_massvector(atoms):
    n_atoms = len(atoms)
    diag_reciprocal = np.zeros(3*n_atoms)
    MASS_INFO = get_mass_information()
    for i in range(0, n_atoms):
        diag_reciprocal[3*i:3*i+3] = 1/(MASS_INFO.loc[atoms[i].symbol.strip(string.digits)])
    return diag_reciprocal

def strip_numbers(string):
    return ''.join([char for char in string if not char.isdigit()])

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    args = parser.parse_args()
    return args


start_time = time.time()

def main():
    # Reading Cartesian Coordinates and Hessian
    args = get_args()
    with open(args.output) as inputfile:
        atoms=molpro_parser.parse_xyz_from_inputfile(inputfile)
        n_atoms = len(atoms) 
    with open(args.output) as inputfile:
        CartesianF_Matrix = molpro_parser.parse_Cartesian_F_Matrix_from_inputfile(inputfile) 
        outputfile = logfile.create_new_filename(inputfile.name)


    # initialize log file
    #TODO: extract
    if os.path.exists(outputfile):
        i = 1
        while True:
            new_outputfile_name = f"{outputfile}_{i}"
            if not os.path.exists(new_outputfile_name):
                os.rename(outputfile, new_outputfile_name)
                break
            i +=1
   
    logging.basicConfig(filename=outputfile, filemode='a', format='%(message)s', level=logging.DEBUG)
    logfile.write_logfile_header() 
    

    # Determining molecular symmetry
    molecule = mg.Molecule([strip_numbers(atom.symbol) for atom in atoms], [atom.coordinates for atom in atoms])
    molecule_pg = PointGroupAnalyzer(molecule)
    point_group_sch = molecule_pg.sch_symbol

    # Generation of all possible bonding and bending internal coordinates
    bonds = icgen.initialize_bonds(atoms)
    angles, linear_angles = icgen.initialize_angles(atoms)

    # Setting specifications for calculation: check if molecule is linear, planar or a general molecule
    specification = dict()
    specification = specifications.calculation_specification(specification, atoms, molecule_pg, bonds, angles, linear_angles)

    # Generation of all possible out-of-plane motions

    if specification["planar"] == "yes":
        out_of_plane = icgen.initialize_oop(atoms)
    elif specification["planar"] == "no" and specification["planar submolecule(s)"] == []:
        out_of_plane = []
    elif specification["planar"] == "no" and not (specification["planar submolecule(s)"] == []):
        out_of_plane = icgen.initialize_oop_planar_subunits(atoms, specification["planar submolecule(s)"])
    else:
        return logging.error("Classification of whether topology is planar or not could not be determined!")
    dihedrals = icgen.initialize_dihedrals(atoms)

    # determine internal degrees of freedom 
    idof = 0
    if specification["linearity"] == "fully linear":
        idof = 3*n_atoms-5
    else:
        idof = 3*n_atoms-6

   # update log file

    logfile.write_logfile_oop_treatment(specification["planar"], specification["planar submolecule(s)"])
    logfile.write_logfile_symmetry_treatment(specification, point_group_sch)

    # Computation of the diagonal mass matrices with 
    # the reciprocal and square root reciprocal masses
    diag_reciprocal_square = reciprocal_square_massvector(atoms)
    reciprocal_square_massmatrix = np.diag(diag_reciprocal_square)
    diag_reciprocal = reciprocal_massvector(atoms)
    reciprocal_massmatrix = np.diag(diag_reciprocal)
    
    # Determination of the Normal Modes and eigenvalues 
    # via the diagonalization of the mass-weighted Cartesian F Matrix
    
    Mass_weighted_CartesianF_Matrix = np.transpose(reciprocal_square_massmatrix) @ CartesianF_Matrix @ reciprocal_square_massmatrix
    
    Cartesian_eigenvalues, L = np.linalg.eigh(Mass_weighted_CartesianF_Matrix)
    #print("Cartesian_eigenvalues (EV from mw hessian):", Cartesian_eigenvalues)

    # Determination of the normal modes of zero and low Frequencies

    rottra = L[:,0:(3*n_atoms-idof)]
    
    logfile.write_logfile_generated_IC(bonds, angles, linear_angles, out_of_plane, dihedrals, idof)

    ic_dict = icsel.get_sets(idof,atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, specification)

    optimal_set = icset_opt.find_optimal_coordinate_set(ic_dict, idof, reciprocal_massmatrix, reciprocal_square_massmatrix, rottra, CartesianF_Matrix, atoms, L)

    '''''
    Final calculation with optimal set
    '''''
    bonds = optimal_set["bonds"]
    angles = optimal_set["angles"]
    linear_angles = optimal_set["linear valence angles"]
    out_of_plane = optimal_set["out of plane angles"]
    dihedrals = optimal_set["dihedrals"]

    n_internals = len(bonds) + len(angles) + len(linear_angles) + len(out_of_plane) + len(dihedrals)
    red = n_internals - idof
    
    # Augmenting the B-Matrix with rottra, calculating 
    # and printing the final B-Matrix

    B = np.concatenate((bmatrix.b_matrix(atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof),
                        np.transpose(rottra)),axis=0)

    # Calculating the G-Matrix

    G = B @ reciprocal_massmatrix @ np.transpose(B)
    e,K = np.linalg.eigh(G)

    # Sorting eigenvalues and eigenvectors (just for the case)
    # Sorting highest eigenvalue/eigenvector to lowest!

    idx = e.argsort()[::-1]   
    e = e[idx]
    K = K[:,idx]

    # if redundancies are present, then approximate the inverse of the G-Matrix
    if red > 0:
        K = np.delete(K, -red, axis=1)
        e = np.delete(e, -red, axis=0)

    e = np.diag(e)
    try:
        G_inv = K @ np.linalg.inv(e) @ np.transpose(K)
    except np.linalg.LinAlgError:
        G_inv = K @ np.linalg.pinv(e) @ np.transpose(K)

    # Calculating the inverse augmented B-Matrix

    B_inv = reciprocal_massmatrix @ np.transpose(B) @ G_inv
    InternalF_Matrix = np.transpose(B_inv) @ CartesianF_Matrix @ B_inv

    logfile.write_logfile_information_results(n_internals, red, bonds, angles, 
            linear_angles, out_of_plane, dihedrals)

    ''''' 
    --------------------------- Main-Calculation ------------------------------
    ''''' 

    # Calculation of the mass-weighted normal modes in Cartesian Coordinates

    l = reciprocal_square_massmatrix @ L

    # Calculation of the mass-weighted normal modes in Internal Coordinates

    D = B @ l
 
    # Calculation of the Vibrational Density Matrices / PED, KED and TED matrices
    
    eigenvalues = np.transpose(D) @ InternalF_Matrix @ D
    eigenvalues = np.diag(eigenvalues)
    # print("eigenvalues (from IC space):", eigenvalues)

    num_rottra = 3*n_atoms - idof

    ''''' 
    ------------------------------- Results --------------------------------------
    ''''' 
    
    P = np.zeros((n_internals-red,n_internals+num_rottra,n_internals+num_rottra))
    T = np.zeros((n_internals-red,n_internals+num_rottra,n_internals+num_rottra)) 
    E = np.zeros((n_internals-red,n_internals+num_rottra,n_internals+num_rottra))

    for i in range(0,n_internals-red):
        for m in range(0,n_internals + num_rottra):
            for n in range(0,n_internals + num_rottra):
                k = i + num_rottra
                P[i][m][n] = D[m][k]*InternalF_Matrix[m][n]*D[n][k] / eigenvalues[k] #PED
                T[i][m][n] = D[m][k]*G_inv[m][n]*D[n][k]  #KED
                E[i][m][n] = 0.5 *(T[i][m][n] + P[i][m][n]) #TED

    # check normalization
    sum_check_PED = np.zeros(n_internals)
    sum_check_KED = np.zeros(n_internals)
    sum_check_TED = np.zeros(n_internals)
    for i in range(0, n_internals - red):
        for m in range(0, n_internals + num_rottra):
            for n in range(0, n_internals + num_rottra):
                sum_check_PED[i] += P[i][m][n] 
                sum_check_KED[i] += T[i][m][n] 
                sum_check_TED[i] += E[i][m][n] 

    # Summarized vibrational energy distribution matrix - can be calculated by either PED/KED/TED
    # rows are ICs, columns are harmonic frequencies!
    sum_check_VED = 0
    ved_matrix = np.zeros((n_internals - red, n_internals + num_rottra))
    for i in range(0,n_internals-red):
        for m in range(0, n_internals + num_rottra):
            for n in range (0,n_internals + num_rottra):
                ved_matrix[i][m] += P[i][m][n]
            sum_check_VED += ved_matrix[i][m]
    
    sum_check_VED = np.around(sum_check_VED / (n_internals-red), 2)
    
    # currently: rows are harmonic modes and columns are ICs ==> need to transpose
    ved_matrix = np.transpose(ved_matrix)
    
    # remove the rottra
    ved_matrix = ved_matrix[0:n_internals, 0:n_internals]

    # compute contribution table
    #TODO: normalize if negativ valuea are present?
    contribution_table = normalize(ved_matrix, axis=0, norm='l1') * 100

    # compute intrinsic frequencies
    nu = np.zeros(n_internals) 
    for n in range(0,n_internals):
        for m in range(0,n_internals):
            for i in range(0,n_internals-red):
                k = i + num_rottra
                nu[n] += D[m][k] * InternalF_Matrix[m][n] * D[n][k]
     
    nu_final = np.sqrt(nu) *  5140.4981

    normal_coord_harmonic_frequencies = np.sqrt(eigenvalues[(3*n_atoms-idof):3*n_atoms]) * 5140.4981
    normal_coord_harmonic_frequencies = np.around(normal_coord_harmonic_frequencies, decimals=2)
    normal_coord_harmonic_frequencies_string = normal_coord_harmonic_frequencies.astype('str')

    all_internals = bonds + angles + linear_angles + out_of_plane + dihedrals

    Results = pd.DataFrame()
    Results['Internal Coordinate'] = all_internals
    Results['Intrinsic Frequencies'] = pd.DataFrame(nu_final).applymap("{0:.2f}".format)
    Results = Results.join(pd.DataFrame(ved_matrix).applymap("{0:.2f}".format))

    ContributionTable = pd.DataFrame()
    ContributionTable['Internal Coordinate'] = all_internals
    ContributionTable['Intrinsic Frequencies'] = pd.DataFrame(nu_final).applymap("{0:.2f}".format)
    ContributionTable = ContributionTable.join(pd.DataFrame(contribution_table).applymap("{0:.2f}".format))
    
    columns = {}
    keys = range(3*n_atoms-((3*n_atoms-idof)))
    for i in keys:
        columns[i] = normal_coord_harmonic_frequencies_string[i]

    Results = Results.rename(columns=columns)
    ContributionTable = ContributionTable.rename(columns=columns)

    logfile.write_logfile_results(Results, ContributionTable, sum_check_VED)

    # here the individual matrices can be computed, one can comment them out
    # if not needed

    all_internals_string = []
    for internal in all_internals:
        all_internals_string.append('(' + ', '.join(internal) + ')')
    
    columns = {}
    keys = range(n_internals)
    for i in keys:
        columns[i] = all_internals_string[i]

    for mode in range(0, len(normal_coord_harmonic_frequencies)):

        PED = pd.DataFrame()
        KED = pd.DataFrame()
        TED = pd.DataFrame()

        PED['Internal Coordinate'] = all_internals
        KED['Internal Coordinate'] = all_internals
        TED['Internal Coordinate'] = all_internals
        PED = PED.join(pd.DataFrame(P[mode][0:n_internals, 0:n_internals]).applymap("{0:.2f}".format))
        KED = KED.join(pd.DataFrame(T[mode][0:n_internals, 0:n_internals]).applymap("{0:.2f}".format))
        TED = TED.join(pd.DataFrame(E[mode][0:n_internals, 0:n_internals]).applymap("{0:.2f}".format))
        PED = PED.rename(columns=columns)
        KED = KED.rename(columns=columns)
        TED = TED.rename(columns=columns)
        
        logfile.write_logfile_extended_results(PED,KED,TED, sum_check_PED[mode], sum_check_KED[mode], sum_check_TED[mode], normal_coord_harmonic_frequencies[mode])

    logfile.call_shutdown()
    
    print("Runtime: %s seconds" % (time.time() - start_time))

if __name__ == '__main__':
    main()


