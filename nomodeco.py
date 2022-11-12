from __future__ import annotations

import itertools
from typing import NamedTuple
from typing import Iterable

import string
import os
import numpy as np
import pandas as pd
import argparse
import logging
import time
from mendeleev.fetch import fetch_table


import icgen
import icsel
import bmatrix
import logfile
import parser
import userinterface

class Atom(NamedTuple):
    symbol: str
    coordinates: tuple


def get_mass_information() -> pd.DataFrame:
    df = fetch_table('elements')
    mass_info = df.loc[:, ['symbol', 'atomic_weight']]
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

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("oop_directive")
    args = parser.parse_args()
    return args

def oop_directive() -> string:
    args = get_args()
    if args.oop_directive == "oop":
        return "oop"
    elif args.oop_directive == "no-oop":
        return "no-oop"
    else:
        return "no specification"

start_time = time.time()

def main():
    args = get_args()
    with open(args.output) as inputfile:
        atoms=parser.parse_xyz_from_inputfile(inputfile)
        n_atoms = len(atoms) 
    with open(args.output) as inputfile:
        CartesianF_Matrix = parser.parse_Cartesian_F_Matrix_from_inputfile(inputfile) 
        outputfile = logfile.create_new_filename(inputfile.name)
    
    # initialize log file
    if os.path.exists(outputfile):
        os.remove(outputfile)
    logging.basicConfig(filename=outputfile, filemode='a', format='%(message)s', level=logging.DEBUG)
    logfile.write_logfile_header()
    logfile.write_logfile_oop_treatment()
    
    # Generation of all possible internal coordinates
    bonds = icgen.initialize_bonds(atoms)
    angles, linear_angles = icgen.initialize_angles(atoms)
    if oop_directive() == "oop":
        out_of_plane = icgen.initialize_oop(atoms)
    elif oop_directive() == "no-oop":
        out_of_plane = []
    else:
        return logging.error("You need to specify the oop or no-oop directive!")
    dihedrals = icgen.initialize_dihedrals(atoms)

    # determine internal degrees of freedom 
    idof = 0
    if (linear_angles and not angles) or (not linear_angles and not angles):
        idof = 3*n_atoms-5
    else:
        idof = 3*n_atoms-6
    
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
    
    # Determination of the normal modes of zero and low Frequencies

    rottra = L[:,0:(3*n_atoms-idof)]
    
    logfile.write_logfile_generated_IC(bonds, angles, linear_angles, out_of_plane, dihedrals, idof)
    userinterface.write_generated_IC(bonds, angles, linear_angles, out_of_plane, dihedrals, idof)

    ic_dict = userinterface.generate_set(n_atoms, idof, bonds, angles, linear_angles, out_of_plane, dihedrals)

    bonds = ic_dict["bonds"]
    angles = ic_dict["angles"]
    linear_angles = ic_dict["linear valence angles"]
    out_of_plane = ic_dict["out of plane angles"]
    dihedrals = ic_dict["dihedrals"]

    n_internals = len(bonds) + len(angles) + len(linear_angles) + len(out_of_plane) + len(dihedrals)
    red = n_internals - idof
    
    # Augmenting the B-Matrix with rottra, calculating 
    # and printing the final B-Matrix

    B = np.concatenate((bmatrix.b_matrix(atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof),
                        np.transpose(rottra)),axis=0)

    # Calculating the G-Matrix

    G = B @ reciprocal_massmatrix @ np.transpose(B)
    G_int = bmatrix.b_matrix(atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof) @ reciprocal_massmatrix @ np.transpose(bmatrix.b_matrix(atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof))
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
    
    #TODO: e inv does not exist - make error handling
    e = np.diag(e)
    G_inv = K @ np.linalg.inv(e) @ np.transpose(K)

    # Calculating the inverse augmented B-Matrix

    B_inv = reciprocal_massmatrix @ np.transpose(B) @ G_inv
    InternalF_Matrix = np.transpose(B_inv) @ CartesianF_Matrix @ B_inv

    logfile.write_logfile_information_results(B, B_inv, CartesianF_Matrix, InternalF_Matrix, n_internals, red, bonds, 
    angles, linear_angles, out_of_plane, dihedrals)
        
    ''''' 
    --------------------------- Main-Calculation ------------------------------
    ''''' 

    # Calculation of the mass-weighted normal modes in Cartesian Coordinates

    l = reciprocal_square_massmatrix @ L

    # Calculation of the mass-weighted normal modes in Internal Coordinates

    D = B @ l

    # Calculation of the Vibrational Density Matrices
    
    eigenvalues = np.transpose(D) @ InternalF_Matrix @ D
    eigenvalues = np.diag(eigenvalues)

    P = np.zeros((n_internals-red,n_internals,n_internals))

    for i in range(0,n_internals-red):
        for m in range(0,n_internals):
            for n in range(0,n_internals):
                k = i + (3*n_atoms-idof)
                P[i][m][n] = D[m][k]*InternalF_Matrix[m][n]*D[n][k] / eigenvalues[k]
                
    ''''' 
    ------------------------------- Results --------------------------------------
    ''''' 

    # Results part 1 
    Diag_elements = np.zeros((n_internals-red,n_internals))
    for i in range(0,n_internals-red):
        for n in range (0,n_internals):
            Diag_elements[i][n] = np.diag(P[i])[n]

    Diag_elements = np.transpose(Diag_elements)


    nu = np.zeros(n_internals) 
    for n in range(0,n_internals):
        for m in range(0,n_internals):
            for i in range(0,n_internals-red):
                k = i + (3*n_atoms-idof)
                nu[n] += D[m][k] * InternalF_Matrix[m][n] * D[n][k]
                
    if np.any(nu < 0) == True:
        logfile.write_logfile_nan_freq()

    nu_final = np.sqrt(nu) *  5140.4981

    normal_coord_harmonic_frequencies = np.sqrt(eigenvalues[(3*n_atoms-idof):3*n_atoms]) * 5140.4981
    normal_coord_harmonic_frequencies = np.around(normal_coord_harmonic_frequencies, decimals=2)
    normal_coord_harmonic_frequencies_string = normal_coord_harmonic_frequencies.astype('str')

    all_internals = bonds + angles + linear_angles + out_of_plane + dihedrals

    Results1 = pd.DataFrame()
    Results1['Internal Coordinate'] = all_internals
    Results1['Intrinsic Frequencies'] = pd.DataFrame(nu_final).applymap("{0:.2f}".format)
    Results1 = Results1.join(pd.DataFrame(Diag_elements).applymap("{0:.2f}".format))

    columns = {}
    keys = range(3*n_atoms-((3*n_atoms-idof)))
    for i in keys:
        columns[i] = normal_coord_harmonic_frequencies_string[i]

    Results1 = Results1.rename(columns=columns)
    
    # Results part 2
    
    nu_perNormalCoordinate = np.zeros((n_internals,n_internals-red)) 
    for n in range(0,n_internals):
        for i in range(0,n_internals-red):
            for m in range(0,n_internals):
                k = i + (3*n_atoms-idof)
                nu_perNormalCoordinate[n][i] += D[m][k] * InternalF_Matrix[m][n] * D[n][k]

    Results2 = pd.DataFrame()
    Results2['Internal Coordinate'] = all_internals
    Results2 = Results2.join(pd.DataFrame(nu_perNormalCoordinate).applymap("{0:.2f}".format))
    Results2 = Results2.rename(columns=columns)

    sum_array1 = np.zeros(n_internals)

    for n in range(0,n_internals):
        for i in range(0, n_internals-red):
            sum_array1[i] += Diag_elements[n][i]

    Contribution_Matrix1 = np.zeros((n_internals, n_internals-red))
    for i in range(0, n_internals-red):
        Contribution_Matrix1[:,i] =((Diag_elements[:,i] / sum_array1[i]) * 100).astype(float)
    Contribution_Table1 = pd.DataFrame()
    Contribution_Table1['Internal Coordinate'] = all_internals
    Contribution_Table1 = Contribution_Table1.join(pd.DataFrame(Contribution_Matrix1).applymap("{0:.2f}".format))
    Contribution_Table1 = Contribution_Table1.rename(columns=columns)

    sum_array2 = np.zeros(n_internals)

    for n in range(0,n_internals):
        for i in range(0,n_internals-red):
            sum_array2[i] += nu_perNormalCoordinate[n][i]

    Contribution_Matrix2 = np.zeros((n_internals,n_internals-red))
    for i in range(0, n_internals-red):
        Contribution_Matrix2[:,i] = (nu_perNormalCoordinate[:,i] / sum_array2[i]) * 100

    Contribution_Table2 = pd.DataFrame()
    Contribution_Table2['Internal Coordinate'] = all_internals
    Contribution_Table2 = Contribution_Table2.join(pd.DataFrame(Contribution_Matrix2).applymap("{0:.2f}".format))
    Contribution_Table2 = Contribution_Table2.rename(columns=columns)
    
    logfile.write_logfile_results(Results1, Results2, Contribution_Table1, Contribution_Table2,G_int)

print("Runtime: %s seconds" % (time.time() - start_time))

if __name__ == '__main__':
    main()
