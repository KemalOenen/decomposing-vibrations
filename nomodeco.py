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

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("specifications")
    args = parser.parse_args()
    return args

def calculation_specification():
    args = get_args()
    specification = dict()
    with open(args.specifications) as inputfile:
        for line in inputfile:
            if line.strip().startswith('Out-of-plane angles specification'):
                specification = {"out-of-plane_treatment" : next(inputfile, '').strip()}
    return specification

start_time = time.time()

def main():
    # Reading Cartesian Coordinates and Hessian
    args = get_args()
    with open(args.output) as inputfile:
        atoms=parser.parse_xyz_from_inputfile(inputfile)
        n_atoms = len(atoms) 
    with open(args.output) as inputfile:
        CartesianF_Matrix = parser.parse_Cartesian_F_Matrix_from_inputfile(inputfile) 
        outputfile = inputfile.name + "_generalized_internal_coordinates.log"
    
    # Setting specifications for calculation
    specification = calculation_specification()


    # initialize log file
    if os.path.exists(outputfile):
        os.remove(outputfile)
    logging.basicConfig(filename=outputfile, filemode='a', format='%(message)s', level=logging.DEBUG)
    logfile.write_logfile_oop_treatment(specification["out-of-plane_treatment"])
    

    # Generation of all possible internal coordinates
    # For delocalized internal coordinates the B-Matrix is formed for all primitive internal coordinates
    bonds = icgen.initialize_bonds(atoms)
    angles, linear_angles = icgen.initialize_angles(atoms)
    if specification["out-of-plane_treatment"] == "oop":
        out_of_plane = icgen.initialize_oop(atoms)
    elif specification["out-of-plane_treatment"] == "no-oop":
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
    
    logfile.write_logfile_generated_IC(bonds, angles, linear_angles, out_of_plane, dihedrals)

    # Computation of the diagonal mass matrices with 
    # the square root reciprocal masses
    diag_reciprocal_square = reciprocal_square_massvector(atoms)
    reciprocal_square_massmatrix = np.diag(diag_reciprocal_square)
    
    # Determination of the Normal Modes and eigenvalues 
    # via the diagonalization of the mass-weighted Cartesian F Matrix
    
    Mass_weighted_CartesianF_Matrix = np.transpose(reciprocal_square_massmatrix) @ CartesianF_Matrix @ reciprocal_square_massmatrix
    
    Cartesian_eigenvalues, L = np.linalg.eigh(Mass_weighted_CartesianF_Matrix)
    
    # Determination of the normal modes of zero and low Frequencies

    rottra = L[:,0:(3*n_atoms-idof)]
    
    # Augmenting the B-Matrix with rottra, calculating 
    # and printing the final B-Matrix

    B = np.concatenate((bmatrix.b_matrix(atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof),
                            np.transpose(rottra)),axis=0)

    """""  Generation of delocalized internal coordinates """""

    G = B @ np.transpose(B)
    full_eigenvalues, full_eigenvectors = np.linalg.eig(G)
   
   # sorting eigenvalues and eigenvectors to seperate redundant and non-redundant space

    idx = full_eigenvalues.argsort()[::-1]
    full_eigenvalues = full_eigenvalues[idx]
    full_eigenvectors = full_eigenvectors[:,idx]

    U = []
    R = []
    
    i = 0
    
    for eigenvalue in full_eigenvalues:
        if np.isclose(eigenvalue, 0):
            R.append(full_eigenvectors[i,:])
        else:
            U.append(full_eigenvectors[i,:])
        i += 1

    R = np.transpose(R)
    U = np.transpose(U)

    logfile.write_logfile_results(full_eigenvalues, full_eigenvectors, U, R)



    

    print("Runtime: %s seconds" % (time.time() - start_time))

if __name__ == '__main__':
    main()
