#general packages
from __future__ import annotations

from typing import NamedTuple
import string
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import time
import pymatgen.core as mg
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from mendeleev.fetch import fetch_table

# for heatmap
mpl.rcParams['backend'] = "Agg"
# do not show messages
logging.getLogger("matplotlib").setLevel(logging.ERROR)

# nomodeco libraries
from modules import icgen
from modules import icsel
from modules import bmatrix
from modules import logfile
from modules import molpro_parser
from modules import specifications
from modules import icset_opt
from modules import arguments
from modules import alt_icgen
from modules import specifications as sp
from modules import dfs_connected
from modules import topology as tp
from modules.nomodeco_classes import Molecule, InternalCoordinates
from modules import icset_opt_multproc

# TODO: isotopes for all elements with command line input?
def get_mass_information() -> pd.DataFrame:
    df = fetch_table('elements')
    mass_info = df.loc[:, ['symbol', 'atomic_weight']]
    deuterium_info = pd.DataFrame({'symbol': ['D'], 'atomic_weight': [2.014102]})
    mass_info = pd.concat([mass_info, deuterium_info])
    mass_info.set_index("symbol", inplace=True)
    return mass_info

def get_bond_information() -> pd.DataFrame:
    df = fetch_table('elements')

    bond_info = df.loc[:, ['symbol', 'covalent_radius_pyykko', 'vdw_radius']]

    bond_info.set_index("symbol", inplace=True)
    bond_info /= 100
    return bond_info

BOND_INFO = get_bond_information()

def reciprocal_square_massvector(atoms):
    n_atoms = len(atoms)
    diag_reciprocal_square = np.zeros(3 * n_atoms)
    MASS_INFO = get_mass_information()
    for i in range(0, n_atoms):
        diag_reciprocal_square[3 * i:3 * i + 3] = 1 / np.sqrt(MASS_INFO.loc[atoms[i].symbol.strip(string.digits)])
    return diag_reciprocal_square


def reciprocal_massvector(atoms):
    n_atoms = len(atoms)
    diag_reciprocal = np.zeros(3 * n_atoms)
    MASS_INFO = get_mass_information()
    for i in range(0, n_atoms):
        diag_reciprocal[3 * i:3 * i + 3] = 1 / (MASS_INFO.loc[atoms[i].symbol.strip(string.digits)])
    return diag_reciprocal


def strip_numbers(string):
    return ''.join([char for char in string if not char.isdigit()])


start_time = time.time() # measure runtime

def main():
    # Reading Cartesian Coordinates and Hessian
    # from user input
    args = arguments.get_args()
    with open(args.output) as inputfile:
        atoms = molpro_parser.parse_xyz_from_inputfile(inputfile)
        n_atoms = len(atoms)
    with open(args.output) as inputfile:
        CartesianF_Matrix = molpro_parser.parse_Cartesian_F_Matrix_from_inputfile(inputfile)
        outputfile = logfile.create_filename_out(inputfile.name)

    # initialize out file
    if os.path.exists(outputfile):
        i = 1
        while True:
            new_outputfile_name = f"{outputfile}_{i}"
            if not os.path.exists(new_outputfile_name):
                os.rename(outputfile, new_outputfile_name)
                break
            i += 1
    out = logfile.setup_logger('outputfile', outputfile)
    logfile.write_logfile_header(out)

    # BUGFIX for Chloride which for some reason gets written as CL by molpro
    updated_atoms = [Atom(symbol='Cl', coordinates=atom.coordinates) if atom.symbol == 'CL' else atom for atom in atoms]
    atoms = updated_atoms
    
    # Determining molecular symmetry
    molecule = mg.Molecule([strip_numbers(atom.symbol) for atom in atoms], [atom.coordinates for atom in atoms])
    molecule_pg = PointGroupAnalyzer(molecule)
    point_group_sch = molecule_pg.sch_symbol   
    

    atoms = Molecule(atoms)
    
    '''
    IC Generation
    '''
    # Intermolecular Bonds:
    # Use Degree of Covalance https://doi.org/10.1002/qua.21049 for hydrogen bond detection
    degofc_table = atoms.degree_of_covalance()    
    cov_bonds = atoms.covalent_bonds(degofc_table)
    sp.covalent_bonds = cov_bonds   
    
    # Detect and generate covalent submolecules
    _, _, cov_submolecules_symbols = atoms.detect_submolecules()

    # Generate Hydrogen Bond and the Acceptor-Donor Coordinate
    
    h_bonds = atoms.intermolecular_h_bond(degofc_table,cov_submolecules_symbols)
    acc_don_bonds = atoms.intermolecular_acceptor_donor(degofc_table,cov_submolecules_symbols)
    
            
    Total_IC_dict = InternalCoordinates()

    #Add covalent ICs
    Total_IC_dict.add_coordinate("cov_bond", cov_bonds)
    Total_IC_dict.add_coordinate("h_bond", h_bonds)
    Total_IC_dict.add_coordinate("acc_don", acc_don_bonds)
    Total_IC_dict.add_coordinate("cov_angles", atoms.generate_angles(cov_bonds)[0])
    Total_IC_dict.add_coordinate("cov_linear_angles", atoms.generate_angles(cov_bonds)[1])
    Total_IC_dict.add_coordinate("cov_dihedrals", atoms.generate_dihedrals(cov_bonds))
    Total_IC_dict.add_coordinate("cov_oop", atoms.generate_out_of_plane(cov_bonds))
     


    # Define Total Bonds
    bonds = list(set(cov_bonds).union(set(h_bonds)))
    bond_acc_don = list(set(cov_bonds).union(set(acc_don_bonds)))
    

    Total_IC_dict.add_coord_diff("h_bond_angles",atoms.generate_angles(bonds)[0],Total_IC_dict["cov_angles"])
    Total_IC_dict.add_coord_diff_linear("h_bond_linear_angles", atoms.generate_angles(bonds)[1],Total_IC_dict["cov_linear_angles"])
    Total_IC_dict.add_coord_diff("h_bond_dihedrals", atoms.generate_dihedrals(bonds),Total_IC_dict["cov_dihedrals"])
    Total_IC_dict.add_coord_diff("h_bond_oop", atoms.generate_out_of_plane(bonds),Total_IC_dict["cov_oop"])
    Total_IC_dict.add_coord_diff("acc_don_angles", atoms.generate_angles(bond_acc_don)[0],Total_IC_dict["cov_angles"])
    Total_IC_dict.add_coord_diff_linear("acc_don_linear_angles",atoms.generate_angles(bond_acc_don)[1],Total_IC_dict["cov_linear_angles"])
    Total_IC_dict.add_coord_diff("acc_don_dihedrals", atoms.generate_dihedrals(bond_acc_don),Total_IC_dict["cov_dihedrals"])
    Total_IC_dict.add_coord_diff("acc_don_oop", atoms.generate_out_of_plane(bond_acc_don), Total_IC_dict["cov_oop"])

    print(Total_IC_dict)  
    # Assing Coordinates from IC Dict to Variables
    angles = Total_IC_dict["cov_angles"] + Total_IC_dict["h_bond_angles"] 
    linear_angles = Total_IC_dict["cov_linear_angles"] + Total_IC_dict["h_bond_linear_angles"]
    dihedrals = Total_IC_dict["cov_dihedrals"] + Total_IC_dict["h_bond_dihedrals"]
    
    # If no valid dihedrals found append acc_don_dihedrals
    if len(dihedrals) == 0:
       dihedrals = Total_IC_dict["cov_dihedrals"] + Total_IC_dict["acc_don_dihedrals"]


    '''
    Generating the Specification for the primary calculation
    '''
   
    # Setting specifications for calculation: check if molecule is linear, planar or a general molecule
    specification = dict()
    specification = specifications.calculation_specification(specification, atoms, molecule_pg, bonds, angles,
                                                                                   linear_angles)
    
    # Generation of all possible out-of-plane motions

    if specification["planar"] == "yes":
         out_of_plane = Total_IC_dict["cov_oop"] + Total_IC_dict["h_bond_oop"]   
    elif specification["planar"] == "no" and specification["planar submolecule(s)"] == []:
        out_of_plane = []
    elif specification["planar"] == "no" and not (specification["planar submolecule(s)"] == []):
        out_of_plane = icgen.initialize_oop_planar_subunits(atoms, specification["planar submolecule(s)"])
        out_of_plane = alt_icgen.initialize_alt_oop_planar_subunits(atoms,bonds,specification["planar submolecule(s)"])
    else:
        return out.error("Classification of whether topology is planar or not could not be determined!")
    
    # determine internal degrees of freedom 
    idof = 0
    if specification["linearity"] == "fully linear":
        idof = atoms.idof_linear()
    else:
        idof = atoms.idof_general()

    # update out file
    
    print(specification)
    logfile.write_logfile_oop_treatment(out, specification["planar"], specification["planar submolecule(s)"])
    logfile.write_logfile_symmetry_treatment(out, specification, point_group_sch)

    '''
    Passing Section for Topology.py
    '''
    tp.atoms_list = atoms
    tp.Total_IC_dict = Total_IC_dict
    
    '''
    Diag Mass Matrix and Reciprocal Square root masses
    '''


    # Computation of the diagonal mass matrices with 
    # the reciprocal and square root reciprocal masses
    diag_reciprocal_square = reciprocal_square_massvector(atoms)
    reciprocal_square_massmatrix = np.diag(diag_reciprocal_square)
    diag_reciprocal = reciprocal_massvector(atoms)
    reciprocal_massmatrix = np.diag(diag_reciprocal)

    # Determination of the Normal Modes and eigenvalues 
    # via the diagonalization of the mass-weighted Cartesian F Matrix

    Mass_weighted_CartesianF_Matrix = np.transpose(
        reciprocal_square_massmatrix) @ CartesianF_Matrix @ reciprocal_square_massmatrix

    Cartesian_eigenvalues, L = np.linalg.eigh(Mass_weighted_CartesianF_Matrix)
    # print("Cartesian_eigenvalues (EV from mw hessian):", Cartesian_eigenvalues)

    # Determination of the normal modes of zero and low Frequencies

    rottra = L[:, 0:(3 * n_atoms - idof)]

    logfile.write_logfile_generated_IC(out, bonds, angles, linear_angles, out_of_plane, dihedrals, idof)

    # get symmetric coordinates
    if args.penalty1 != 0:
        symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
        symmetric_angles = icsel.get_symm_angles(angles, specification)
        symmetric_dihedrals = icsel.get_symm_angles(dihedrals, specification)
        symmetric_coordinates = {**symmetric_bonds, **symmetric_angles, **symmetric_dihedrals}
    else:
        symmetric_coordinates = dict()

    ic_dict = icsel.get_sets(idof, out, atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, specification)
    
    # Here also a multiprocessing variant is implemented to use comment out!

    optimal_set = icset_opt_multproc.find_optimal_coordinate_set(ic_dict, args, idof, reciprocal_massmatrix,
                                                        reciprocal_square_massmatrix, rottra, CartesianF_Matrix, atoms,
                                                        symmetric_coordinates, L, args.penalty1, args.penalty2) 
    #optimal_set = icset_opt.find_optimal_coordinate_set(ic_dict, args, idof, reciprocal_massmatrix,
    #                                                    reciprocal_square_massmatrix, rottra, CartesianF_Matrix, atoms,
    #                                                    symmetric_coordinates, L, args.penalty1, args.penalty2) 
    


    '''
    Baker Description of Internal Coordinates
    '''
    # This whole section is only for testing
    

    import sympy

    B_baker  = bmatrix.b_matrix(atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof)

    # in order to evaluate which internal coordinates are linearly independent we first transpose the b matrix
    # row space --> internal coordinates
    B_baker_transpose = np.transpose(B_baker)
    _, inds = sympy.Matrix(B_baker).T.rref()
    G_baker = np.matmul(B_baker, B_baker_transpose)
    G_eigenvalues, G_eigenvectors = np.linalg.eig(G_baker)
     
    G_eigensum = np.sum(G_eigenvalues)
    G_eigenvalues_norm = G_eigenvalues / G_eigensum
    # In order to find the redundant vectors we specify a tolerance and evaluate the position where eigenvalues occure
    tol = 0.1
    index_arr = []
    for i , eigval in enumerate(G_eigenvalues):
        if eigval < tol:
           index_arr.append(i)
    U_eigenvalues= np.delete(G_eigenvalues, index_arr)
    U_eigenvector = np.delete(G_eigenvectors, index_arr, axis=1)
    # We now normalize this Eigenvalues in order to get the contribution of each internal coordinate
    eigen_sum = np.sum(U_eigenvalues)
    U_eigenvalues_norm = U_eigenvalues / eigen_sum
    D_matrix = np.diag(U_eigenvalues)
     
    weights = np.linalg.lstsq(G_eigenvectors, B_baker_transpose[2], rcond=None)[0]

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
                        np.transpose(rottra)), axis=0)
    
   


    
    # Calculating the G-Matrix

    G = B @ reciprocal_massmatrix @ np.transpose(B)
    e, K = np.linalg.eigh(G)

    # Sorting eigenvalues and eigenvectors (just for the case)
    # Sorting highest eigenvalue/eigenvector to lowest!

    idx = e.argsort()[::-1]
    e = e[idx]
    K = K[:, idx]

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

    logfile.write_logfile_information_results(out, n_internals, red, bonds, angles,
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

    num_rottra = 3 * n_atoms - idof

    ''''' 
    ------------------------------- Results --------------------------------------
    '''''

    P = np.zeros((n_internals - red, n_internals + num_rottra, n_internals + num_rottra))
    T = np.zeros((n_internals - red, n_internals + num_rottra, n_internals + num_rottra))
    E = np.zeros((n_internals - red, n_internals + num_rottra, n_internals + num_rottra))

    for i in range(0, n_internals - red):
        for m in range(0, n_internals + num_rottra):
            for n in range(0, n_internals + num_rottra):
                k = i + num_rottra
                P[i][m][n] = D[m][k] * InternalF_Matrix[m][n] * D[n][k] / eigenvalues[k]  # PED
                T[i][m][n] = D[m][k] * G_inv[m][n] * D[n][k]  # KED
                E[i][m][n] = 0.5 * (T[i][m][n] + P[i][m][n])  # TED

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
    for i in range(0, n_internals - red):
        for m in range(0, n_internals + num_rottra):
            for n in range(0, n_internals + num_rottra):
                ved_matrix[i][m] += P[i][m][n]
            sum_check_VED += ved_matrix[i][m]

    sum_check_VED = np.around(sum_check_VED / (n_internals - red), 2)

    # currently: rows are harmonic modes and columns are ICs ==> need to transpose
    ved_matrix = np.transpose(ved_matrix)

    # remove the rottra
    ved_matrix = ved_matrix[0:n_internals, 0:n_internals]

    # compute diagonal elements of PED matrix

    Diag_elements = np.zeros((n_internals - red, n_internals))
    for i in range(0, n_internals - red):
        for n in range(0, n_internals):
            Diag_elements[i][n] = np.diag(P[i])[n]

    Diag_elements = np.transpose(Diag_elements)

    # compute contribution matrix
    sum_diag = np.zeros(n_internals)

    for n in range(0, n_internals):
        for i in range(0, n_internals - red):
            sum_diag[i] += Diag_elements[n][i]

    contribution_matrix = np.zeros((n_internals, n_internals - red))
    for i in range(0, n_internals - red):
        contribution_matrix[:, i] = ((Diag_elements[:, i] / sum_diag[i]) * 100).astype(float)

    # compute intrinsic frequencies
    nu = np.zeros(n_internals)
    for n in range(0, n_internals):
        for m in range(0, n_internals):
            for i in range(0, n_internals - red):
                k = i + num_rottra
                nu[n] += D[m][k] * InternalF_Matrix[m][n] * D[n][k]

    nu_final = np.sqrt(nu) * 5140.4981

    normal_coord_harmonic_frequencies = np.sqrt(eigenvalues[(3 * n_atoms - idof):3 * n_atoms]) * 5140.4981
    normal_coord_harmonic_frequencies = np.around(normal_coord_harmonic_frequencies).astype(int)
    normal_coord_harmonic_frequencies_string = normal_coord_harmonic_frequencies.astype('str')

    all_internals = bonds + angles + linear_angles + out_of_plane + dihedrals

    all_internals_string = []
    for internal in all_internals:
        all_internals_string.append('(' + ', '.join(internal) + ')')

    Results = pd.DataFrame()
    Results['Internal Coordinate'] = all_internals_string
    Results['Intrinsic Frequencies'] = pd.DataFrame(nu_final).applymap("{0:.2f}".format)
    Results = Results.join(pd.DataFrame(ved_matrix).applymap("{0:.2f}".format))

    DiagonalElementsPED = pd.DataFrame()
    DiagonalElementsPED['Internal Coordinate'] = all_internals_string
    DiagonalElementsPED['Intrinsic Frequencies'] = pd.DataFrame(nu_final).applymap("{0:.2f}".format)
    DiagonalElementsPED = DiagonalElementsPED.join(pd.DataFrame(Diag_elements).applymap("{0:.2f}".format))

    ContributionTable = pd.DataFrame()
    ContributionTable['Internal Coordinate'] = all_internals_string
    ContributionTable['Intrinsic Frequencies'] = pd.DataFrame(nu_final).applymap("{0:.2f}".format)
    ContributionTable = ContributionTable.join(pd.DataFrame(contribution_matrix).applymap("{0:.2f}".format))

    columns = {}
    keys = range(3 * n_atoms - ((3 * n_atoms - idof)))
    for i in keys:
        columns[i] = normal_coord_harmonic_frequencies_string[i]

    Results = Results.rename(columns=columns)
    DiagonalElementsPED = DiagonalElementsPED.rename(columns=columns)
    ContributionTable = ContributionTable.rename(columns=columns)
    ContributionTable_Index = ContributionTable.set_index("Internal Coordinate")
  
    Contribution_Sliced = []
    for freq in ContributionTable_Index.columns[1:]:
        row = {"Frequency" : freq}
        details = []
        for i,coord in enumerate(ContributionTable["Internal Coordinate"]):
            percentage = ContributionTable.loc[i,freq]
            # weird quickfix maybe repair this some time
            if isinstance(percentage, pd.Series):
               percentage = float(percentage[1])
            elif isinstance(percentage, str) and percentage != ".":
               percentage = float(percentage)
            if percentage > 10:
               details.append(f"{coord}: {percentage:.2f}\%")
        row["Details"] = "; ".join(details)
        Contribution_Sliced.append(row)
    annotated_df = pd.DataFrame(Contribution_Sliced)
    
    ContributionTable_T = ContributionTable.transpose()
   
    latex_table = annotated_df.to_latex(column_format="c|c",
                                        index = False,
                                        multicolumn = True,
                                        header=["Frequencies","Contributions"],
                                        escape = False 
                                        )
    # Adjustable table width
    table_width = ""
    latex_table_adjusted = '\\begin{adjustbox}{width=' + str(table_width) + '\\textwidth}\n' + latex_table + '\\end{adjustbox}' 
    print(latex_table_adjusted) 
 
   # new_data = []
   # for freq in ContributionTable_T.columns[2:]:  # Ensuring we're getting the right frequency columns
   #     row = {"Frequency": freq}
   #     details = []
   #     for i, coord in enumerate(ContributionTable_T["Internal Coordinate"]):
   ##         percentage = ContributionTable_T.loc[i, freq]
   #         if percentage > 0:
   #            details.append(f"{coord}: {percentage:.2f}%")
   #     row["Details"] = "; ".join(details)
   #     new_data.append(row)


    # TODO:  line breaks in output file

    logfile.write_logfile_results(out, Results, DiagonalElementsPED, ContributionTable, sum_check_VED)

    # heat map results
    # TODO: clean up
    if args.heatmap:
        columns = {}
        keys = range(3 * n_atoms - ((3 * n_atoms - idof)))
        for i in keys:
            columns[i] = normal_coord_harmonic_frequencies[i]

        for matrix_type in args.heatmap:
            if matrix_type == "ved":
                rows, cols = ved_matrix.shape
                figsize = (cols, rows)
                plt.figure(figsize=figsize)
                
                heatmap_df = pd.DataFrame(ved_matrix).applymap("{0:.2f}".format)
                heatmap_df = heatmap_df.rename(columns=columns)
                heatmap_df = heatmap_df.astype('float')
                heatmap_df.index = all_internals_string

                heatmap = sns.heatmap(heatmap_df, cmap="Blues", annot=True, square=True, annot_kws={"size": 35 / np.sqrt(len(ved_matrix))})
                heatmap.figure.savefig("heatmap_ved_matrix.png", bbox_inches="tight", dpi=500)
                plt.close(heatmap.figure)
            if matrix_type == "diag":
                rows, cols = Diag_elements.shape
                figsize = (cols, rows)
                plt.figure(figsize=figsize)
                
                heatmap_df = pd.DataFrame(Diag_elements).applymap("{0:.2f}".format)
                heatmap_df = heatmap_df.rename(columns=columns)
                heatmap_df = heatmap_df.astype('float')
                heatmap_df.index = all_internals_string

                heatmap = sns.heatmap(heatmap_df, cmap="Blues", annot=True, square=True, annot_kws={"size": 35 / np.sqrt(len(Diag_elements))})
                heatmap.figure.savefig("heatmap_ped_diagonal.png", bbox_inches="tight", dpi=500)
                plt.close(heatmap.figure)
            if matrix_type == "contr":
                #rows, cols = contribution_matrix.shape
                #figsize = (cols, rows)
                #plt.figure(figsize=figsize)

                heatmap_df = pd.DataFrame(contribution_matrix).applymap("{0:.2f}".format)
                heatmap_df = heatmap_df.rename(columns=columns)
                heatmap_df = heatmap_df.astype('float')
                heatmap_df.index = all_internals_string
                heatmap_df.to_csv("ped_contribution_raw", sep='\t')

                
                #TODO Maybe Restructure the Heatmap for bigger Molecules
                sns.set(font_scale=0.35)
                heatmap = sns.heatmap(heatmap_df, cmap="Blues", annot=True, fmt='.1f', cbar_kws={"label": "Contribution %"}) #, annot_kws={"size": 20 / np.sqrt(len(contribution_matrix))})
                heatmap.figure.savefig("heatmap_contribution_table.png", bbox_inches="tight", dpi=600)
                plt.close(heatmap.figure)

    if args.csv:
        for matrix_type in args.csv:
            if matrix_type == "ved":
                Results.to_csv("ved_matrix.csv")
            if matrix_type == "diag":
                DiagonalElementsPED.to_csv("ped_diagonal.csv")
            if matrix_type == "contr":
                ContributionTable.to_csv("contribution_table.csv")


    # here the individual matrices can be computed, one can comment them out
    # if not needed

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

        logfile.write_logfile_extended_results(out, PED, KED, TED, sum_check_PED[mode], sum_check_KED[mode],
                                               sum_check_TED[mode], normal_coord_harmonic_frequencies[mode])

    logfile.call_shutdown()

    print("Runtime: %s seconds" % (time.time() - start_time))


if __name__ == '__main__':
    main()

