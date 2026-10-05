# general packages

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
import re
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from mendeleev.fetch import fetch_table
import inquirer
import pyfiglet
import pubchempy as pcp  # import pubchemmpy --> pubchempy as a databank
import plotly.graph_objects as go
import sys
import random

# for heatmap
mpl.rcParams["backend"] = "Agg"
# do not show messages
logging.getLogger("matplotlib").setLevel(logging.ERROR)


"""
Nomodeco Modules
"""

from nomodeco.libraries import icsel
from nomodeco.libraries import bmatrix
from nomodeco.libraries import logfile
from nomodeco.libraries import molpro_parser
from nomodeco.libraries import specifications
from nomodeco.libraries import icset_opt
from nomodeco.libraries import linear_bends
from nomodeco.libraries import arguments
from nomodeco.libraries import specifications as sp
from nomodeco.libraries import topology as tp
from nomodeco.libraries.molecule_class import Molecule
from nomodeco.libraries.ic_class import InternalCoordinates
from nomodeco.libraries import gaussian_parser
from nomodeco.libraries import orca_parser
from nomodeco.libraries import zmat
from nomodeco.libraries import pymolpro as pymolpro_lib


def _generate_heatmap(
    matrix, columns_map, row_labels, filename, cbar_label=None, dpi=500
    ):
    """
    Helper function to process data and make heatmap
    """
    df = pd.DataFrame(matrix).map("{0:.2f}".format).astype(float)
    df = df.rename(columns=columns_map)
    df.index = row_labels

    rows, cols = df.shape

    cell_size = 0.85 if max(rows, cols) <= 20 else 0.5
    fig_width = max(10, cols * cell_size)
    fig_height = max(8, rows * cell_size)
    
    max_dim = max(rows, cols)
    font_scale = max(0.5, 1.6 - (max_dim * 0.015))
    annot_size = max(9, 45/np.sqrt(max_dim))

    sns.set_theme(font_scale=font_scale)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    cbar_kws = {"label": cbar_label} if cbar_label else {}

    heatmap = sns.heatmap(
        df,
        cmap = "Blues",
        annot=True,
        fmt=".2f" if cbar_label is None else ".1f",
        ax=ax,
        cbar_kws=cbar_kws,
        annot_kws={"size": annot_size}
    )

    # Save and reset
    fig.savefig(filename, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    sns.reset_defaults()


def get_mass_information() -> pd.DataFrame:
    """
    Returns atomic mass for all elements using mendeleev
    """
    df = fetch_table("elements")
    mass_info = df.loc[:, ["symbol", "atomic_weight"]]
    deuterium_info = pd.DataFrame({"symbol": ["D"], "atomic_weight": [2.014102]})
    mass_info = pd.concat([mass_info, deuterium_info])
    mass_info.set_index("symbol", inplace=True)
    return mass_info


def get_bond_information() -> pd.DataFrame:
    """
    Get bond information for all elements

    Returns:
        A pd.dataframe with symbol, covalent radius and vdw radius
    """

    df = fetch_table("elements")

    bond_info = df.loc[:, ["symbol", "covalent_radius_pyykko", "vdw_radius"]]

    bond_info.set_index("symbol", inplace=True)
    bond_info /= 100
    return bond_info


BOND_INFO = get_bond_information()


def reciprocal_square_massvector(atoms) -> pd.DataFrame:
    """
    Get the reciprocal square massvector using atoms.symbols

    Parameters:
        atoms: a list of atoms of the atom class

    Returns:
        A dataframe with the reciprocal square massvector
    """
    n_atoms = len(atoms)
    diag_reciprocal_square = np.zeros(3 * n_atoms)
    MASS_INFO = get_mass_information()
    for i in range(0, n_atoms):
        diag_reciprocal_square[3 * i : 3 * i + 3] = 1 / np.sqrt(
            MASS_INFO.loc[atoms[i].symbol.strip(string.digits)]
        )
    return diag_reciprocal_square


def reciprocal_massvector(atoms) -> pd.DataFrame:
    """
    Returns the reciprocal massvector out of a list of atoms

    Parameters:
        atoms: a list of atoms of the atom class

    Returns:
        A dataframe with the reciprocal massvector
    """
    n_atoms = len(atoms)
    diag_reciprocal = np.zeros(3 * n_atoms)
    MASS_INFO = get_mass_information()
    for i in range(0, n_atoms):
        diag_reciprocal[3 * i : 3 * i + 3] = 1 / (
            MASS_INFO.loc[atoms[i].symbol.strip(string.digits)]
        )
    return diag_reciprocal


# change to deuterium


def strip_numbers(string) -> str:
    """
    Strips the numbers of a string, e.q. H2 -> H
    """
    return "".join([char for char in string if not char.isdigit()])


# welcome Message
def print_welcome_message():
    """
    Prints welcoming message in the command line
    """
    header = pyfiglet.figlet_format(
        "NOMODECO.PY", font="starwars", width=110, justify="center"
    )
    print(header)


start_time = time.time()  # measure runtime


def main():
    os.system("clear")
    # get args
    args = arguments.get_args()
    # print welcome message
    print_welcome_message()

    if args.pymolpro:
        atoms, n_atoms, CartesianF_Matrix, outputfile = pymolpro_lib.run_pymolpro_workflow()

    # Gaussian Parser
    if not args.gv == None:

        # The first file for the gaussian importer specifies the atoms and coordinates
        with open(args.gv[0]) as inputfile:
            atoms = gaussian_parser.parse_xyz(inputfile)
            n_atoms = len(atoms)
        with open(args.gv[0]) as inputfile:
            CartesianF_Matrix = gaussian_parser.parse_cartesian_force_constants(
                inputfile, n_atoms
            )
            outputfile = logfile.create_filename_out(inputfile.name)

    # Molpro Parser Intitialization
    if not args.molpro == None:
        with open(args.molpro[0]) as inputfile:
            atoms = molpro_parser.parse_xyz_from_inputfile(inputfile)
            n_atoms = len(atoms)
        with open(args.molpro[0]) as inputfile:
            CartesianF_Matrix = molpro_parser.parse_Cartesian_F_Matrix_from_inputfile(
                inputfile
            )
            outputfile = logfile.create_filename_out(inputfile.name)

    # Orca Parser
    if not args.orca == None:
        with open(args.orca[0]) as inputfile:
            atoms = orca_parser.parse_xyz_from_inputfile(inputfile)
            n_atoms = len(atoms)
        with open(args.orca[0]) as inputfile:
            CartesianF_Matrix = orca_parser.parse_cartesian_force_constants(
                inputfile, n_atoms
            )
            outputfile = logfile.create_filename_out(inputfile.name)

    if os.path.exists(outputfile):
        i = 1
        while True:
            new_outputfile_name = f"{outputfile}_{i}"
            if not os.path.exists(new_outputfile_name):
                os.rename(outputfile, new_outputfile_name)
                break
            i += 1
    out = logfile.setup_logger("outputfile", outputfile)
    logfile.write_logfile_header(out)

    # BUGFIX for Chloride which for some reason gets written as CL by molpro
    updated_atoms = [
        (
            Molecule.Atom(symbol="Cl", coordinates=atom.coordinates)
            if atom.symbol == "CL"
            else atom
        )
        for atom in atoms
    ]
    atoms = updated_atoms

    atoms = Molecule(atoms)

    molecule = mg.Molecule(
        [strip_numbers(atom.symbol) for atom in atoms],
        [atom.coordinates for atom in atoms],
    )

    molecule_pg = PointGroupAnalyzer(molecule)
    point_group_sch = molecule_pg.sch_symbol


    if args.nomodeco_coords == None:

        """
        IC Generation
        """
        # Intermolecular Bonds:
        # Use Degree of Covalance https://doi.org/10.1002/qua.21049 for hydrogen bond detection
        print("Generating intra- and intermolecular internal coordinates...")
        degofc_table = atoms.degree_of_covalance()
        cov_bonds = atoms.covalent_bonds(degofc_table)

        # Pass covalent bonds to specification
        sp.covalent_bonds = cov_bonds

        # Detect and generate covalent submolecules

        _, _, cov_submolecules_symbols = atoms.detect_submolecules()

        # Generate Hydrogen Bond and the Acceptor-Donor Coordinate

        h_bonds = atoms.intermolecular_h_bond(degofc_table, cov_submolecules_symbols)
        acc_don_bonds = atoms.intermolecular_acceptor_donor(
            degofc_table, cov_submolecules_symbols
        )

        Total_IC_dict = InternalCoordinates()

        # Add covalent ICs
        Total_IC_dict.add_coordinate("cov_bond", cov_bonds)
        Total_IC_dict.add_coordinate("h_bond", h_bonds)
        Total_IC_dict.add_coordinate("acc_don", acc_don_bonds)
        Total_IC_dict.add_coordinate("cov_angles", atoms.generate_angles(cov_bonds)[0])
        Total_IC_dict.add_coordinate(
            "cov_linear_angles", atoms.generate_angles(cov_bonds)[1]
        )
        Total_IC_dict.add_coordinate(
            "cov_dihedrals", atoms.generate_dihedrals(cov_bonds)
        )

        # Define Total Bonds
        bonds = list(set(cov_bonds).union(set(h_bonds)))
        bond_acc_don = list(set(cov_bonds).union(set(acc_don_bonds)))

        Total_IC_dict.add_coord_diff(
            "h_bond_angles",
            atoms.generate_angles(bonds)[0],
            Total_IC_dict["cov_angles"],
        )
        Total_IC_dict.add_coord_diff_linear(
            "h_bond_linear_angles",
            atoms.generate_angles(bonds)[1],
            Total_IC_dict["cov_linear_angles"],
        )
        Total_IC_dict.add_coord_diff(
            "h_bond_dihedrals",
            atoms.generate_dihedrals(bonds),
            Total_IC_dict["cov_dihedrals"],
        )
        Total_IC_dict.add_coord_diff(
            "acc_don_angles",
            atoms.generate_angles(bond_acc_don)[0],
            Total_IC_dict["cov_angles"],
        )
        Total_IC_dict.add_coord_diff_linear(
            "acc_don_linear_angles",
            atoms.generate_angles(bond_acc_don)[1],
            Total_IC_dict["cov_linear_angles"],
        )
        Total_IC_dict.add_coord_diff(
            "acc_don_dihedrals",
            atoms.generate_dihedrals(bond_acc_don),
            Total_IC_dict["cov_dihedrals"],
        )

        # Assing Coordinates from IC Dict to Variables
        angles = Total_IC_dict["cov_angles"] + Total_IC_dict["h_bond_angles"]
        linear_angles = (
            Total_IC_dict["cov_linear_angles"] + Total_IC_dict["h_bond_linear_angles"]
        )
        dihedrals = Total_IC_dict["cov_dihedrals"] + Total_IC_dict["h_bond_dihedrals"]

        # If no valid dihedrals found append acc_don_dihedrals
        if len(dihedrals) == 0:
            dihedrals = (
                Total_IC_dict["cov_dihedrals"] + Total_IC_dict["acc_don_dihedrals"]
            )
        print("Initial Generation of Internal Coordinates finished.")
        print(f"RAM usage of Total IC dictionary: {sys.getsizeof(Total_IC_dict) / 1000} MB")


        """
       Generating the Specification for the primary calculation
       """

        # Setting specifications for calculation: check if molecule is linear, planar or a general molecule
        specification = dict()
        specification = specifications.calculation_specification(
            specification, atoms, molecule_pg, bonds, angles, linear_angles
        )

        # Generation of all possible out-of-plane motions

        # Keys in Total IC_dict get initialized but not filled if system is not planar
        Total_IC_dict.add_coordinate("cov_oop", [])
        Total_IC_dict.add_coordinate("h_bond_oop", [])
        Total_IC_dict.add_coordinate("acc_don_oop", [])

        if specification["planar"] == "yes":
            Total_IC_dict.add_coordinate(
                "cov_oop", atoms.generate_out_of_plane(cov_bonds)
            )
            Total_IC_dict.add_coord_diff(
                "h_bond_oop",
                atoms.generate_out_of_plane(bonds),
                Total_IC_dict["cov_oop"],
            )
            Total_IC_dict.add_coord_diff(
                "acc_don_oop",
                atoms.generate_out_of_plane(bond_acc_don),
                Total_IC_dict["cov_oop"],
            )
            out_of_plane = Total_IC_dict["cov_oop"] + Total_IC_dict["h_bond_oop"]
        elif (
            specification["planar"] == "no"
            and specification["planar submolecule(s)"] == []
        ):
            out_of_plane = []
        elif specification["planar"] == "no" and not (
            specification["planar submolecule(s)"] == []
        ):
            out_of_plane = atoms.generate_oop_planar_subunits(
                bonds, specification["planar submolecule(s)"]
            )
        else:
            return out.error(
                "Classification of whether topology is planar or not could not be determined!"
            )

        # determine internal degrees of freedom
        idof = 0
        if specification["linearity"] == "fully linear":
            idof = atoms.idof_linear()
        else:
            idof = atoms.idof_general()

        # update out file

        logfile.write_logfile_oop_treatment(
            out, specification["planar"], specification["planar submolecule(s)"]
        )
        logfile.write_logfile_symmetry_treatment(out, specification, point_group_sch)

        # print the specification to out
        print("Specification of the system for the primary calculation:")
        for key, value in specification.items():
            print(f"{key}: {value}")

        """
       Passing Section for Topology.py
       """

        tp.atoms_list = atoms
        tp.Total_IC_dict = Total_IC_dict

        """
       Diag Mass Matrix and Reciprocal Square root masses
       """

        # Computation of the diagonal mass matrices with
        # the reciprocal and square root reciprocal masses
        diag_reciprocal_square = reciprocal_square_massvector(atoms)
        reciprocal_square_massmatrix = np.diag(diag_reciprocal_square)
        diag_reciprocal = reciprocal_massvector(atoms)
        reciprocal_massmatrix = np.diag(diag_reciprocal)

        # Determination of the Normal Modes and eigenvalues
        # via the diagonalization of the mass-weighted Cartesian F Matrix
        Mass_weighted_CartesianF_Matrix = (
            np.transpose(reciprocal_square_massmatrix)
            @ CartesianF_Matrix
            @ reciprocal_square_massmatrix
        )

        Cartesian_eigenvalues, L = np.linalg.eigh(Mass_weighted_CartesianF_Matrix)
        # print("Cartesian_eigenvalues (EV from mw hessian):", Cartesian_eigenvalues)

        # Align the linear bends with the normal modes: rotates degenerate pairs in L and gives
        # the reference vectors used by every b_matrix call below
        L, bend_refs = linear_bends.linear_bend_references(
            atoms, linear_angles, L, Cartesian_eigenvalues, diag_reciprocal_square, idof
        )

        # Determination of the normal modes of zero and low Frequencies

        rottra = L[:, 0 : (3 * n_atoms - idof)]

        logfile.write_logfile_generated_IC(
            out, bonds, angles, linear_angles, out_of_plane, dihedrals, idof
        )

        # Print Mass Weighted F Matrix to out
        logfile.write_mass_weighted_f_matrix(
            out, Mass_weighted_CartesianF_Matrix, atoms
        )

        # Print L matrix and Cartesian Eigenvalues

        # logfile.write_cartesian_eigenvalues(out, Cartesian_eigenvalues)

        logfile.write_l_matrix(out, L, Cartesian_eigenvalues, atoms)

        """
       Logging Intermolecular IC Information
       """
        # Information only gets logged if hydrogen bonds where found in the structure
        if len(Total_IC_dict["h_bond"]) != 0:
            logfile.write_hydrogen_bond_information(
                out,
                Total_IC_dict["h_bond"],
                Total_IC_dict["acc_don"],
                Total_IC_dict["h_bond_angles"],
                Total_IC_dict["acc_don_angles"],
                Total_IC_dict["h_bond_linear_angles"],
                Total_IC_dict["acc_don_linear_angles"],
                Total_IC_dict["h_bond_dihedrals"],
                Total_IC_dict["acc_don_dihedrals"],
                Total_IC_dict["h_bond_oop"],
                Total_IC_dict["acc_don_oop"],
            )

        # get symmetric coordinates
        if args.penalty1 != 0:
            symmetric_bonds = icsel.get_symm_bonds(bonds, specification)
            symmetric_angles = icsel.get_symm_angles(angles, specification)
            symmetric_dihedrals = icsel.get_symm_dihedrals(dihedrals, specification)
            symmetric_coordinates = {
                **symmetric_bonds,
                **symmetric_angles,
                **symmetric_dihedrals,
            }
        else:
            symmetric_coordinates = dict()

        if not args.graph_ic:
            print("Generating IC sets based on Topology and Specification...")
            ic_dict = icsel.get_sets(
                idof,
                out,
                atoms,
                bonds,
                angles,
                linear_angles,
                out_of_plane,
                dihedrals,
                specification,
            )
            print("Initial Generation of IC sets finished.")
        if args.graph_ic:
            ic_dict = {}
            with open(args.graph_ic[0]) as set_file:
                lines = set_file.readlines()
                for i in range(len(lines)):
                    ic_dict[i] = {
                        "bonds": [],
                        "angles": [],
                        "linear valence angles": [],
                        "out of plane angles": [],
                        "dihedrals": [],
                    }

                # now we fill up this dictionary
                for line in lines:
                    if ":" in line:
                        key, value = line.split(":", 1)
                        value = eval(value.strip())
                        ic_dict[int(key)]["bonds"] = value["bonds"]
                        ic_dict[int(key)]["angles"] = value["angles"]
                        ic_dict[int(key)]["linear valence angles"] = value[
                            "linear valence angles"
                        ]
                        ic_dict[int(key)]["out of plane angles"] = value[
                            "out of plane angles"
                        ]
                        ic_dict[int(key)]["dihedrals"] = value["dihedrals"]

            print("Length of Imported IC set", len(ic_dict))

        _MAX_IC_SETS = 100_000
        _total_sets = len(ic_dict)
        if _total_sets > _MAX_IC_SETS:
            print(
                f"Large search space detected: {_total_sets:,} sets. "
                f"Randomly sampling {_MAX_IC_SETS:,} for evaluation..."
            )
            _sampled_keys = random.sample(range(_total_sets), _MAX_IC_SETS)
            ic_dict = {i: ic_dict[k] for i, k in enumerate(_sampled_keys)}

        print(f"{len(ic_dict):,} IC sets will be evaluated. "
              f"Memory footprint: {sys.getsizeof(ic_dict) / 1e6:.1f} MB")

        result = icset_opt.find_optimal_coordinate_set(
            ic_dict,
            args,
            idof,
            reciprocal_massmatrix,
            reciprocal_square_massmatrix,
            rottra,
            CartesianF_Matrix,
            atoms,
            symmetric_coordinates,
            L,
            args.penalty1,
            args.penalty2,
            bend_refs,
        )
        optimal_set = result["set"]
        if optimal_set is None:
            out.error("No valid internal coordinate set was found (0 complete sets evaluated).")
            sys.exit("Nomodeco: no valid internal coordinate set was found. See the .out file.")

    if not args.nomodeco_coords == None:

        ic_dict = {}
        with open(args.nomodeco_coords[0]) as set_file:

            # TODO we need another for loop
            ic_set = {}
            lines = set_file.readlines()
            for line in lines:
                line = line.strip()

                # now split with = as a seperator
                if "=" in line:
                    key, value = line.split("=", 1)

                    # remove whitespaces in keys
                    key = key.strip()

                    # remove newlines or other stuff from value
                    value = value.strip()

                    # remove the trailing commas and the brackets
                    value = value.rstrip("]\n").lstrip("[").strip()

                    # after this point we have tuples
                    if value:
                        # convert string into into actual list
                        value = eval("[" + value + "]")

                    else:
                        value = []

                    ic_set[key] = value

        ic_dict[0] = ic_set

        idof = 3 * n_atoms - 6

        # Computation of the diagonal mass matrices with
        # the reciprocal and square root reciprocal masses
        diag_reciprocal_square = reciprocal_square_massvector(atoms)
        reciprocal_square_massmatrix = np.diag(diag_reciprocal_square)
        diag_reciprocal = reciprocal_massvector(atoms)
        reciprocal_massmatrix = np.diag(diag_reciprocal)

        # Determination of the Normal Modes and eigenvalues
        # via the diagonalization of the mass-weighted Cartesian F Matrix
        Mass_weighted_CartesianF_Matrix = (
            np.transpose(reciprocal_square_massmatrix)
            @ CartesianF_Matrix
            @ reciprocal_square_massmatrix
        )
        # Write the mass_weighted F matrix to out

        Cartesian_eigenvalues, L = np.linalg.eigh(Mass_weighted_CartesianF_Matrix)
        # print("Cartesian_eigenvalues (EV from mw hessian):", Cartesian_eigenvalues)

        L, bend_refs = linear_bends.linear_bend_references(
            atoms, ic_set.get("linear valence angles", []), L, Cartesian_eigenvalues,
            diag_reciprocal_square, idof,
        )

        # Determination of the normal modes of zero and low Frequencies

        rottra = L[:, 0 : (3 * n_atoms - idof)]

        # TODO: Implement argspenalty for user specified IC sets
        symmetric_coordinates = dict()

        result = icset_opt.find_optimal_coordinate_set(
            ic_dict,
            args,
            idof,
            reciprocal_massmatrix,
            reciprocal_square_massmatrix,
            rottra,
            CartesianF_Matrix,
            atoms,
            symmetric_coordinates,
            L,
            args.penalty1,
            args.penalty2,
            bend_refs,
        )
        optimal_set = result["set"]
        if optimal_set is None:
            out.error("The user-defined internal coordinate set is not complete or gives imaginary intrinsic frequencies.")
            sys.exit("Nomodeco: the user-defined internal coordinate set is not valid. See the .out file.")

    """''
    Final calculation with optimal set
    """ ""
    bonds = optimal_set["bonds"]
    angles = optimal_set["angles"]
    linear_angles = optimal_set["linear valence angles"]
    out_of_plane = optimal_set["out of plane angles"]
    dihedrals = optimal_set["dihedrals"]

    n_internals = (
        len(bonds)
        + len(angles)
        + len(linear_angles)
        + len(out_of_plane)
        + len(dihedrals)
    )
    red = n_internals - idof

    # Augmenting the B-Matrix with rottra, calculating

    B = np.concatenate(
        (
            bmatrix.b_matrix(
                atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof, bend_refs
            ),
            np.transpose(rottra),
        ),
        axis=0,
    )
    # log the argumented b matrix

    logfile.write_b_matrix(
        out, B, rottra, atoms, bonds, angles, linear_angles, out_of_plane, dihedrals
    )

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

    logfile.write_logfile_information_results(
        out, n_internals, red, bonds, angles, linear_angles, out_of_plane, dihedrals
    )

    """
    Logging the final intermolecular_ICs
    """

    # TODO: Maybe find another way to build a switch her
    if args.nomodeco_coords == None:

        if len(Total_IC_dict["h_bond"]) != 0:
            logfile.write_logfile_final_intermolecular_ics(
                out,
                Total_IC_dict.common_coordinate("h_bond", bonds),
                Total_IC_dict.common_coordinate("acc_don", bonds),
                Total_IC_dict.common_coordinate("h_bond_angles", angles),
                Total_IC_dict.common_coordinate("acc_don_angles", angles),
                Total_IC_dict.common_coordinate("h_bond_linear_angles", linear_angles),
                Total_IC_dict.common_coordinate("acc_don_linear_angles", linear_angles),
                Total_IC_dict.common_coordinate("h_bond_dihedrals", dihedrals),
                Total_IC_dict.common_coordinate("acc_don_dihedrals", dihedrals),
                Total_IC_dict.common_coordinate("h_bond_oop", out_of_plane),
                Total_IC_dict.common_coordinate("acc_don_oop", out_of_plane),
            )
    

    # printing the B matrix without rottra

    b_without_rottra = bmatrix.b_matrix(
        atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof, bend_refs
    )

    logfile.write_b_matrix_raw(
        out,
        b_without_rottra,
        atoms,
        bonds,
        angles,
        linear_angles,
        out_of_plane,
        dihedrals,
    )

    """'' 
    --------------------------- Main-Calculation ------------------------------
    """ ""

    # Calculation of the mass-weighted normal modes in Cartesian Coordinates

    l = reciprocal_square_massmatrix @ L

    # Calculation of the mass-weighted normal modes in Internal Coordinates

    D = B @ l

    logfile.write_d_matrix(
        out, D, atoms, rottra, bonds, angles, linear_angles, out_of_plane, dihedrals
    )

    # Calculation of the Vibrational Density Matrices / PED, KED and TED matrices

    eigenvalues = np.transpose(D) @ InternalF_Matrix @ D
    eigenvalues = np.diag(eigenvalues)
    # print("eigenvalues (from IC space):", eigenvalues)

    num_rottra = 3 * n_atoms - idof

    """'' 
    ------------------------------- Results --------------------------------------
    """ ""

    P = np.zeros(
        (n_internals - red, n_internals + num_rottra, n_internals + num_rottra)
    )
    T = np.zeros(
        (n_internals - red, n_internals + num_rottra, n_internals + num_rottra)
    )
    E = np.zeros(
        (n_internals - red, n_internals + num_rottra, n_internals + num_rottra)
    )

    for i in range(0, n_internals - red):
        for m in range(0, n_internals + num_rottra):
            for n in range(0, n_internals + num_rottra):
                k = i + num_rottra
                P[i][m][n] = (
                    D[m][k] * InternalF_Matrix[m][n] * D[n][k] / eigenvalues[k]
                )  # PED
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
        contribution_matrix[:, i] = ((Diag_elements[:, i] / sum_diag[i]) * 100).astype(
            float
        )

    # compute intrinsic frequencies
    nu = np.zeros(n_internals)
    for n in range(0, n_internals):
        for m in range(0, n_internals):
            for i in range(0, n_internals - red):
                k = i + num_rottra
                nu[n] += D[m][k] * InternalF_Matrix[m][n] * D[n][k]

    nu_final = np.sqrt(nu) * 5140.4981

    normal_coord_harmonic_frequencies = (
        np.sqrt(eigenvalues[(3 * n_atoms - idof) : 3 * n_atoms]) * 5140.4981
    )
    normal_coord_harmonic_frequencies = np.around(
        normal_coord_harmonic_frequencies
    ).astype(int)
    normal_coord_harmonic_frequencies_string = normal_coord_harmonic_frequencies.astype(
        "str"
    )

    all_internals = bonds + angles + linear_angles + out_of_plane + dihedrals

    all_internals_string = []
    for internal in all_internals:
        all_internals_string.append("(" + ", ".join(internal) + ")")

    Results = pd.DataFrame()
    Results["Internal Coordinate"] = all_internals_string
    Results["Intrinsic Frequencies"] = pd.DataFrame(nu_final).map("{0:.2f}".format)
    Results = Results.join(pd.DataFrame(ved_matrix).map("{0:.2f}".format))

    DiagonalElementsPED = pd.DataFrame()
    DiagonalElementsPED["Internal Coordinate"] = all_internals_string
    DiagonalElementsPED["Intrinsic Frequencies"] = pd.DataFrame(nu_final).map(
        "{0:.2f}".format
    )
    DiagonalElementsPED = DiagonalElementsPED.join(
        pd.DataFrame(Diag_elements).map("{0:.2f}".format)
    )

    ContributionTable = pd.DataFrame()
    ContributionTable["Internal Coordinate"] = all_internals_string
    ContributionTable["Intrinsic Frequencies"] = pd.DataFrame(nu_final).map(
        "{0:.2f}".format
    )
    ContributionTable = ContributionTable.join(
        pd.DataFrame(contribution_matrix).map("{0:.2f}".format)
    )

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
        row = {"Frequency": freq}
        details = []
        for i, coord in enumerate(ContributionTable["Internal Coordinate"]):
            percentage = ContributionTable.loc[i, freq]
            # weird quickfix maybe repair this some time
            if isinstance(percentage, pd.Series):
                percentage = float(percentage[1])
            elif isinstance(percentage, str) and percentage != ".":
                percentage = float(percentage)
            if percentage > 10:
                details.append(f"{coord}: {percentage:.2f}%")
        row["Details"] = "; ".join(details)
        Contribution_Sliced.append(row)
    annotated_df = pd.DataFrame(Contribution_Sliced)

    ContributionTable_T = ContributionTable.transpose()

    # TODO:  line breaks in output file

    logfile.write_logfile_results(
        out, Results, DiagonalElementsPED, ContributionTable, sum_check_VED
    )

    if args.barplot:
        from matplotlib.patches import Patch
        from collections import defaultdict

        # Pre-build frozensets for O(1) IC-type lookup
        _bond_set  = {tuple(b)  for b  in bonds}
        _angle_set = {tuple(a)  for a  in angles}
        _oop_set   = {tuple(o)  for o  in out_of_plane}
        _dih_set   = {tuple(d)  for d  in dihedrals}
        _la_set    = {tuple(la) for la in linear_angles}

        def _ic_type(coord_str):
            tup = tuple(c.strip() for c in coord_str.strip("()").split(","))
            rev = tup[::-1]
            if tup in _bond_set  or rev in _bond_set:  return "Bond"
            if tup in _angle_set or rev in _angle_set: return "Angle"
            if tup in _oop_set:                         return "Out-of-plane"
            if tup in _dih_set   or rev in _dih_set:   return "Dihedral"
            if tup in _la_set    or rev in _la_set:    return "Linear angle"
            return "Other"

        _type_cmaps = {
            "Bond":          plt.cm.Blues,
            "Angle":         plt.cm.Oranges,
            "Out-of-plane":  plt.cm.Greens,
            "Dihedral":      plt.cm.RdPu,
            "Linear angle":  plt.cm.Purples,
            "Other":         plt.cm.Greys,
        }
        _type_order = list(_type_cmaps.keys())

        threshold = 20.0
        value_cols = [c for c in ContributionTable.columns
                      if c not in ["Internal Coordinate", "Intrinsic Frequencies"]]

        long_df = ContributionTable.melt(
            id_vars=["Internal Coordinate", "Intrinsic Frequencies"],
            value_vars=value_cols,
            var_name="Frequency",
            value_name="Contribution",
        )
        long_df["Contribution"] = pd.to_numeric(long_df["Contribution"], errors="coerce").fillna(0.0)
        long_df = long_df[long_df["Contribution"] > threshold]
        long_df["IC Type"] = long_df["Internal Coordinate"].apply(_ic_type)

        # Individual-IC pivot — each IC gets its own segment
        pivot_df = long_df.pivot_table(
            index="Frequency",
            columns="Internal Coordinate",
            values="Contribution",
            aggfunc="sum",
            fill_value=0.0,
        )
        pivot_df = pivot_df[pivot_df.sum(axis=1) > 0]

        # Sort x-axis from lowest to highest harmonic frequency
        pivot_df.index = pd.to_numeric(pivot_df.index, errors="coerce")
        pivot_df = pivot_df.sort_index()
        pivot_df.index = pivot_df.index.map(lambda f: f"{f:.1f}")

        # Sort columns: group by type, then alphabetically within type
        ic_type_map = {ic: _ic_type(ic) for ic in pivot_df.columns}
        pivot_df = pivot_df[sorted(pivot_df.columns,
            key=lambda ic: (_type_order.index(ic_type_map.get(ic, "Other")), ic))]

        # Group ICs by type (in column order so indices stay consistent)
        type_ics: dict = defaultdict(list)
        for ic in pivot_df.columns:
            type_ics[ic_type_map[ic]].append(ic)

        # Assign color + hatch per IC
        ic_colors:       dict = {}
        type_repr_color: dict = {}
        for ic_type, ics in type_ics.items():
            shades = _type_cmaps[ic_type](np.linspace(0.15, 0.95, max(len(ics), 1)))
            type_repr_color[ic_type] = shades[len(shades) // 2]
            for ic, color in zip(ics, shades):
                ic_colors[ic] = color

        bar_colors = [ic_colors[ic] for ic in pivot_df.columns]

        fig, ax = plt.subplots(figsize=(max(14, len(pivot_df) * 0.55), 6))
        pivot_df.plot(
            kind="bar", stacked=True, ax=ax,
            color=bar_colors, width=0.8, legend=False,
        )

        ax.set_xlabel("Harmonic Frequency (cm⁻¹)", fontsize=11)
        ax.set_ylabel("Contribution (%)", fontsize=11)
        ax.set_title("Contributions of Internal Coordinates to Vibrational Modes", fontsize=12)
        ax.tick_params(axis="x", rotation=45, labelsize=8)

        # Legend: type header + each IC with its exact color AND hatch
        legend_handles = []
        for ic_type in _type_order:
            ics = type_ics.get(ic_type, [])
            if not ics:
                continue
            legend_handles.append(
                Patch(facecolor=type_repr_color[ic_type], edgecolor="black",
                      linewidth=1.2, label=f"── {ic_type} ──")
            )
            for ic in ics:
                legend_handles.append(
                    Patch(facecolor=ic_colors[ic], edgecolor="none", label=f"   {ic}")
                )

        ax.legend(
            handles=legend_handles,
            loc="upper left", fontsize=7,
            bbox_to_anchor=(1.01, 1),
            borderaxespad=0, framealpha=0.9,
            title="Internal Coordinates", title_fontsize=8,
        )
        plt.tight_layout()
        plt.savefig("contribution_barplot.png", dpi=300, bbox_inches="tight")
        plt.close()


    if args.sankey_plot != 0:
        # If given always create a sankey diagramm 

        # ano
        # Helper Function to determine the type of coordinate
        def get_coordinate_type(coord):
            # Check if its already a tuple
            if not isinstance(coord, tuple):
                coord_tuple = tuple(coord.strip("()").split(","))
                coord_tuple = tuple(element.strip() for element in coord_tuple)
            else:
                coord_tuple = coord
            
            if coord_tuple in bonds:
                return "bond"
            elif coord_tuple in angles:
                return "angle"
            elif coord_tuple in linear_angles:
                return "linear_angle"
            elif coord_tuple in dihedrals:
                return "dihedral"
            elif coord_tuple in out_of_plane:
                return "out-of-plane"

        # Define a color mapping for the different coordinate types
        color_map = {
            "bond": "rgba(255, 99, 132, 0.8)",  # Red
            "angle": "rgba(54, 162, 235, 0.8)",  # Blue
            "dihedral": "rgba(75, 192, 192, 0.8)",  # Green
            "out-of-plane": "rgba(255, 206, 86, 0.8)",  # Yellow
            "linear_angle": "rgba(153, 102, 255, 0.8)",  # Purple
        }

        # Prepare Sankey Diagramm
        labels = []
        source = []
        target = []
        value = []
        link_labels = [] # Store Contribution Percentages
        link_colors = [] # Store colors for links
        coord_types = []

        # Create Nodes first all internal coordinates
        int_coords_label = ContributionTable["Internal Coordinate"].tolist()
        modes = [str(col) for col in ContributionTable.columns if col not in ["Internal Coordinate", "Intrinsic Frequencies"]] 


        


        
        labels = int_coords_label + modes


        # Create Mappings from names to indices
        coord_indices = {coord: idx for idx, coord in enumerate(int_coords_label)}
        mode_indices = {mode: idx + len(int_coords_label) for idx, mode in enumerate(modes)}

        # Build up the links
        for _, row in ContributionTable.iterrows():
            coord = row["Internal Coordinate"] 
            # Here we decide which type of coordinates we take
            if args.sankey_plot == 2:
                coord_check = tuple(coord.strip("()").split(","))
                coord_check = tuple(element.strip() for element in coord_check)
                print(coord_check)
                # Skipt the loop element if the coordinate is not contained in the h_bond, h_bond_angles ..
                if not (
                    coord_check in Total_IC_dict["h_bond"]
                    or coord_check in Total_IC_dict["h_bond_angles"]
                    or coord_check in Total_IC_dict["h_bond_linear_angles"]
                    or coord_check in Total_IC_dict["h_bond_dihedrals"]
                    or coord_check in Total_IC_dict["h_bond_oop"]
                    or coord_check in Total_IC_dict["acc_don"]
                    or coord_check in Total_IC_dict["acc_don_angles"]
                    or coord_check in Total_IC_dict["acc_don_linear_angles"]
                    or coord_check in Total_IC_dict["acc_don_dihedrals"]
                    or coord_check in Total_IC_dict["acc_don_oop"]
                ):
                    continue
                
                

            
            coord_type = get_coordinate_type(coord)

            coord_types.append(coord_type)
            for mode in modes:
                contribution = float(row[mode])
                if contribution > args.min_contr_sankey: 
                    source.append(coord_indices[coord])
                    target.append(mode_indices[mode])
                    value.append(contribution)
                    link_labels.append(f"{coord} to {mode}: {contribution:.1f}%")
                    link_colors.append(color_map[coord_type])  # Assign color based on type
        
        legend_trace = []
        for coord_type, color in color_map.items():
            legend_trace.append(
                go.Scatter(
                    x=[None], y=[None],
                    mode = 'markers',
                    marker=dict(color=color, size=10),
                    name=coord_type,
                    hoverinfo="none"
                )
            ) 

        # Create Sankey Diagram
        fig = go.Figure(
            data=[go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=0.5),
                    label=labels,
                    color="blue"
                ),
                link=dict(
                    source=source,
                    target=target,
                    value=value,
                    label=link_labels,
                    color=link_colors,
                    hovertemplate="%{label}<extra></extra>",
                )
            ),
            *legend_trace
            ]
        )

        fig.update_layout(
            title_text = "Sankey Diagram of Internal Coordinates and Frequencies",
            font_size = 10,
            height = 800,
            showlegend=True,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
        )
        )
        fig.show()

    # heat map results
    # TODO: clean up
    if args.heatmap:
        # 3*n_atoms - (3*n_atoms - idof) == idof
        heatmap_columns = {
            i: normal_coord_harmonic_frequencies[i] for i in range(int(idof))
        }

        for matrix_type in args.heatmap:
            if matrix_type == "ved":
                _generate_heatmap(ved_matrix, columns_map=heatmap_columns, row_labels=all_internals_string, filename="heatmap_ved_matrix.png")
            if matrix_type == "diag":
                _generate_heatmap(Diag_elements, columns_map=heatmap_columns, row_labels=all_internals_string, filename="heatmap_diag_ped.png", cbar_label="Diagonal PED")
            if matrix_type == "contr":
                _generate_heatmap(contribution_matrix, columns_map=heatmap_columns, row_labels=all_internals_string, filename="heatmap_contribution_matrix.png", cbar_label="Contribution (%)")

    if args.csv:
        for matrix_type in args.csv:
            if matrix_type == "ved":
                Results.to_csv("ved_matrix.csv")
            if matrix_type == "diag":
                DiagonalElementsPED.to_csv("ped_diagonal.csv")
            if matrix_type == "contr":
                ContributionTable.to_csv("contribution_table.csv")
    # Generate a Latex Table containing the results of the contribution table
    if args.latex_tab:
        latex_table = annotated_df.to_latex(
            column_format="c|c",
            index=False,
            multicolumn=True,
            header=["Frequencies", "Contributions"],
            escape=False,
        )
        # Adjustable table width
        table_width = ""
        latex_table_adjusted = (
            "Latex Table with Contributions over 10 %\n"
            + "\\begin{adjustbox}{width="
            + str(table_width)
            + "\\textwidth}\n"
            + latex_table
            + "\\end{adjustbox}"
        )

        if os.path.isfile("./contribution_table_latex.txt"):
            f = open("contribution_table_latex.txt", "w")
            f.write(latex_table_adjusted)
            f.close()
        else:
            f = open("contribution_table_latex.txt", "x")
            f.write(latex_table_adjusted)
            f.close()

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

        PED["Internal Coordinate"] = all_internals
        KED["Internal Coordinate"] = all_internals
        TED["Internal Coordinate"] = all_internals
        PED = PED.join(
            pd.DataFrame(P[mode][0:n_internals, 0:n_internals]).map("{0:.2f}".format)
        )
        KED = KED.join(
            pd.DataFrame(T[mode][0:n_internals, 0:n_internals]).map("{0:.2f}".format)
        )
        TED = TED.join(
            pd.DataFrame(E[mode][0:n_internals, 0:n_internals]).map("{0:.2f}".format)
        )
        PED = PED.rename(columns=columns)
        KED = KED.rename(columns=columns)
        TED = TED.rename(columns=columns)

        logfile.write_logfile_extended_results(
            out,
            PED,
            KED,
            TED,
            sum_check_PED[mode],
            sum_check_KED[mode],
            sum_check_TED[mode],
            normal_coord_harmonic_frequencies[mode],
        )

    logfile.call_shutdown()

    print("Runtime: %s seconds" % (time.time() - start_time))


if __name__ == "__main__":
    main()
