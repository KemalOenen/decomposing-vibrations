import argparse
import numpy as np
import pyfiglet
import string
import logging
import icsel

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("oop_directive")
    args = parser.parse_args()
    return args

def create_new_filename(old_filename):
    pieces = old_filename.split(".")
    return "_".join([pieces[0], "nomodeco"]) +  ".log"

def write_logfile_header():
    header=pyfiglet.figlet_format("NOMODECO.PY", font="starwars", width = 110, justify="center")
    logging.info(header)
    logging.info("")
    logging.info("*".center(110, "*"))
    logging.info("Kemal Önen, Dennis F. Dinu, Klaus R. Liedl".center(110))
    logging.info("Liedl Lab, University of Innsbruck".center(110))
    logging.info("*".center(110,"*"))
    logging.info("")

def write_logfile_oop_treatment(oop_directive):
    if oop_directive == "oop":
        logging.info("You have included the consideration of out-of-plane angles. Note that out-of-plane angles tend to perform worse")
        logging.info("than the other internal coordinates in terms of computational cost and in the decomposition. However they can be")
        logging.info("useful for planar systems.")
    elif oop_directive == "no-oop":
        logging.info("You did NOT include the consideration of out-of-plane angles. Note that out-of-plane angles tend to perform worse")
        logging.info("than the other internal coordinates in terms of computational cost and in the decomposition. However they can be")
        logging.info("useful for planar systems.")
    logging.info("")

def write_logfile_symmetry_treatment(specification):
    if specification["angle_symmetry"] == "angle_symmetry":
        logging.info("You have included symmetry considerations for angles. This can be useful for larger systems")
        logging.info("to reduce the number of generally not useful possibilites. Note that for some systems, ")
        logging.info("it is still advisable to not include symmetry considerations for angles with the keyword 'no-angle_symmetry'.")
    if specification["angle_symmetry"] == "no-angle_symmetry":
        logging.info("You have NOT included symmetry considerations for angles. Note that for some systems, ")
        logging.info("symmetry considerations shrink down the computational runtime considerably, while only removing chemically")
        logging.info("uninteristing coordinate sets. If you want to use angle symmetry, use the keyword 'angle_symmetry'.")
    logging.info("")
    if specification["out-of-plane angle_symmetry"] == "out-of-plane_symmetry":
        logging.info("You have included symmetry considerations for out-of-plane angles. This can be useful for larger systems")
        logging.info("to reduce the number of generally not useful possibilites. Note that for some systems, it is still advisable")
        logging.info("to not include symmetry considerations for out-of-plane angles with the keyword 'no-out-of-plane_symmetry'.")
    if specification["out-of-plane angle_symmetry"] == "no-out-of-plane_symmetry":
        logging.info("You have NOT included symmetry considerations for out-of-plane angles. Note that for some systems, ")
        logging.info("symmetry considerations shrink down the computational runtime considerably, while only removing chemically")
        logging.info("uninteristing coordinate sets. If you want to use out-of-plane angle symmetry, use the keyword 'out-of-plane_symmetry'.")
    logging.info("")
    if specification["dihedral_reduction"] == "no-dihedral_reduction":
        logging.info("You have NOT included the sampling of symmetric dihedrals. Note that for some systems, the reduction of dihedrals")
        logging.info("shrinks down the computational runtime considerably, while only removing chemically uninteristing coordinate sets.")
        logging.info("If you want to reduce dihedral dimensions, use the keyword 'dihedral_reduction = n' with n being the number of")
        logging.info("dihedral representatives per symmetric dihedral.")
    if specification["dihedral_reduction"][0] == "dihedral_reduction":
        logging.info("You have included the sampling of symmetric dihedrals (%s dihedrals per symmetry). This can be useful for larger systems",specification["dihedral_reduction"][1])
        logging.info("to reduce the number of generally not useful possibilites. Note that for some systems, it is still advisable")
        logging.info("to not sample symmetric dihedrals with the keyword 'no-dihedral_reduction'.")
    logging.info("")

def write_logfile_generated_IC(bonds, angles, linear_angles, out_of_plane, dihedrals, idof):
    logging.info("The following primitive internals (bonds, (in-plane) angles, linear valence angles, out-of-plane angles and torsions) were generated:".center(110))
    logging.info("bonds: %s", bonds)
    logging.info("in-plane angles: %s", angles)
    logging.info("linear valence-angles: %s", linear_angles)
    logging.info("out-of-plane angles: %s", out_of_plane)
    logging.info("dihedrals: %s", dihedrals)
    logging.info("")
    logging.info("%s internal coordinates are at least needed for the decomposition scheme.".center(50), idof)
    logging.info("")

def write_logfile_nan_freq():
    logging.info("")
    logging.error('Imaginary frequency values were detected for the intrinsic frequencies.')
    logging.info("The computation will be stopped")
    logging.info("")

def write_logfile_information_results(B, B_inv, CartesianF_Matrix, InternalF_Matrix, n_internals, red, bonds, angles, linear_angles, out_of_plane, dihedrals):
    logging.info("")
    logging.info(" Initialization of an internal coordinate set ".center(110, "-"))
    logging.info("")
    logging.info("The following %s Internal Coordinates are used in your analysis:".center(50), n_internals)
    logging.info("bonds: %s", bonds)
    logging.info("in-plane angles: %s", angles)
    logging.info("linear valence-angles: %s", linear_angles)
    logging.info("out-of-plane angles: %s", out_of_plane)
    logging.info("dihedrals: %s", dihedrals)
    logging.info("")
    if red == 1:
        logging.info('There is %s redundant internal coordinates used.', red)
    else:
        logging.info('There are %s redundant internal coordinates used.', red)
    logging.info("")
    logging.info('The condition number of the B-Matrix is given by:')
    logging.info('Sum-Norm of B: %s', icsel.matrix_norm(B, B_inv, 1))
    logging.info('Eucledian-Norm of B: %s', icsel.matrix_norm(B, B_inv, 2))
    logging.info('Maximum-Norm of B: %s', icsel.matrix_norm(B, B_inv, np.inf))
    logging.info("")
    logging.info("")
    logging.info("Testing if the Internal Coordinate Set is complete ...")
    if icsel.test_completeness(CartesianF_Matrix, B, B_inv, InternalF_Matrix) != True:
        logging.error('No! The double transformed f-matrix is NOT the same as in the input.')
        logging.info("")
        logging.error("This set is not complete and will not be computed!")
        logging.info("")
    else:
        logging.info('Yes! The double transformed f-matrix is the same as in the input.')
        logging.info("")

def write_logfile_results(Results1, Results2, Contribution_Table1, Contribution_Table2):
    logging.info(" Results of the Decomposition Scheme ".center(90, "-"))
    logging.info("")
    logging.info("1.) Intrinsic Frequencies for all normal-coordinate frequencies and")
    logging.info("    diagonal elements of the vibrational density matrices at the following")
    logging.info("    harmonic frequencies (all frequency values in cm-1):")
    logging.info("")
    logging.info(Results1.to_string(index=False))
    
    logging.info("")
    logging.info("Contribution Table generated from the Vibrational Density Matrices (values in percent)")
    logging.info("")
    logging.info(Contribution_Table1.to_string(index=False))
    
    logging.info("")
    logging.info("-".center(90, "-"))
    logging.info("")
    
    logging.info("")
    logging.info("2.) Eigenvalues for each internal coordinate per normal-coordinate ")
    logging.info("")
    logging.info(Results2.to_string(index=False))
    
    logging.info("")
    logging.info("Contribution Table generated from the eigenvalues per coordinate (values in percent)")
    logging.info("")
    logging.info(Contribution_Table2.to_string(index=False))
    
    logging.shutdown()
