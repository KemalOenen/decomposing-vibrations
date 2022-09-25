import argparse
import numpy as np
import pyfiglet
import string
import logging

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("oop_directive")
    args = parser.parse_args()
    return args

def create_new_filename(old_filename):
    pieces = old_filename.split(".")
    return "_".join([pieces[0], "nomodeco"]) +  ".log"

def oop_directive() -> string:
    args = get_args()
    if args.oop_directive == "oop":
        return "oop"
    elif args.oop_directive == "no-oop":
        return "no-oop"
    else:
        return "no specification"

def write_logfile_header():
    header=pyfiglet.figlet_format("NOMODECO.PY", font="starwars", width = 110, justify="center")
    logging.info(header)
    logging.info("")
    logging.info("*".center(110, "*"))
    logging.info("Kemal Önen, Dennis F. Dinu, Klaus R. Liedl".center(110))
    logging.info("Liedl Lab, University of Innsbruck".center(110))
    logging.info("*".center(110,"*"))
    logging.info("")

def write_logfile_oop_treatment():
    if oop_directive() == "oop":
        logging.info("You have included the consideration of out-of-plane angles. Note that out-of-plane angles tend to perform worse")
        logging.info("than the other internal coordinates in terms of computational cost and in the decomposition. However they can be")
        logging.info("useful for planar systems.")
    elif oop_directive() == "no-oop":
        logging.info("You did NOT include the consideration of out-of-plane angles. Note that out-of-plane angles tend to perform worse")
        logging.info("than the other internal coordinates in terms of computational cost and in the decomposition. However they can be")
        logging.info("useful for planar systems.")
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
    l,c = B.shape
    logging.info('The condition numbers of the B-Matrix are given by:')
    if l == c:
        B_Norm_1 = np.linalg.cond(B, 1)
        B_Norm_2 = np.linalg.cond(B, 2)
        B_Norm_inf = np.linalg.cond(B, np.inf)
        logging.info('Sum-Norm of B: %s', B_Norm_1)
        logging.info('Eucledian-Norm of B: %s', B_Norm_2)
        logging.info('Infinity-Norm of B: %s', B_Norm_inf)
    else:
        logging.warning("The condition number can not be calculated as B is not quadratic!")
    logging.info("")
    logging.info("")
    logging.info("Testing if the Internal Coordinate Set is complete ...")

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
