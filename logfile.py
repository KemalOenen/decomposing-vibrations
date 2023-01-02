import argparse
import numpy as np
import pandas as pd
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


def write_logfile_generated_IC(bonds, angles, linear_angles, out_of_plane, dihedrals):
    logging.info("The following primitive internals (bonds, (in-plane) angles, linear valence angles, out-of-plane angles and torsions) were generated:".center(110))
    logging.info("bonds: %s", bonds)
    logging.info("in-plane angles: %s", angles)
    logging.info("linear valence-angles: %s", linear_angles)
    logging.info("out-of-plane angles: %s", out_of_plane)
    logging.info("dihedrals: %s", dihedrals)
    logging.info("")

def write_logfile_results(full_eigenvalues, full_eigenvectors, U, R):
    logging.info("")
    logging.info('-----------Generating Delocalized Internal Coordinates-----------')
    logging.info("")
    logging.info("Eigenvalues of B*B(t)")
    logging.info("")
    for eigenvalue in full_eigenvalues:
        logging.info(eigenvalue)
    logging.info("")
    logging.info("Eigenvectors of B*B(t)")
    logging.info("")
    for eigenvector in full_eigenvectors:
        logging.info(pd.DataFrame(eigenvector))
    logging.info("")
    logging.info("Coordinates of the redundant space:")
    logging.info(pd.DataFrame(R))
    logging.info("")
    logging.info("Coordinates of the non-redundant space:")
    logging.info(pd.DataFrame(U))
    logging.shutdown()



