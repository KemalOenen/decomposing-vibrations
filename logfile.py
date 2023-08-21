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

def write_logfile_oop_treatment(oop_directive, planar_subunit_list):
    if oop_directive == "yes":
        logging.info("Planar system detected. Therefore out-of-plane angles are used in your analysis. Note that out-of-plane angles tend to perform worse")
        logging.info("than the other internal coordinates in terms of computational cost and in the decomposition. However they can be")
        logging.info("useful for planar systems.")
    elif oop_directive == "no" and len(planar_subunit_list) == 0:
        logging.info("Non-planar system detected. Therefore out-of-plane angles are NOT used in your analysis. Note that out-of-plane angles tend to perform worse")
        logging.info("than the other internal coordinates in terms of computational cost and in the decomposition. However they can be")
        logging.info("useful for planar systems.")
    elif oop_directive == "no" and len(planar_subunit_list) != 0:
            logging.info("Planar submolecule(s) detected. Out-of-plane angles will be used for the central atoms of the planar subunits ")
    logging.info("")

def write_logfile_symmetry_treatment(specification, point_group_sch):
    logging.info("The detected point group is %s", point_group_sch)
    logging.info("The following equivalent atoms were found: %s ", specification["equivalent_atoms"])
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

def write_logfile_nan_freq_debug():
    logging.info("")
    logging.error('Imaginary frequency values were detected for the intrinsic frequencies.')
    logging.info("In debug mode intrinsic frequencies are set to 0")

def write_logfile_information_results(n_internals, red, bonds, angles, linear_angles, out_of_plane, dihedrals):
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

def write_logfile_updatedICs_cyclic(bonds, angles, linear_angles, out_of_plane, dihedrals, specification):
    logging.info("")
    logging.info("Due to having a cyclic molecule, %s bond(s) is/are removed in the analysis", specification["mu"]) 
    logging.info("")
    logging.info("The following updated Internal Coordinates are used in your analysis:".center(50))
    logging.info("bonds: %s", bonds)
    logging.info("in-plane angles: %s", angles)
    logging.info("linear valence-angles: %s", linear_angles)
    logging.info("out-of-plane angles: %s", out_of_plane)
    logging.info("dihedrals: %s", dihedrals)
    logging.info("")

def write_logfile_updatedICs_linunit(out_of_plane, dihedrals):
    logging.info("")
    logging.info("Due to having a linear subunit, out-of-plane / dihedral angles were removed") 
    logging.info("")
    logging.info("The following updated Internal Coordinates are used in your analysis:".center(50))
    logging.info("out-of-plane angles: %s", out_of_plane)
    logging.info("dihedrals: %s", dihedrals)
    logging.info("")

def write_logfile_results(Results, sum_check_VED):
    logging.info(" Results of the Decomposition Scheme ".center(90, "-"))
    logging.info("")
    logging.info("1.) Intrinsic Frequencies for all normal-coordinate frequencies and")
    logging.info("    vibrational distribution matrix (exact approximation) at the following")
    logging.info("    harmonic frequencies (all frequency values in cm-1):")
    logging.info("")
    logging.info(Results.to_string(index=False))
    logging.info("")
    logging.info("Check if matrix is normalized: %s", sum_check_VED)
    logging.info("")

def write_logfile_extended_results(PED,KED,TED, sum_check_PED, sum_check_KED, sum_check_TED, harmonic_frequency):
    logging.info(" Extended Results: Exact PED, KED and TED analysis per mode ".center(70, "-"))
    logging.info("")
    logging.info("Normal mode at the followin harmonic frequency: %s cm-1", harmonic_frequency)
    logging.info("")
    logging.info("Potential energy distribution matrix (normalization check: %s)", np.round(sum_check_PED,2))
    logging.info("")
    logging.info(PED.to_string(index=False))
    logging.info("")
    #TODO: explain in output, why KED does not perfectly sum up to 1
    logging.info("Kinetic energy distribution matrix (normalization check: %s)", np.round(sum_check_KED,2))
    logging.info("")
    logging.info(KED.to_string(index=False))
    logging.info("")
    logging.info("Total energy distribution matrix (normalization check: %s)", np.round(sum_check_TED,2))
    logging.info("")
    logging.info(TED.to_string(index=False))
    logging.info("")

def call_shutdown():
    print("Done!")
    logging.shutdown()
