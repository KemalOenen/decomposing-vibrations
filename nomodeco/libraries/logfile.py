import numpy as np
import pyfiglet
import logging

formatter = logging.Formatter('%(message)s')
from pprint import pformat



# Helper Function for submolecule treatment

def getList(dictionary):
    list_keys = []
    for key in dictionary.keys():
        list_keys.append(key)
    return list_keys

def setup_logger(name, log_file, level=logging.DEBUG):
    """To setup multiple loggers"""

    handler = logging.FileHandler(log_file, mode='a')
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    # do not print in stdout
    logger.propagate = False
    return logger


def create_filename_log(old_filename):
    pieces = old_filename.split(".")
    return "_".join([pieces[0], "nomodeco"]) + ".log"


def create_filename_out(old_filename):
    pieces = old_filename.split(".")
    return "_".join([pieces[0], "nomodeco"]) + ".out"


def write_logfile_header(logger):
    header = pyfiglet.figlet_format("NOMODECO.PY", font="starwars", width=110, justify="center")
    logger.info(header)
    logger.info("")
    logger.info("*".center(110, "*"))
    logger.info("Kemal Önen, Lukas Meinschad, Dennis F. Dinu, Klaus R. Liedl".center(110))
    logger.info("Liedl Lab, University of Innsbruck".center(110))
    logger.info("*".center(110, "*"))
    logger.info("")


def write_logfile_oop_treatment(logger, oop_directive, planar_subunit_list):
    if oop_directive == "yes":
        logger.info(
            "Planar system detected. Therefore out-of-plane angles are used in your analysis. Note that out-of-plane "
            "angles tend to perform worse")
        logger.info(
            "than the other internal coordinates in terms of computational cost and in the decomposition. However "
            "they can be")
        logger.info("useful for planar systems.")
    elif oop_directive == "no" and len(planar_subunit_list) == 0:
        logger.info(
            "Non-planar system detected. Therefore out-of-plane angles are NOT used in your analysis. Note that "
            "out-of-plane angles tend to perform worse")
        logger.info(
            "than the other internal coordinates in terms of computational cost and in the decomposition. However "
            "they can be")
        logger.info("useful for planar systems.")
    elif oop_directive == "no" and len(planar_subunit_list) != 0:
        logger.info(
            "Planar submolecule(s) detected. Out-of-plane angles will be used for the central atoms of the planar "
            "subunits ")
    logger.info("")


def write_logfile_symmetry_treatment(logger, specification, point_group_sch):
    logger.info("The detected point group is %s", point_group_sch)
    logger.info("The following equivalent atoms were found: %s ", pformat(specification["equivalent_atoms"], width=110, compact=True))
    logger.info("")

def write_hydrogen_bond_information(logger, h_bonds, acc_don_bonds, h_bond_angles, acc_don_angles,
        h_bond_linear_angles, acc_don_linear_angles, h_bond_dihedrals, acc_don_dihedrals, h_bond_oop, acc_don_oop):
    logger.info("Generation of intermolecular internal coordinates".center(110,"-"))
    logger.info("")
    logger.info("The following intermolecular internal coordinates (bonds, (in-plane) angles," +  
    " linear valance angles, out-of-plane angles and torsions) were generated:")
    logger.info(f"Hydrogen bonds: {h_bonds}")
    logger.info(f"Acceptor donor bonds: {acc_don_bonds}")
    logger.info(f"Hydrogen bond in-plane angles: {h_bond_angles}")
    logger.info(f"Acceptor donor in-plane angles: {acc_don_angles}")
    logger.info(f"Hydrogen bond bond linear valence-angles: {h_bond_linear_angles}")
    logger.info(f"Acceptor donor bond linear valence-angles: {acc_don_linear_angles}")
    logger.info(f"Hydrogen bond dihedrals: {h_bond_dihedrals}")
    logger.info(f"Acceptor donor bond dihedrals: {acc_don_dihedrals}")
    logger.info(f"Hydrogen bond out-of-plane angles: {h_bond_oop}")
    logger.info(f"Acceptor donor bond out-of-plane angles: {acc_don_oop}")

def write_logfile_submolecule_treatment(logger,submolecule_atom, connected_components):
    logger.info("Generation of internal coordinate sets for submolecules".center(110,"-"))
    logger.info("")
    logger.info(f"The following covalent submolecules where found in the structure: {connected_components}")

    for i, subdict in enumerate(submolecule_atom):
        logger.info(f"Generated IC sets for submolecule {i}: {max(getList(subdict)) + 1}")


def write_logfile_generated_IC(logger, bonds, angles, linear_angles, out_of_plane, dihedrals, idof):
    logger.info(
        "The following primitive internals (bonds, (in-plane) angles, linear valence angles, out-of-plane angles and "
        "torsions) were generated:".center(
            110))
    logger.info("bonds: %s", pformat(bonds, width=110, compact=True))
    logger.info("in-plane angles: %s", pformat(angles, width=110, compact=True))
    logger.info("linear valence-angles: %s",  pformat(linear_angles, width=110, compact=True))
    logger.info("out-of-plane angles: %s", pformat(out_of_plane, width=110, compact=True))
    logger.info("dihedrals: %s", pformat(dihedrals, width=110, compact=True))
    logger.info("")
    logger.info("%s internal coordinates are at least needed for the decomposition scheme.".center(50), idof)
    logger.info("")

def write_logfile_nan_freq(logger):
    logger.info("")
    logger.error('Imaginary frequency values were detected for the intrinsic frequencies.')
    logger.info("The computation will be stopped")
    logger.info("")


def write_logfile_not_complete_sets(logger):
    logger.info("")
    logger.error('No! The double transformed f-matrix is NOT the same as in the input.')
    logger.info("")
    logger.error("This set is not complete and will not be computed!")
    logger.info("")


def write_logfile_nan_freq_debug(logger):
    logger.info("")
    logger.error('Imaginary frequency values were detected for the intrinsic frequencies.')
    logger.info("In debug mode intrinsic frequencies are set to 0")


def write_logfile_information_results(logger, n_internals, red, bonds, angles, linear_angles, out_of_plane, dihedrals):
    logger.info("")
    logger.info(" Initialization of an internal coordinate set ".center(110, "-"))
    logger.info("")
    logger.info("The following %s Internal Coordinates are used in your analysis:".center(50), n_internals)
    logger.info("bonds: %s", pformat(bonds, width=110, compact=True))
    logger.info("in-plane angles: %s", pformat(angles, width=110, compact=True))
    logger.info("linear valence-angles: %s", pformat(linear_angles, width=110, compact=True))
    logger.info("out-of-plane angles: %s", pformat(out_of_plane, width=110, compact=True))
    logger.info("dihedrals: %s", pformat(dihedrals, width=110, compact=True))
    logger.info("")
    if red == 1:
        logging.info('There is %s redundant internal coordinates used.', red)
    else:
        logging.info('There are %s redundant internal coordinates used.', red)
    logger.info("")

def write_logfile_final_intermolecular_ics(logger,h_bonds,acc_don, h_bond_angles,
        acc_don_angles, h_bond_linear_angles, acc_don_linear_angles,h_bond_dihedrals, acc_don_dihedrals, h_bond_oop,
        acc_don_oop):
    logger.info("")
    logger.info("The final set contains the following intermolecular internal coordinates:")
    logger.info(f"Hydrogen bonds: {h_bonds}")
    logger.info(f"Acceptor donor bonds: {acc_don}")
    logger.info(f"Hydrogen bond in-plane angles: {h_bond_angles}")
    logger.info(f"Acceptor donor in-plane angles: {acc_don_angles}")
    logger.info(f"Hydrogen bond linear valence-angles: {h_bond_linear_angles}")
    logger.info(f"Acceptor bond linear valence-angles: {acc_don_linear_angles}")
    logger.info(f"Hydrogen bond dihedrals: {h_bond_dihedrals}")
    logger.info(f"Acceptor bond dihedrals: {acc_don_dihedrals}")
    logger.info(f"Hydrogen bond out-of-plane angles: {h_bond_oop}")
    logger.info(f"Acceptor bond out-of-plane angles: {acc_don_oop}")

def write_logfile_updatedICs_cyclic(logger, bonds, angles, linear_angles, out_of_plane, dihedrals, specification):
    logger.info("")
    logger.info("Due to having a cyclic molecule, %s bond(s) is/are removed in the analysis", specification["mu"])
    logger.info("")
    logger.info("The following updated Internal Coordinates are used in your analysis:".center(50))
    logger.info("bonds: %s", pformat(bonds, width=110, compact=True))
    logger.info("in-plane angles: %s", pformat(angles, width=110, compact=True))
    logger.info("linear valence-angles: %s", pformat(linear_angles, width=110, compact=True))
    logger.info("out-of-plane angles: %s", pformat(out_of_plane, width=110, compact=True))
    logger.info("dihedrals: %s", pformat(dihedrals, width=110, compact=True))
    logger.info("")


def write_logfile_updatedICs_linunit(logger, out_of_plane, dihedrals):
    logger.info("")
    logger.info("Due to having a linear subunit, out-of-plane / dihedral angles were removed")
    logger.info("")
    logger.info("The following updated Internal Coordinates are used in your analysis:".center(50))
    logger.info("out-of-plane angles: %s", pformat(out_of_plane, width=110, compact=True))
    logger.info("dihedrals: %s", pformat(dihedrals, width=110, compact=True))
    logger.info("")


def write_logfile_results(logger, Results, DiagonalElementsPED, ContributionTable, sum_check_VED):
    logger.info(" Results of the Decomposition Scheme ".center(90, "-"))
    logger.info("")
    logger.info("1.) Intrinsic Frequencies for all normal-coordinate frequencies and")
    logger.info("    vibrational distribution matrix (exact approximation) at the following")
    logger.info("    harmonic frequencies (all frequency values in cm-1):")
    logger.info("")
    logger.info(Results.to_string(index=False))
    logger.info("")
    logger.info("Check if matrix is normalized: %s", sum_check_VED)
    logger.info("")
    logger.info("2.) Intrinsic Frequencies for all normal-coordinate frequencies and")
    logger.info("    diagonal elements of the PED matrix at the following")
    logger.info("    harmonic frequencies (all frequency values in cm-1):")
    logger.info("")
    logger.info(DiagonalElementsPED.to_string(index=False))
    logger.info("")
    logger.info("3.) Intrinsic Frequencies for all normal-coordinate frequencies and")
    logger.info("    contribution table calculated from the diagonal elements of the PED at the following")
    logger.info("    harmonic frequencies (all frequency values in cm-1):")
    logger.info("")
    logger.info(ContributionTable.to_string(index=False))
    logger.info("")


def write_logfile_extended_results(logger, PED, KED, TED, sum_check_PED, sum_check_KED, sum_check_TED,
                                   harmonic_frequency):
    logger.info(" Extended Results: Exact PED, KED and TED analysis per mode ".center(70, "-"))
    logger.info("")
    logger.info("Normal mode at the following harmonic frequency: %s cm-1", harmonic_frequency)
    logger.info("")
    logger.info("Potential energy distribution matrix (normalization check: %s)", np.round(sum_check_PED, 2))
    logger.info("")
    logger.info(PED.to_string(index=False))
    logger.info("")
    # TODO: explain in output, why KED does not perfectly sum up to 1
    logger.info("Kinetic energy distribution matrix (normalization check: %s)", np.round(sum_check_KED, 2))
    logger.info("")
    logger.info(KED.to_string(index=False))
    logger.info("")
    logger.info("Total energy distribution matrix (normalization check: %s)", np.round(sum_check_TED, 2))
    logger.info("")
    logger.info(TED.to_string(index=False))
    logger.info("")


def call_shutdown():
    print("Done!")
    logging.shutdown()
