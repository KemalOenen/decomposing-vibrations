
"""
This module contains the get_sets method as well as symmetry functions for selecting internal coordinates
"""

import itertools
import logging
from collections import Counter
import numpy as np
from nomodeco.libraries import logfile
from nomodeco.libraries import topology


def Kemalian_metric(matrix, Diag_elements, counter, intfreq_penalty, intfc_penalty, args) -> float:
    """
    Given a matrix calculates the metric implemented in Nomodeco
    """
    # Exit Check if any value is negative < -1 we return 0 immediately
    if np.any(matrix < -1):
        return 0

    # Axis=1 when maximum of each row, 
    max_values = np.max(matrix, axis=1)

    # penalty for high fc values
    penalty1 = intfreq_penalty * counter

    # penalty for high fc values
    max_values_int_fc = np.max(Diag_elements, axis=1)
    excess_fc = max_values_int_fc - 1
    penality2 = np.sum(np.where(excess_fc > 0, excess_fc * 10 * intfc_penalty, 0))

    return float(np.mean(max_values) - penalty1 - penality2)


def Kemalian_metric_log(matrix, Diag_elements, counter, intfreq_penalty, intfc_penalty, log) -> float:
    """
    Optimized metric implemented in Nomodeco with logging using vectorized operations.
    """
    # Early exit check: Global minimum check is faster than row-by-row
    if np.any(matrix < -1):
        log.info("Negative values in energy distribution matrix found!")
        return 0.0

    # Find the maximum normal mode contribution (columns) for each IC (rows)
    max_values = np.max(matrix, axis=1)

    # Penalty 1: Static calculation
    penalty1 = intfreq_penalty * counter

    # Penalty 2: Vectorized replacement of the Python loop
    max_values_int_fc = np.max(Diag_elements, axis=1)
    excess_fc = max_values_int_fc - 1
    
    # Vectorized conditional math: ((val - 1) / 0.1) is equivalent to (excess * 10)
    penalty2 = float(np.sum(np.where(excess_fc > 0, excess_fc * 10 * intfc_penalty, 0.0)))

    # Calculate the base metric value
    mean_max = float(np.mean(max_values))

    # Log results with clean formatting
    log.info(
        "The following diagonalization parameter, penalty for asymmetric intrinsic frequencies "
        "and penalty for unphysical contributions has been determined: %s, %s, %s",
        np.around(mean_max, 2), penalty1, penalty2
    )

    return mean_max - penalty1 - penalty2


def are_two_elements_same(tup1, tup2) -> bool:
    """
    For two given tuples checks if they are the same and return True/False
    """
    return sum(x ==y for x, y in zip(tup1, tup2)) >=2

def get_different_elements(tup1, tup2) -> list:
    """
    Given two tuples calculates the difference in elements and outputs as list
    
    Conversion to list then using the symmetric difference operator ^
    """
    return list(set(tup1) ^ set(tup2))

def avoid_double_oop(test_oop, used_out_of_plane) -> bool:
    """  
    Returns True if test_oop is not a subset of any element in used_out_of_plane
    """
    test_set = set(test_oop)
    return not any(test_set.issubset(set(oop)) for oop in used_out_of_plane)

def remove_enumeration(atom_list) -> list:
    """
    Removes numeric characters from all string elements inside the tuples
    """
    remove_digits = str.maketrans('', '', '0123456789')
    return [
        tuple(element.translate(remove_digits) for element in tup)
        for tup in atom_list
    ]
    

def remove_enumeration_tuple(atom_tuple) -> tuple:
    """ 
    Removes numeric characters from all string elements inside a single tuple
    """
    remove_digits = str.maketrans('', '', '0123456789')
    return tuple(element.translate(remove_digits) for element in atom_tuple)
    


def check_in_nested_list(check_list, nested_list):
    """  
    Returns True if check_list is a subset of any element in nested_list
    Short Circuit instantly if a match is found
    """
    check_set = set(check_list)
    return any(check_set.issubset(set(single_list)) for single_list in nested_list)
    


def _make_equiv_sets(nested_equivalent_atoms):
    return [frozenset(s) for s in nested_equivalent_atoms]


def _atom_pair_match(a, b, equiv_sets):
    if a == b:
        return True
    pair = frozenset((a, b))
    return any(pair <= s for s in equiv_sets)


def _positions_match(test_ic, key_ic, equiv_sets):
    return all(_atom_pair_match(t, k, equiv_sets) for t, k in zip(test_ic, key_ic))


def all_atoms_can_be_superimposed_bond(test_bond, key_bond, nested_equivalent_atoms):
    return _positions_match(test_bond, key_bond, _make_equiv_sets(nested_equivalent_atoms))


def all_atoms_can_be_superimposed(test_angle, key_angle, nested_equivalent_atoms):
    return _positions_match(test_angle, key_angle, _make_equiv_sets(nested_equivalent_atoms))


def all_atoms_can_be_superimposed_dihedral(test_dihedral, key_dihedral, nested_equivalent_atoms):
    return _positions_match(test_dihedral, key_dihedral, _make_equiv_sets(nested_equivalent_atoms))


def get_symm_bonds(bonds, specification):
    symmetric_bonds = {key: [] for key in Counter(bonds)}
    equiv_sets = _make_equiv_sets(specification["equivalent_atoms"])
    for bond in bonds:
        for key in symmetric_bonds:
            if _positions_match(bond, key, equiv_sets) or _positions_match((bond[1], bond[0]), key, equiv_sets):
                symmetric_bonds[key].append(bond)
    return symmetric_bonds


def get_bond_subsets(symmetric_bonds) -> list:
    seen, result = set(), []
    for val in symmetric_bonds.values():
        key = frozenset(val)
        if key not in seen:
            seen.add(key)
            result.append(val)
    return result


def get_symm_angles(angles, specification, equiv_sets=None):
    symmetric_angles = {key: [] for key in Counter(angles)}
    if equiv_sets is None:
        equiv_sets = _make_equiv_sets(specification["equivalent_atoms"])
    for ang in angles:
        for key in symmetric_angles:
            if _positions_match(ang, key, equiv_sets) or _positions_match((ang[2], ang[1], ang[0]), key, equiv_sets):
                symmetric_angles[key].append(ang)
    return symmetric_angles


def _unique_groups(symmetric_dict):
    """Deduplicate symmetric dict values in O(n) using frozenset identity."""
    seen, groups = set(), []
    for val in symmetric_dict.values():
        key = frozenset(val)
        if key not in seen:
            seen.add(key)
            groups.append(val)
    return groups


# Maximum number of subsets returned from a single _flat_subsets_of_size call.
# Without this, no-symmetry molecules hit C(N, k) which can be billions.
_MAX_SUBSETS = 10_000


def _flat_subsets_of_size(groups, target):
    """Return flattened combinations of groups whose total length == target.

    Key optimisation: tighten the outer loop bounds from the group sizes so
    that i values that can never sum to target are skipped entirely.

    For no-symmetry molecules every group has size 1, so min_i = max_i = target
    and the function jumps directly to the one valid i, avoiding the
    ~sum(C(N,i) for i<target) wasted iterations that caused the hang.
    """
    # Choosing zero coordinates is one valid (empty) selection, e.g. no dihedrals in H2O
    if target == 0:
        return [[]]
    sizes = [len(g) for g in groups]
    if not sizes:
        return []
    min_size = min(sizes)
    max_size = max(sizes)
    # Smallest number of groups that could sum to target: ceil(target / max_size)
    min_i = -(-target // max_size)          # ceiling division without math.ceil
    # Largest number of groups that could sum to target: floor(target / min_size)
    max_i = min(target // min_size, len(groups))
    result = []
    for i in range(min_i, max_i + 1):
        for idx in itertools.combinations(range(len(groups)), i):
            if sum(sizes[j] for j in idx) == target:
                result.append([item for j in idx for item in groups[j]])
                if len(result) >= _MAX_SUBSETS:
                    logfile.search_log.warning(
                        "subset enumeration stopped at %s subsets of size %s (from %s groups)",
                        f"{_MAX_SUBSETS:,}", target, len(groups),
                    )
                    return result
    return result


def get_angle_subsets(symmetric_angles, num_bonds, num_angles, idof, n_phi) -> list:
    groups = _unique_groups(symmetric_angles)
    redundancy_msgs = [
        None,
        "In order to obtain symmetry in the angles and hence intrinsic frequencies, inclusion of 1 redundant angle coordinate will be attempted",
        "In order to obtain symmetry in the angles and hence intrinsic frequencies, inclusion of 2 redundant angle coordinates will be attempted",
    ]
    for extra in range(3):
        if extra:
            logging.info(redundancy_msgs[extra])
        angles = _flat_subsets_of_size(groups, n_phi + extra)
        if angles:
            return angles
    return []


def get_symm_dihedrals(dihedrals, specification, equiv_sets=None):
    symmetric_dihedrals = {key: [] for key in Counter(dihedrals)}
    if equiv_sets is None:
        equiv_sets = _make_equiv_sets(specification["equivalent_atoms"])
    for dihedral in dihedrals:
        rev = (dihedral[3], dihedral[2], dihedral[1], dihedral[0])
        for key in symmetric_dihedrals:
            if _positions_match(dihedral, key, equiv_sets) or _positions_match(rev, key, equiv_sets):
                symmetric_dihedrals[key].append(dihedral)
    return symmetric_dihedrals


def get_oop_subsets(out_of_plane, n_gamma):
    if n_gamma == 0:
        return [[]]
    # Group by central atom so we never pick two OOPs with the same centre.
    # Then choose n_gamma distinct central-atom groups and one OOP from each.
    by_central = {}
    for oop in out_of_plane:
        by_central.setdefault(oop[0], []).append(oop)
    groups = list(by_central.values())
    if len(groups) < n_gamma:
        return []
    result = []
    for idx_combo in itertools.combinations(range(len(groups)), n_gamma):
        for oops in itertools.product(*[groups[i] for i in idx_combo]):
            result.append(list(oops))
    return result


def get_dihedral_subsets(symmetric_dihedrals, num_bonds, num_angles, idof, n_tau) -> list:
    groups = _unique_groups(symmetric_dihedrals)
    for extra, msg in enumerate([
        None,
        "In order to obtain symmetry in the dihedrals and hence intrinsic frequencies, inclusion of 1 redundant dihedral coordinate will be attempted",
    ]):
        if extra:
            logging.info(msg)
        dihedrals = _flat_subsets_of_size(groups, n_tau + extra)
        if dihedrals:
            return dihedrals
    return []


def test_completeness(CartesianF_Matrix, B, B_inv, InternalF_Matrix) -> bool:
    return bool(np.allclose(np.transpose(B) @ InternalF_Matrix @ B, CartesianF_Matrix))


def check_evalue_f_matrix(reciprocal_square_massmatrix, B, B_inv, InternalF_Matrix):
    CartesianF_Matrix_check = np.transpose(B) @ InternalF_Matrix @ B
    evalue, evect = np.linalg.eigh(
        np.transpose(reciprocal_square_massmatrix) @ CartesianF_Matrix_check @ reciprocal_square_massmatrix)
    return evalue


def number_terminal_bonds(mult_list):
    return sum(1 for _, mult in mult_list if mult == 1)


def not_same_central_atom(list_oop_angles) -> bool:
    central_atoms = set()
    not_same_central_atom = True
    for oop_angle in list_oop_angles:
        if oop_angle[0] in central_atoms:
            not_same_central_atom = False
            break
        else:
            central_atoms.add(oop_angle[0])
    return not_same_central_atom


def matrix_norm(matrix, matrix_inv, p):
    return np.linalg.norm(matrix, p) * np.linalg.norm(matrix_inv, p)

def get_sets(idof, out, atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, specification) -> dict:
    """
    get_sets is the decicion tree of nomodeco where for a given molecular topology all possible sets are generated

    Attributes:
        idof:
            a integer with the vibrational degrees of freedom
        out:
            the output file of nomodeco
        atoms:
            a object of the molecule class
        bonds:
            a list of tuples containing bonds
        angles:
            a list of tuples containing angles
        linear_angles:
            a list of tuples containing linear angles
        out_of_plane:
            a list of tuples containing oop's
        dihedrals:
            a list of tuples containing dihedrals
        specficiation:
            the specification dictionary of nomodeco
    
    Returns:
        a dictionary of all the possible IC sets this is generated using the Topology module
    """ 
    ic_dict = dict()
    num_bonds = len(bonds)
    num_atoms = len(atoms)

    num_of_red = 6 * specification["mu"]
 
    # @decision tree: linear
    if specification["linearity"] == "fully linear" and specification["intermolecular"] == "no":
        ic_dict = topology.fully_linear_molecule(ic_dict, bonds, angles, linear_angles, out_of_plane, dihedrals)

    # @decision tree: planar, acyclic and no linear submolecules 
    if specification["planar"] == "yes" and specification["linearity"] == "not linear" and (
            num_of_red == 0) and specification["intermolecular"] == "no":
        ic_dict = topology.planar_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                             dihedrals, num_bonds, num_atoms,
                                                             number_terminal_bonds(specification["multiplicity"]),
                                                             specification)

    # @decision tree: planar, cyclic and no linear submolecules 
    if specification["planar"] == "yes" and specification["linearity"] == "not linear" and (
            num_of_red != 0) and specification["intermolecular"] == "no":
        ic_dict = topology.planar_cyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification)

    # @decision tree: general molecule, acyclic and no linear submolecules
    if specification["planar"] == "no" and specification["linearity"] == "not linear" and (
            num_of_red == 0) and specification["intermolecular"] == "no":
        ic_dict = topology.general_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                              dihedrals, num_bonds, num_atoms,
                                                              number_terminal_bonds(specification["multiplicity"]),
                                                              specification)

    # @decision tree: general molecule, cyclic and no linear submolecules
    if specification["planar"] == "no" and specification["linearity"] == "not linear" and (
            num_of_red != 0) and specification["intermolecular"] == "no":
        ic_dict = topology.general_cyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                             dihedrals, num_bonds, num_atoms, num_of_red,
                                                             number_terminal_bonds(specification["multiplicity"]),
                                                             specification)

    # @decision tree: planar, acyclic molecules with linear submolecules
    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (
            num_of_red == 0) and specification["intermolecular"] == "no":
        ic_dict = topology.planar_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                           dihedrals, num_bonds, num_atoms,
                                                           number_terminal_bonds(specification["multiplicity"]),
                                                           specification["length of linear submolecule(s) l"],
                                                           specification)

    # @decision tree: planar, cyclic molecules with linear submolecules
    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "no":
        ic_dict = topology.planar_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                          dihedrals, num_bonds, num_atoms,
                                                          number_terminal_bonds(specification["multiplicity"]),
                                                          specification["length of linear submolecule(s) l"],
                                                          specification)

        # @decision tree: general, acyclic molecule with linear submolecules
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found" and (
            num_of_red == 0) and specification["intermolecular"] == "no":
        ic_dict = topology.general_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification["length of linear submolecule(s) l"],
                                                            specification)

    # @decision tree: general, cyclic molecule with linear submolecules
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "no":
        ic_dict = topology.general_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                           dihedrals, num_bonds, num_atoms, num_of_red,
                                                           number_terminal_bonds(specification["multiplicity"]),
                                                           specification["length of linear submolecule(s) l"],
                                                           specification)


# For intermolecular complexes through determination of connectivity c we can manipulate our specification
# Therefore all the cases on the top get duplicated for the intermolecular systems

# General Molecules:
   
    # This is already done and working
    if specification["planar"] == "no" and specification["linearity"] == "not linear" and (num_of_red != 0) and specification["intermolecular"] == "yes":
       ic_dict = topology.intermolecular_general_cyclic_nolinsub(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms, num_of_red,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification)

    # This is already done and working 
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found" and (num_of_red == 0) and specification["intermolecular"] == "yes":
       ic_dict = topology.intermolecular_general_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification["length of linear submolecule(s) l"],
                                                            specification)
    # This is already done and woring
    if specification["planar"] == "no" and specification["linearity"] == "not linear" and (
            num_of_red == 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_general_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                              dihedrals, num_bonds, num_atoms,
                                                              number_terminal_bonds(specification["multiplicity"]),
                                                              specification)
    
    if specification["planar"] == "no" and specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_general_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                           dihedrals, num_bonds, num_atoms, num_of_red,
                                                           number_terminal_bonds(specification["multiplicity"]),
                                                           specification["length of linear submolecule(s) l"],
                                                           specification)
# Linear Molecules:
    
    # basically done just needs a example
    if specification["linearity"] == "fully linear" and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_fully_linear_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                              dihedrals, num_bonds, num_atoms,
                                                              number_terminal_bonds(specification["multiplicity"]),
                                                              specification)

# Planar Molecules


    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (
            num_of_red != 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_planar_cyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                          dihedrals, num_bonds, num_atoms,
                                                          number_terminal_bonds(specification["multiplicity"]),
                                                          specification["length of linear submolecule(s) l"],
                                                          specification)
 
    if specification["planar"] == "yes" and specification["linearity"] == "not linear" and (
            num_of_red != 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_planar_cyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                            dihedrals, num_bonds, num_atoms,
                                                            number_terminal_bonds(specification["multiplicity"]),
                                                            specification)
    # Here we need the specification of x+y = l-1 with kemal
    if specification["planar"] == "yes" and specification["linearity"] == "linear submolecules found" and (num_of_red == 0) and specification["intermolecular"] == "yes":
       ic_dict = topology.intermolecular_planar_acyclic_linunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                           dihedrals, num_bonds, num_atoms,
                                                           number_terminal_bonds(specification["multiplicity"]),
                                                           specification["length of linear submolecule(s) l"],
                                                           specification)


    if specification["planar"] == "yes" and specification["linearity"] == "not linear" and (num_of_red == 0) and specification["intermolecular"] == "yes":
        ic_dict = topology.intermolecular_planar_acyclic_nolinunit_molecule(ic_dict, out, idof, bonds, angles, linear_angles, out_of_plane,
                                                             dihedrals, num_bonds, num_atoms,
                                                             number_terminal_bonds(specification["multiplicity"]),
                                                             specification)

 


    print(len(ic_dict), "internal coordinate sets were generated.")
    print("The optimal coordinate set will be determined...")
    return ic_dict

