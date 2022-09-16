import itertools
import numpy as np
import nomodeco

"""""
step 1: generate all IC sets that have 3N-6 ICs and have all bonds included - NOT DONE
step 2: do a calculation run for every possibiliy - NOT DONE
step 3: compute the B-Matrix for all possibilites and check for completeness, if not complete remove - NOT DONE

step 4: order the decomposition according to new metric - NOT DONE
"""""

def generate_all_possible_sets(n_atoms, idof, bonds, angles, linear_angles, out_of_plane, dihedrals):
    num_bonds = len(bonds)
    num_angles = len(angles)
    num_linear_angles = len(linear_angles)
    num_out_of_plane = len(out_of_plane)
    num_dihedrals = len(dihedrals)

    all_non_bond_coordinates = angles + linear_angles + out_of_plane + dihedrals
    ic_sets = [] + bonds
    for ic_subsets in itertools.combinations(all_non_bond_coordinates, idof - num_bonds):
        ic_sets.append(ic_subsets)
        print(ic_sets, type(ic_sets))