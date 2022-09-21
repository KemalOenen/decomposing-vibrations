import itertools
import numpy as np
import pprint
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

    ic_dict = dict()

    # TODO: approach to update angles, linear_angles etc. and then feed to the dictionary
    for i in range(0, (3*n_atoms) - idof):
        k = 0
        for ic_subset in itertools.combinations(angles + linear_angles + out_of_plane + dihedrals, idof - num_bonds + i):
            ic_subset = list(ic_subset)

            ic_dict[k] = {
                "bonds" : bonds,
                "angles" : angles,
                "linear valence angles" : linear_angles,
                "out of plane angles" : out_of_plane,
                "dihedrals" : dihedrals
            }
            k +=1
    pprint.pprint(ic_dict)
    return ic_dict