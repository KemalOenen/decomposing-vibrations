import itertools
import numpy as np
import pprint
import nomodeco

"""""
step 1: generate all IC sets that have 3N-6 ICs and have all bonds included - DONE
step 2: do a calculation run for every possibiliy - DONE
step 3: compute the B-Matrix for all possibilites and check for completeness, if not complete remove - DONE

step 4: order the decomposition according to new metric - NOT DONE
step 5: think of more methods to reduce the high dimensionaliy of new IC sets - NOT DONEs
"""""

def test_completeness(CartesianF_Matrix, B, B_inv, InternalF_Matrix) -> bool:
    CartesianF_Matrix_check = np.transpose(B) @ InternalF_Matrix @ B
    if (np.allclose(CartesianF_Matrix_check, CartesianF_Matrix)) == True:
        return True
    else:
        return False

def generate_all_possible_sets(n_atoms, idof, bonds, angles, linear_angles, out_of_plane, dihedrals):
    num_bonds = len(bonds)
    num_angles = len(angles)
    num_linear_angles = len(linear_angles)
    num_out_of_plane = len(out_of_plane)
    num_dihedrals = len(dihedrals)

    ic_dict = dict()

    # TODO: make approach more elegant?
    # TODO: remove calculations that have NAN intrinsic frequencies and also "double" out of plane angles
    # TODO: remove calculations that do not include useful angle combinations
    # TODO: remove calculations that have multiple dihedrals
    k = 0
    for i in range(0, (3*n_atoms) - idof):
        for ic_subset in itertools.combinations(angles + linear_angles + out_of_plane + dihedrals, idof - num_bonds + i):
            ic_subset = list(ic_subset)
            used_angles, used_linear_angles, used_out_of_plane, used_dihedrals = [], [], [], []
            for i in range(0, len(ic_subset)):
                if ic_subset[i] in angles:
                    used_angles.append(ic_subset[i])
                if ic_subset[i] in linear_angles:
                    used_linear_angles.append(ic_subset[i])
                if ic_subset[i] in out_of_plane:
                    used_out_of_plane.append(ic_subset[i])
                if ic_subset[i] in dihedrals:
                    used_dihedrals.append(ic_subset[i])
            ic_dict[k] = {
                "bonds" : bonds,
                "angles" : used_angles,
                "linear valence angles" : used_linear_angles,
                "out of plane angles" : used_out_of_plane,
                "dihedrals" : used_dihedrals
            }
            k +=1
            used_angles, used_linear_angles, used_out_of_plane, used_dihedrals = [], [], [], []
    pprint.pprint(ic_dict)
    return ic_dict