import numpy as np
import pandas as pd
from . import icsel
from . import bmatrix
from . import logfile
import os
import itertools

# import concurrent futures for multiprocessing

import concurrent.futures

def process_ic_set(num_of_set,ic_dict, args, idof, reciprocal_massmatrix, reciprocal_square_massmatrix, rottra,
                                CartesianF_Matrix, atoms, symmetric_coordinates, L, intfreq_penalty, intfc_penalty):

        if args.log:
            with open(args.output) as inputfile:
                 outputfile = logfile.create_filename_log(inputfile.name)
        #if os.path.exists(outputfile):
        #    os.remove(outputfile)
        log = logfile.setup_logger('logfile', outputfile)
        logfile.write_logfile_header(log)

        
        bonds = ic_dict[num_of_set]["bonds"]
        angles = ic_dict[num_of_set]["angles"]
        linear_angles = ic_dict[num_of_set]["linear valence angles"]
        out_of_plane = ic_dict[num_of_set]["out of plane angles"]
        dihedrals = ic_dict[num_of_set]["dihedrals"]

        n_internals = len(bonds) + len(angles) + len(linear_angles) + len(out_of_plane) + len(dihedrals)
        red = n_internals - idof

        # Augmenting the B-Matrix with rottra, calculating 
        # and printing the final B-Matrix

        B = np.concatenate((bmatrix.b_matrix(atoms, bonds, angles, linear_angles, out_of_plane, dihedrals, idof),
                            np.transpose(rottra)), axis=0)

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
 

        # Calculation of the mass-weighted normal modes in Cartesian Coordinates

        l = reciprocal_square_massmatrix @ L

        # Calculation of the mass-weighted normal modes in Internal Coordinates

        D = B @ l

        eigenvalues = np.transpose(D) @ InternalF_Matrix @ D
        eigenvalues = np.diag(eigenvalues)

        num_rottra = 3 * len(atoms) - idof

        nu = np.zeros(n_internals)
        for n in range(0, n_internals):
            for m in range(0, n_internals):
                for i in range(0, n_internals - red):
                    k = i + num_rottra
                    nu[n] += D[m][k] * InternalF_Matrix[m][n] * D[n][k]

        # Calculation of the Vibrational Density Matrices / PED, KED and TED matrices
        P = np.zeros((n_internals - red, n_internals + num_rottra, n_internals + num_rottra))
        T = np.zeros((n_internals - red, n_internals + num_rottra, n_internals + num_rottra))
        E = np.zeros((n_internals - red, n_internals + num_rottra, n_internals + num_rottra))

        for i in range(0, n_internals - red):
            for m in range(0, n_internals + num_rottra):
                for n in range(0, n_internals + num_rottra):
                    k = i + num_rottra
                    P[i][m][n] = D[m][k] * InternalF_Matrix[m][n] * D[n][k] / eigenvalues[k]  # PED
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

        sum_check_VED = 0
        ved_matrix = np.zeros((n_internals - red, n_internals + num_rottra))
        for i in range(0, n_internals - red):
            for m in range(0, n_internals + num_rottra):
                for n in range(0, n_internals + num_rottra):
                    ved_matrix[i][m] += P[i][m][n]
                sum_check_VED += ved_matrix[i][m]

        sum_check_VED = np.around(sum_check_VED / (n_internals - red), 2)
        ved_matrix = np.transpose(ved_matrix)
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
            contribution_matrix[:, i] = ((Diag_elements[:, i] / sum_diag[i]) * 100).astype(float)

        nu = np.zeros(n_internals)
        for n in range(0, n_internals):
            for m in range(0, n_internals):
                for i in range(0, n_internals - red):
                    k = i + num_rottra
                    nu[n] += D[m][k] * InternalF_Matrix[m][n] * D[n][k]

        nu_final = np.sqrt(nu) * 5140.4981

        if intfreq_penalty != 0:

            all_internals = bonds + angles + linear_angles + out_of_plane + dihedrals

            # check how often the intrinsic frequencies are the same for symmetric counterparts
            # get counter for asymmetry

            nu_dict = dict()
            for n in range(0, n_internals):
                nu_dict[all_internals[n]] = nu[n]

            counter_same_intrinsic_frequencies = 0
            counter_expected_symmetric_coordinates = 0
            for key1, value1 in nu_dict.items():
                for key2, value2 in nu_dict.items():
                    if key1 != key2 and np.isclose(value1, value2):
                        counter_same_intrinsic_frequencies += 1
            for key in nu_dict:
                if len(symmetric_coordinates[key]) > 1:
                    counter_expected_symmetric_coordinates += 1

            # no double counting
            counter_same_intrinsic_frequencies = (counter_same_intrinsic_frequencies // 2)
            counter_expected_symmetric_coordinates = (counter_expected_symmetric_coordinates // 2)

            counter = np.abs(counter_expected_symmetric_coordinates - counter_same_intrinsic_frequencies)
        else:
            counter = 0

        matrix = contribution_matrix

        metric_analysis_entry = icsel.Kemalian_metric(matrix,Diag_elements,counter,intfreq_penalty,
                                                       intfc_penalty,log)
        # return metric analysis and the num of the set for appending
        return num_of_set, metric_analysis_entry


def find_optimal_coordinate_set(ic_dict, args, idof, reciprocal_massmatrix, reciprocal_square_massmatrix, rottra,
                                CartesianF_Matrix, atoms, symmetric_coordinates, L, intfreq_penalty, intfc_penalty):
    metric_analysis = np.zeros(len(ic_dict))
    

    # Prepare the iterable of all arguments of the intial function
    all_args = ((num_of_set, ic_dict, args, idof, reciprocal_massmatrix, reciprocal_square_massmatrix, rottra,
                 CartesianF_Matrix, atoms, symmetric_coordinates, L, intfreq_penalty, intfc_penalty) 
                 for num_of_set in ic_dict.keys())

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_ic_set, *args_tuple) for args_tuple in all_args]        
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    
    for num_of_set, metric_analysis_entry in results:
        metric_analysis[num_of_set] = metric_analysis_entry
    
    # return the optimal set  with the metric valaue
    print("Optimal coordinate set has the following assigned metric value:",metric_analysis[np.argmax(metric_analysis)])
    return ic_dict[np.argmax(metric_analysis)]
