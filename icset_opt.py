import numpy as np
import pandas as pd
import icsel
import bmatrix
import logging
import logfile
import os


def find_optimal_coordinate_set(ic_dict, args, idof, reciprocal_massmatrix, reciprocal_square_massmatrix, rottra,
                                CartesianF_Matrix, atoms, symmetric_coordinates, L, intfreq_penalty, intfc_penalty):
    metric_analysis = np.zeros(len(ic_dict))
    if args.log:
        with open(args.output) as inputfile:
            outputfile = logfile.create_filename_log(inputfile.name)
        if os.path.exists(outputfile):
            os.remove(outputfile)
        log = logfile.setup_logger('logfile', outputfile)
        logfile.write_logfile_header(log)

    for num_of_set in ic_dict.keys():
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

        if args.log:
            logfile.write_logfile_information_results(log, n_internals, red, bonds, angles,
                                                      linear_angles, out_of_plane, dihedrals)

        # remove not complete sets here
        # if you want the information where not completeness does occur
        # you can first call logfile.write_logfile_information_results
        if not icsel.test_completeness(CartesianF_Matrix, B, B_inv, InternalF_Matrix):
            if args.log:
                logfile.write_logfile_not_complete_sets(log)
            continue

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

        # if you want the information where imaginary freq. occur
        # uncomment below
        if np.any(nu < 0) == True:
            logfile.write_logfile_nan_freq(log)
            continue

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
        if args.matrix_opt == "diag":
            matrix = Diag_elements
        if args.matrix_opt == "ved":
            matrix = ved_matrix

        if args.log:
            n_atoms = len(atoms)
            normal_coord_harmonic_frequencies = np.sqrt(eigenvalues[(3 * n_atoms - idof):3 * n_atoms]) * 5140.4981
            normal_coord_harmonic_frequencies = np.around(normal_coord_harmonic_frequencies, decimals=2)
            normal_coord_harmonic_frequencies_string = normal_coord_harmonic_frequencies.astype('str')

            all_internals = bonds + angles + linear_angles + out_of_plane + dihedrals

            all_internals_string = []
            for internal in all_internals:
                all_internals_string.append('(' + ', '.join(internal) + ')')

            Results = pd.DataFrame()
            Results['Internal Coordinate'] = all_internals_string
            Results['Intrinsic Frequencies'] = pd.DataFrame(nu_final).applymap("{0:.2f}".format)
            Results = Results.join(pd.DataFrame(ved_matrix).applymap("{0:.2f}".format))

            DiagonalElementsPED = pd.DataFrame()
            DiagonalElementsPED['Internal Coordinate'] = all_internals_string
            DiagonalElementsPED['Intrinsic Frequencies'] = pd.DataFrame(nu_final).applymap("{0:.2f}".format)
            DiagonalElementsPED = DiagonalElementsPED.join(pd.DataFrame(Diag_elements).applymap("{0:.2f}".format))

            ContributionTable = pd.DataFrame()
            ContributionTable['Internal Coordinate'] = all_internals_string
            ContributionTable['Intrinsic Frequencies'] = pd.DataFrame(nu_final).applymap("{0:.2f}".format)
            ContributionTable = ContributionTable.join(pd.DataFrame(contribution_matrix).applymap("{0:.2f}".format))

            columns = {}
            keys = range(3 * n_atoms - ((3 * n_atoms - idof)))
            for i in keys:
                columns[i] = normal_coord_harmonic_frequencies_string[i]

            Results = Results.rename(columns=columns)
            DiagonalElementsPED = DiagonalElementsPED.rename(columns=columns)
            ContributionTable = ContributionTable.rename(columns=columns)
            logfile.write_logfile_results(log, Results, DiagonalElementsPED, ContributionTable, sum_check_VED)

        if args.log:
            metric_analysis[num_of_set] = icsel.Kemalian_metric_log(matrix, Diag_elements, counter,
                                                                intfreq_penalty, intfc_penalty, log)
        else:
            metric_analysis[num_of_set] = icsel.Kemalian_metric(matrix, Diag_elements, counter,
                                                                intfreq_penalty, intfc_penalty, args)

    print("Optimal coordinate set has the following assigned metric value:",
          metric_analysis[np.argmax(metric_analysis)])
    return ic_dict[np.argmax(metric_analysis)]
