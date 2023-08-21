import numpy as np
import pandas as pd
import icsel
import bmatrix

def find_optimal_coordinate_set(ic_dict, idof, reciprocal_massmatrix, reciprocal_square_massmatrix, rottra, CartesianF_Matrix, atoms, L):
    metric_analysis = np.zeros(len(ic_dict))

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
                            np.transpose(rottra)),axis=0)

        # Calculating the G-Matrix

        G = B @ reciprocal_massmatrix @ np.transpose(B)
        e,K = np.linalg.eigh(G)

        # Sorting eigenvalues and eigenvectors (just for the case)
        # Sorting highest eigenvalue/eigenvector to lowest!

        idx = e.argsort()[::-1]   
        e = e[idx]
        K = K[:,idx]

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

        # remove not complete sets here
        # if you want the information where not completeness does occur
        # you can first call logfile.write_logfile_information_results
        if not icsel.test_completeness(CartesianF_Matrix, B, B_inv, InternalF_Matrix):
             continue
        
        # Calculation of the mass-weighted normal modes in Cartesian Coordinates

        l = reciprocal_square_massmatrix @ L

        # Calculation of the mass-weighted normal modes in Internal Coordinates

        D = B @ l
 
        eigenvalues = np.transpose(D) @ InternalF_Matrix @ D
        eigenvalues = np.diag(eigenvalues)
        # print("eigenvalues (from IC space):", eigenvalues)

        num_rottra = 3*len(atoms) - idof
 
        nu = np.zeros(n_internals) 
        for n in range(0,n_internals):
            for m in range(0,n_internals):
                for i in range(0,n_internals-red):
                    k = i + num_rottra
                    nu[n] += D[m][k] * InternalF_Matrix[m][n] * D[n][k]
        
        # if you want the information where imaginary freq. occur
        # uncomment below
        if np.any(nu < 0) == True:
            #logfile.write_logfile_nan_freq()
            continue
        
        # Calculation of the Vibrational Density Matrices / PED, KED and TED matrices
        P = np.zeros((n_internals-red,n_internals+num_rottra,n_internals+num_rottra))
        T = np.zeros((n_internals-red,n_internals+num_rottra,n_internals+num_rottra)) 
        E = np.zeros((n_internals-red,n_internals+num_rottra,n_internals+num_rottra))

        for i in range(0,n_internals-red):
            for m in range(0,n_internals + num_rottra):
                for n in range(0,n_internals + num_rottra):
                    k = i + num_rottra
                    P[i][m][n] = D[m][k]*InternalF_Matrix[m][n]*D[n][k] / eigenvalues[k] #PED
                    T[i][m][n] = D[m][k]*G_inv[m][n]*D[n][k]  #KED
                    E[i][m][n] = 0.5 *(T[i][m][n] + P[i][m][n]) #TED

        # check normalization
        sum_check_PED = np.zeros(n_internals)
        sum_check_KED = np.zeros(n_internals)
        sum_check_TED = np.zeros(n_internals)
        for i in range(0, n_internals):
            for m in range(0, n_internals + num_rottra):
                for n in range(0, n_internals + num_rottra):
                    sum_check_PED[i] += P[i][m][n] 
                    sum_check_KED[i] += T[i][m][n] 
                    sum_check_TED[i] += E[i][m][n] 

        # Summarized vibrational energy distribution matrix - can be calculated by either PED/KED/TED
        # rows are ICs, columns are harmonic frequencies!
        sum_check_VED = 0
        ved_matrix = np.zeros((n_internals - red, n_internals + num_rottra))
        for i in range(0,n_internals-red):
            for m in range(0, n_internals + num_rottra):
                for n in range (0,n_internals + num_rottra):
                    ved_matrix[i][m] += P[i][m][n]
                sum_check_VED += ved_matrix[i][m]
        
        sum_check_VED = np.around(sum_check_VED / (n_internals-red), 2)
        
        # currently: rows are harmonic modes and columns are ICs ==> need to transpose
        ved_matrix = np.transpose(ved_matrix)
        
        # remove the rottra
        ved_matrix = ved_matrix[0:n_internals, 0:n_internals]

        metric_analysis[num_of_set] = icsel.Kemalian_metric(ved_matrix)
    
    return ic_dict[np.argmax(metric_analysis)] 
