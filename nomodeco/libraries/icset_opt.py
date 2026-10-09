import numpy as np
import pandas as pd
from nomodeco.libraries import icsel
from nomodeco.libraries import bmatrix
from nomodeco.libraries import logfile
from nomodeco.libraries import metric
import os
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor

# a singular value of B counts towards the rank if it is above RANK_RTOL * sigma_max
# (the same criterion as bmatrix_metrics.condition_metrics)
RANK_RTOL = 1e-8
_RANK_BATCH = 1024


def _ic_universe(ic_dict):
    """
    All distinct ICs over all sets, in first-seen order, and their row in the B matrix built from
    them. A linear angle appears twice per set (both bends), so its rows are keyed by occurrence.
    Returns (universe, ic_to_row, la_occ_to_row), universe = (bonds, angles, linear angles,
    oops, dihedrals) for bmatrix.b_matrix.
    """
    bond_univ, angle_univ, la_univ, oop_univ, dih_univ = [], [], [], [], []
    seen_bonds, seen_angles, seen_oops, seen_dihs = set(), set(), set(), set()
    seen_la_keys: set = set()

    for s in ic_dict.values():
        for ic in s["bonds"]:
            ic = tuple(ic)
            if ic not in seen_bonds:
                seen_bonds.add(ic); bond_univ.append(ic)
        for ic in s["angles"]:
            ic = tuple(ic)
            if ic not in seen_angles:
                seen_angles.add(ic); angle_univ.append(ic)
        la_local: dict = {}
        for ic in s["linear valence angles"]:
            ic = tuple(ic)
            cnt = la_local.get(ic, 0)
            key = (ic, cnt)
            if key not in seen_la_keys:
                seen_la_keys.add(key); la_univ.append(ic)
            la_local[ic] = cnt + 1
        for ic in s["out of plane angles"]:
            ic = tuple(ic)
            if ic not in seen_oops:
                seen_oops.add(ic); oop_univ.append(ic)
        for ic in s["dihedrals"]:
            ic = tuple(ic)
            if ic not in seen_dihs:
                seen_dihs.add(ic); dih_univ.append(ic)

    ic_to_row: dict = {}
    _r = 0
    for ic in bond_univ:  ic_to_row[('b',   ic)] = _r; _r += 1
    for ic in angle_univ: ic_to_row[('a',   ic)] = _r; _r += 1
    la_occ_to_row: dict = {}
    la_occ_cnt: dict = {}
    for ic in la_univ:
        cnt = la_occ_cnt.get(ic, 0)
        la_occ_to_row[(ic, cnt)] = _r
        la_occ_cnt[ic] = cnt + 1
        _r += 1
    for ic in oop_univ: ic_to_row[('oop', ic)] = _r; _r += 1
    for ic in dih_univ: ic_to_row[('d',   ic)] = _r; _r += 1

    universe = (bond_univ, angle_univ, la_univ, oop_univ, dih_univ)
    return universe, ic_to_row, la_occ_to_row


def _row_index(ic_set, ic_to_row, la_occ_to_row) -> np.ndarray:
    """Rows of the universe B matrix for one set, in bonds, angles, linear, oop, dihedral order."""
    row_idx = [ic_to_row[('b', tuple(b))] for b in ic_set["bonds"]]
    row_idx += [ic_to_row[('a', tuple(a))] for a in ic_set["angles"]]
    _la_cnt: dict = {}
    for la in ic_set["linear valence angles"]:
        la = tuple(la)
        cnt = _la_cnt.get(la, 0)
        row_idx.append(la_occ_to_row[(la, cnt)])
        _la_cnt[la] = cnt + 1
    row_idx += [ic_to_row[('oop', tuple(o))] for o in ic_set["out of plane angles"]]
    row_idx += [ic_to_row[('d', tuple(d))] for d in ic_set["dihedrals"]]
    return np.array(row_idx, dtype=np.intp)


def _complete_keys(keys_and_rows, B_master, idof, workers=None) -> list:
    """
    Keys whose rows of B_master have rank idof, i.e. the set describes every internal motion.
    Sets with the same number of ICs are stacked into batches; each batch is one batched SVD,
    and the batches run in a thread pool (numpy's SVD releases the GIL).
    """
    workers = workers or min(os.cpu_count() or 1, 16)
    complete = []
    pending = defaultdict(list)
    in_flight = deque()

    def complete_in(batch):
        s = np.linalg.svd(B_master[np.stack([rows for _, rows in batch])], compute_uv=False)
        rank = (s > RANK_RTOL * s[:, :1]).sum(axis=1)
        return [key for (key, _), r in zip(batch, rank) if r >= idof]

    with ThreadPoolExecutor(workers) as pool:
        def submit(n):
            batch = pending.pop(n)
            if n < idof:
                return  # fewer ICs than internal degrees of freedom: never complete
            in_flight.append(pool.submit(complete_in, batch))
            # bounded: the stacked batches of a large search would not fit in memory at once
            while len(in_flight) > 2 * workers:
                complete.extend(in_flight.popleft().result())

        for key, rows in keys_and_rows:
            pending[len(rows)].append((key, rows))
            if len(pending[len(rows)]) >= _RANK_BATCH:
                submit(len(rows))
        for n in list(pending):
            submit(n)
        while in_flight:
            complete.extend(in_flight.popleft().result())
    return sorted(complete)


def complete_set_keys(ic_dict, atoms, idof, bend_refs=None) -> list:
    """Keys of the IC sets whose B matrix has full rank idof (sorted, generation order)."""
    if len(ic_dict) == 0:
        return []
    universe, ic_to_row, la_occ_to_row = _ic_universe(ic_dict)
    B_master = bmatrix.b_matrix(atoms, *universe, idof, bend_refs)
    return _complete_keys(
        ((k, _row_index(ic_dict[k], ic_to_row, la_occ_to_row)) for k in ic_dict.keys()),
        B_master, idof,
    )


def find_optimal_coordinate_set(ic_dict, args, idof, reciprocal_massmatrix, reciprocal_square_massmatrix, rottra,
                                CartesianF_Matrix, atoms, symmetric_coordinates, L, intfreq_penalty, intfc_penalty,
                                bend_refs=None) -> dict:
    """
    Returns a dictionary with the optimal coordinate set. For each entry in the ic_dict, the metric of Nomodeco gets calculated, then the set with the highest metric gets selected.

    Key optimisation: B-matrix rows, G-matrix blocks, H-matrix blocks, and D-matrix blocks are
    precomputed once for the full IC universe and then assembled per set via cheap row-indexing,
    replacing the per-set  B @ M_inv @ B^T  and  B_inv^T @ F @ B_inv  matrix multiplications.
    """
    metric_analysis = {}

    if len(ic_dict) == 0:
        return {"best_key": None, "metric": 0, "set": None}

    if args.log:
        if not args.gv == None:
            with open(args.gv[0]) as inputfile:
                outputfile = logfile.create_filename_log(inputfile.name)
        if not args.molpro == None:
            with open(args.molpro[0]) as inputfile:
                outputfile = logfile.create_filename_log(inputfile.name)
        if not args.orca == None:
            with open(args.orca[0]) as inputfile:
                outputfile = logfile.create_filename_log(inputfile.name)
        if args.pymolpro:
            with open(os.getenv('OUT_FILE_LINK')) as inputfile:
                outputfile = logfile.create_filename_log(inputfile.name)
        if os.path.exists(outputfile):
            os.remove(outputfile)
        log = logfile.setup_logger('logfile', outputfile)
        logfile.write_logfile_header(log)

    # ------------------------------------------------------------------
    # Build IC universe: collect all unique ICs across every set so we
    # can compute one B_master covering them all.
    # Linear angles must preserve occurrence order (same tuple appears
    # twice: once for first-plane, once for second-plane bending).
    # ------------------------------------------------------------------
    universe, ic_to_row, la_occ_to_row = _ic_universe(ic_dict)

    # One B-matrix call for all ICs (replaces N per-set calls)
    B_master = bmatrix.b_matrix(atoms, *universe, idof, bend_refs)

    # Rank check first: an incomplete set (rank B < idof) can not be the optimal set, so it is
    # dropped before any G/H/PED work
    set_rows = {k: _row_index(ic_dict[k], ic_to_row, la_occ_to_row) for k in ic_dict.keys()}
    complete_keys = _complete_keys(set_rows.items(), B_master, idof)
    n_incomplete = len(set_rows) - len(complete_keys)
    if n_incomplete:
        print(f"{n_incomplete:,} of {len(set_rows):,} IC sets are incomplete (rank B < {idof}) and are skipped")
    if args.log and n_incomplete:
        log.info("%s of %s IC sets are incomplete (rank of B < %s) and were skipped",
                 n_incomplete, len(set_rows), idof)

    # ------------------------------------------------------------------
    # Precompute shared blocks (computed once, indexed per set)
    #
    #   G = B_aug @ M_inv @ B_aug^T
    #   H = B_aug_Minv @ F_cart @ B_aug_Minv^T   (used as InternalF = G_inv @ H @ G_inv)
    #   D = B_aug @ l
    #
    # Subscripts: _u = universe ICs, _r = rottra rows
    # ------------------------------------------------------------------
    diag_m    = np.diag(reciprocal_massmatrix)   # (3N,)
    B_rottra  = rottra.T                          # (n_rottra, 3N)
    B_u_Minv  = B_master * diag_m                # (n_univ, 3N)
    B_r_Minv  = B_rottra * diag_m                # (n_rottra, 3N)

    G_uu = B_u_Minv @ B_master.T                 # (n_univ, n_univ)
    G_ur = B_u_Minv @ rottra                     # (n_univ, n_rottra)
    G_rr = B_r_Minv @ B_rottra.T                 # (n_rottra, n_rottra)

    _tmp_u = B_u_Minv @ CartesianF_Matrix        # (n_univ, 3N)
    H_uu   = _tmp_u @ B_u_Minv.T                 # (n_univ, n_univ)
    H_ur   = _tmp_u @ B_r_Minv.T                 # (n_univ, n_rottra)
    H_rr   = B_r_Minv @ CartesianF_Matrix @ B_r_Minv.T  # (n_rottra, n_rottra)

    l    = reciprocal_square_massmatrix @ L
    D_u  = B_master @ l                          # (n_univ, 3N)
    D_r  = B_rottra @ l                          # (n_rottra, 3N)

    n_rottra = rottra.shape[1]

    # ------------------------------------------------------------------
    # Main loop — per-set work is now cheap: index + eigh + G_inv
    # ------------------------------------------------------------------
    for num_of_set in complete_keys:
        bonds         = [tuple(b)   for b   in ic_dict[num_of_set]["bonds"]]
        angles        = [tuple(a)   for a   in ic_dict[num_of_set]["angles"]]
        linear_angles = [tuple(la)  for la  in ic_dict[num_of_set]["linear valence angles"]]
        out_of_plane  = [tuple(oop) for oop in ic_dict[num_of_set]["out of plane angles"]]
        dihedrals     = [tuple(d)   for d   in ic_dict[num_of_set]["dihedrals"]]

        n_internals = len(bonds) + len(angles) + len(linear_angles) + len(out_of_plane) + len(dihedrals)
        red = n_internals - idof

        row_idx = set_rows[num_of_set]

        # Assemble G_aug from precomputed blocks (no B @ M_inv @ B^T per set)
        G_11  = G_uu[np.ix_(row_idx, row_idx)]
        G_12  = G_ur[row_idx, :]
        G_aug = np.block([[G_11, G_12], [G_12.T, G_rr]])

        e, K = np.linalg.eigh(G_aug)
        idx_sort = e.argsort()[::-1]
        e = e[idx_sort]; K = K[:, idx_sort]

        # drop the red smallest eigenvalues (the last red after the descending sort);
        # np.delete(K, -red) removed only one column, wrong for red >= 2
        if red > 0:
            K = K[:, :-red]
            e = e[:-red]

        e = np.diag(e)
        try:
            G_inv = K @ np.linalg.inv(e) @ K.T
        except np.linalg.LinAlgError:
            G_inv = K @ np.linalg.pinv(e) @ K.T

        # InternalF = G_inv @ H_aug @ G_inv  (no B_inv^T @ F @ B_inv per set)
        H_11       = H_uu[np.ix_(row_idx, row_idx)]
        H_12       = H_ur[row_idx, :]
        H_aug      = np.block([[H_11, H_12], [H_12.T, H_rr]])
        InternalF_Matrix = G_inv @ H_aug @ G_inv

        # B_aug and B_inv only needed for the completeness check (cheap row slice)
        B_aug      = np.concatenate([B_master[row_idx], B_rottra], axis=0)
        B_Minv_aug = np.concatenate([B_u_Minv[row_idx], B_r_Minv], axis=0)
        B_inv      = B_Minv_aug.T @ G_inv

        if args.log:
            logfile.write_logfile_information_results(log, n_internals, red, bonds, angles,
                                                      linear_angles, out_of_plane, dihedrals)

        if not icsel.test_completeness(CartesianF_Matrix, B_aug, B_inv, InternalF_Matrix):
            if args.log:
                logfile.write_logfile_not_complete_sets(log)
            continue

        D = np.concatenate([D_u[row_idx], D_r], axis=0)

        eigenvalues = np.diag(D.T @ InternalF_Matrix @ D)

        num_rottra = n_rottra
        n_vib      = n_internals - red
        D_ni = D[:n_internals, num_rottra:num_rottra + n_vib]
        nu   = np.einsum('mi,mn,ni->n', D_ni, InternalF_Matrix[:n_internals, :n_internals], D_ni)
        if np.any(nu < 0):
            if args.log:
                logfile.write_logfile_nan_freq(log)
            continue

        ev_vib = eigenvalues[num_rottra:num_rottra + n_vib]
        D_vib  = D[:, num_rottra:num_rottra + n_vib].T
        outer  = D_vib[:, :, None] * D_vib[:, None, :]
        P = outer * InternalF_Matrix[None, :, :] / ev_vib[:, None, None]
        T = outer * G_inv[None, :, :]
        E = 0.5 * (P + T)

        ved_matrix     = P.sum(axis=2)
        sum_check_VED  = np.around(ved_matrix.sum() / n_vib, 2)
        ved_matrix     = ved_matrix.T[:n_internals, :n_internals]

        Diag_elements       = np.diagonal(P, axis1=1, axis2=2)[:, :n_internals].T
        sum_diag            = Diag_elements.sum(axis=0)
        contribution_matrix = Diag_elements / sum_diag * 100

        nu_final = np.sqrt(nu) * 5140.4981

        if intfreq_penalty != 0:
            all_internals = bonds + angles + linear_angles + out_of_plane + dihedrals
            nu_dict = {all_internals[n]: nu[n] for n in range(n_internals)}

            counter_same_intrinsic_frequencies = 0
            counter_expected_symmetric_coordinates = 0
            for key1, value1 in nu_dict.items():
                for key2, value2 in nu_dict.items():
                    if key1 != key2 and np.isclose(value1, value2):
                        counter_same_intrinsic_frequencies += 1
            for key in nu_dict:
                if len(symmetric_coordinates.get(key, ())) > 1:
                    counter_expected_symmetric_coordinates += 1

            counter_same_intrinsic_frequencies    = counter_same_intrinsic_frequencies // 2
            counter_expected_symmetric_coordinates = counter_expected_symmetric_coordinates // 2
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
            all_internals_string = ['(' + ', '.join(internal) + ')' for internal in all_internals]

            Results = pd.DataFrame()
            Results['Internal Coordinate'] = all_internals_string
            Results['Intrinsic Frequencies'] = pd.DataFrame(nu_final).map("{0:.2f}".format)
            Results = Results.join(pd.DataFrame(ved_matrix).map("{0:.2f}".format))

            DiagonalElementsPED = pd.DataFrame()
            DiagonalElementsPED['Internal Coordinate'] = all_internals_string
            DiagonalElementsPED['Intrinsic Frequencies'] = pd.DataFrame(nu_final).map("{0:.2f}".format)
            DiagonalElementsPED = DiagonalElementsPED.join(pd.DataFrame(Diag_elements).map("{0:.2f}".format))

            ContributionTable = pd.DataFrame()
            ContributionTable['Internal Coordinate'] = all_internals_string
            ContributionTable['Intrinsic Frequencies'] = pd.DataFrame(nu_final).map("{0:.2f}".format)
            ContributionTable = ContributionTable.join(pd.DataFrame(contribution_matrix).map("{0:.2f}".format))

            columns = {}
            for i in range(3 * n_atoms - (3 * n_atoms - idof)):
                columns[i] = normal_coord_harmonic_frequencies_string[i]

            Results = Results.rename(columns=columns)
            DiagonalElementsPED = DiagonalElementsPED.rename(columns=columns)
            ContributionTable = ContributionTable.rename(columns=columns)
            logfile.write_logfile_results(log, Results, DiagonalElementsPED, ContributionTable, sum_check_VED)

        if args.metric == "kemalian":
            # with the --matrix_opt choice and the penalties of --penalty1/--penalty2
            if args.log:
                metric_analysis[num_of_set] = icsel.Kemalian_metric_log(matrix, Diag_elements, counter,
                                                                         intfreq_penalty, intfc_penalty, log)
            else:
                metric_analysis[num_of_set] = icsel.Kemalian_metric(matrix, Diag_elements, counter,
                                                                     intfreq_penalty, intfc_penalty, args)
        else:
            nu_harmonic = np.sqrt(ev_vib) * 5140.4981
            metric_fn = metric.METRICS[args.metric][0]
            metric_analysis[num_of_set] = float(metric_fn(contribution_matrix,
                                                          P[:, :n_internals, :n_internals],
                                                          nu_final, nu_harmonic))
            if args.log:
                log.info("%s metric: %s", args.metric, np.around(metric_analysis[num_of_set], 4))

    if not metric_analysis:
        return {"best_key": None, "metric": 0, "set": None}

    best_key = max(metric_analysis, key=metric_analysis.get)
    best_metric_value = metric_analysis[best_key]
    print(f"Optimal coordinate set has the following assigned {args.metric} metric value:",
          metric_analysis[best_key])

    return {
        "best_key": best_key,
        "metric": best_metric_value,
        "set": ic_dict[best_key]
    }
