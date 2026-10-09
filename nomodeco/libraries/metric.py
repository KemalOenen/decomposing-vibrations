"""
Quality metrics for an internal coordinate set.

Every metric gets the same inputs and returns one number. Convention: higher is better
(metrics that should be small are returned with a minus sign), so the optimizer always
takes the maximum.

Inputs:
    contribution:  (n_ic, n_vib) contribution matrix in %, every column sums to 100
    ped:           (n_vib, n_ic, n_ic) PED matrix of every mode
    nu_intrinsic:  (n_ic,) intrinsic frequencies in cm-1
    nu_harmonic:   (n_vib,) harmonic frequencies in cm-1
""" 



import numpy as np
from scipy.optimize import linear_sum_assignment

def best_assignment(contribution):
    """  
    One-to-one assignment of ICs to vibrational modes with the largest total contribution
    happens via the Hungarian algorithm. Returns two index arrays: ics, modes
    """
    ics, modes = linear_sum_assignment(-contribution)
    return ics, modes


def kemalian(contribution, ped, nu_intrinsic, nu_harmonic):
    """ 
    Mean over all ICs of the largest contribution to each IC (without penalties)
    """
    return np.mean(np.max(contribution, axis=1))

def mode_centric(contribution, ped, nu_intrinsic, nu_harmonic):
    """  
    Mean over all modes of the largest IC contribution to each mode
    """
    return np.mean(np.max(contribution, axis=0))

def assignment(contribution, ped, nu_intrinsic, nu_harmonic):
    """  
    Mean over all ICs of the contribution of the assigned mode to each IC
    """
    ics, modes = best_assignment(contribution)
    n_vib = contribution.shape[1]
    return np.sum(contribution[ics, modes]) / n_vib

def freq_match(contribution, ped, nu_intrinsic, nu_harmonic):
    """  
    Minus the mean relative difference between intrinsic and harmonic freuqenciy (best is 0)
    """
    ics, modes = best_assignment(contribution)
    difference = np.abs(nu_intrinsic[ics] - nu_harmonic[modes]) / nu_harmonic[modes]
    return -np.mean(difference)

METRICS = {
    "kemalian": (kemalian, "Kemalian metric: mean of the largest contribution to each IC"),
    "mode_centric": (mode_centric, "Mode-centric metric: mean of the largest contribution to each mode"),
    "assignment": (assignment, "Assignment metric: mean of the contribution of the assigned mode to each IC"),
    "freq_match": (freq_match, "Frequency match metric: mean relative difference between intrinsic and harmonic frequencies"),
}