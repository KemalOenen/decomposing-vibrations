"""
Track C=O stretching bands across ORCA frequency calculations.

For every normal mode k, the C=O stretch amplitude is the change of the C-O distance,
    s_k = (u_O - u_C) . e_CO
with Cartesian (not mass-weighted) displacements u and the unit bond vector e_CO.
The mode with the largest share s_k^2 / sum(s^2) is reported as the C=O band.

Usage:
    python track_co.py file1.property.txt [file2.property.txt ...]
    (without arguments all *.property.txt files in this directory are used)
"""
import sys
import re
import glob
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../.."))
from nomodeco.libraries import orca_parser

# average atomic masses, as ORCA uses by default
MASSES = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
          "F": 18.998, "S": 32.06, "Cl": 35.45}
FREQ_CONV = 5140.4981  # sqrt(Eh / (bohr^2 amu)) -> cm^-1, as in nomodeco
CO_MAX_DIST = 1.30  # Angstrom; carbonyl C=O is ~1.21, C-O single ~1.34


def load(path):
    with open(path) as f:
        atoms = orca_parser.parse_xyz_from_inputfile(f)
    with open(path) as f:
        hessian = orca_parser.parse_cartesian_force_constants(f, len(atoms))
    symbols = [re.sub(r"\d+", "", atom.symbol) for atom in atoms]
    names = [atom.symbol for atom in atoms]
    coords = np.array([atom.coordinates for atom in atoms])
    return symbols, names, coords, hessian


def normal_modes(symbols, hessian):
    masses = np.repeat([MASSES[s] for s in symbols], 3)
    inv_sqrt_m = 1 / np.sqrt(masses)
    evals, evecs = np.linalg.eigh(hessian * np.outer(inv_sqrt_m, inv_sqrt_m))
    # drop the 6 translations/rotations (smallest |eigenvalue|), assumes a non-linear minimum
    vib = np.sort(np.argsort(np.abs(evals))[6:])
    freqs = np.sign(evals[vib]) * np.sqrt(np.abs(evals[vib])) * FREQ_CONV
    displacements = (evecs[:, vib] * inv_sqrt_m[:, None]).T  # Cartesian, one mode per row
    displacements /= np.linalg.norm(displacements, axis=1)[:, None]
    return freqs, displacements.reshape(len(vib), -1, 3)


def carbonyl_pairs(symbols, coords):
    pairs = []
    for o, s_o in enumerate(symbols):
        if s_o != "O":
            continue
        dists = np.linalg.norm(coords - coords[o], axis=1)
        neighbours = [i for i in range(len(symbols)) if i != o and dists[i] < 1.6]
        # carbonyl O: exactly one neighbour, a carbon at double-bond distance
        if len(neighbours) == 1 and symbols[neighbours[0]] == "C" and dists[neighbours[0]] < CO_MAX_DIST:
            pairs.append((neighbours[0], o))
    return pairs


def main(paths):
    for path in paths:
        symbols, names, coords, hessian = load(path)
        freqs, modes = normal_modes(symbols, hessian)
        print(f"\n{os.path.basename(path)}")
        for c, o in carbonyl_pairs(symbols, coords):
            e_co = (coords[o] - coords[c]) / np.linalg.norm(coords[o] - coords[c])
            s = np.einsum("kx,x->k", modes[:, o] - modes[:, c], e_co)
            share = s**2 / np.sum(s**2)
            order = np.argsort(share)[::-1]
            print(f"  {names[c]}={names[o]}  r = {np.linalg.norm(coords[o] - coords[c]):.4f} A")
            for k in order[:3]:
                # ORCA numbering: modes 0-5 are translations/rotations
                print(f"    mode {k + 6:3d}  {freqs[k]:8.1f} cm-1  C=O share {100 * share[k]:5.1f} %")


if __name__ == "__main__":
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "*.property.txt")))
    main(files)
