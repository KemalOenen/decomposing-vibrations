import numpy as np
from typing import NamedTuple
from scipy import constants

class Atom(NamedTuple):
    symbol: str
    coordinates: tuple

#TODO: make a parse module
#TODO: adaption for other output formats like Gaussian, Orca
def parse_xyz_from_inputfile(inputfile):
    bohr = constants.value(u'Bohr radius')
    angstrom = constants.value(u'Angstrom star')
    BOHR_PER_ANGSTROM = angstrom/bohr 
    for line in inputfile:
        if line.strip().startswith('ATOMIC COORDINATES'):
            break
    next(inputfile)
    next(inputfile)
    next(inputfile)
    names = []
    for line in inputfile:
        entries = line.strip().split()
        if len(entries) == 0:
            break
        names.append(entries[1])
    for line in inputfile:
        if line.strip().startswith('FREQUENCIES * CALCULATION OF NORMAL MODES'):
            break
    for _ in range(6):
        next(inputfile)
    coordinates = []
    for line in inputfile:
        entries = line.strip().split()
        if len(entries) == 0:
            break
        xyz = [float(f) / BOHR_PER_ANGSTROM for f in entries[3:]]
        coordinates.append(xyz)

    return [Atom(name, tuple(coordinate)) for name, coordinate in zip(names, coordinates)]

#TODO: general parsing

def can_be_float(string):
    try:
        float(string)
        return True
    except:
        return False
    
def parse_Cartesian_F_Matrix_from_start_of_matrix(file):
    all_columns = []
    all_rows = []
    matrix = dict()
    for line in file:
        entries = line.split()
        if len(entries) == 0:
            break
        #print(line)
        if all(not can_be_float(e) for e in entries):
            columns = list(entries)
            all_columns.extend(columns)
        else:
            row = entries[0]
            if row not in all_rows:
                all_rows.append(row)
            for entry, col in zip(entries[1:], columns):
                matrix[row, col] = float(entry)
                matrix[col, row] = float(entry)
    out = np.array([
        [matrix[row, col] for col in all_columns]
        for row in all_rows
    ])
    return out
 
def parse_Cartesian_F_Matrix_from_inputfile(inputfile):
    for line in inputfile:
        if line.strip().startswith('Force Constants (Second Derivatives of the Energy) in [a.u.]'):
            break
    return parse_Cartesian_F_Matrix_from_start_of_matrix(inputfile)
