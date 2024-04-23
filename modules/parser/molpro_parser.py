import numpy as np
from typing import NamedTuple
from scipy import constants
from collections import Counter

class Atom(NamedTuple):
    symbol: str
    coordinates: tuple

def numerate_strings(string_list):
    string_counts = Counter(string_list)
    numeration = {}
    for string, count in string_counts.items():
        if count > 1:
            numeration[string] = 1
        else:
            numeration[string] = 0
    numerated_strings = []
    for string in string_list:
        if numeration[string] > 0:
            numerated_strings.append(f"{string}{numeration[string]}")
            numeration[string] += 1
        else:
            numerated_strings.append(string)
    return numerated_strings

#TODO: make a parse module
#TODO: adaption for other output formats like Gaussian, Orca
def parse_xyz_from_inputfile(inputfile):
    bohr = constants.value(u'Bohr radius')
    angstrom = constants.value(u'Angstrom star')
    BOHR_PER_ANGSTROM = angstrom/bohr 
    for line in inputfile:
        if line.strip().startswith('FREQUENCIES * CALCULATION OF NORMAL MODES'):
            break
    for _ in range(6):
        next(inputfile)
    names = []
    coordinates = []
    for line in inputfile:
        entries = line.strip().split()
        if len(entries) == 0:
            break
        names.append(entries[1])
        xyz = [float(f) / BOHR_PER_ANGSTROM for f in entries[3:]]
        coordinates.append(xyz)
    
    names = numerate_strings(names)
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


