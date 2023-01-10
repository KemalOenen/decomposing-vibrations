import pandas as pd
import numpy as np
import re
import pprint
from typing import NamedTuple
from collections import Counter
import argparse

class Atom(NamedTuple):
    symbol: str
    coordinates: tuple

#TODO: once integration in nomodeco.py -> remove get_args()

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("log")
    parser.add_argument("fchk")
    args = parser.parse_args()
    return args

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

def parse_xyz(inputfile):
    # parse names:
    names = []
    in_name_section = False
    for line in inputfile.split('\n'):
        if 'Final structure in terms of initial Z-matrix:' in line:
            in_name_section = True
            continue
        elif in_name_section and 'Variables:' in line:
            in_name_section = False
            continue
        elif in_name_section and line.strip():
            names.append(line[1])
    names = numerate_strings(names)
    # parse cartesian coordinates:
    coord_pattern = re.compile(r'\s*\d+\s+\d+\s+\w+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
    coords = []
    in_coord_section = False 
    skip_lines = 0
    for line in inputfile.split('\n'):
        if 'Standard orientation:' in line:
            in_coord_section = True
            coords = []
            skip_lines = 4
            continue
        elif skip_lines > 0:
            skip_lines -= 1
            continue
        elif in_coord_section and '---' in line:
            in_coord_section = False
            continue
        elif in_coord_section and line.strip():
            m = coord_pattern.match(line)
            if m:
                coords.append([m.group(1), float(m.group(2)), float(m.group(3))])
    return [Atom(name, tuple(coordinate)) for name, coordinate in zip(names, coords)]

#TODO: hessian information
def parse_cartesian_force_constants(inputfile, n_atoms):
    force_constants = []
    matrix_size = 3*n_atoms
    print((matrix_size**2)/5)
    in_force_const_section = False
    for line in inputfile.split('\n'):
        print("This is one line")
        if 'Cartesian Force Constants' in line:
            in_force_const_section = True
            continue
        elif in_force_const_section and 'Nonadiabatic coupling' in line:
            in_force_const_section = False
            continue
        elif in_force_const_section and line.strip():
            print("I am here!")
            force_constants = np.fromfile(inputfile,np.float64, matrix_size**2).reshape((matrix_size, matrix_size))
        return force_constants

args = get_args()
with open(args.log) as f:
    inputfile = f.read()
    atoms = parse_xyz(inputfile)
    n_atoms = len(atoms)

with open(args.fchk, 'rb') as f:
    inputfile = f.read()
    Cartesian_F_Matrix = parse_cartesian_force_constants(inputfile, n_atoms)

pprint.pprint(atoms)
print("n_atoms = ", n_atoms)
pprint.pprint(Cartesian_F_Matrix)
