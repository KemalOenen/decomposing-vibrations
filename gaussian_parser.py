import pandas as pd
import numpy as np
import re
import pprint
from typing import NamedTuple

class Atom(NamedTuple):
    symbol: str
    coordinates: tuple


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

with open('/media/storage_4_240T/ddi/Kemal/C5H10/envelope/c5h10.log', 'r') as f:
    inputfile = f.read()

atoms = parse_xyz(inputfile)
#pprint.pprint(atoms)

#TODO: hessian information
