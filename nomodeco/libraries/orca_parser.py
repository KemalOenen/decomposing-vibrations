import pandas as pd
import numpy as np
import re
import pprint
from typing import NamedTuple
from collections import Counter
import argparse
from scipy import constants
import re
from mendeleev import element

from nomodeco.libraries.nomodeco_classes import Molecule


'''
Parser for Orca
'''


def numerate_strings(string_list) -> list:
    """
    For a given list of string enumerates all the strings based on occurence
    """
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

def parse_xyz_from_inputfile(inputfile) -> Molecule:
    """
    Parses the xyz coordinates out of the orca.property.txt file
    """

    # Read in lines
    lines = inputfile.readlines()
    elements_list = []
    bohr = constants.value(u'Bohr radius')
    angstrom = constants.value(u'Angstrom star')
    BOHR_PER_ANGSTROM = angstrom/bohr
     
    sections = []
    capturing = False
    current_section = []
    indices = []

    for i, line in enumerate(lines):
        # start collecting when CartesianCoordinates appears
        if line.strip().startswith("&CartesianCoordinates"):
            if capturing and current_section:
                sections.append((start_index,current_section))
            capturing = True
            start_index = i
            current_section = []
            indices.append(start_index)

            continue
        
        if line.strip().startswith("$End"):
            capturing = False
            if current_section:
                sections.append((start_index,current_section))
            current_section = []


        if capturing:
           current_section.append(line)

    # we now select the largest index --> which denotes the last geom opt
    max_index = max(indices)
    coordinates_raw = []
    for sec in sections:
        if sec[0] == max_index:
            for i , el in enumerate(sec):
                if i == 0:
                    continue # skip first element
                coordinates_raw.extend(el)

    # now we generate a atom list and a list of xyz_coordinates in Bohr
    atom_list = []
    xyz_coordinates = []
    for atom in coordinates_raw:
        atom_data = atom.split()
        # the first element is the atom
        atom_list.append(atom_data[0])
        xyz = []
        xyz.append(float(atom_data[1])/BOHR_PER_ANGSTROM)
        xyz.append(float(atom_data[2])/BOHR_PER_ANGSTROM)
        xyz.append(float(atom_data[3])/BOHR_PER_ANGSTROM)
        xyz_coordinates.append(xyz)
    
    names = numerate_strings(atom_list) # automatically numerate atoms if they occure more than once
    return [Molecule.Atom(name, tuple(coordinate)) for name, coordinate in zip(names,xyz_coordinates)]
    
def parse_cartesian_force_constants(inputfile,n_atoms) -> np.array:
    """
    Parses the second derivative matrix out of the orca.property.txt file
    """
    lines = inputfile.readlines()
    sections = []
    capturing = False
    current_section = []
    indices = []

    for line in lines:
        # start collecting when &Hessian Appears
        if line.strip().startswith("&HESSIAN"):
           capturing = True
           current_section = []
           
        if line.strip().startswith("&MODES"):
            capturing = False
            if current_section:
                sections.extend(current_section)
            current_section = []

        if capturing:
            current_section.append(line)
    
    # the hessian is given as a full matrix in orca property.txt
    # therefore we first split again into blocks
    block_dict = {}
    block_index = 0
    for line in sections:
        # if the line containes more then 6 spaces it determines a column index line
        if re.match(r'^ {6,}',line):
            block_index += 1
            block_dict.setdefault(block_index)
            block_dict[block_index] = []
        else:
            if block_index in block_dict:
                block_element = line.strip().split()
                block_element = block_element[1:] # slice the row index away
                # TODO QUICKFIX:
                if block_element:
                    block_dict[block_index].append(block_element)
    
    # now we generate our final hessian rows

    hessian_rows = []
     
    # we need to calculate the dimension
    n = 3 * n_atoms
    idx = 1
    for i, entries in enumerate(block_dict[idx]):
        if len(entries) == n: # a whole row is present
            hessian_rows.append(entries)
        elif len(entries) < n:
            idx_next = idx + 1
            # take the next block with index i
            entries.extend(block_dict[idx_next][i])
            if len(entries) < n:
                while len(entries) < n:
                   # now we jump again to the next block
                   idx_next += 1
                   # we aggain append our elements
                   entries.extend(block_dict[idx_next][i])
                hessian_rows.append(entries)
            else:
                hessian_rows.append(entries)
    
    
    # initialize hessian
    # because we already have full matrix we just return as a numpy array with floats
    return np.array(hessian_rows, dtype=float)
    
    
    
 #   for i in range(n):
 #       for j in range(i + 1):
 #           f_matrix[i,j] = float(hessian_rows[i][j])
 #   print(f_matrix)
 #   return f_matrix
