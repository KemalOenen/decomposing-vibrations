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

class Atom(NamedTuple):
    symbol: str
    coordinates: tuple


def numerate_strings(string_list) -> list:
    """
    Enumerates all strings in a list based on occurence
    
    Attributes:
        a list of strings representing atom symbols ["H","O","H"]
    
    Return:
        Enumerated list based on occurence
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

def parse_xyz(inputfile)-> Molecule():
    """
    Parses the xyz coordinates from the gaussian log file and converts into a object of molecule class

    Attributes:
        inputfile:
            the gaussian.log file as an inputfile 
    
    Return:
        A object of the Molecule class
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
        # Start collecting when **** Axes restored to original set ***** appears
        if line.strip().startswith("Input orientation"):
            if capturing and current_section:
                sections.append((start_index,current_section))
            capturing = True
            start_index = i
            current_section = []
            indices.append(start_index)

            continue
        
        if " Distance matrix" in line:
            capturing = False
            if current_section:
                sections.append((start_index,current_section))
            current_section = []


        if capturing:
           current_section.append(line)
    # now select the maximum index out of the indices list and build up the correct list
    
    max_index = max(indices)
    coordinates_raw = []
    for sec in sections:
        if sec[0] == max_index:
            for i, el in enumerate(sec):
                if i == 0:
                    continue # skip the first element
                coordinates_raw.extend(el)
    # remove the first four elements
    coordinates_raw = coordinates_raw[4:]
    # also remove the last element
    coordinates_raw = coordinates_raw[:-1]

    
    # now each line in this list corresponds to a coordinate and element
    # now the first element corresponds to the number of the atom, the second to an 
    atom_list =[]
    xyz_coordinates = []
    for atom in coordinates_raw:
        atom_data = atom.split()
        # we now build the atom list: take the second element retrieve it in the number dict and add the
        # first element of the atom_data list
        element_symbol = element(int(atom_data[1])).symbol
        atom_list.append(element_symbol)
        # the xyz data is in angstrom here
        xyz = []
        xyz.append(float(atom_data[3])) # x
        xyz.append(float(atom_data[4])) # y
        xyz.append(float(atom_data[5])) # z
        xyz_coordinates.append(xyz)

    atom_list = numerate_strings(atom_list)

    return [Molecule.Atom(name,tuple(coordinate)) for name, coordinate in zip(atom_list,xyz_coordinates)]



def parse_cartesian_force_constants(inputfile, n_atoms) -> np.array:
    """
    Parses the second derivative matrix out of the gaussian.log file. In Gaussian the second derivative matrix is given as lower triangular, Nomodeco.py uses
    the full matrix for the calculation

    Attributes:
        inputfile:
            the gaussian.log inputfile
        n_atoms: int
            a integer with the number of atoms in the structure
    """
    lines = inputfile.readlines()
    sections = []
    capturing = False
    current_section = []
    indices = []

    for i, line in enumerate(lines):
        # start collecting when Force constants in cartesian coordinates begin 
        if " Force constants in Cartesian " in line: 
            if capturing and current_section:
                sections.append((start_index,current_section))
            capturing = True
            start_index = i
            current_section = []
            indices.append(start_index)

            continue
        
        if " FormGI is" in line:
            capturing = False
            if current_section:
                sections.append((start_index,current_section))
            current_section = []


        if capturing:
           current_section.append(line)
    
    max_index = max(indices)
    f_matrix_raw = []
    for sec in sections:
        if sec[0] == max_index:
            for i, el in enumerate(sec):
                if i == 0:
                    continue # skip the first element
                f_matrix_raw.extend(el)
    
    rows = []
    current_row = []
    block_dict = {}
    block_index = 0
    for line in f_matrix_raw:
        # Find all numbers and convert to flow
        numbers = re.findall(r'[-+]?\d*\.\d+E[-+]?\d+', line.replace('D','E'))
        numbers = list(map(float,numbers)) # convs
        
        rows.append(numbers) 
        
    # in this rows for each [] we generate a new block
    for row in rows:
        if len(row) == 0:
            block_index += 1
            block_dict.setdefault(block_index)
            block_dict[block_index] = []
        else:
            if block_index in block_dict:
                block_dict[block_index].append(row)
    
    
    # now we generate our final hessian rows
    hessian_rows = []
    row_index = 1
    # random idx
    idx = 1
    for entries in block_dict[idx]:
        if len(entries) == row_index:
            hessian_rows.append(entries)
        if len(entries) < row_index:
                # now we need to add elements out of the second block
           idx_next_dict = idx + 1
           # we need for example extend the first element out of the 
           # the absolute value of the difference gives us the elements we need to add -1 because index starts = 0
           index_elements_difference = abs(row_index - len(entries)) - 1
           entries.extend(block_dict[idx_next_dict][index_elements_difference])
           if len(entries) < row_index:
               while  len(entries) < row_index:
                   # now we jump into the next block
                   idx_next_dict += 1
                   # for example jump from 2 to three and so on
                   index_elements_difference = abs(row_index - len(entries)) - 1 
                   entries.extend(block_dict[idx_next_dict][index_elements_difference])
               hessian_rows.append(entries)
           else: 
               hessian_rows.append(entries)
        row_index += 1
    # now initialize our hessian matrix
    n = n_atoms * 3
    f_matrix = np.zeros((n,n))
    
    for i in range(n):
        for j in range(i + 1):
            f_matrix[i,j] = hessian_rows[i][j]
    # the f matrix is now lower triangular
    symmetric_f_matrix = np.copy(f_matrix)
    for i in range(n):
        for j in range(i + 1):
            symmetric_f_matrix[j,i] = f_matrix[i,j]
    return symmetric_f_matrix
