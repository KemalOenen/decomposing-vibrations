#general packages
from __future__ import annotations

from typing import NamedTuple
import string
import os
import numpy as np
import pandas as pd
import seaborn as sns
import pymatgen.core as mg
import itertools
import networkx as nx


from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from mendeleev.fetch import fetch_table

# TODO: isotopes for all elements with command line input?
def get_mass_information():
    df = fetch_table('elements')
    mass_info = df.loc[:, ['symbol', 'atomic_weight']]
    deuterium_info = pd.DataFrame({'symbol': ['D'], 'atomic_weight': [2.014102]})
    mass_info = pd.concat([mass_info, deuterium_info])
    mass_info.set_index("symbol", inplace=True)
    return mass_info

def get_bond_information():
    df = fetch_table('elements')
    bond_info = df.loc[:, ['symbol', 'covalent_radius_pyykko', 'vdw_radius']]
    
    bond_info.set_index("symbol", inplace=True)
    bond_info /= 100
    # manually add deuterium
    bond_info.loc['D'] = {'covalent_radius_pyykko': 0.25, 'vdw_radius': 120}
    return bond_info

BOND_INFO = get_bond_information()



'''
Molecule Class in Nomodeco
'''

class Molecule(list):
    """
    Nomodeco molecule class

    Attributes:

        a list of of atoms
    """ 
    
    # Sublass Atoms Definition
    class Atom:
        """
        The atom class in nomodeco
        
        Attributes:
        
            symbol:
                Enumerated atomic symbol of the molecule
            coordinates:
                3D coordinates of a atom in a tuple
        """
        
        def __init__(self,symbol, coordinates):
            """
            Constructor of the Atom class

            Attributes:
                symbol:
                    Enumerated atomic symbol of the molecule
                coordinates:
                    3D coordinates of atom in a tuple
            """
            self.symbol = symbol
            self.coordinates = coordinates
        #symbol: str
        #coordinates: tuple
        
        def __repr__(self): 
            """
            Printable representation of the Atom class

            Returns:
                str 
            """
            # Returns a string of the coordinates and the symbol
            return f"{self.symbol} at {self.coordinates}"

        # Subclass Variables:
        donor_atoms = ['N','O','F','Cl'] 

        # Useful Methods for Hydrogen Bonding
        def is_donor(self) -> bool:
            """
            Checks if a given atom is a donor molecule

            Attributes:
                self: a object of the Atom class
            
            Returns:
                True if a atom is a donor atom
            """
            if self.symbol.strip(string.digits) in self.donor_atoms:
                return True
        def is_hydrogen(self) -> bool:
            """
            Checks if a given atom is a hydrogen atom

            Attributes:
                self: a object of the Atom class
            
            Returns:
                True if the atom is hydrogen
            """
            if self.symbol.strip(string.digits) == "H":
                return True

        # Swaping the symbol in a atom    
        def swap_deuterium(self) -> self:
            """
            Swaps a given hydrogen atom with a deuterium atom

            Attributes:
                self: a object of the Atom class

            Returns:
                a object of the atom class with H swapped to D
            """
            self.symbol = self.symbol.replace('H','D',1)
    
    # retrieve class variables
    def get_donor_atoms(self) -> list:
        """
        Returns the class variable donor_atom of the Atom class
        """
        return self.Atom.donor_atoms



    # Dunder Methods for Debugging
    def __repr__(self):
        """
        Printable representation of the Molecule class

        Returns a string with string of the type

            Molecule[(atom1 at coordinates), ...]
        """
        atom_reprs = ', '.join([repr(atom) for atom in self])
        # here with the Molecule Class Keyword
        return f"Molecule([{atom_reprs}])"

    # Printing of the Molecule Type can be adjusted here
    def __str__(self):
        """
        Representation of the molecule when print() in python gets used

        Returns a string of the type Molecule with atoms: atoms at coordinates
        """
        atom_strs = ", ".join([str(atom) for atom in self])
        return f"Molecule with atoms: {atom_strs}"

    def num_atoms(self) -> int:
        """
        Returns the number of atoms given of the given Molecule
        """
        num_of_atoms = 0
        for atom in self:
            num_of_atoms +=1
        return num_of_atoms
    
    # List of symbols
    def list_of_atom_symbols(self) -> list:
        """
        Returns a list of just the symbols in the given Molecule

            e.q ['H1','O2','H2']
        """
        symbol_list = []
        for atom in self:
            symbol_list.append(atom.symbol)
        return symbol_list 

    # IDOF Functions

    def idof_linear(self) -> int:
        """
        Returns the degrees of freedom in the linear case 3N-5
        """
        num_of_atoms = self.num_atoms()
        return 3*num_of_atoms - 5

    def idof_general(self) -> int:
        """
        Returns the degrees of freedom in the general case 3N-6
        """
        num_of_atoms = self.num_atoms()
        return 3*num_of_atoms -6

    def append(self, value):
        """
        Appends a given atom to the Molecule Class
        
        Attributes:
            self: a object of the Molecule Class
            value: a object of the Atom Class

        """
        if not isinstance(value, self.Atom):
            raise TypeError(f"only instances of Atoms can be added not {type(value).__name__}")
        super().append(value) # super imports the method of the list type


    # covalent_length
    def theoretical_length(self,symbol1,symbol2) -> float:
        """
        Returns the theoretical covalent bond length between two atoms

        Attributes:
            self: 
                a object of the Molecule Class
            symbol1: str
                a string of the first atom in question
            symbol2: str
                a string of the second atom in question
        Returns:
            a floating point number with the covalent bond length
        """
        element1 = next((el for el in self if el.symbol == symbol1),None)
        element2 = next((el for el in self if el.symbol == symbol2), None)

        # case symbols not found
        if element1 is None or element2 is None:
            raise ValueError("One of the Elements not found in the Molecule")
        
        rad_a = BOND_INFO.loc[element1.symbol.strip(string.digits)]
        rad_b = BOND_INFO.loc[element2.symbol.strip(string.digits)]
        return rad_a[0] + rad_b[0]
    
    def theoretical_length_vdw(self,symbol1,symbol2) -> float:
        """
        Returns the theoretical vdw length between two atoms

        Attributes:
            self: 
                a object of the Molecule Class
            symbol1: str
                a string of the first atom in question
            symbol2: str
                a string of the second atom in question
        
        Returns:
            a floating point number with the vdw bond length
        """
        element1 = next((el for el in self if el.symbol == symbol1),None)
        element2 = next((el for el in self if el.symbol == symbol2), None)

        # case symbols not found
        if element1 is None or element2 is None:
            raise ValueError("One of the Elements not found in the Molecule")
        
        rad_a = BOND_INFO.loc[element1.symbol.strip(string.digits)]
        rad_b = BOND_INFO.loc[element2.symbol.strip(string.digits)]
        return rad_a[1] + rad_b[1]
    

    def actual_length(self,symbol1,symbol2) -> float:
        """
        Calculates the actual bond length between two atoms

        Attributes:
            self: 
                a object of the Molecule Class
            symbol1: str
                a string of the first atom in question
            symbol2: str
                a string of the second atom in question
        
        Returns:
            a floating point number with the actual bond length
        """
        element1 = next((el for el in self if el.symbol == symbol1),None)
        element2 = next((el for el in self if el.symbol == symbol2), None)

        # case symbols not found
        if element1 is None or element2 is None:
            raise ValueError("One of the Elements not found in the Molecule")
        return abs(np.linalg.norm(np.array(element1.coordinates) - np.array(element2.coordinates)))
    
    def bond_angle(self, symbol1, symbol2, symbol3) -> float:
        """
        Calculates the bond angle between three atoms in the Molecule object

        Attributes:
            self:
                a object of the Molecule class
            symbol1: str
                a string of the atom
            symbol2: str
                a string of the atom
            symbol3: str
                a string of the atom
        Returns:
            The angle between the three atoms (ABC)
        """
        element1 = next((el for el in self if el.symbol == symbol1),None)
        element2 = next((el for el in self if el.symbol == symbol2), None)
        element3 = next((el for el in self if el.symbol == symbol3),None)

        coord_a = np.array(element1.coordinates)
        coord_b = np.array(element2.coordinates)
        coord_c = np.array(element3.coordinates)
        cosine_angle = np.clip((np.inner((coord_a - coord_b), (coord_c - coord_b))) / (
            self.actual_length(element1.symbol, element2.symbol) * (self.actual_length(element3.symbol, element2.symbol))), -1.0, 1.0)
        return np.arccos(cosine_angle) 
    
    '''
    Degree of Covalance
    '''


    def degree_of_covalance(self) -> dict:
        """
        Calculates the degree of covalance between all the combinations of two atoms
        
        Reference  https://doi.org/10.1002/qua.21049

        Returns:
            a python dictionary with the two atoms as a tuple and the value of the degree of covalance
        """
        degofc = {} # initialize array
        for atom_a,atom_b in itertools.combinations(self,2):
            value = np.exp(-(abs(self.actual_length(atom_a.symbol,atom_b.symbol))/self.theoretical_length(atom_a.symbol,atom_b.symbol) - 1))
            degofc.update({(atom_a.symbol,atom_b.symbol) : value})
        return degofc


    '''
    Covalent_bond_detection
    '''
    # Gives a list of all covalent bonds using the degofc_table

    def covalent_bonds(self,degofc_table) -> list:
        """
        Searches for covalent bonds between two atoms. A covalent bond is found when the degree of covalence is greater than 0.75

        Attributes:
            self:
                a object of the Molecule class
            degofc_table: dict
                a dictionary containing the degrees of covalence
            
        Returns:
            a list of possible bonds e.q [(O1,H2), ...]
        """
        combinations = [(atom_a.symbol,atom_b.symbol) for atom_a,atom_b in itertools.combinations(self,2)]
        hits = []
        for key,value in degofc_table.items():
            if key in combinations and value > 0.75:
                hits.append(key)
        return hits
    
        
    '''
    Graph Theory Functions
    '''
    
    def detect_submolecules(self):
        """
        A graph theory oriented function to detect the covalent submolecules

        Attributes:
            self:
                a object of the molecule class
        Returns:
            connected_components: more information see networkx.connected_components
            submolecules: a dictionary containing the subgraphs of the molecule
            submolecule_symbols: a list of the atomic symbols of the submolecules

        """
        bonds = self.covalent_bonds(self.degree_of_covalance())
        molecular_graph = nx.Graph()
        molecular_graph.add_edges_from(bonds)
        connected_components = list(nx.connected_components(molecular_graph))
        # initialize the submolecule list
        submolecules = []
        for component in connected_components:
             subgraph = molecular_graph.subgraph(component)
             submolecule_bonds = list(subgraph.edges)
             submolecules.append(submolecule_bonds)
        # also generate just the symbols because its convenient
        submolecule_symbols = {}
        i = 0
        symbols = set()
        for submolecule in submolecules:
            for bond in submolecule: # Use the individual bonds
                symbols.update(bond)
            submolecule_symbols.update({i: symbols})
            symbols = set() # clear the symbols every loop
            i = i +1
        return connected_components, submolecules, submolecule_symbols
    
    def graph_rep(self):
        """
        Generates the covalent bonds and using the covalent bonds a graph representation of a molecule

        Attributes:
            self:
                a object of the Molecule class
        
        Returns:
            a dictionary where the graph representation of the molecule is stored
        """
        graph = {}
        bonds = self.covalent_bonds(self.degree_of_covalance())
        for bond in bonds:
            a,b = bond
            graph.setdefault(a,[]).append(b)
            graph.setdefault(b,[]).append(a)
        return graph

    
    @staticmethod
    def dfs(graph,start,visited):
        """
        A static method that applies dfs algorithm on a given graph

        Attributes:
            graph:
                a molecular graph created by Molecule.graph_rep()
            start:
                start node for the dfs algorithm
            visitied:
                visited list for the dfs algorithm
        """
        visited.add(start)
        for neighbor in graph.get(start,[]):
            if neighbor not in visited:
                Molecule.dfs(graph,neighbor,visited) # Recursion

    @staticmethod
    def is_connected(graph):
        """
        Uses the dfs method and a molecular graph to check if a graph is connected

        Attributes:
            graph:
                a molecular graph created by Molecule.graph_rep()
        
        Returns:
            True if the length of the visited list is equal in length of the total graph
        """
        if not graph:
            return True
        visited = set()
        start_node = next(iter(graph))
        Molecule.dfs(graph,start_node,visited)
        return len(visited) == len(graph)
    
    @staticmethod
    def count_connected_components(graph) -> int:
        """
        Using the dfs algorithm counts the fully covalently connected components inside a molecule

        Attributes:
            graph:
                a molecular graph created by Molecule.graph_rep()
        
        Returns:
            a integer which represents the total count of connected components ranges from 1..N
        """
        if not graph:
            return 0
        visited = set()
        count = 0
        for node in graph:
            if node not in visited:
                count += 1
                Molecule.dfs(graph,node,visited)
        return count
   
    # Is used in Topology.py
    def bond_dict(self,bonds) -> dict:
        """
        Created a dictionary where for each atom the corresponding "bonded" atoms are stored

        Attributes:
            self:
                a object of the Molecule class
            bonds:
                a list of tuples containing the bonds which are used for the dict creation
        
        Returns:
            A dictionary containing the bonded atom for each atom e.q {H1: [O1,O2]}
        """
        bond_dict = {}
        for atom_1,atom_2 in bonds:
            if atom_1 not in bond_dict:
                bond_dict[atom_1] = []
            if atom_2 not in  bond_dict:
                bond_dict[atom_2] = []
            bond_dict[atom_1].append(atom_2)
            bond_dict[atom_2].append(atom_1)
        return bond_dict


    '''
    Retrieving Coordinates
    '''

    def get_atom_coords_by_symbol(self,symbol):
        """
        Retries a atom by specifying a atomic symbol

        Attributes:
            self:
                a object of the Molecule class
            symbol:
                a atomic symbol as a string
        """
        for atom in self:
            if atom.symbol == symbol:
                return atom.coordinates
        raise ValueError(f"Atom with symbol {symbol} not found")
    def retrieve_index_in_list(self,symbol):
        for i, atom in enumerate(self):
            if atom.symbol == symbol:
                return i

    '''
    Hydrogen_bond_dection
    '''
    def intermolecular_h_bond(self,degofc_table,submolecule_symbols):
        combinations = [(atom_a.symbol,atom_b.symbol) for atom_a,atom_b in itertools.combinations(self,2)]
        possible_h_bonds = []
        donor_atoms = self.get_donor_atoms()
        for key,value in degofc_table.items():
            if key in combinations and value > 0.16 and value < 0.7:
                possible_h_bonds.append(key)
        h_bonds = []
        index = range(len(submolecule_symbols))
        for index_a, index_b in itertools.combinations(index,2):
            for bond in possible_h_bonds:
                # Long statment because we need symmetric thing also
                if (bond[0] in submolecule_symbols[index_a] and bond[1] in submolecule_symbols[index_b]) or (bond[1] in submolecule_symbols[index_a] and bond[0] in submolecule_symbols[index_b]):
                    if ((bond[0].strip(string.digits) == "H" or bond[1].strip(string.digits) == "H") and
                    (bond[0].strip(string.digits) in donor_atoms or bond[1].strip(string.digits) in donor_atoms)):
                        h_bonds.append(bond)
        return list([tuple(sorted(i)) for i in h_bonds])

    def intermolecular_acceptor_donor(self, degofc_table,submolecule_symbols):
        combinations = [(atom_a.symbol,atom_b.symbol) for atom_a,atom_b in itertools.combinations(self,2)]
        possible_acc_don_bond = []
        for key,value in degofc_table.items():
            if key in combinations and value > 0.25 and value < 0.7:
                possible_acc_don_bond.append(key)
        acc_don_bonds = []
        index = range(len(submolecule_symbols))
        for index_a, index_b in itertools.combinations(index, 2):
            for bond in possible_acc_don_bond:
                if bond[0] in submolecule_symbols[index_a] and bond[1] in submolecule_symbols[index_b] and bond[0].strip(string.digits) in self.Atom.donor_atoms and bond[1].strip(string.digits) in self.Atom.donor_atoms:
                    acc_don_bonds.append(bond)
        return acc_don_bonds
   
    '''
    Adjacency Matrix
    '''


    # Generate adjacency matrix of a molecule
    def covalent_adjacency_matrix(self):
        atom_symbols = self.list_of_atom_symbols()
        bonds = self.covalent_bonds(self.degree_of_covalance())
        molecular_graph = nx.Graph()
        molecular_graph.add_nodes_from(atom_symbols)
        molecular_graph.add_edges_from(bonds)
        return nx.to_numpy_array(molecular_graph)

    # Generate hydrogen bond adjacency matrix
    def hydrogen_adjacency_matrix(self):
        _, _, submolecule_symbols = self.detect_submolecules() 
        degofc_table = self.degree_of_covalance()
        atom_symbols = self.list_of_atom_symbols()
        h_bonds = self.intermolecular_h_bond(degofc_table,submolecule_symbols)
        molecular_graph = nx.Graph()
        molecular_graph.add_nodes_from(atom_symbols)
        molecular_graph.add_edges_from(h_bonds)
        return nx.to_numpy_array(molecular_graph)

    '''
    Mu and Beta Calculation
    '''
    def mu(self):
        # evaluate degofc
        # mu = total_bonds + atoms - 1
        degofc = self.degree_of_covalance()
        num_of_atoms = self.num_atoms()
        cov_bonds = self.covalent_bonds(degofc)
        connectivity_c = self.count_connected_components(self.graph_rep())
        _,_,submolecule_symbols = self.detect_submolecules()
        len_bonds = len(cov_bonds) + len(self.intermolecular_h_bond(degofc,submolecule_symbols))
        return len_bonds - num_of_atoms + 1
    def beta(self):
        # evaluate doegofc
        degofc = self.degree_of_covalance()
        num_of_atoms = self.num_atoms()
        cov_bonds = self.covalent_bonds(degofc)
        connectivity_c = self.count_connected_components(self.graph_rep())
        return len(cov_bonds) - num_of_atoms + connectivity_c

    '''
    Generation of primitive ICs
    '''
    #TODO  Similar function defined above maybe delete sometime :)
    @staticmethod
    def generate_connectivity(bonds):
        connectivity_dict = {}
        for atom1,atom2 in bonds:
            if atom1 not in connectivity_dict:
               connectivity_dict[atom1] = [] 
            if atom2 not in connectivity_dict:
               connectivity_dict[atom2] = []
            connectivity_dict[atom1].append(atom2)
            connectivity_dict[atom2].append(atom1)
        return connectivity_dict

    
    def generate_angles(self,bonds):
        connectivity_dict = self.generate_connectivity(bonds)
        possible_angles = []
        for atom,bonded_atoms in connectivity_dict.items():
            bonded_list = list(bonded_atoms)
            for i in range(len(bonded_list)):
                for j in range(i+1,len(bonded_list)):
                    possible_angles.append((bonded_list[i],atom, bonded_list[j]))
        angles = []
        linear_angles = []
        for element in possible_angles:
            if self.bond_angle(element[0],element[1],element[2])*(180/np.pi) > 10 and self.bond_angle(element[0],element[1],element[2])*(180/np.pi) < 169:
                angles.append(element)
            elif self.bond_angle(element[0],element[1],element[2])*(180/np.pi) >= 169:
                linear_angles.append(element) # 2x degenerate
                linear_angles.append(element)
        return angles,linear_angles
    
    def generate_dihedrals(self,bonds):
        connectivity_dict = self.generate_connectivity(bonds)
        possible_dihedrals = []
        for atom_b,bonded_atoms_b in connectivity_dict.items():
            for atom_a in bonded_atoms_b:
                for atom_c in bonded_atoms_b:
                    if atom_a == atom_c:
                        continue # use continue to end the current iteration
                    bonded_atoms_c = connectivity_dict[atom_c]
                    for atom_d in bonded_atoms_c:
                        if atom_d != atom_b and atom_d != atom_a:
                            # bugfix for symmetric dupletts
                           if (atom_a,atom_b,atom_c,atom_d) not in possible_dihedrals and (atom_d,atom_c,atom_b,atom_a) not in possible_dihedrals:
                            possible_dihedrals.append((atom_a,atom_b,atom_c,atom_d))
        # bugfix still symmetric dupletts ?!
        unique_dihedrals = set()
        for dihedral in possible_dihedrals:

            reversed_dihedral = tuple(reversed(dihedral))
            unique_dihedrals.add(min(dihedral,reversed_dihedral))

        return list(unique_dihedrals)

    def generate_out_of_plane(self,bonds):
        connectivity_dict = self.generate_connectivity(bonds)
        out_of_plane = []
        for atom, bonded_atoms in connectivity_dict.items():
            bonded_list = list(bonded_atoms)
            if len(bonded_list) >=3: # the only possible way we can generate out of plane angles
                for i in range(len(bonded_list)):
                    for j in range(i+1, len(bonded_list)):
                        for k in range(j + 1, len(bonded_list)):
                            # Append all three combinations
                            out_of_plane.append((atom,bonded_list[i],bonded_list[j],bonded_list[k]))
                            out_of_plane.append((atom,bonded_list[j],bonded_list[i],bonded_list[k]))
                            out_of_plane.append((atom,bonded_list[k],bonded_list[i],bonded_list[j]))
        return out_of_plane
    def generate_oop_planar_subunits(self,bonds, central_atoms_list):
        hits = []
        connectivity_dict = self.generate_connectivity(bonds)
        for central_atom in central_atoms_list:
            for atom, bonded_atoms in connectivity_dict.items():
                if atom == central_atom[0]:
                    bonded_list = list(bonded_atoms)
                    if len(bonded_list) >=3: # the only possible way we can generate out of plane angles
                       for i in range(len(bonded_list)):
                           for j in range(i+1, len(bonded_list)):
                               for k in range(j + 1, len(bonded_list)):
                               # Append all three combinations
                                   hits.append((atom,bonded_list[i],bonded_list[j],bonded_list[k]))
                                   hits.append((atom,bonded_list[j],bonded_list[i],bonded_list[k]))
                                   hits.append((atom,bonded_list[k],bonded_list[i],bonded_list[j]))
        return hits


'''
Internal Coordinate Class
'''
# This class provides a way to safe store and manipulate the primitive ICs
# primitive ICs are safed in a Dictionary as a list of tuples with the coordinate type as key

class InternalCoordinates:
    def __init__(self):
        # Initialize an empty dictionary to store internal coordinates
        self.coordinates = {}

    def add_coordinate(self,key,coordinate_list):
        
        self.coordinates[key] = coordinate_list
    
    def get_coordinate(self,key):
        return self.coordinates.get(key,None)

    def common_coordinate(self,key,coordinate_list):
        result = [i for i in coordinate_list if i in self[key]]
        return result

    def add_coord_diff(self,key,ic_list1,ic_list2):
        # here we eliminate the symmetry problem:
        # if the length of a tuple is ==3 we also check the symmetric options
        if len(ic_list1) != 0:
           if len(ic_list1[0]) == 3:
               sym_ic_list_2 = [(c,b,a) for a,b,c in ic_list2]
               diff_ic_list = list(set(ic_list1).difference(set(ic_list2)))
               diff_ic_list = list(set(diff_ic_list).difference(set(sym_ic_list_2)))
               self.coordinates[key] = diff_ic_list
           else:
               diff_ic_list = list(set(ic_list1).difference(set(ic_list2)))
               self.coordinates[key] = diff_ic_list
        else:
            # If no coordinates to generate just return a empty list
            self.coordinates[key] = []

    def add_coord_diff_linear(self,key,ic_list1,ic_list2):
        sym_ic_list_2 = [(c,b,a) for a,b,c in ic_list2]
        diff_ic_list = list(set(ic_list1).difference(set(ic_list2)))
        diff_ic_list = list(set(diff_ic_list).difference(set(sym_ic_list_2)))
        dup_diff_ic_list = diff_ic_list * 2
        self.coordinates[key] = dup_diff_ic_list


    # add dunder for string rep
    def __str__(self):

        coords_str = ", ".join(f"{name}: {value}" for name,value in self.coordinates.items())
        return f"InternalCoordinates({coords_str})"
    def __getitem__(self,index):
        return self.coordinates[index]


