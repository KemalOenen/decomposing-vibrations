import numpy as np
import pandas as pd
from nomodeco.libraries.nomodeco_classes import Molecule





def get_linear_bonds(linear_angles) -> list:
    """
    Uses the linear angles (i.e lists of the form ["H","O","H"]) to determine linear bonds and safe them in a list of tuples
    """
    linear_bonds = []
    for i in range(len(linear_angles)):
        # extract elements from current tuple
        first = linear_angles[i][0]
        middle = linear_angles[i][1]
        last = linear_angles[i][2]

        # create new tuples by including adjacent elements
        linear_bonds.append((middle, first))
        linear_bonds.append((middle,last))
    unique_linear_bonds = list(set(frozenset(bond) for bond in linear_bonds))
    linear_bonds = [tuple(unique_bond) for unique_bond in unique_linear_bonds]
    return linear_bonds

def is_string_in_tuples(string, list_of_tuples) -> bool:
    """
    checks if a given string is contained in a list of tuples
    """
    for tuple_ in list_of_tuples:
        if string in tuple_:
            return True
    return False 


def is_system_planar(coordinates, tolerance=1e-1) ->bool:
    """
    Calculates if a given system is planar, here three non linear atoms are selected, a normal vector perpenticular to the plane of the three atoms gets calculates.
    Then one calculates the dot product and compares it with a given tolerance

    Attributes:
        coordinates:
            a list of coordinates for each atom
        tolerance = 1e-1
            the given tolerance to determine planarity
    """
    # convert tuples first to arrays
    if len(coordinates) == 0:
        return True
    # Also if the length of the molecule is 2 this is always gtrue
    if len(coordinates) == 2:
       return True
    coordinates = [np.array(coord) for coord in coordinates]
    # select three non-linearly aligned atoms (already pre-filtered)
    atom1, atom2, atom3 = coordinates[:3]

    # calculate the vector perpendicular to the plane:
    vec1 = atom2 - atom1
    vec2 = atom3 - atom1
    normal_vector = np.cross(vec1, vec2)
    
    # Iterate through remaining atoms and calculate dot products
    for atom in coordinates[3:]:
        vec3 = atom - atom1
        dot_product = np.dot(normal_vector, vec3)
        
        #if connectivity_c > 1:
        #   tolerance = 1e-1
        # Check if dot product is close to zero within the tolerance
        if abs(dot_product) > tolerance:
            return False

    return True

def bound_to_atom(query_atom, bonds, central_atom) -> bool:
    """
    given a query_atom and a central atom and a bond checks if the query atom is bound to the central atom
    """
    for bond in bonds:
        if query_atom in bond and central_atom in bond:
            return True
    return False

def calculation_specification(specification, atoms, molecule_pg, bonds, angles, linear_angles) -> dict:
    """
    Calculates the specification dictionary, for example the numbers mu, beta, the planarity of a system, the multiplicity of the atoms, and a list of linear bond

    Attributes:
        atoms:
            a object of the Molecule class
        molecule_pg: 
            the point group of a given molecular structure
        bonds:
            the bonds given as a list of tuples
        angles:
            the angles given as a list of tuples
        linear_angles:
            the linear angles given as a list of tuples
    """
    # map every atom on their multiplicity 
    atoms_multiplicity = dict()
    atoms_multiplicity_list = []
    for bond_pair in bonds:
        for atom in bond_pair:
            if atom in atoms_multiplicity:
                atoms_multiplicity[atom] += 1
            else:
                atoms_multiplicity[atom] = 1
    for atom, multiplicity in atoms_multiplicity.items():
        atoms_multiplicity_list.append((atom, multiplicity))

    specification = {"multiplicity": atoms_multiplicity_list}
    
    # check if molecule is planar or general
    
    all_coordinates = []
    for atom in atoms:
        # atoms should not be colinear, otherwise we can not form a useful cross product
        if is_string_in_tuples(atom.symbol, angles) or not is_string_in_tuples(atom.symbol, linear_angles):
            all_coordinates.append(atom.coordinates)
    if is_system_planar(all_coordinates):
        specification.update({"planar": "yes"})
        specification.update({"planar submolecule(s)": [] })
    else:
        specification.update({"planar": "no"})
        # check if there are planar subunits, if so, check which atoms
        valid_atoms = []
        for tup in atoms_multiplicity_list:
            if tup[1] > 2:
                valid_atoms.append(tup)

        central_atom_planar_submolecule = []
        for valid_atom in valid_atoms:
            check_coordinates = []
            for atom in atoms:
                if atom.symbol in valid_atom[0] or bound_to_atom(atom.symbol, bonds, valid_atom[0]):
                    check_coordinates.append(atom.coordinates)
            if is_system_planar(check_coordinates):
                central_atom_planar_submolecule.append(valid_atom)

        specification.update({"planar submolecule(s)": central_atom_planar_submolecule})

    # check if molecule is linear or has a linear submolecule
    # if the molecule has a linear submolecule, then determine the number of linear bonds
    if (linear_angles and not angles) or (not linear_angles and not angles):
        specification.update({"linearity": "fully linear"})
    if (linear_angles and angles):
        specification.update({"linearity": "linear submolecules found"})
        linear_bonds = get_linear_bonds(linear_angles)
        specification.update({"length of linear submolecule(s) l": len(linear_bonds)})
    if (not linear_angles and angles):
        specification.update({"linearity": "not linear"})
    

    atoms = Molecule(atoms)
    connectivity_c = atoms.count_connected_components(atoms.graph_rep())
    mu = atoms.mu()
    beta = atoms.beta()
        
    if connectivity_c >= 2:
       specification.update({"intermolecular": "yes"})
       if mu > 0 and beta == 0:
          specification.update({"cyclic": "yes"})
          specification.update({"beta": beta})
          specification.update({"mu": mu})
          specification.update({"intermolecular ring": "yes"})
       elif mu > 0 and beta == mu:
          specification.update({"beta" : beta})
          specification.update({"mu" : mu})
          specification.update({"intermolecular ring": "no"})
       elif mu <= 0:
          specification.update({"cyclic": "no"})
          specification.update({"mu": mu})
    else:   
      # check if molecule is acyclic or cyclic
       if mu > 0:
        specification.update({"cyclic": "yes"})
        specification.update({"mu": mu})
        specification.update({"intermolecular": "no"})
       elif mu <= 0:
        specification.update({"cyclic": "no"})
        specification.update({"mu": 0})
        specification.update({"intermolecular": "no"})
       else:
        specification.update({"cyclic": "no"})
        specification.update({"mu": mu})
        specification.update({"intermolecular": "no"})

    # define which atoms are equal or not, based on group theory
    equivalent_atoms = molecule_pg.get_equivalent_atoms()
    atom_names = [atom.symbol for atom in atoms]
    equivalent_atoms_list = []
    
    for atom_number_set in equivalent_atoms["eq_sets"].values():
        equivalent_atoms_list.append(list(atom_number_set))


    for i, sublist in enumerate(equivalent_atoms_list):
        for j, element in enumerate(sublist):
            equivalent_atoms_list[i][j] = atom_names[element]
   
   
    specification.update({
        "equivalent_atoms": equivalent_atoms_list 
        })
    return specification
