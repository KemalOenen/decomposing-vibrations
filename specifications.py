import numpy as np
import pandas as pd

def get_linear_bonds(linear_angles):
    linear_bonds = []
    for i in range(len(linear_angles)):
        # extract elements from current tuple
        first = linear_angles[i][0]
        middle = linear_angles[i][1]
        last = linear_angles[i][2]

        # create new tuples by including adjacent elements
        linear_bonds.append((middle, first))
        linear_bonds.append((middle,last))
    linear_bonds = list(set(linear_bonds))
    return linear_bonds

def is_string_in_tuples(string, list_of_tuples):
    for tuple_ in list_of_tuples:
        if string in tuple_:
            return True
    return False 

def is_system_planar(coordinates, tolerance=1e-3):
    # convert tuples first to arrays
    if len(coordinates) == 0:
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
        # Check if dot product is close to zero within the tolerance
        if abs(dot_product) > tolerance:
            return False

    return True

def bound_to_atom(query_atom, bonds, central_atom):
    for bond in bonds:
        if query_atom in bond and central_atom in bond:
            return True
    return False

def calculation_specification(specification, atoms, molecule_pg, bonds, angles, linear_angles):
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

    # check if molecule is acyclic or cyclic 
    mu = len(bonds) - len(atoms) + 1
    if mu > 0:
        specification.update({"cyclic": "yes"})
        specification.update({"mu": mu})
    else:
        specification.update({"cyclic": "no"})
        specification.update({"mu": mu})


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
