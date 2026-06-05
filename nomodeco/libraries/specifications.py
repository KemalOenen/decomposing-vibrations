import numpy as np
from collections import Counter
import time
from nomodeco.libraries.ic_class import Molecule


def get_linear_bonds(linear_angles) -> list:
    bonds = set()
    for a1, mid, a3 in linear_angles:
        bonds.add(frozenset((mid, a1)))
        bonds.add(frozenset((mid, a3)))
    return [tuple(b) for b in bonds]


def is_string_in_tuples(string, list_of_tuples) -> bool:
    return any(string in t for t in list_of_tuples)


def is_system_planar(coordinates, tolerance=1e-1) -> bool:
    if len(coordinates) < 3:
        return True
    coords = np.array(coordinates)
    normal = np.cross(coords[1] - coords[0], coords[2] - coords[0])
    if np.linalg.norm(normal) == 0:
        return True
    return bool(np.all(np.abs((coords[3:] - coords[0]) @ normal) <= tolerance))


def bound_to_atom(query_atom, bonds, central_atom) -> bool:
    return any(query_atom in bond and central_atom in bond for bond in bonds)


def calculation_specification(specification, atoms, molecule_pg, bonds, angles, linear_angles) -> dict:
    multiplicity = Counter(atom for bond in bonds for atom in bond)
    atoms_multiplicity_list = list(multiplicity.items())
    specification = {"multiplicity": atoms_multiplicity_list}

    # Planarity
    all_coordinates = [
        atom.coordinates for atom in atoms
        if is_string_in_tuples(atom.symbol, angles) or not is_string_in_tuples(atom.symbol, linear_angles)
    ]
    if is_system_planar(all_coordinates):
        specification.update({"planar": "yes", "planar submolecule(s)": []})
    else:
        specification["planar"] = "no"
        specification["planar submolecule(s)"] = [
            (sym, mult) for sym, mult in atoms_multiplicity_list
            if mult > 2 and is_system_planar([
                a.coordinates for a in atoms
                if a.symbol == sym or bound_to_atom(a.symbol, bonds, sym)
            ])
        ]

    # Linearity
    if not angles:
        specification["linearity"] = "fully linear"
    elif linear_angles:
        linear_bonds = get_linear_bonds(linear_angles)
        specification.update({
            "linearity": "linear submolecules found",
            "length of linear submolecule(s) l": len(linear_bonds),
        })
    else:
        specification["linearity"] = "not linear"

    # Topology
    atoms_mol = Molecule(atoms)
    connectivity_c = atoms_mol.count_connected_components(atoms_mol.graph_rep())
    mu = atoms_mol.mu()
    beta = atoms_mol.beta()

    if connectivity_c >= 2:
        specification["intermolecular"] = "yes"
        if mu > 0 and beta == 0:
            specification.update({"cyclic": "yes", "beta": beta, "mu": mu, "intermolecular ring": "yes"})
        elif mu > 0 and beta == mu:
            specification.update({"beta": beta, "mu": mu, "intermolecular ring": "no"})
        else:
            specification.update({"cyclic": "no", "mu": mu})
    else:
        specification["intermolecular"] = "no"
        specification.update({"cyclic": "yes" if mu > 0 else "no", "mu": mu if mu > 0 else 0})

    # Equivalent atoms from point group
    atom_names = [atom.symbol for atom in atoms_mol]
    specification["equivalent_atoms"] = [
        [atom_names[idx] for idx in group]
        for group in molecule_pg.get_equivalent_atoms()["eq_sets"].values()
    ]
    return specification


if __name__ == "__main__":
    # Check runtime of planarity and linearity checks
    N = [10, 20, 30, 40, 50, 100, 200, 300, 400, 500]
    # generate random coordinates and check planarity
    for n in N:
        coords = np.random.rand(n, 3)
        start_time = time.time()
        planar = is_system_planar(coords)
        print(f"Planarity check for {n} atoms: {planar}, runtime: {time.time() - start_time:.4f} seconds")