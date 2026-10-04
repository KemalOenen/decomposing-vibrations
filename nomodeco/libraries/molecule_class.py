from __future__ import annotations

import string
import itertools

import numpy as np
import networkx as nx
from numba import njit, prange

from mendeleev.fetch import fetch_table


def get_bond_information():
    df = fetch_table("elements")
    bond_info = df.loc[:, ["symbol", "covalent_radius_pyykko", "vdw_radius"]]
    bond_info.set_index("symbol", inplace=True)
    bond_info /= 100
    bond_info.loc["D"] = {"covalent_radius_pyykko": 0.25, "vdw_radius": 120}
    return bond_info





BOND_INFO = get_bond_information()

# Bond angles at or above this are treated as linear (two linear-bend coordinates); a dihedral
# through such an angle is undefined, so generate_dihedrals uses the same cutoff
LINEAR_ANGLE_DEG = 169.0


class Molecule(list):
    """
    Nomodeco molecule class

    Attributes:

        a list of of atoms
    """

    class Atom:
        """
        The atom class in nomodeco

        Attributes:

            symbol:
                Enumerated atomic symbol of the molecule
            coordinates:
                3D coordinates of a atom in a tuple
        """

        def __init__(self, symbol, coordinates):
            self.symbol = symbol
            self.coordinates = coordinates

        def __repr__(self):
            return f"{self.symbol} at {self.coordinates}"

        donor_atoms = ["N", "O", "F", "Cl"]

        def is_donor(self) -> bool:
            return self.symbol.strip(string.digits) in self.donor_atoms

        def is_hydrogen(self) -> bool:
            return self.symbol.strip(string.digits) == "H"

        def swap_deuterium(self) -> None:
            self.symbol = self.symbol.replace("H", "D", 1)

    def get_donor_atoms(self) -> list:
        return self.Atom.donor_atoms

    def __repr__(self):
        atom_reprs = ", ".join([repr(atom) for atom in self])
        return f"Molecule([{atom_reprs}])"

    def __str__(self):
        atom_strs = ", ".join([str(atom) for atom in self])
        return f"Molecule with atoms: {atom_strs}"

    def num_atoms(self) -> int:
        return len(self)

    def list_of_atom_symbols(self) -> list:
        return [atom.symbol for atom in self]

    def list_of_atom_symbols_xyz(self) -> list:
        result = []
        for atom in self:
            result.extend([f"{atom.symbol}x", f"{atom.symbol}y", f"{atom.symbol}z"])
        return result

    def idof_linear(self) -> int:
        return 3 * self.num_atoms() - 5

    def idof_general(self) -> int:
        return 3 * self.num_atoms() - 6

    @staticmethod
    def from_xyz_file(xyz_file) -> Molecule:
        atoms = []
        with open(xyz_file, "r") as f:
            lines = f.readlines()
            for line in lines[2:]:  # skip first two lines
                parts = line.split()
                if len(parts) >= 4:
                    symbol = parts[0]
                    coords = tuple(map(float, parts[1:4]))
                    atoms.append(Molecule.Atom(symbol, coords))
        return Molecule(atoms)

    def append(self, value):
        if not isinstance(value, self.Atom):
            raise TypeError(
                f"only instances of Atoms can be added not {type(value).__name__}"
            )
        super().append(value)

    def theoretical_length(self, symbol1, symbol2) -> float:
        return float(
            BOND_INFO.loc[symbol1.strip(string.digits)].iloc[0]
            + BOND_INFO.loc[symbol2.strip(string.digits)].iloc[0]
        )

    def theoretical_length_vdw(self, symbol1, symbol2) -> float:
        return float(
            BOND_INFO.loc[symbol1.strip(string.digits)].iloc[1]
            + BOND_INFO.loc[symbol2.strip(string.digits)].iloc[1]
        )

    def actual_length(self, symbol1, symbol2) -> float:
        sym_coords = {atom.symbol: atom.coordinates for atom in self}
        if symbol1 not in sym_coords or symbol2 not in sym_coords:
            raise ValueError("One of the Elements not found in the Molecule")
        return float(np.linalg.norm(
            np.array(sym_coords[symbol1]) - np.array(sym_coords[symbol2])
        ))

    def bond_angle(self, symbol1, symbol2, symbol3) -> float:
        sym_coords = {atom.symbol: np.array(atom.coordinates) for atom in self}
        ba = sym_coords[symbol1] - sym_coords[symbol2]
        bc = sym_coords[symbol3] - sym_coords[symbol2]
        cosine_angle = np.clip(
            np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc)),
            -1.0, 1.0,
        )
        return np.arccos(cosine_angle)

    def degree_of_covalance(self) -> dict:
        """
        Reference https://doi.org/10.1002/qua.21049
        """
        atoms = list(self)
        n = len(atoms)
        coords = np.array([a.coordinates for a in atoms])  # (n, 3)

        rad_cache = {}
        for a in atoms:
            s = a.symbol.strip(string.digits)
            if s not in rad_cache:
                rad_cache[s] = float(BOND_INFO.loc[s].iloc[0])
        radii = np.array([rad_cache[a.symbol.strip(string.digits)] for a in atoms])

        diff = coords[:, np.newaxis] - coords[np.newaxis]       # (n, n, 3)
        dist = np.sqrt((diff ** 2).sum(axis=2))                 # (n, n)
        theo = radii[:, np.newaxis] + radii[np.newaxis]         # (n, n)
        doc  = np.exp(-(dist / theo - 1))                       # (n, n)

        symbols = [a.symbol for a in atoms]
        return {(symbols[i], symbols[j]): doc[i, j]
                for i in range(n) for j in range(i + 1, n)}

    def covalent_bonds(self, degofc_table) -> list:
        return [key for key, value in degofc_table.items() if value > 0.75]

    def detect_submolecules(self, degofc_table=None):
        if degofc_table is None:
            degofc_table = self.degree_of_covalance()
        bonds = self.covalent_bonds(degofc_table)
        molecular_graph = nx.Graph()
        molecular_graph.add_edges_from(bonds)
        connected_components = list(nx.connected_components(molecular_graph))
        submolecules = []
        for component in connected_components:
            subgraph = molecular_graph.subgraph(component)
            submolecules.append(list(subgraph.edges))
        submolecule_symbols = {}
        for i, submolecule in enumerate(submolecules):
            symbols = set()
            for bond in submolecule:
                symbols.update(bond)
            submolecule_symbols[i] = symbols
        return connected_components, submolecules, submolecule_symbols

    def graph_rep(self, bonds=None):
        if bonds is None:
            bonds = self.covalent_bonds(self.degree_of_covalance())
        graph = {}
        for a, b in bonds:
            graph.setdefault(a, []).append(b)
            graph.setdefault(b, []).append(a)
        return graph

    @staticmethod
    def dfs(graph, start, visited):
        stack = [start]
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                stack.extend(n for n in graph.get(node, []) if n not in visited)

    @staticmethod
    def is_connected(graph):
        if not graph:
            return True
        visited = set()
        Molecule.dfs(graph, next(iter(graph)), visited)
        return len(visited) == len(graph)

    @staticmethod
    def count_connected_components(graph) -> int:
        if not graph:
            return 0
        visited = set()
        count = 0
        for node in graph:
            if node not in visited:
                count += 1
                Molecule.dfs(graph, node, visited)
        return count
    
    

    def bond_dict(self, bonds) -> dict:
        return self.generate_connectivity(bonds)

    def get_atom_coords_by_symbol(self, symbol):
        for atom in self:
            if atom.symbol == symbol:
                return atom.coordinates
        raise ValueError(f"Atom with symbol {symbol} not found")

    def retrieve_index_in_list(self, symbol):
        for i, atom in enumerate(self):
            if atom.symbol == symbol:
                return i

    def intermolecular_h_bond(self, degofc_table, submolecule_symbols):
        possible_h_bonds = [
            key for key, value in degofc_table.items()
            if 0.27 < value < 0.7
        ]
        donor_atoms = self.get_donor_atoms()
        h_bonds = []
        index = range(len(submolecule_symbols))
        for index_a, index_b in itertools.combinations(index, 2):
            for bond in possible_h_bonds:
                if (
                    bond[0] in submolecule_symbols[index_a]
                    and bond[1] in submolecule_symbols[index_b]
                ) or (
                    bond[1] in submolecule_symbols[index_a]
                    and bond[0] in submolecule_symbols[index_b]
                ):
                    if (
                        bond[0].strip(string.digits) in ("H", "D")
                        or bond[1].strip(string.digits) in ("H", "D")
                    ) and (
                        bond[0].strip(string.digits) in donor_atoms
                        or bond[1].strip(string.digits) in donor_atoms
                    ):
                        h_bonds.append(bond)
        return list([tuple(sorted(i)) for i in h_bonds])

    def intermolecular_acceptor_donor(self, degofc_table, submolecule_symbols):
        possible_acc_don_bond = [
            key for key, value in degofc_table.items()
            if 0.25 < value < 0.7
        ]
        acc_don_bonds = []
        index = range(len(submolecule_symbols))
        for index_a, index_b in itertools.combinations(index, 2):
            for bond in possible_acc_don_bond:
                if (
                    bond[0] in submolecule_symbols[index_a]
                    and bond[1] in submolecule_symbols[index_b]
                    and bond[0].strip(string.digits) in self.Atom.donor_atoms
                    and bond[1].strip(string.digits) in self.Atom.donor_atoms
                ):
                    acc_don_bonds.append(bond)
        return acc_don_bonds

    def covalent_adjacency_matrix(self):
        atom_symbols = self.list_of_atom_symbols()
        bonds = self.covalent_bonds(self.degree_of_covalance())
        molecular_graph = nx.Graph()
        molecular_graph.add_nodes_from(atom_symbols)
        molecular_graph.add_edges_from(bonds)
        return nx.to_numpy_array(molecular_graph)

    def hydrogen_adjacency_matrix(self):
        degofc_table = self.degree_of_covalance()
        _, _, submolecule_symbols = self.detect_submolecules(degofc_table)
        atom_symbols = self.list_of_atom_symbols()
        h_bonds = self.intermolecular_h_bond(degofc_table, submolecule_symbols)
        molecular_graph = nx.Graph()
        molecular_graph.add_nodes_from(atom_symbols)
        molecular_graph.add_edges_from(h_bonds)
        return nx.to_numpy_array(molecular_graph)

    def mu(self):
        degofc = self.degree_of_covalance()
        cov_bonds = self.covalent_bonds(degofc)
        _, _, submolecule_symbols = self.detect_submolecules(degofc)
        len_bonds = len(cov_bonds) + len(
            self.intermolecular_h_bond(degofc, submolecule_symbols)
        )
        return len_bonds - len(self) + 1

    def beta(self):
        degofc = self.degree_of_covalance()
        cov_bonds = self.covalent_bonds(degofc)
        connectivity_c = self.count_connected_components(self.graph_rep(cov_bonds))
        return len(cov_bonds) - len(self) + connectivity_c

    @staticmethod
    def generate_connectivity(bonds):
        connectivity_dict = {}
        for atom1, atom2 in bonds:
            connectivity_dict.setdefault(atom1, []).append(atom2)
            connectivity_dict.setdefault(atom2, []).append(atom1)
        return connectivity_dict

    def generate_angles(self, bonds):
        connectivity_dict = self.generate_connectivity(bonds)
        sym_coords = {atom.symbol: np.array(atom.coordinates) for atom in self}
        angles = []
        linear_angles = []
        for atom, bonded_atoms in connectivity_dict.items():
            bonded_list = list(bonded_atoms)
            cb = sym_coords[atom]
            for i in range(len(bonded_list)):
                for j in range(i + 1, len(bonded_list)):
                    a1, a3 = bonded_list[i], bonded_list[j]
                    ba = sym_coords[a1] - cb
                    bc = sym_coords[a3] - cb
                    cos_a = np.clip(np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc)), -1, 1)
                    angle_deg = np.arccos(cos_a) * 180 / np.pi
                    triple = (a1, atom, a3)
                    if 10 < angle_deg < LINEAR_ANGLE_DEG:
                        angles.append(triple)
                    elif angle_deg >= LINEAR_ANGLE_DEG:
                        linear_angles.append(triple)
                        linear_angles.append(triple)
        return angles, linear_angles
     
    # Make static method for numba compilation of dihedral filtering, we will call this from the main method to generate dihedrals
    @staticmethod
    @njit
    def _filter_dihedrals_numba(coords_array, dihedral_indices, rad_min, rad_linear):
        """
        Numba-compiled angle filtering for dihedrals: keeps a dihedral only if both of its
        bond angles lie in (rad_min, rad_linear); through a linear angle a torsion is undefined

        coords_array: (n, 3) array of atomic coordinates
        dihedral_indices: list of tuples (a, b, c, d) with atom indices for dihedrals
        rad_min: lower angle bound in radians (15 degrees)
        rad_linear: linear-angle cutoff in radians (LINEAR_ANGLE_DEG)
        """
        result_mask = np.zeros(len(dihedral_indices), dtype=np.bool_)

        for i in range(len(dihedral_indices)):
            a_idx = dihedral_indices[i, 0]
            b_idx = dihedral_indices[i, 1]
            c_idx = dihedral_indices[i, 2]
            d_idx = dihedral_indices[i, 3]

            ba = coords_array[a_idx] - coords_array[b_idx]
            bc = coords_array[c_idx] - coords_array[b_idx]
            cb = -bc
            ce = coords_array[d_idx] - coords_array[c_idx]

            cos1 = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
            cos1 = max(-1.0, min(1.0, cos1))

            cos2 = np.dot(cb, ce) / (np.linalg.norm(cb) * np.linalg.norm(ce))
            cos2 = max(-1.0, min(1.0, cos2))

            angle1 = np.arccos(cos1)
            angle2 = np.arccos(cos2)

            if rad_min < angle1 < rad_linear and rad_min < angle2 < rad_linear:
                result_mask[i] = True
        return result_mask


    def generate_dihedrals(self, bonds):
        connectivity_dict = self.generate_connectivity(bonds)
        sym_coords = {atom.symbol: np.array(atom.coordinates) for atom in self}

        # Collect unique dihedrals directly into a set — O(1) membership vs O(k) list check
        seen = set()
        for atom_b, bonded_b in connectivity_dict.items():
            for atom_a in bonded_b:
                for atom_c in bonded_b:
                    if atom_a == atom_c:
                        continue
                    for atom_d in connectivity_dict.get(atom_c, []):
                        if atom_d != atom_b and atom_d != atom_a:
                            d = (atom_a, atom_b, atom_c, atom_d)
                            seen.add(min(d, d[::-1]))
        # Convert to array for numba filtering
        atom_list = list(self.list_of_atom_symbols())
        coords_array = np.array([atoms.coordinates for atoms in self])

        dihedral_indices = []
        dihedral_tuples = []
        for d in seen:
            indices = tuple(atom_list.index(atom) for atom in d)
            dihedral_indices.append(indices)
            dihedral_tuples.append(d)

        dihedral_indices = np.array(dihedral_indices, dtype=np.int32).reshape(-1, 4)
        mask = self._filter_dihedrals_numba(
            coords_array, dihedral_indices, np.radians(15.0), np.radians(LINEAR_ANGLE_DEG)
        )

        return [d for d, keep in zip(dihedral_tuples, mask) if keep]

    def generate_out_of_plane(self, bonds):
        connectivity_dict = self.generate_connectivity(bonds)
        out_of_plane = []
        for atom, bonded_atoms in connectivity_dict.items():
            bonded_list = list(bonded_atoms)
            if len(bonded_list) >= 3:
                for i in range(len(bonded_list)):
                    for j in range(i + 1, len(bonded_list)):
                        for k in range(j + 1, len(bonded_list)):
                            out_of_plane.append((atom, bonded_list[i], bonded_list[j], bonded_list[k]))
                            out_of_plane.append((atom, bonded_list[j], bonded_list[i], bonded_list[k]))
                            out_of_plane.append((atom, bonded_list[k], bonded_list[i], bonded_list[j]))
        return out_of_plane

    def generate_oop_planar_subunits(self, bonds, central_atoms_list):
        hits = []
        connectivity_dict = self.generate_connectivity(bonds)
        for central_atom in central_atoms_list:
            for atom, bonded_atoms in connectivity_dict.items():
                if atom == central_atom[0]:
                    bonded_list = list(bonded_atoms)
                    if len(bonded_list) >= 3:
                        for i in range(len(bonded_list)):
                            for j in range(i + 1, len(bonded_list)):
                                for k in range(j + 1, len(bonded_list)):
                                    hits.append((atom, bonded_list[i], bonded_list[j], bonded_list[k]))
                                    hits.append((atom, bonded_list[j], bonded_list[i], bonded_list[k]))
                                    hits.append((atom, bonded_list[k], bonded_list[i], bonded_list[j]))
        return hits

    def interatomic_distance_matrix(self):
        coords = np.array([a.coordinates for a in self])     # (n, 3)
        diff = coords[:, np.newaxis] - coords[np.newaxis]    # (n, n, 3)
        return np.sqrt((diff ** 2).sum(axis=2))              # (n, n)
        

# ----------------------------
# Profiling and Memory Usage
# ----------------------------

if __name__ == "__main__":
    import time
    import numpy as _np

    # build molecule
    def build_random_molecule(n):
        symbols = ["H", "C", "O", "N", "Cl", "S"]
        mol = Molecule()
        for i in range(n):
            sym = symbols[i % len(symbols)] + str(i + 1)
            coords = tuple(_np.random.rand(3))
            mol.append(Molecule.Atom(sym, coords))
        return mol
    
    def run_benchmarks(mol):
        results = {}
        t0 = time.time(); mol.degree_of_covalance(); results["degree_of_covalance"] = time.time() - t0
        t0 = time.time(); mol.covalent_bonds(mol.degree_of_covalance()); results["covalent_bonds"] = time.time() - t0
        t0 = time.time(); mol.detect_submolecules(); results["detect_submolecules"] = time.time() - t0
        t0 = time.time(); mol.generate_angles(mol.covalent_bonds(mol.degree_of_covalance())); results["generate_angles"] = time.time() - t0
        t0 = time.time(); mol.generate_dihedrals(mol.covalent_bonds(mol.degree_of_covalance())); results["generate_dihedrals"] = time.time() - t0
        return results
    
    molecule_sizes = [10, 15, 20, 30]
    for size in molecule_sizes:
        print(f"Building molecule with {size} atoms...")
        mol = build_random_molecule(size)
        print(f"Running benchmarks for molecule with {size} atoms...")
        benchmark_results = run_benchmarks(mol)
        print(f"Results for {size} atoms: {benchmark_results}")

    # Make a test function that the right dihedrals are generated
    methanol_test = "/home/lme/decomposing-vibrations/test_calculations/xyz_tests/methanol.xyz"

    # Generate mol
    mol = Molecule.from_xyz_file(methanol_test)
    # Generate bonds
    degofc = mol.degree_of_covalance()
    bonds = mol.covalent_bonds(degofc)
    print("Bonds:", bonds)
    # Generate angles
    angles, linear_angles = mol.generate_angles(bonds)
    print("Angles:", angles)
    print("Linear Angles:", linear_angles)
    # Generate dihedrals
    dihedrals = mol.generate_dihedrals(bonds)
    print("Dihedrals:", dihedrals)
    # Generate out of plane
    oop = mol.generate_out_of_plane(bonds)
    print("Out of Plane:", oop) 


