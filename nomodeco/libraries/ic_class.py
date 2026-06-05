from __future__ import annotations

import pandas as pd

from mendeleev.fetch import fetch_table

from nomodeco.libraries.molecule_class import Molecule

__all__ = ["Molecule", "InternalCoordinates", "get_mass_information"]


# TODO: isotopes for all elements with command line input?
def get_mass_information():
    df = fetch_table("elements")
    mass_info = df.loc[:, ["symbol", "atomic_weight"]]
    deuterium_info = pd.DataFrame({"symbol": ["D"], "atomic_weight": [2.014102]})
    mass_info = pd.concat([mass_info, deuterium_info])
    mass_info.set_index("symbol", inplace=True)
    return mass_info


class InternalCoordinates:
    def __init__(self):
        self.coordinates = {}

    def add_coordinate(self, key, coordinate_list):
        self.coordinates[key] = coordinate_list

    def get_coordinate(self, key):
        return self.coordinates.get(key, None)

    def common_coordinate(self, key, coordinate_list):
        stored_set = set(self[key])
        return [i for i in coordinate_list if i in stored_set]

    def add_coord_diff(self, key, ic_list1, ic_list2):
        if not ic_list1:
            self.coordinates[key] = []
            return

        set1 = set(ic_list1)
        set2 = set(ic_list2)

        if len(ic_list1[0]) == 3:
            sym_ic_set_2 = {(c, b, a) for a, b, c in ic_list2}
            self.coordinates[key] = list(set1 - set2 - sym_ic_set_2)
        else:
            self.coordinates[key] = list(set1 - set2)

    def add_coord_diff_linear(self, key, ic_list1, ic_list2):
        sym_ic_set_2 = {(c, b, a) for a, b, c in ic_list2}
        set1 = set(ic_list1)
        set2 = set(ic_list2)
        diff = list(set1 - set2 - sym_ic_set_2)
        self.coordinates[key] = diff * 2

    def __str__(self):
        coords_str = ", ".join(
            f"{name}: {value}" for name, value in self.coordinates.items()
        )
        return f"InternalCoordinates({coords_str})"

    def __getitem__(self, index):
        return self.coordinates[index]
    
if __name__ == "__main__":
    # Import Tracemalloc and test memory usage 
    import tracemalloc
    tracemalloc.start()
    # Test the internal coordinates class
    ic = InternalCoordinates()
    # Add big numbers of ICs to  test memory usage
    ic.add_coordinate("bonds", [(i, i+1) for i in range(10**6)])
    print(ic)
    # Check memory usage    
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage: {current / 10**6} MB; Peak memory usage: {peak / 10**6} MB")
    tracemalloc.stop()
