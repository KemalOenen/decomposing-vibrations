from __future__ import annotations

from typing import Iterable
import string
import numpy as np
from mendeleev.fetch import fetch_table

def get_bond_information() -> pd.DataFrame:
    df = fetch_table('elements')

    bond_info = df.loc[:, ['symbol', 'covalent_radius_pyykko', 'covalent_radius_pyykko_double',
        'covalent_radius_pyykko_triple']]

    bond_info.set_index("symbol", inplace=True)
    bond_info /= 100
    return bond_info

BOND_INFO = get_bond_information()

def actual_length(coord_a: Iterable, coord_b: Iterable) -> float:
    return np.linalg.norm(np.array(coord_b) - np.array(coord_a))

#TODO: remove this horrible if statement
def theoretical_length(symbol_a: str, symbol_b: str) -> pd.Series:
    if symbol_a.strip(string.digits) == 'D':
        rad_a = BOND_INFO.loc['H']
        rad_b = BOND_INFO.loc[symbol_b.strip(string.digits)]
    elif symbol_b.strip(string.digits) == 'D':
        rad_b = BOND_INFO.loc['H']
        rad_a = BOND_INFO.loc[symbol_a.strip(string.digits)]
    else:
        rad_a = BOND_INFO.loc[symbol_a.strip(string.digits)]
        rad_b = BOND_INFO.loc[symbol_b.strip(string.digits)]
    return  rad_a + rad_b 


def is_valid(atom_a: Atom, atom_b: Atom, thr: float = 0.3) -> bool:
    actual = actual_length(atom_a.coordinates, atom_b.coordinates)
    theoretical = theoretical_length(atom_a.symbol, atom_b.symbol)

    return any(abs(theoretical - actual) < thr)
