from __future__ import annotations

import itertools
from typing import NamedTuple
import pprint
import numpy as np
from scipy import constants
import argparse
from mendeleev.fetch import fetch_table
import string

from . import bond
from . import angle
from . import oop
from . import dihedral
from . import molpro_parser
from . import degofc
from .nomodeco_classes import Molecule


# we first generate a connectivity dictionary (which we will later put inside the class)



def eliminate_symmetric_tuples(list_tuples):
    new_list = []
    seen_tuples= set()
    for tp1 in list_tuples:
        reversed_tp1 = tuple(reversed(tp1))
        if reversed_tp1 not in seen_tuples:
           new_list.append(tp1)
           seen_tuples.add(tp1)
    return new_list


def initialize_angles(atoms,bonds):
    hits = []
    linear_hits = []
    for a,b,c in itertools.permutations(atoms,3):
        if (a.symbol,b.symbol) in bonds and (c.symbol,b.symbol) in bonds and not (c.symbol,b.symbol,a.symbol) in hits:
           if atoms.bond_angle(a.symbol,b.symbol,c.symbol)*180/np.pi > 169:
              linear_hits.append((a.symbol,b.symbol,c.symbol)) # 2x degenerate angles
              linear_hits.append((a.symbol,b.symbol,c.symbol))
           else:
              hits.append((a.symbol,b.symbol,c.symbol))
        elif (a.symbol,b.symbol) in bonds and (a.symbol,c.symbol) in bonds and not (b.symbol,a.symbol,c.symbol) in hits:
           if atoms.bond_angle(c.symbol,a.symbol,b.symbol)*180/np.pi > 169:
              linear_hits.append((c.symbol,a.symbol,b.symbol))
              linear_hits.append((c.symbol,a.symbol,b.symbol))
           else:
             hits.append((c.symbol,a.symbol,b.symbol))
        elif (a.symbol,c.symbol) in bonds and (c.symbol,b.symbol) in bonds and not (a.symbol,c.symbol,b.symbol) in hits:
           if atoms.bond_angle(a.symbol,c.symbol,b.symbol)*180/np.pi > 169:
              linear_hits.append((a.symbol,c.symbol,b.symbol))
              linear_hits.append((a.symbol,c.symbol,b.symbol))
           else:
             hits.append((a.symbol,c.symbol,b.symbol)) 
    return eliminate_symmetric_tuples(hits), eliminate_symmetric_tuples(linear_hits)




def initialize_dihedrals(atoms,ic_bonds):
    hits = []
    inproper_dihedrals = []
    for a,b,c,d in itertools.permutations(atoms,4):
        if (a.symbol,b.symbol) in ic_bonds and (b.symbol,c.symbol) in ic_bonds and (c.symbol,d.symbol) in ic_bonds and not (d.symbol,c.symbol,b.symbol,a.symbol) in hits:
          if atoms.bond_angle(a.symbol,b.symbol,c.symbol)*180/np.pi > 169 or  atoms.bond_angle(b.symbol,c.symbol,d.symbol)*180/np.pi > 169:
             inproper_dihedrals.append((a.symbol,b.symbol,c.symbol,d.symbol))
          else:
             if atoms.bond_angle(a.symbol,b.symbol,c.symbol)*180/np.pi > 20 and atoms.bond_angle(b.symbol,c.symbol,d.symbol)*180/np.pi > 20:
                hits.append((a.symbol,b.symbol,c.symbol,d.symbol))
        elif (a.symbol,b.symbol) in ic_bonds and (c.symbol,b.symbol) in ic_bonds and (c.symbol,d.symbol) in ic_bonds and not (d.symbol,c.symbol,b.symbol,a.symbol) in hits:
          if atoms.bond_angle(a.symbol,b.symbol,c.symbol)*180/np.pi > 169 or  atoms.bond_angle(b.symbol,c.symbol,d.symbol)*180/np.pi > 169:
             inproper_dihedrals.append((a.symbol,b.symbol,c.symbol,d.symbol))
          else:
             if atoms.bond_angle(a.symbol,b.symbol,c.symbol)*180/np.pi > 20 and atoms.bond_angle(b.symbol,c.symbol,d.symbol)*180/np.pi > 20:
                hits.append((a.symbol,b.symbol,c.symbol,d.symbol))
        elif (b.symbol,a.symbol) in ic_bonds and (c.symbol,b.symbol) in ic_bonds and (c.symbol,d.symbol) in ic_bonds and not (d.symbol,c.symbol,b.symbol,a.symbol) in hits:
          if atoms.bond_angle(a.symbol,b.symbol,c.symbol)*180/np.pi > 169 or  atoms.bond_angle(b.symbol,c.symbol,d.symbol)*180/np.pi > 169:
             inproper_dihedrals.append((a.symbol,b.symbol,c.symbol,d.symbol))
          else: 
             if atoms.bond_angle(a.symbol,b.symbol,c.symbol)*180/np.pi > 20 and atoms.bond_angle(b.symbol,c.symbol,d.symbol)*180/np.pi > 20:
                hits.append((a.symbol,b.symbol,c.symbol,d.symbol))
        elif (a.symbol,b.symbol) in ic_bonds and (c.symbol,b.symbol) in ic_bonds and (d.symbol,c.symbol) in ic_bonds and not (d.symbol,c.symbol,b.symbol,a.symbol) in hits:
          if atoms.bond_angle(a.symbol,b.symbol,c.symbol)*180/np.pi > 169 or  atoms.bond_angle(b.symbol,c.symbol,d.symbol)*180/np.pi > 169:
             inproper_dihedrals.append((a.symbol,b.symbol,c.symbol,d.symbol))
          else:
             hits.append((a.symbol,b.symbol,c.symbol,d.symbol))
    hits = list(set(hits))
    inproper_dihedrals = (list(set(inproper_dihedrals))) 
    return hits, inproper_dihedrals


def alt_oop(atoms,bonds):
    hits = []
    for a,b,c,d in itertools.permutations(atoms,4):
        if (a.symbol,b.symbol) in bonds and (a.symbol,c.symbol) in bonds and (a.symbol,d.symbol) in bonds and not (a.symbol,b.symbol,d.symbol,c.symbol) in hits:
           hits.append((a.symbol,b.symbol,c.symbol,d.symbol))
        elif (b.symbol, a.symbol) in bonds and (a.symbol,c.symbol) in bonds and (a.symbol,d.symbol) in bonds and not (a.symbol,b.symbol,d.symbol,c.symbol) in hits:
           hits.append((b.symbol,a.symbol,c.symbol,d.symbol))
        elif (b.symbol, a.symbol) in bonds and (c.symbol,a.symbol) in bonds and (a.symbol,d.symbol) in bonds and not (a.symbol,b.symbol,d.symbol,c.symbol) in hits:
           hits.append((a.symbol,b.symbol,c.symbol,d.symbol))             
        elif (b.symbol,a.symbol) in bonds and (c.symbol,a.symbol) in bonds and (d.symbol,a.symbol) in bonds and not (a.symbol,b.symbol,d.symbol,c.symbol) in hits:
           hits.append((a.symbol,b.symbol,c.symbol,d.symbol))
        elif (b.symbol,a.symbol) in bonds and (a.symbol,c.symbol) in bonds and (d.symbol,a.symbol) in bonds and not (a.symbol,c.symbol,b.symbol,d.symbol) in hits:
           hits.append((a.symbol,c.symbol,d.symbol,b.symbol))
    return list(set(hits))


def initialize_alt_oop_planar_subunits(atoms,bonds,central_atoms_list):
    hits = []
    for central_atom in central_atoms_list:
        for a,b,c,d in itertools.permutations(atoms, 4):
            if a.symbol == central_atom[0] and (a.symbol,b.symbol) in bonds and (a.symbol,c.symbol) in bonds and (a.symbol,d.symbol) in bonds  and not (a.symbol, b.symbol, d.symbol, c.symbol) in hits:
                hits.append((a.symbol, b.symbol, c.symbol, d.symbol))
            elif a.symbol == central_atom[0] and (b.symbol,a.symbol) in bonds and (a.symbol,c.symbol) in bonds and (a.symbol,d.symbol) in bonds and not (a.symbol,b.symbol,d.symbol,c.symbol) in hits:
                hits.append((a.symbol, b.symbol, c.symbol, d.symbol))
    return hits

