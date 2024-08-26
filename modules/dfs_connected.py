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

class Atom(NamedTuple):
    symbol: str
    coordinates: tuple

# For this algorithm we perform DFS seaerch in order to optain the number c the connectivity of the underlying atoms


# This function gives us the Graph as a dictionary for example:
# {'O1': ['H4','H1','H2']} gives the atom as a key and the "bonds" as values of the given key

def graph_rep(bonds):
    graph = {}
    for bond in bonds:
        a, b = bond 
        graph.setdefault(a,[]).append(b)
        graph.setdefault(b,[]).append(a)
    return graph

# Small description of the DFS Algorithm:
# We start at a arbitrary and then explore each branch
# visit a point and put its adjecent nodes into a stack
# now we go for example to the next node and throw all connected notes into the stack


def dfs(graph,start,visited):
    visited.add(start)
    for neighbor in graph.get(start,[]):
        if neighbor not in visited:
           dfs(graph,neighbor,visited) # Recursive Algorithm

def is_connected(graph):
    if not graph:
       return True # Not really necessary
    visited = set()
    start_node = next(iter(graph)) #in this case next specifices the first element of the dictionary iter just returns the iterator
    dfs(graph,start_node,visited)
    return len(visited) == len(graph)


# with a similar function we can now evaluate the number of connected components c

def count_connected_components(graph):
    if not graph:
       return 0
    visited = set()
    count = 0
    for node in graph:
        if node not in visited:
           count += 1
           dfs(graph,node,visited)
    return count
