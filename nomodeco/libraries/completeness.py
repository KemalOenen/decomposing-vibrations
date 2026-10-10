""" 
Completeness check for internal coordinates


An IC set is complete if its B matrix has
rank equal to idof (it can describe ever internal motion)

Goal of this module is to collect several ways to test this


1. is_complete_svd  -- computes rank via SVD decomposition of B matrix
2. is_complete_gram -- eigenvalues of the idof x idof Gram matrix
3. complete_mask_product -- all sets of a LazyIcDict at once, using its product structure
"""
import time 

import numpy as np
from nomodeco.libraries import bmatrix

KINDS = ["bonds", "angles", "linear valence angles", "out of plane angles", "dihedrals"]

#---------
# Helper Methods
#---------

def set_b_matrix(atoms, ic_set, bend_refs=None):
    """B matrix (n_ic x 3N) of one IC set."""
    return bmatrix.b_matrix(atoms, *[ic_set[k] for k in KINDS], 0, bend_refs)

def numerical_rank(singular_values, tol=1e-8):
    """  
    Number of Singular Values above rtol
    """
    if len(singular_values) == 0 or singular_values[0] == 0:
        return 0
    return int(np.sum(singular_values > tol * singular_values[0]))

