import numpy as np
from nomodeco.libraries.nomodeco_classes import Molecule


def replace_vars(vlist, variables):
    """ Replaces a list of variable names (vlist) with their values
        from a dictionary (variables).
    """
    for i, v in enumerate(vlist):
        if v in variables:
            vlist[i] = variables[v]
        else:
            try:
                # assume the "variable" is a number
                vlist[i] = float(v)
            except:
                print("Problem with entry " + str(v))

def read_zmatrix(file_path):
    """
    Reads in a z-matrix in standard format, returning a list of atoms and coordinates
    """
    zmatf = open(file_path, 'r')
    atomnames = []
    rconnect = []  # bond connectivity
    rlist = []     # list of bond length values
    aconnect = []  # angle connectivity
    alist = []     # list of bond angle values
    dconnect = []  # dihedral connectivity
    dlist = []     # list of dihedral values
    variables = {} # dictionary of named variables

    for line in zmatf:
        words = line.split()
        eqwords = line.split("=")
        if len(eqwords) > 1:
            # these are the named variables
            varname = str(eqwords[0]).strip()
            # Make some exeption handeling
            try:
                varval = float(eqwords[1])
                variables[varname] = varval
            except:
                print("Invalid variable definition: " + line)
        else:
            # no variable, just number
            if len(words) > 0:
                atomnames.append(words[0])
            if len(words) > 1:
                rconnect.append(int(words[1]))
            if len(words) > 2:
                rlist.append(words[2])
            if len(words) > 3:
                aconnect.append(int(words[3]))
            if len(words) > 4:
                alist.append(words[4])
            if len(words) > 5:
                dconnect.append(int(words[5]))
            if len(words) > 6:
                dlist.append(words[6])
    # replace all the variables
    replace_vars(rlist,variables)
    replace_vars(alist,variables)
    replace_vars(dlist,variables)
    return (atomnames,rconnect,rlist,aconnect,alist, dconnect,dlist)

def write_xyz(atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist):
    """Prints out an xyz file from a decomposed z-matrix"""
    npart = len(atomnames)
        
    # put the first atom at the origin
    xyzarr = np.zeros([npart, 3])
    if (npart > 1):
        # second atom at [r01, 0, 0]
        xyzarr[1] = [rlist[0], 0.0, 0.0]

    if (npart > 2):
        # third atom in the xy-plane
        # such that the angle a012 is correct 
        i = rconnect[1] - 1
        j = aconnect[0] - 1
        r = rlist[1]
        theta = alist[0] * np.pi / 180.0
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        a_i = xyzarr[i]
        b_ij = xyzarr[j] - xyzarr[i]
        if (b_ij[0] < 0):
            x = a_i[0] - x
            y = a_i[1] - y
        else:
            x = a_i[0] + x
            y = a_i[1] + y
        xyzarr[2] = [x, y, 0.0]

    for n in range(3, npart):
        # back-compute the xyz coordinates
        # from the positions of the last three atoms
        r = rlist[n-1]
        theta = alist[n-2] * np.pi / 180.0
        phi = dlist[n-3] * np.pi / 180.0
        
        sinTheta = np.sin(theta)
        cosTheta = np.cos(theta)
        sinPhi = np.sin(phi)
        cosPhi = np.cos(phi)

        x = r * cosTheta
        y = r * cosPhi * sinTheta
        z = r * sinPhi * sinTheta
        
        i = rconnect[n-1] - 1
        j = aconnect[n-2] - 1
        k = dconnect[n-3] - 1
        a = xyzarr[k]
        b = xyzarr[j]
        c = xyzarr[i]
        
        ab = b - a
        bc = c - b
        bc = bc / np.linalg.norm(bc)
        nv = np.cross(ab, bc)
        nv = nv / np.linalg.norm(nv)
        ncbc = np.cross(nv, bc)
        
        new_x = c[0] - bc[0] * x + ncbc[0] * y + nv[0] * z
        new_y = c[1] - bc[1] * x + ncbc[1] * y + nv[1] * z
        new_z = c[2] - bc[2] * x + ncbc[2] * y + nv[2] * z
        xyzarr[n] = [new_x, new_y, new_z]
    
    molecule = []

    for i in range(npart):
        atom = Molecule.Atom(atomnames[i], tuple(xyzarr[i]))
        molecule.append(atom)
    return molecule

