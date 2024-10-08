from cmath import isclose
import numpy as np
import scipy
import jax
import jax.numpy as jnp
from jax import grad
from jax import jacfwd

def bond_vector_norm(r_a,r_b):
    bond_length = jnp.linalg.norm(r_b - r_a)

    return (r_b - r_a)/bond_length # bond vector


def calculate_bond_vector(bond,atom_list):
    bond_row = np.zeros((len(atom_list)*3))
    atom1_idx, atom2_idx = bond
    
    i1 = atom_list.retrieve_index_in_list(atom1_idx)*3
    i2 = atom_list.retrieve_index_in_list(atom2_idx)*3

    atom1_coords = jnp.array(atom_list.get_atom_coords_by_symbol(atom1_idx))
    atom2_coords = jnp.array(atom_list.get_atom_coords_by_symbol(atom2_idx))
    bond_length = jnp.linalg.norm(atom1_coords - atom2_coords)
    
    # calculate bond vector
    bond_vector = (atom2_coords - atom1_coords) / bond_length
    
    # for a --> negative
    neg_bond_vector = bond_vector * -1

    bond_row[i2:(i2+3)] = bond_vector
    bond_row[i1:(i1+3)] = neg_bond_vector
    
    return bond_row
    

def bond_rows(atoms, bonds):
    bond_rows = []
    for bond in bonds:
        bond_row = calculate_bond_vector(bond,atoms)
        bond_rows.append(bond_row)
    return jnp.array(bond_rows)

def angle_func(a,b,c):
    return jnp.arccos(((a[0]-b[0])*(c[0]-b[0]) + (a[1]-b[1])*(c[1]-b[1]) + (a[2]-b[2])*(c[2]-b[2])) /
            (jnp.sqrt((a[0] - b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) * 
                jnp.sqrt((c[0]-b[0])**2 + (c[1]-b[1])**2 + (c[2]-b[2])**2)))


def calculate_angle_vector(angle, atom_list):
    atom1_idx, atom2_idx, atom3_idx = angle

    a_coords = jnp.array(atom_list.get_atom_coords_by_symbol(atom1_idx))
    b_coords = jnp.array(atom_list.get_atom_coords_by_symbol(atom2_idx))
    c_coords = jnp.array(atom_list.get_atom_coords_by_symbol(atom3_idx))
    
    # Calculate_Gradient
    grad_a = grad(angle_func,argnums=0)(a_coords,b_coords,c_coords)
    grad_b = grad(angle_func,argnums=1)(a_coords,b_coords,c_coords)
    grad_c = grad(angle_func,argnums=2)(a_coords,b_coords,c_coords)
    
    # ahh so nice

    bond_row = np.zeros(len(atom_list*3))

    i1 = atom_list.retrieve_index_in_list(atom1_idx)*3
    i2 = atom_list.retrieve_index_in_list(atom2_idx)*3
    i3 = atom_list.retrieve_index_in_list(atom3_idx)*3

    bond_row[i1:(i1+3)] = grad_a
    bond_row[i2:(i2+3)] = grad_b
    bond_row[i3:(i3+3)] = grad_c
    return bond_row

def angle_rows(atoms,angles):
    angle_rows = []
    for angle in angles:
        angle_row = calculate_angle_vector(angle, atoms)
        angle_rows.append(angle_row)
        
    return jnp.array(angle_rows)

def dihedral_func(a,b,c,d):
    angle = jnp.dot(jnp.cross((b-a),(b-c)),jnp.cross((c-d),(b-c))) / (jnp.linalg.norm(
                jnp.cross((b-a),(b-c))) * jnp.linalg.norm(jnp.cross((c-d),(b-c))))
    print(angle)
    return jnp.arccos(jnp.dot(jnp.cross((b-a),(b-c)),jnp.cross((c-d),(b-c))) / (jnp.linalg.norm(
                jnp.cross((b-a),(b-c))) * jnp.linalg.norm(jnp.cross((c-d),(b-c))))) 
#jnp.sign(jnp.dot((b-c), jnp.cross(jnp.cross((b-a),(c-b)),jnp.cross((b-c),(d-c))))) * 
    #    return jnp.arccos( ( (((b[1]-a[1])*(b[2]-c[2]) - (b[2]-a[2])*(b[1]-c[1])) * ((b[1]-c[1])*(d[2]-c[2]) - (b[2]-c[2])*(d[1]-c[1]))) +
#           (((b[2]-a[2])*(b[0]-c[0]) - (b[0]-a[0])*(b[2]-c[2])) * ((b[2]-c[2])*(d[0]-c[0]) - (b[0] - c[0])*(d[2]-c[2]))) + 
#           (((b[0]-a[0])*(b[1]-c[1]) - (b[1]-a[1])*(b[0]-c[0])) * ((b[0]-c[0])*(d[1]-c[1]) - (b[1] - c[1])*(d[0]-c[0])))) / (jnp.sqrt( ( (b[1]-a[1])*(b[2]-c[2]) 
#               - (b[2]-a[2])*(b[1]-c[1]) )**2 + ( (b[2]-a[2])*(b[0]-c[0]) - (b[0]-a[0])*(b[2]-c[2]))**2 +
#               ( (b[0]-a[0])*(b[1]-c[1]) - (b[1]-a[1])*(b[0]-c[0]))**2) * 
#               jnp.sqrt( ((b[1]-c[1])*(d[2]-c[2]) - (b[2]-c[2])*(d[1]-c[1]))**2 + ((b[2]-c[2])*(d[0]-c[0]) - (b[0] - c[0])*(d[2]-c[2]))**2 +
#                        ((b[0]-c[0])*(d[1]-c[1]) - (b[1] - c[1])*(d[0]-c[0]))**2)))

def calculate_dihedral_row(dihedral,atom_list):
    atom1_idx, atom2_idx, atom3_idx, atom4_idx = dihedral
    
    dihedral_row = np.zeros((len(atom_list)*3))
    print(dihedral)
    a_coords = jnp.array(atom_list.get_atom_coords_by_symbol(atom1_idx))
    b_coords = jnp.array(atom_list.get_atom_coords_by_symbol(atom2_idx))
    c_coords = jnp.array(atom_list.get_atom_coords_by_symbol(atom3_idx))
    d_coords = jnp.array(atom_list.get_atom_coords_by_symbol(atom4_idx))
   
    a = dihedral_func(a_coords,b_coords,c_coords,d_coords)
#    grad_a = grad(dihedral_func, argnums=0)(a_coords,b_coords,c_coords,d_coords)
#    grad_b = grad(dihedral_func, argnums=1)(a_coords,b_coords,c_coords,d_coords)
#    grad_c = grad(dihedral_func, argnums=2)(a_coords,b_coords,c_coords,d_coords)
#    grad_d = grad(dihedral_func, argnums=3)(a_coords,b_coords,c_coords,d_coords)
#
#    i1 = atom_list.retrieve_index_in_list(atom1_idx)*3
#    i2 = atom_list.retrieve_index_in_list(atom2_idx)*3
#    i3 = atom_list.retrieve_index_in_list(atom3_idx)*3
#    
#    dihedral_row[i1:(i1+3)] = grad_a
#    dihedral_row[i2:(i2+3)] = grad_b
#    dihedral_row[i3:(i3+3)] = grad_c
#    print(dihedral_row)

def dihedral_rows(atoms,dihedrals):
    dihedral_rows = []
    for dihedral in dihedrals:
        dihedral_row = calculate_dihedral_row(dihedral,atoms)


def b_matrix(atoms,bonds,angles,dihedrals):
    bond_mat = bond_rows(atoms, bonds)
    angle_mat = angle_rows(atoms,angles)
    dihedral_mat = dihedral_rows(atoms,dihedrals)
    
    return jnp.vstack((bond_mat,angle_mat))
    


