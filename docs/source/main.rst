Introduction
"""""""""""
Welcome to the Nomodeco.py documentation!
=========================================

This is the online documentation for the python package Nomoceco.py.

Nomodeco.py enables the automatic determination of an optimal coordinate set for a given molecular structure. Using the atomic coordinates of the given molecule Nomodeco constructs all possible primitive internal coordinates:

* bond stretching
* in-plane angle bending
* linear valance bending
* torsion
* out-of-plane angle bending

The transformation between cartesian and internal coordinates is given by the *Wilson B Matrix*, Nomodeco initializes all possible internal coordinates sets and selects the optimal set based on topolgy and symmetry considerations.

The advantage of Nomodeco is that molecular vibrations can be represented using primitive internal coordinates, further the contribution of each internal coordinate to a particular normal mode of vibration can be visualized.



For further information on the algorithm an theory: `Determining internal coordinate sets for optimal representation of molecular vibration`_

.. _Determining internal coordinate sets for optimal representation of molecular vibration: https://pubs.aip.org/aip/jcp/article/160/1/014104/2932467/Determining-internal-coordinate-sets-for-optimal




Special Thanks
--------------

I would like to extend my deepest thanks to:

- Kemal Önen for initializing the project and guiding the future development 
- Klaus R. Liedl for the supervision of the project
- Dennis F. Dinu for important contributions of theory aspects
- Leonardo Pedri for the sfinal setup and structuring of the python package