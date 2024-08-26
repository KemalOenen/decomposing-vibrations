# nomodeco.py - a normal mode decomposition tool
A Decomposition scheme for the force constants and frequencies of molecular vibrations, which can be used to quantify the vibrational problem for larger systems. 
Internal coordinate sets are automatically generated and optimized for optimal seperation of the normal modes in terms of internal coordinates.

## Installation

Clone the Github repository

```
git clone https://github.com/KemalOenen/decomposing-vibrations.git
```
Install the requirements for python

```
pip install -r requirements.txt
```
Generate an alias on you machine to run repository

```
alias nomodeco = 'python -W ignore /path/to/local/repo/clone'
```

## Usage

To run Nomodeco.py a Molpro output file containing the analytical hessian matrix and the atom coordinates is needed, the file template_molpro.inp contains a basic input script for a Molpro calculation that writes the output file needed for the nomodeco algorithm. 

With the Molpro.out file one can subsequently using the alias run Nomodeco.py

```
nomodeco --args molpro.out
```
For a list of all arguments use

```
nomodeco --help
```


## Publication for further information:
Kemal Oenen, Dennis F. Dinu, Klaus R. Liedl; Determining internal coordinate sets for optimal representation of molecular vibration. J. Chem. Phys. 7 January 2024; 160 (1): 014104. https://doi.org/10.1063/5.0180657

