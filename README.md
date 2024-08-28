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

## Example

Here as an example the Normal Mode Decomposition for the Water Dimer is presented.

```
molpro -n6 template_molpro.inp
```
Using the Molpro.out file one can now carry out the automatic normal mode decompositon and generate a contribution table

```
nomodeco --log --heatmap contr --latex_tab template_molpro.out
```


## Command Line Arguments:

In this section the general usage and the arguments for tthe command line are further explained:

```
nomodeco.py [-h] [--log] [--matrix_opt [matrix]] [--penalty1 [INTFREQ-PENALTY]] [--penalty2 [INTFC-PENALTY]] [--heatmap matrix [matrix ...]] [--csv matrix [matrix...]] [--latex_tab]
``` 
+ [-h] --help shows the help message for the arguments of Nomodeco.py
+ [--log] displays an additional .log file that containes results for every tested internal coordinate set
+ [--matrix_opt [matrix]] chooses which matrix is used for the optimization: i.e VED matrix (keyword: ved) diagonal elements of PED (keyword: diag) and Contribution Table (keyword: contr) (default: contr)
+ [--penalty1 [INTFREQ-PENALTLY]] penalty value for asymmetric intrinsic frequency values, which can be helpful for cyclic systems (default: 0)
+ [--penalty2 [INTFC-PENALTY]] penalty value for unphysical contributons per internal coordinate (default: 0)
+ [--heatmap matrix [matrix ...]] return a heatmap for the specified matrix: i.e, VED matrix (keyword: ved), Diagonal elements of PED matrix (keyword: diag) and Contribution Table (keyword: contr)
+ [--csv matrix [matrix...]] return a csv file for the specified matrix: i.e VED matrix (keyword: ved), Diagonal elements of PED matrix (keyword: diag) and Contribution Table (keyword: contr)
+ [--latextab] generate additional latex table file which displayes the Contribution Table




## Publication for further information:
Kemal Oenen, Dennis F. Dinu, Klaus R. Liedl; Determining internal coordinate sets for optimal representation of molecular vibration. J. Chem. Phys. 7 January 2024; 160 (1): 014104. https://doi.org/10.1063/5.0180657

