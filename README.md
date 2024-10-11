# Nomodeco.py a normal mode decomposition tool


Nomodeco.py enables the automatic determination of an optimal internal coordinate set for a given molecular structure. Using the atomic coordinates of the given molecule Nomodeco.py constructs all possible primitive internal coordinates.

## Documentation

The following links provide further information and documentation about the package

[decomposing vibrations-docs](https://kemaloenen.github.io/decomposing-vibrations/)


## How to Use:

Nomodeco can be installed using pip:

```
pip install nomodeco
```

Make sure to use a python version >=3.12.5!

For further information and a version history see:

[Nomodeco on PyPI](https://pypi.org/project/nomodeco/)

Alternatively, Nomodeco can be installed using the [github-repository](https://github.com/KemalOenen/decomposing-vibrations). In the repository a enviroment.yml file is included which can be used to create a conda enviroment with python 3.12.5 and poetry installed.

```
conda env create -f enviroment.yml
```

Then using poetry the additional packages can be installed:

```
poetry update
```

With the usage of the enviroment.yml file also pymolpro is installed which runs jobs on the local molpro installation.



## Publication for further information:
Kemal Oenen, Dennis F. Dinu, Klaus R. Liedl; Determining internal coordinate sets for optimal representation of molecular vibration. J. Chem. Phys. 7 January 2024; 160 (1): 014104. https://doi.org/10.1063/5.0180657

