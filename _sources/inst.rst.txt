Installation
=============

.. note::
    Nomodeco.py is now on `pypi`_ and can be installed using pip 


.. _pypi: https://pypi.org/project/nomodeco/ 


In your conda or virtual enviroment make sure to use **python 3.12.5** or higher. Then Nomodeco can be installed using::
    
    pip install nomodeco

All the packages are automatically installed in your virtual enviroment now nomodeco can be run using the alias **nomodeco**::

    nomodeco --molpro h2o.out --log --heatmap contr



**Github Installation**

Nomodeco can be manually installed from the github repository::

    git clone https://github.com/KemalOenen/decomposing-vibrations.git
    cd decomposing-vibrations
    pip install .

**Conda installation**

For the usage of the pymolpro integration one needs to generate a conda enviroment using the *eviroment.yml* file of the nomodeco package.

Again this can be done via cloning of the github repository::

    git clone https://github.com/KemalOenen/decomposing-vibrations.git
    cd decomposing-vibrations
    conda env create -f enviroment.yml
    poetry install

Now the package pymolpro should be loaded into the conda enviroment this can be checked via::

    conda list

Now also the conda-enviroment knows the alias and a Nomodeco calculation can be run, for example::

    nomodeco --gv h2o.log --heatmap contr

