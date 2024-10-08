Installation
==========

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


.. note::
    Hopefully nomedeco will in future be available on pip