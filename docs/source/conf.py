# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
from pathlib import Path

# Autodoc Inport
#current_dir = os.path.abspath(os.path.dirname(__file__))
#parent_dir = os.path.abspath(os.path.join(current_dir, '../..'))
#module_path = os.path.join(parent_dir,"nomodeco","libraries")
#nomodeco_path = os.path.join(parent_dir,"nomdeco")
#sys.path.append(module_path)
##sys.path.append(nomodeco_path)
#sys.path.append(str((Path(__file__).parents[2].absolute() / "nomodeco")))


# Set the path relative to the location of conf.py
current_dir = Path(__file__).parent.resolve()  # Directory where conf.py is located
module_path = current_dir / '../../nomodeco/libraries'
nomodeco_path = current_dir / '../../nomodeco'

# Convert Path objects to strings and append them to sys.path
sys.path.append(str(module_path.resolve()))
sys.path.append(str(nomodeco_path.resolve()))


project = 'nomodeco'
copyright = '2024, Lukas'
author = 'Lukas'
release = '0.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc"] 

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
