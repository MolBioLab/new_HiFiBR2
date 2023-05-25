#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 10:07:41 2023

@author: mobilab5
"""

# LIBRARY
import sys
import subprocess
import pkg_resources
from shutil import which

# Check if programs already installed on your environment 
required_tools = {'fastqc', 'samtools', 'bamtools', 'trimmomatic', 'flash2', 'geneious'}
missing_tools = {tool for tool in required_tools if which(tool) is None}
 
with open("missing_packages.txt", "w") as f:
    if missing_tools:
        f.write('***MISSING TOOLS***\n')
        f.write("\n".join(missing_tools) + "\n\n")
    else:
        f.write("All required tools are already installed" + "\n\n")
    
# Check if programs already installed on your Python  
required_pyPackages = {'pysam', 'pandas', 'numpy', 'scikit-learn', 'scipy', 
                       'diff-match-patch', 'tqdm', 'xlsxwriter'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required_pyPackages - installed

with open("missing_packages.txt", "+a") as f:
    if missing: 
        f.write('***MISSING PYTHON PACKAGES***\n')
        f.write("\n".join(missing))
    else:
        f.write("All required Python packages are already installed")
    