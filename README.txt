# Empirical Relationship between Chemical Structure and Redox Properties: Mathematical Expressions Connecting Structural Features to Energies of Frontier Orbitals and Redox Potentials for Organic Molecules

This repository contains python utility for orbital energies and redox potential prediction. 

Requirements: 
python3.6
numpy

Usage:

The utility can be accessed using a simple python script: 

'''
import sys
sys.path.append('./code/')
from properties import properties
properties = properties(smiles, homo=None, lumo=None)
'''

The utility admits SMILES representation of a molecule as input. HOMO and LUMO values can be provided as optional arguments. 
If optional arguments are not provided, the utility calculates HOMO, LUMO, Oxidation potential and Reduction potential from SMILES input.  
If HOMO/LUMO values are provided, these values are used by the utility to calculate Oxidation potential and Reduction potential. 
Output is a dictionary with predicted properties. The predicted properties can be accessed using
'''
homo = properties['homo'] 
lumo = properties['lumo'] 
oxidation = properties['oxidation']
reduction = properties['reduction'] 
'''




