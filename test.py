# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 12:45:43 2018

@author: p.tagade
"""

import sys
sys.path.append('./code/')
from properties import properties
smiles = 'C1=CC=CC=C1'
properties = properties(smiles, homo=None, lumo=None)

homo = properties['homo'] 
lumo = properties['lumo'] 
oxidation = properties['oxidation']
reduction = properties['reduction']