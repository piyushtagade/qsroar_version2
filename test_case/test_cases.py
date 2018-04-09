# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 12:45:29 2018

@author: p.tagade
"""
import sys
sys.path.append('../code/')

import numpy as np 
import matplotlib.pyplot as plt 
from properties import properties
import json 
# -----------------------------------------------------------------------------
fjson = 'PAH_data.jsn' 
with open(fjson, 'r') as fid:  
    mol_properties = json.load(fid) 
# -----------------------------------------------------------------------------    
smiles = mol_properties['PAH']['SMILES']
oxid = mol_properties['PAH']['Oxid']
red = mol_properties['PAH']['Red']
homo = mol_properties['PAH']['Homo']
lumo = mol_properties['PAH']['Lumo']
# -----------------------------------------------------------------------------

cases = []

homo_pred = np.zeros(len(smiles))
lumo_pred = np.zeros(len(smiles))
oxidation_pred = np.zeros(len(smiles))
reduction_pred = np.zeros(len(smiles))

# Test case 1: Predicting only from SMILES
for i in range(0, len(smiles)):
    mol = smiles[i]
    prediction = properties(mol)
    lumo_pred[i] = prediction['lumo']
    homo_pred[i] = prediction['homo']
    oxidation_pred[i] = prediction['oxidation']
    reduction_pred[i] = prediction['reduction']
    
cases.append([lumo_pred, homo_pred, oxidation_pred, reduction_pred])

homo_pred = np.zeros(len(smiles))
lumo_pred = np.zeros(len(smiles))
oxidation_pred = np.zeros(len(smiles))
reduction_pred = np.zeros(len(smiles))

# Test case 2: SMILES and lumo as input
for i in range(0, len(smiles)):
    mol = smiles[i]
    prediction = properties(mol, homo = homo[i])
    lumo_pred[i] = prediction['lumo']
    homo_pred[i] = prediction['homo']
    oxidation_pred[i] = prediction['oxidation']
    reduction_pred[i] = prediction['reduction']

cases.append([lumo_pred, homo_pred, oxidation_pred, reduction_pred]) 

homo_pred = np.zeros(len(smiles))
lumo_pred = np.zeros(len(smiles))
oxidation_pred = np.zeros(len(smiles))
reduction_pred = np.zeros(len(smiles))
# Test case 3: SMILES and homo-lumo as input 
for i in range(0, len(smiles)):
    mol = smiles[i]
    prediction = properties(mol, lumo = lumo[i], homo = homo[i])
    lumo_pred[i] = prediction['lumo']
    homo_pred[i] = prediction['homo']
    oxidation_pred[i] = prediction['oxidation']
    reduction_pred[i] = prediction['reduction']

cases.append([lumo_pred, homo_pred, oxidation_pred, reduction_pred]) 

true = [np.array(lumo), np.array(homo), np.array(oxidation_pred), np.array(reduction_pred)]


fig = plt.figure() 

# For test case 1
for icase in range(0, 4):
    fig.add_subplot(3, 4, icase+1)
    h1 = plt.plot(true[icase], cases[0][icase], 'b.')
    h2 = plt.plot(true[icase], true[icase], 'k-')

# For test case 2
for icase in range(1, 4):
    fig.add_subplot(3, 4, 4 + icase+1)
    h1 = plt.plot(true[icase], cases[1][icase], 'b.')
    h2 = plt.plot(true[icase], true[icase], 'k-')

# For test case 3
for icase in range(2, 4):
    fig.add_subplot(3, 4, 8 + icase+1)
    h1 = plt.plot(true[icase], cases[2][icase], 'b.')
    h2 = plt.plot(true[icase], true[icase], 'k-')

