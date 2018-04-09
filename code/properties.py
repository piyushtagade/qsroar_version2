# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 11:47:15 2018

@author: p.tagade
"""
import os
import numpy as np
import json
from extract_features import * 
import platform

# -----------------------------------------------------------------------------
def properties(smiles, lumo=None, homo=None):
# -----------------------------------------------------------------------------    
    ''' Function for orbital energies and redox potentials prediction. 
    The function accepts smiles of an organic molecule as input and provides 
    orbital energies and redox potential as output. 
    
    Inputs: 
            smiles - smiles representation of organic molecule
    Optional inputs:
            lumo - Lowest Unoccopied Molecular Orbital Energy in eV.
            homo - Highest Occuped Molecular Orbital Energy in eV.
    Output:
            output_properties - A dictionary of output properties         
    '''
# -----------------------------------------------------------------------------    
    property_limits = {}
    property_limits['oxidation'] = [-0.485437, 3.68406]
    property_limits['reduction'] = [-4.89283, 0.201063]
    property_limits['bandgap']   = [-9.06482453, -1.06047362]
    property_limits['lumo']      = [-3.91109452, 0.70178201]
    property_limits['homo']      = [-8.495, -4.625]
# Output is sent in a dictionary
    output_properties = {}
# -----------------------------------------------------------------------------        
    functional_groups = create_functional_groups()    
# -----------------------------------------------------------------------------
# Extracting features from SMILES representation
# -----------------------------------------------------------------------------    
    smiles_features = extract_features(smiles, functional_groups)    
# ----------------------------------------------------------------------------- 
# Loading coefficients of the correlation
# -----------------------------------------------------------------------------
    filepath = os.path.join(os.path.dirname(__file__), 'coefficients.jsn')

 
    with open(filepath, 'r') as fp:
        coefficients = json.load(fp)            
# -----------------------------------------------------------------------------    
# Predicting lumo energy
# -----------------------------------------------------------------------------
    lim_lumo = np.array(property_limits['lumo'])
    if lumo == None:        
        lumo_coefs = coefficients['lumo']
        lumo = predict(smiles_features, lumo_coefs)                
        lumo = lumo*(lim_lumo[1] - lim_lumo[0]) + lim_lumo[0]
    output_properties['lumo'] = lumo     
    smiles_features['lumo'] = -2.0 + ((4.0)/(lim_lumo[1] - lim_lumo[0])) * (lumo - lim_lumo[0]);
# -----------------------------------------------------------------------------    
# Predicting homo energy
# -----------------------------------------------------------------------------                           
    lim_bandgap = np.array(property_limits['bandgap']) 
        
    if homo == None:        
        homo_coefs = coefficients['homo']
        bandgap = predict(smiles_features, homo_coefs)            
        bandgap = bandgap*(lim_bandgap[1] - lim_bandgap[0]) + lim_bandgap[0]
        homo = bandgap + np.array(lumo)
    output_properties['homo'] = homo     
    lim_homo = np.array(property_limits['homo'])                  
    smiles_features['homo'] = -2.0 + ((4.0)/(lim_homo[1] - lim_homo[0])) * (homo - lim_homo[0]);    
# -----------------------------------------------------------------------------    
# Redox potentials                    
# -----------------------------------------------------------------------------
    for key in ['reduction', 'oxidation']:
        prop = predict(smiles_features, coefficients[key])
        lim = np.array(property_limits[key])
        output_properties[key] = prop*(lim[1] - lim[0]) + lim[0] 

    return output_properties
    
# -----------------------------------------------------------------------------
# Predicting the properties using correlation. 
# -----------------------------------------------------------------------------            
def predict(input_features, coefs):
# -----------------------------------------------------------------------------    
    ''' Function for predicting the value of the correlation. 
    
    Inputs: 
            input_features - fingerprints of the organic molecule
            coefs - correlation coefficients

    Output:
            prop - prediction from the correlation
    '''
# -----------------------------------------------------------------------------        
    bias = np.array(coefs['bias'])
    prop = bias 
    for key in coefs.keys():
        if key != 'bias':
            x = np.array(input_features[key])            
            c = np.array(coefs[key])            
            prop = prop + c[0]*x + c[1]*hyperbolic_tan(x) + c[2]*sigmoid(x)
    return prop
# -----------------------------------------------------------------------------
# Sigmoid function
# -----------------------------------------------------------------------------
def sigmoid(x):
# -----------------------------------------------------------------------------    
    ''' Function for predicting the value of a sigmoid. 
    
    Inputs: 
            x - input value at which sigmoid is evaluated            
    Output:
            sigmoid(x)
    '''
# -----------------------------------------------------------------------------
    return 1.0/(1.0 + np.exp(-x))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Hyperbolic tangent function
# -----------------------------------------------------------------------------
def hyperbolic_tan(x):
# -----------------------------------------------------------------------------    
    ''' Function for predicting the value of a hyperbolic tangent. 
    
    Inputs: 
            x - input value at which hyperbolic tangent is evaluated            
    Output:
            tanh(x)
    '''
# -----------------------------------------------------------------------------
    return np.tanh(x)
# -----------------------------------------------------------------------------
    
if __name__ == '__main__':
    properties = properties('C1=CC=CC=C1', lumo = -0.513207004)
    print(properties)
        
    