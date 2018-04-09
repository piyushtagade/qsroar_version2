# -*- coding: utf-8 -*-
"""
This module defines a class for creating functional groups. 
"""

class functional_group:
     ''' A class for creating functional groups. 
         This class has following four attributes:
             symbol         - SMILES substring used to identify the functional group
             Identifying_fn - User defined function for identifying the number of functional groups in the molecule
             min_val        - Minimum possible number of the functional group present in the molecule
             max_val        - Maximum possible number of the functional group present in the molecule 
     '''
     symbol = None
     Identifying_fn = None

     def __init__(self,symbol,fn, min_val = 0, max_val = 1):
         ''' Method to initialize the class '''
         self.symbol = symbol
         self.Identifying_fn = fn
         self.min_val = min_val 
         self.max_val = max_val
