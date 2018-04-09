"""
This module provide functions for fingerprinting of an organic molecules. 
"""
# -----------------------------------------------------------------------------
import numpy as np
import re
from functional_group import functional_group as fg
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Defining functions for fingerprinting 
# -----------------------------------------------------------------------------
def default_function(symbol,smile):
    ''' Default function for molecular fingerprinting. 
    Input:
        symbol : substring used to identify the functional group.
        smile  : SMILES representation of the molecule
    Output:
        count  : Number of the functional group present in the SMILES representation
    '''     
    count = smile.count(symbol)
    if symbol[-1] == 'C' and count > 0:
        return count - smile.count('Cl')
    else:
        return count

def number_of_rings(symbol,smile):
    ''' Function to identify number of rings in the molecule. 
        Input:
        symbol : substring used to identify the functional group.
        smile  : SMILES representation of the molecule
    Output:
        Number of rings in the molecule''' 
    return int(len(re.findall('\d+', smile))/2) 

def number_of_ring_type(symbol,smile):
    ''' Function to identify number of particular type of rings in the molecule. 
        Input:
        symbol : substring used to identify the functional group.
        smile  : SMILES representation of the molecule
    Output:
        A dictionary containing a count of number of particular type of rings in the molecule.
    '''         
    num = int(len(re.findall('\d+', smile))/2)
    if num == 0:
       return 0
    nb = 0; npz = 0; nep = 0; npentn = 0; npentnp = 0; 
    npyrid = 0; npyrid2 = 0; ring_length = []   
    for ir in range(1, num+1):
         st = smile.partition(str(ir))[-1].rpartition(str(ir))[0]            
         ring_length.append(len(st))
         if st.count('C') == 5 and st.count('=') == 3:
              nb = nb + 1
              if st.count('C') == 4 and st.count('=') == 0:
                  npz = npz + 1
              if st.count('C') == 1 and st.count('O') == 1:
                  nep = nep + 1
              if st.count('C') == 4 and st.count('N') == 1 and st.count('=') == 3:
                  npyrid = npyrid + 1
              if st.count('N') == 2:
                  npyrid2 = npyrid2 + 1                

    return {
       'bz'           : nb,
       'pz'           : npz,
       'ep'           : nep,
       'pentn'        : npentn,
       'pentnp'       : npentnp,
       'pyrid'        : npyrid,
       'pyrid2'       : npyrid2, 
       'maxringlength': np.max(np.array(ring_length))        
	 }.get(symbol,0)
	   

def chain_endo(symbol,smile):
    ''' Function to identify number of chains with terminal oxygen. 
        Input:
        symbol : substring used to identify the functional group.
        smile  : SMILES representation of the molecule
    Output:
        Number of chains with terminal oxygen'''     
    return int(smile[-1]=='O') + smile.count(symbol)

def chain_endc(symbol,smile):
    ''' Function to identify number of chains with terminal carbon. 
        Input:
        symbol : substring used to identify the functional group.
        smile  : SMILES representation of the molecule
    Output:
        Number of chains with terminal carbon'''    
    return int(smile[-1]=='C') + smile.count(symbol)

def chain_endn(symbol,smile):
    ''' Function to identify number of chains with terminal nitrogen. 
        Input:
        symbol : substring used to identify the functional group.
        smile  : SMILES representation of the molecule
    Output:
        Number of chains with terminal nitrogen'''        
    return int(smile[-1]=='N') + smile.count(symbol)
 
def number_of_haloformyl(symbol,smile):
    ''' Function to identify number of haloformyl groups. 
        Input:
        symbol : substring used to identify the functional group.
        smile  : SMILES representation of the molecule
    Output:
        num    : Number of Haloformyl groups in the molecule. 
    '''            
    num = smile.count('COF') + smile.count('COCl') + smile.count('COBr') + smile.count('COI') + smile.count('COAt')
    num = num + smile.count('C(=O)F') + smile.count('C(=O)Cl') + smile.count('C(=O)Br') + smile.count('C(=O)I') + smile.count('C(=O)At')
    return num

# -----------------------------------------------------------------------------
def create_functional_groups(): 
    ''' This function creates user-defined functional groups 
    Output:
        FG : A dictionary functional groups 
    ''' 
    
    FG = {}
    FG['carbon']        =   fg( 'C'            ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 46) #lumo
    FG['oxygen']        =   fg( 'O'            ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 21) #lumo
    FG['nitro']         =   fg( 'N'            ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 15) #lumo
    FG['dbond']         =   fg( '='            ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 24) #lumo
    FG['tbond']         =   fg( '#'            ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 4.0)
    FG['metals']        =   fg( '['            ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 15.0)
    FG['branch']        =   fg( '('            ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 23.0)
    FG['rings']         =   fg( '\d+'          ,  lambda x,y:number_of_rings(x,y))
    FG['cdb']           =   fg( 'C='           ,  lambda x,y:default_function(x,y))
    FG['cdbcc']         =   fg( 'C=CC'         ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 6.0)
    FG['cdbccdbcc']     =   fg( 'C=CC=CC'      ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 2.0) #lumo
    FG['cdbccdbccdbcc'] =   fg( 'C=CC=CC=CC'   ,  lambda x,y:default_function(x,y))
    FG['co']            =   fg( 'CO'           ,  lambda x,y:default_function(x,y))
    FG['cn']            =   fg( 'CN'           ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 8.0)
    FG['cno']           =   fg( 'CNO'          ,  lambda x,y:default_function(x,y))
    FG['cc']            =   fg( 'CC'           ,  lambda x,y:default_function(x,y))
    FG['coo']           =   fg( 'COO'          ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 5.0)
    FG['ccc']           =   fg( 'CCC'          ,  lambda x,y:default_function(x,y))
    FG['cccc']          =   fg( 'CCCC'         ,  lambda x,y:default_function(x,y))
    FG['dbo']           =   fg( '=O'           ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 9) #lumo
    FG['ndbo']          =   fg( 'N=O'          ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 3) #lumo
    FG['sulphur']       =   fg( 'S'            ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 9)
    FG['fluorine']      =   fg( 'F'            ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 23) #lumo
    FG['chlorine']      =   fg( 'Cl'           ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 12) #lumo
    FG['alkene']        =   fg( 'C=C'          ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 14) #lumo
    FG['alkyne']        =   fg( 'C#C'          ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 4.0)
    FG['ether']         =   fg( 'COC'          ,  lambda x,y:default_function(x,y))
    FG['alde']          =   fg( 'CC=O'         ,  lambda x,y:default_function(x,y))
    FG['ket']           =   fg( 'C(=O)'        ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 8) #lumo
    FG['carbox']        =   fg( 'C(=O)O'       ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 6) #lumo
    FG['anhy']          =   fg( 'CC(=O)OC(=O)C',  lambda x,y:default_function(x,y))
    FG['ester']         =   fg( 'C(=O)O'       ,  lambda x,y:default_function(x,y))
    FG['amide']         =   fg( 'C(=O)N'       ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 6) #lumo
    FG['nitrile']       =   fg( 'C#N'          ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 4.0)
    FG['imine']         =   fg( 'CC(=NC)C'     ,  lambda x,y:default_function(x,y))
    FG['isocyanate']    =   fg( 'N=C=O'        ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 3.0) #lumo
    FG['azo']           =   fg( 'N=N'          ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 3.0) #lumo
    FG['thiol']         =   fg( 'CS'           ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 8.0) #lumo
    FG['achalide_f']    =   fg( 'CC(=O)F'      ,  lambda x,y:default_function(x,y))
    FG['achalide_cl']   =   fg( 'CC(=O)Cl'     ,  lambda x,y:default_function(x,y))
    FG['achalide_br']   =   fg( 'CC(=O)Br'     ,  lambda x,y:default_function(x,y))
    FG['achalide_i']    =   fg( 'CC(=O)I'      ,  lambda x,y:default_function(x,y))
    FG['achalide_at']   =   fg( 'CC(=O)At'     ,  lambda x,y:default_function(x,y))
    FG['plus']          =   fg( '+'            ,  lambda x,y:default_function(x,y))
    FG['minus']         =   fg( '-'            ,  lambda x,y:default_function(x,y))
    FG['br']            =   fg( 'Br'           ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 15) #lumo
    FG['p']             =   fg( 'P'            ,  lambda x,y:default_function(x,y))
    #FG['diff_pm']       =   fg( '_'            ,  lambda x,y:diff_fn(x,y))
    FG['at']            =   fg( '@'            ,  lambda x,y:default_function(x,y))
    FG['np']            =   fg( 'N+'           ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 6) #lumo
    FG['nm']            =   fg( 'N-'           ,  lambda x,y:default_function(x,y))
    FG['om']            =   fg( 'O-'           ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 6) #lumo
    FG['branch_dbo']    =   fg( '(=O)'         ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 8) #lumo
    FG['branch_dbc']    =   fg( '(=C'          ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 8) #lumo
    FG['carb_ester']    =   fg( 'OC(=O)OC'     ,  lambda x,y:default_function(x,y))
    FG['sec_amine']     =   fg( 'CNC'          ,  lambda x,y:default_function(x,y))
    FG['tert_amine']    =   fg( 'CN(C'         ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 5) #lumo
    FG['pri_ketamine']  =   fg( 'C(=N'         ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 6) #lumo
    FG['sec_ketamine']  =   fg( 'C(=NC'        ,  lambda x,y:default_function(x,y))
    FG['pri_aldimine']  =   fg( 'C(=N'         ,  lambda x,y:default_function(x,y))
    FG['sec_aldimine']  =   fg( 'C(=NC'        ,  lambda x,y:default_function(x,y))
    FG['imide']         =   fg( 'C(=O)NC(=O)'  ,  lambda x,y:default_function(x,y))
    FG['azide']         =   fg( 'N=[N+]=[N-]'  ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 3.0)
    FG['cyanate']       =   fg( 'OC#N'         ,  lambda x,y:default_function(x,y))
    FG['nitrate']       =   fg( 'O[N+](=O)[O-]',  lambda x,y:default_function(x,y))
    FG['nitrite']       =   fg( 'ON=O'         ,  lambda x,y:default_function(x,y))
    FG['disulphide']    =   fg( 'SS'           ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 4) #lumo
    FG['sulfinyl']      =   fg( 'SO'           ,  lambda x,y:default_function(x,y))
    FG['sulfo']         =   fg( 'S(=O)(=O)O'   ,  lambda x,y:default_function(x,y))
    FG['sulfonyl']      =   fg( 'S(=O)(=O)'    ,  lambda x,y:default_function(x,y))
    FG['thiocynate']    =   fg( 'SC#N'         ,  lambda x,y:default_function(x,y))
    FG['isothiocynate'] =   fg( 'N=C=S'        ,  lambda x,y:default_function(x,y))
    FG['phosphono']     =   fg( 'P(=O)(O)'     ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 5.0)
    FG['phosphate']     =   fg( 'OP(=O)(O)'    ,  lambda x,y:default_function(x,y))
    FG['thial']         =   fg( 'C(=S'         ,  lambda x,y:default_function(x,y), min_val = 0, max_val = 2.0)
    FG['endo']          =   fg( 'O)'           ,  lambda x,y:chain_endo(x,y), min_val = 0, max_val = 15) #lumo
    FG['endc']          =   fg( 'C)'           ,  lambda x,y:chain_endc(x,y), min_val = 0, max_val = 20)
    FG['endn']          =   fg( 'N)'           ,  lambda x,y:chain_endn(x,y), min_val = 0, max_val = 6) #lumo
    FG['haloformyl']    =   fg( 'X'            ,  lambda x,y:number_of_haloformyl(x,y))
    FG['carbonyl']      =   fg( 'C(=O)'        ,  lambda x,y:default_function(x,y))
    FG['isonitrile']    =   fg( 'N#C'          ,  lambda x,y:default_function(x,y))
    FG['benz']          =   fg( 'bz'           ,  lambda x,y:number_of_ring_type(x,y), min_val = 0, max_val = 7) #lumo
    FG['penz']          =   fg( 'pz'           ,  lambda x,y:number_of_ring_type(x,y))
    FG['epox']          =   fg( 'ep'           ,  lambda x,y:number_of_ring_type(x,y))
    FG['pentn']         =   fg( 'pentn'        ,  lambda x,y:number_of_ring_type(x,y), min_val = 0, max_val = 4.0)
    FG['pentnp']        =   fg( 'pentnp'       ,  lambda x,y:number_of_ring_type(x,y))
    FG['pyrid']         =   fg( 'pyrid'        ,  lambda x,y:number_of_ring_type(x,y), min_val = 0, max_val = 4) #lumo
    FG['pyrid2']        =   fg( 'pyrid2'       ,  lambda x,y:number_of_ring_type(x,y), min_val = 0, max_val = 4) #lumo
    FG['maxringlength'] =   fg( 'maxringlength',  lambda x,y:number_of_ring_type(x,y), min_val = 0, max_val = 91) #lumo

    return FG 
# -----------------------------------------------------------------------------
def extract_features(molecule, functional_groups):     
    ''' Function for molecular fingerprinting. 
    Inputs:
        molecule : SMILES representation of molecule
        functional_groups: Functional group object
    Output:
        features: Dictionary with value of each feature
    ''' 
    features = {} 

    for keys in functional_groups:
        count = functional_groups[keys].Identifying_fn(functional_groups[keys].symbol,molecule)
        min_val = functional_groups[keys].min_val
        max_val = functional_groups[keys].max_val
                                   
        count = -2.0 + ((4.0)/(max_val - min_val)) * (count - min_val);  
        features[keys] = count
         
    return features     
# -----------------------------------------------------------------------------         