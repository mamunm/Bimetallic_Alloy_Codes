#!/bin/env python
# coding: utf-8
"""
Created on Oct 18, 2017
"""
__author__ = "Osman Mamun"
__version__ = "0.1"
__maintainer__ = "Osman Mamun"
__email__ = "mamunm@stanford.edu"
__date__ = "Oct 18, 2017"

from fireworks import LaunchPad, Firework, Workflow, PyTask 
import os 
import sys
from shutil import copyfile
import numpy as np
from ase.db import connect
from numpy.linalg import norm
from qescripts.fwio import atoms_to_encode                                      
from qescripts.fwio import encode_to_atoms                                      
from qescripts.fwespresso import get_potential_energy  

# Read credentials from a secure location                                       
host = 'suncatls2.slac.stanford.edu'                                            
#username, name, password = netrc().authenticators(host)                        
username = 'mamun'                                                              
name = 'mamunm_db'                                                              
password = 'postdoc'                                                            
                                                                                
launchpad = LaunchPad(                                                          
    host=host,                                                                  
    name=name,                                                                  
    username=username,                                                          
    password=password) 

db = connect('slab.db')

def get_chemical_symbol(M):
    M = M[-5:-1]
    N1, N2 = set(M)
    n = 0
    for m in M:
        if m == N1:
            n += 1
    if n == 3:
        return N1 + '3' + N2
    else:
        return N2 + '3' + N1

def get_ads_symbol(M):
    if M == 1:
        return 'H'
    elif M == 6:
        return 'C'
    elif M == 7:
        return 'N'
    elif M == 8:
        return 'O'
    else:
        return 'S'

for i, d in enumerate(db.select(['site=hollow'])):
    
    #if i != 0:
    #    continue
    slab  = d.toatoms()
    kvp = d.key_value_pairs
    param = d.data
    slab.info = param
    
    if kvp['SBS_symbol'] == 'L12':
        kvp['symbol'] = get_chemical_symbol(slab.get_chemical_symbols())
        if kvp['symbol'] == 'Ag3Cd':
            # Encode the atoms                                               
            encoding = atoms_to_encode(slab)                      

            # Define some searching keys  
            kvp_sk = {'symbol': kvp['symbol']}
            kvp_sk['adsorbate'] = get_ads_symbol(kvp['adsorbate'])
            search_keys = {'calc': kvp_sk}                                     
                                                                                
            # Two steps - write the input structure to an input file, 
            #then relax      
            t0 = PyTask(                                                    
                        func='qescripts.fwio.encode_to_atoms',              
                        args=[encoding])                                
            t1 = PyTask(                                                    
                        func='qescripts.fwespresso.get_potential_energy',   
                        stored_data_varname='trajectory')                   
                                                                                
            # Package the tasks into a firework, the fireworks into a workflow,
            # and submit the workflow to the launchpad                      
            firework = Firework([t0, t1], spec={'_priority': 1}, 
                                name=search_keys)     
            workflow = Workflow([firework])                         
            launchpad.add_wf(workflow)    

