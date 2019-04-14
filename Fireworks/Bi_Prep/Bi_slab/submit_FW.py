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

db = connect('A_slab.db')

for i, d in enumerate(db.select()):
    
    #if i != 0:
    #    continue
    slab  = d.toatoms()
    kvp = d.key_value_pairs
    kvp['metal'] = slab.get_chemical_symbols()[0]
    param = d.data
    slab.info = param

    # Encode the atoms                                                          
    encoding = atoms_to_encode(slab)                                           
                                                                                
    # Define some searching keys  
    search_keys = {'calc': kvp}                                           
                                                                                
    # Two steps - write the input structure to an input file, then relax      
    t0 = PyTask(                                                              
        func='qescripts.fwio.encode_to_atoms',                                 
        args=[encoding])                                                     
    t1 = PyTask(                                                              
        func='qescripts.fwespresso.get_potential_energy',                  
        stored_data_varname='trajectory')                              
                                                                                
    # Package the tasks into a firework, the fireworks into a workflow,        
    # and submit the workflow to the launchpad                                  
    firework = Firework([t0, t1], spec={'_priority': 1}, name=search_keys)     
    workflow = Workflow([firework])                                             
    launchpad.add_wf(workflow)    
