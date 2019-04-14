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
from qescripts.fwespresso import get_neb
from ase.visualize import view

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

db = connect('NEB-pairings-fixedbase.db')

slab = []

for i, d in enumerate(db.select()):
    
    if i > 4:
        break
    slab += [d.toatoms()]
    kvp = d.key_value_pairs
    kvp['SBS'] = 'A1'
    param = d.data
    param['properties']['kpts'] = param['properties']['kpts'].tolist()
    del param['connectivity']
    del param['symmetry_index']
    del param['site_positions']
    param = param['properties']
    param['calcstress'] = False
    slab[i].info = param
# Encode the atoms                                                          
encoding = atoms_to_encode(slab)                                           
                                                                                
# Define some searching keys  
search_keys = {'calc': kvp}                             
search_keys['des'] = 'NEB_TEST'
                                                                                
# Two steps - write the input structure to an input file, then relax      
t0 = PyTask(                                                              
        func='qescripts.fwio.encode_to_atoms',                                 
        args=[encoding])                                                     
t1 = PyTask(                                                              
        func='qescripts.fwespresso.get_neb',                  
        stored_data_varname='trajectory')                              
                                                                                
# Package the tasks into a firework, the fireworks into a workflow,        
# and submit the workflow to the launchpad                                  
firework = Firework([t0, t1], spec={'_priority': 1000000}, 
                        name=search_keys)     
workflow = Workflow([firework])                                             
launchpad.add_wf(workflow)   

