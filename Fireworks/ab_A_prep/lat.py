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

import os 
import sys
from shutil import copyfile
import numpy as np
from ase.build import bulk
from ase.db import connect
from ase import Atoms
from ase.build import fcc111
from ase.data import covalent_radii as cradii
from ase.data import atomic_numbers as an
from ase.constraints import FixAtoms
from ase.build import add_adsorbate


db_A = connect('pure_A_surface.db')
db_B = connect('~/DATABASE/Bulk.db')

def get_a(M):                                                                   
    DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Bulk/A/'      
    RC = False                                                                  
                                                                                
    with open(DB_dir + 'EOS_DATA/eos_results.txt') as f:                        
                                                                                
        for line in f:                                                          
            search_st = 'New lattice constant found for {0}:'.format(M)         
            if search_st in line:                                               
                return line.split()[-2]                                         
            else:                                                               
                RC = True                                                       
                                                                                
    if RC:                                                                      
        with open(DB_dir + M + '/bulk.out') as f:                               
            for line in f:                                                      
                if 'Lattice constant:' in line:                                 
                    return line.split()[-2]          

data_lat = {}
for i in db_A.select():
    
    sl = i.toatoms()
    k = sl.get_chemical_symbols()[0]
    
    v_right = get_a(k)
    data_lat[k] = v_right
    for j in db_B.select():
        if j['unique_formula'] == k:
            v_wrong = j['lat_a']
    print('metal: {0}\n lattice_constant_right: {1}\n lattice_constant_wrong: {2}\n difference_percentage: {3}'.format(k, v_right,
                v_wrong, (100 * (float(v_right)-float(v_wrong))/float(v_right))))


print(data_lat)
np.save('data_lat.npy', data_lat)    






