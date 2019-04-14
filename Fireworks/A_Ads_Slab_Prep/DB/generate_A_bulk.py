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


metal = ['Cd', 'Co', 'Cu', 'Fe', 'Ga', 'La', 'Mn', 'Ni', 'Sn', 'Tl', 'Zn']
metal = sorted(metal)
data = {m : get_a(m) for m in metal}

db = connect('A_Bulk.db')

images = []
for k, v in data.items():
    atoms = bulk(k, 'fcc', a=v, cubic=True)
    keys = {'lattice_constant': float(v)}

    db.write(atoms, key_value_pairs=keys)






