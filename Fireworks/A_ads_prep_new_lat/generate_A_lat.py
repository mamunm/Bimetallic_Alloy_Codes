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


def get_height(K, A, S):
    if S == 'ontop':
        return 0.85 * (cradii[an[A]] + cradii[an[K]])
    elif S == 'bridge':
        return 0.75 * (cradii[an[A]] + cradii[an[K]])
    else:
        return 0.65 * (cradii[an[A]] + cradii[an[K]])

# Change any of the desired keywords here.
parameters = {
    'mode': 'relax',
    'opt_algorithm': 'bfgs',
    'xc': 'BEEF-vdW',
    'kpts': (6, 6, 1),
    'pw': 500,
    'dw': 5000,
    'sigma': 0.15,
    'nbands': -10,
    'fmax': 0.05,
    'beefensemble': True,
    'nosym': True,
    'outdir': '.',
    'calcstress': True,
    'dipole': {'status': True},
    'output': {
        'removesave': True,
        'avoidio': False,
        'removewf': True,
        'wf_collect': False,
    },
    'convergence': {
            'energy': 1e-5 * 13.6,
            'mixing': 0.3,
            'nmix': 10,
            'mix': 4,
            'maxsteps': 2000,
            'diag': 'david'
    },
}

keys = {
    'facet': '(1, 1, 1)',
    'layer': 3,
    'vacuum': 10,
}

sites = ['ontop', 'bridge', 'fcc', 'hcp']
adsorbates = ['C', 'H', 'N', 'O', 'S']
pol_metal = ['Fe', 'Co', 'Ni', 'Mn']
surf = ['Al', 'Co', 'Cu', 'Cr', 'Sc', 'Tc', 'V', 'Zr', 
        'Zn', 'Y', 'W', 'Ta', 'Pd', 
        'Nb', 'La', 'Hf', 'Ga', 'Fe', 'Au', 'Ag']

db_A = connect('A_ads_new.db')



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

for i in surf:

    v = get_a(i)
    slab = fcc111(i, (2, 2, 3), float(v), 10) 
    c = FixAtoms(indices=[a.index for a in slab if a.tag > 1])
    slab.set_constraint(c)
    slab.wrap()
    parameters['spinpol'] = False
            
            
    if i in pol_metal:
        parameters['spinpol'] = True
        slab.set_initial_magnetic_moments(
                        [2.0 if atom.symbol in pol_metal else 0 
                            for atom in slab])
            
    keys['metal'] = i
    keys['description'] = 'New lattice pure metal ads stuffs'
    
    for j in adsorbates:
        for k in sites:
            ads_slab = slab.copy()
            add_adsorbate(ads_slab, j, get_height(i, j, k), k)

            db_A.write(ads_slab, key_value_pairs=keys, data=parameters)

 
   





