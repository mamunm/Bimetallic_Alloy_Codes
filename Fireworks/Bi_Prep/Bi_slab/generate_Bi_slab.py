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
from ase.build import fcc111
from ase.data import covalent_radii as cradii
from ase.data import atomic_numbers as an
from ase.constraints import FixAtoms
from ase.build import add_adsorbate

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


metal = ['Bi']
metal = sorted(metal)
data = {m : get_a(m) for m in metal}
print(data)

db = connect('A_slab.db')


for k, v in data.items():

    slab = fcc111(k, (2, 2, 3), float(v), 10)
    c = FixAtoms(indices=[a.index for a in slab if a.tag > 1])
    slab.set_constraint(c)
    
    slab.wrap()
    parameters['spinpol'] = False
            
    keys['metal'] = k

    db.write(slab, key_value_pairs=keys, data=parameters)

 
   





