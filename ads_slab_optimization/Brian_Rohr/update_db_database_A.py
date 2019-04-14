#!/usr/bin/env python
# coding: utf-8

from fireworks import LaunchPad
from qescripts.fwio import encode_to_atoms
import ase
import os
from ase.db import connect
from ase.io import Trajectory, read
from ase.visualize import view
from ase import Atoms

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
    'crystal': 'fcc'
}

db = connect('Brian_ads_slab.db')

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/A'
list_A = os.listdir(DB_dir)
list_A = sorted(list_A)
print(list_A)
bad_lat = ['Al', 'Cu', 'Co', 'Cr', 'Sc', 'Tc','V', 'Zr', 'Zn', 'Y', 
           'W', 'Ta', 'Pd', 'Nb', 'La', 'Hf', 'Ga', 'Fe', 'Au', 'Ag']
print(bad_lat)
for i in bad_lat:
    if i in list_A:
        list_A.remove(i)

for a in list_A:
    s_files = os.listdir(DB_dir + '/' + a + '/Calculations')
    for b in s_files:
        
        print('Processing file:{0}'.format(a + '_' + b))
        traj_file = (DB_dir + '/' + a + '/Calculations/' + b + '/qn.traj')
        try:
            atoms = read(traj_file, ":")
        except ValueError:
            continue
        keys['metal'] = a
        keys['unique_slab_formula'] = keys['metal']
        ads, st = b.split('_')
        keys['adsorbate'] = ads
        keys['site'] = st
        keys['SBS_symbol'] = 'A1'
        
        for ato in atoms:
            db.write(ato, key_value_pairs=keys, data=parameters) 


