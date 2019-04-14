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

db = connect('All_ads_slab.db')

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/A3B'
list_A3B_hol = os.listdir(DB_dir + '/Hollow')
list_A3B_hol = sorted(list_A3B_hol)

for a in list_A3B_hol:
    s_files = os.listdir(DB_dir + '/Hollow/' + a + '/Calculations')
    for b in s_files:
        
        print('Processing file:{0}'.format(a + '_' + b))
        traj_file = (DB_dir + '/Hollow/' + a + '/Calculations/' + b + '/qn.traj')
        
        try:
            atoms = read(traj_file, ":")
        except ValueError:
            continue
        keys['metal'] = a
        keys['unique_slab_formula'] = keys['metal']
        ads, st, _, _ = b.split('_')
        keys['adsorbate'] = ads
        keys['site'] = st
        keys['SBS_symbol'] = 'L12'

        for ato in atoms:
            db.write(ato, key_value_pairs=keys, data=parameters) 


DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/A3B'
list_A3B_top = os.listdir(DB_dir + '/Top')
list_A3B_top = sorted(list_A3B_top)

for a in list_A3B_top:
    s_files = os.listdir(DB_dir + '/Top/' + a + '/Calculations')
    for b in s_files:
        
        print('Processing file:{0}'.format(a + '_' + b))
        traj_file = (DB_dir + '/Top/' + a + '/Calculations/' + b + '/qn.traj')
        try:
            atoms = read(traj_file, ":")
        except ValueError:
            continue

        keys['metal'] = a
        keys['unique_slab_formula'] = keys['metal']
        ads, st, _, _ = b.split('_')
        keys['adsorbate'] = ads
        keys['site'] = st
        keys['SBS_symbol'] = 'L12'

        for ato in atoms:
            db.write(ato, key_value_pairs=keys, data=parameters) 

DB_dir = '/scratch/users/mamunm/SHERLOCK/Binary_Alloy/Adsorbate-Surface/A3B/'
DB_dir += 'Calculations'
list_A3B_scratch = os.listdir(DB_dir)
list_A3B_scrtch = sorted(list_A3B_scratch)

for a in list_A3B_scratch:
    try:
        s_files = os.listdir(DB_dir + '/' + a + '/Calculations')
    except OSError:
        continue
    
    for b in s_files:
        c, d = a.split('3')
        if os.path.exists(fin_traj_path):
            print('Processing file:{0}'.format(a + '_' + b))
            traj_file = (DB_dir + '/' + a + '/Calculations/' + b + '/qn.traj')
            try:
                atoms = read(traj_file, ":")
            except ValueError:
                continue
            keys['metal'] = a
            ads, st, _, _ = b.split('_')
            keys['unique_slab_formula'] = keys['metal']
            keys['adsorbate'] = ads
            keys['site'] = st
            keys['SBS_symbol'] = 'L12'

            for ato in atoms:
                 db.write(ato, key_value_pairs=keys, data=parameters) 

