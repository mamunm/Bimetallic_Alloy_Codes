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

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/AB'
list_AB_top = os.listdir(DB_dir)
list_AB_top = sorted(list_AB_top)

for a in list_AB_top:
    s_files = os.listdir(DB_dir + '/' + a + '/Calculations')
    for b in s_files:
        
        print('Processing file:{0}'.format(a + '_' + b))
        traj_file = (DB_dir + '/' + a + '/Calculations/' + b + '/qn.traj')
        try:
            atoms = read(traj_file, ":")
        except ValueError:
            continue
        keys['metal'] = a
        m1, m2, _ = a.split('2')
        if m1 < m2:
            AB_name = m1 + '2' + m2 + '2'
        else:
            AB_name = m2 + '2' + m1 + '2'
        
        keys['unique_slab_formula'] = AB_name
        ads, st, _, _ = b.split('_')
        keys['adsorbate'] = ads
        keys['site'] = st
        keys['SBS_symbol'] = 'L10'

        for ato in atoms:
            db.write(ato, key_value_pairs=keys, data=parameters) 

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/'
DB_dir += 'AB_from_SLAC_SUNCAT/Calculations'
list_AB_scratch = os.listdir(DB_dir)
list_AB_scrtch = sorted(list_AB_scratch)

for a in list_AB_scratch:
    try:
        s_files = os.listdir(DB_dir + '/' + a + '/Calculations')
    except OSError:
        continue
    
    for b in s_files:
        c, d, _ = a.split('2')
        fin_traj_path = (DB_dir + '/' + a + '/Calculations/' + b + '/final_'
                + b + '_' + c + '6' + d + '6.traj')
        if os.path.exists(fin_traj_path):
            print('Processing file:{0}'.format(a + '_' + b))
            in_traj_file = (DB_dir + '/' + a + '/Calculations/' + b + 
                    '/Initial_slab_with_adsorbate.traj')
            traj_file = (DB_dir + '/' + a + '/Calculations/' + b + '/qn.traj')
            try:
                atoms = read(traj_file, ":")
            except ValueError:
                continue
            keys['metal'] = a
            m1, m2, _ = a.split('2')
            
            if m1 < m2:
                AB_name = m1 + '2' + m2 + '2'
            else:
                AB_name = m2 + '2' + m1 + '2'
        
            keys['unique_slab_formula'] = AB_name
            ads, st, _, _ = b.split('_')
            keys['adsorbate'] = ads
            keys['site'] = st
            keys['SBS_symbol'] = 'L10'

            for ato in atoms:
                db.write(ato, key_value_pairs=keys, data=parameters) 


