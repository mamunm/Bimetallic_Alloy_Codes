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
}

# Read credentials from a secure location
host = 'suncatls2.slac.stanford.edu'
username, name, password = 'mamun', 'mamunm_db', 'postdoc'

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)

# Select the ID of the first completed calcualtion
IDS = launchpad.get_fw_ids(query={'state': 'COMPLETED'}) 
IDS = [i for i in IDS if i < 1000 or i > 24200]
IDS = sorted(IDS)
print(IDS)

if os.path.exists('A_ads_slab.db'):
    os.remove('A_ads_slab.db')

db = connect('A_ads_slab.db')

metal = ['Bi', 'Cd', 'Co', 'Cu', 'Fe', 'Ga', 'La', 'Mn', 'Ni', 'Sn', 'Tl', 'Zn']

import json
for ID in IDS:
    print('Processing id: {0}'.format(ID))
    launch = launchpad.get_fw_dict_by_id(ID)
    
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
    atoms = encode_to_atoms(encoding)
    
    try: 
        keys['metal'] = launch['name']['calc']['metal']
        keys['adsorbate'] = launch['name']['calc']['adsorbate']
        keys['site'] = launch['name']['calc']['site']
    
    except KeyError:
        continue
        

    if keys['metal'] in metal:
        keys['unique_name'] = keys['metal'] + keys['adsorbate'] + keys['site'] 
        if db.count('unique_name={}'.format(ID)) > 0:
            print('ALLREADY IN DB???')
            continue

        db.write(atoms[-1], key_value_pairs=keys, data=parameters)

if os.path.exists('input.traj'):
    os.remove('input.traj')

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/A'
list_A = os.listdir(DB_dir)
list_A = sorted(list_A)

for a in list_A:
    s_files = os.listdir(DB_dir + '/' + a + '/Calculations')
    for b in s_files:
        
        print('Processing file:{0}'.format(a + '_' + b))
        traj_file = (DB_dir + '/' + a + '/Calculations/' + b + '/qn.traj')
        atoms = read(traj_file)
        #atoms = atoms[-1]
        keys['metal'] = a
        ads, st = b.split('_')
        keys['adsorbate'] = ads
        keys['site'] = st

        db.write(atoms, key_value_pairs=keys, data=parameters) 




