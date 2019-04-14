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
import numpy as np


data = np.load('dict.npy', encoding='latin1')[()]

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
IDS_A = [i for i in IDS if i < 700]
IDS_A = sorted(IDS_A)

if os.path.exists('All_ads_slab.db'):
    os.remove('All_ads_slab.db')

db = connect('All_ads_slab.db')

#Processing A calculations
for ID in IDS_A:
    print('Processing id: {0}'.format(ID))
    launch = launchpad.get_fw_dict_by_id(ID)
    launch2 = launchpad.get_fw_by_id(ID)
    
    atoms0 = encode_to_atoms(launch2.spec['_tasks'][0]['args'][0])[0]
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
    atoms = encode_to_atoms(encoding)
    
    try: 
        keys['metal'] = launch['name']['calc']['metal']
        keys['unique_slab_formula'] = keys['metal']
        keys['adsorbate'] = launch['name']['calc']['adsorbate']
        keys['site'] = launch['name']['calc']['site']
        keys['SBS_symbol'] = 'A1'
    
    except KeyError:
        continue

    db.write(atoms0, key_value_pairs=keys, data=parameters)
    db.write(atoms[-1], key_value_pairs=keys, data=parameters)

IDS_A_New = [i for i in range(30688, 31046)]

for ID in IDS_A_New:
    print('Processing id: {0}'.format(ID))
    launch = launchpad.get_fw_dict_by_id(ID)
    launch2 = launchpad.get_fw_by_id(ID)
    
    atoms0 = encode_to_atoms(launch2.spec['_tasks'][0]['args'][0])[0]
    try:
        encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
        atoms = encode_to_atoms(encoding)
    except (TypeError, KeyError) as e:
        continue
    
    try: 
        a = launch['name']['calc']['description']
        if a == 'New lattice pure metal ads stuffs':
            keys['metal'] = launch['name']['calc']['metal']
            keys['unique_slab_formula'] = keys['metal']
            keys['adsorbate'] = data[ID]['adsorbate']
            keys['site'] = data[ID]['site']
            keys['SBS_symbol'] = 'A1'
    
    except KeyError:
        continue

    db.write(atoms0, key_value_pairs=keys, data=parameters)
    db.write(atoms[-1], key_value_pairs=keys, data=parameters)


IDS_Bi = [i for i in range(24219, 24239)]
for ID in IDS_Bi:
    print('Processing id: {0}'.format(ID))
    launch = launchpad.get_fw_dict_by_id(ID)
    launch2 = launchpad.get_fw_by_id(ID)
    
    atoms0 = encode_to_atoms(launch2.spec['_tasks'][0]['args'][0])[0]
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
    atoms = encode_to_atoms(encoding)
    
    try: 
        keys['metal'] = launch['name']['calc']['metal']
        keys['unique_slab_formula'] = keys['metal']
        keys['adsorbate'] = launch['name']['calc']['adsorbate']
        keys['site'] = launch['name']['calc']['site']
        keys['SBS_symbol'] = 'A1'
    
    except KeyError:
        continue

    db.write(atoms0, key_value_pairs=keys, data=parameters)
    db.write(atoms[-1], key_value_pairs=keys, data=parameters)



#Processing A3B calculations
IDS = sorted(IDS)

for ID in IDS:
    print('Processing id: {0}'.format(ID))
    launch = launchpad.get_fw_dict_by_id(ID)
    launch2 = launchpad.get_fw_by_id(ID)
    
    atoms0 = encode_to_atoms(launch2.spec['_tasks'][0]['args'][0])[0]
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
    atoms = encode_to_atoms(encoding)
    
    try: 
        keys['metal'] = launch['name']['calc']['symbol']
        keys['unique_slab_formula'] = keys['metal']
        keys['adsorbate'] = launch['name']['calc']['adsorbate']
        keys['site'] = 'hollow'
        keys['SBS_symbol'] = 'L12'
    
    except KeyError:
        continue

    db.write(atoms0, key_value_pairs=keys, data=parameters)
    db.write(atoms[-1], key_value_pairs=keys, data=parameters)


##Note:: include Bi slab calculations from fireworks


if os.path.exists('input.traj'):
    os.remove('input.traj')
