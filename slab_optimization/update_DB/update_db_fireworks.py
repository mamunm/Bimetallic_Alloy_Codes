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

def get_u_f(A):
    if len(set(A)) == 1:
        return A[0]
    else:
        m1, m2 = set(A)
        n1 = sum([1 for i in A if i == m1])
        if n1 == 6:
            if m1 < m2:
                return m1 + '2' + m2 + '2'
            else:
                return m2 + '2' + m1 + '2'
        else:
            if n1 == 9:
                return m1 + '3' + m2
            else:
                return m2 + '3' + m1

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
IDS = sorted(launchpad.get_fw_ids(query={'state': 'COMPLETED'}))
db = connect('Surface.db')
omit = ['Fe2Sc2', 'Fe2Ti2', 'Fe2V2', 'Mn2Sc2', 'Mn2Ti2', 'Mn2V2', 'Ni2Sc2',
        'Ni2Ti2', 'Ni2V2']
#Processing A calculations
for ID in IDS:
    print('Processing id: {0}'.format(ID))
    launch = launchpad.get_fw_dict_by_id(ID)
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
    atoms = encode_to_atoms(encoding)
    if not len(atoms[0]) == 12:
        continue
    
    try: 
        usf = get_u_f(atoms[0].get_chemical_symbols())
    except KeyError:
        continue
    
    if usf in omit:
        print('Found the culprit!')
        continue
    keys['unique_slab_formula'] = usf
    keys_ini = keys.copy()
    keys_ini['PE'] = atoms[0].get_potential_energy()
    keys['PE'] = atoms[-1].get_potential_energy()
    db.write(atoms[0], key_value_pairs=keys_ini, data=parameters)
    db.write(atoms[-1], key_value_pairs=keys, data=parameters)
