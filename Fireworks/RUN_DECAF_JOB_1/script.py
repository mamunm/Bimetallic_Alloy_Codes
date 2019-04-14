#!/usr/bin/env python
# -*-coding: utf-8 -*-

#script.py
#Osman Mamun
#DATE CREATED: 01-07-2019

from fireworks import LaunchPad
import ase
from ase.io import Trajectory, read
from ase.visualize import view
from ase import Atoms
import json
import os
import numpy as np
import datetime
import sys
from catkit.flow.fwio import atoms_to_encode, encode_to_atoms

host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
        host=host,
        name=name,
        username=username,
        password=password)

IDS = np.load('fwid.npy')[()]
IDS = IDS[:-1]
pol_metal = ['Fe', 'Ni', 'Mn', 'Co']

for i in IDS:
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_dict_by_id(i)
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
    atoms = encode_to_atoms(encoding)[-1]

    atoms.info['calculator_parameters'] =  {
               'calculator_name': 'decaf.Espresso',
               'calculation'  : 'relax',
               'input_dft': 'BEEF-vdW',
               'ecutwfc': 500,
               'ecutrho': 5000,
               'kpts': [6, 6, 1],
               'nbnd': -10,
               'mixing_beta': 0.3,
               'dipfield': True,
               'degauss': 0.15,
               'beefensemble': True,
               'spinpol': False,
               'mixing_beta': 0.2,
               'electron_maxstep': 2000}


    formula = set(atoms.get_chemical_symbols())
    if np.array([a in pol_metal for a in formula]).any():
        atoms.set_initial_magnetic_moments([3 if atom.symbol == 'Fe'
                                       else 2 if atom.symbol == 'Co'
                                       else 1 if atom.symbol == 'Ni'
                                       else 3 if atom.symbol == 'Mn'
                                       else 0 for atom in atoms])
        atoms.info['calculator_parameters'].update({'spinpol': True,
                                                    'mixing_beta': 0.1,
                                                   'electron_maxstep': 4000})
    print(atoms.info)
    atoms._calc = None

    encoding = atoms_to_encode(atoms)
    print('Submitting job ....')
    launchpad.rerun_fw(i)
    launchpad.update_spec([i],
                spec_document={'_tasks.0.args.0': encoding,
                '_tasks.0.func': 'catkit.flow.fwio.encode_to_atoms',
                '_tasks.1.args': ['decaf.Espresso'],
                '_tasks.1.func':
                'catkit.flow.fwase.get_potential_energy'})

if os.path.exists('input.traj'):
    os.remove('input.traj')
