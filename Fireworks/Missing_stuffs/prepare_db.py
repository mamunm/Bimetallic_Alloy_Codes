#!/usr/bin/env python
# -*-coding: utf-8 -*-

#prepare_db.py
#Osman Mamun
#LAST UPDATED: 02-27-2018

import numpy as np
from ase.build import bulk
from catkit.surface import SlabGenerator
from ase.db import connect
from ase import Atoms
from catkit import utils

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

M_calc = ['Mn3Cr', 'Mn3Ir', 'Fe3Os', 'Bi3Co', 'Ti3Fe', 'Mn3Bi', 'Sc3Mn', 
          'Fe3Tc', 'Co3Bi', 'Zr3Bi', 'Co2Nb2', 'Fe2Mo2', 'Co2Tc2', 'Ga2Rh2', 
          'Co2Mn2', 'Fe2W2', 'Mn2Bi2', 'Co2Bi2', 'Mn2Tc2', 'Mn2Os2', 'Co2Ta2', 
          'Mn2Hf2', 'Mn2Fe2', 'Fe3Cr', 'Fe3Mn']

pol_metal = ['Fe', 'Ni', 'Co', 'Mn']

dbb = connect('Bulk.db')
db = connect('missing.db')

for m in M_calc:
    for d in dbb.select():
        if d['unique_formula'] == m:
            bulk = d.toatoms()
    gen = SlabGenerator(
        bulk,
        miller_index=(1, 1, 1),
        layers=3,
        fixed=2,
        vacuum=10)
    atoms = gen.get_slab()

    keys['metal'] = m
    atoms.wrap()
    parameters['spinpol'] = False
    
    if '3' in m: 
        metals = m.split('3')
    else:
        metals = m.split('2')

    if metals[0] in pol_metal or metals[1] in pol_metal:
        parameters['spinpol'] = True

        atoms.set_initial_magnetic_moments(
        [2.0 if atom.symbol in pol_metal else 0 for atom in atoms]
                    )
    db.write(atoms, key_value_pairs=keys, data=parameters)

