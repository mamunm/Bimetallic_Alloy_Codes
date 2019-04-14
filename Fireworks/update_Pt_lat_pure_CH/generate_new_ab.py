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

adsorbates = ['C-H', 'O-H', 'N-H', 'S-H', 'C-H2', 'C-H3', 'O-H2']
sites = ['ontop', 'bridge', 'fcc', 'hcp']
pol_metal = ['Fe', 'Co', 'Ni', 'Mn']

db_ab_A = connect('new_ab_A.db')


def get_height(K, A, S):
    if S == 'ontop':
        return 0.85 * (cradii[an[A]] + cradii[an[K]])
    elif S == 'bridge':
        return 0.75 * (cradii[an[A]] + cradii[an[K]])
    else:
        return 0.65 * (cradii[an[A]] + cradii[an[K]])


metals = ['Ir', 'Mo', 'Os', 'Pt', 'Re', 'Ti', 'Ir', 'Mo', 'Rh', 'Ru']

a_dict = {}

db_bulk = connect('bulk-mamunm.db')

for d in db_bulk.select():
    if d['unique_formula'] in metals:
        a_dict[d['unique_formula']] = d['lat_a']

print(a_dict)

for i in metals:

    v = a_dict[i]
    slab = fcc111(i, (2, 2, 3), float(v), 10)
    c = FixAtoms(indices=[a.index for a in slab if a.tag > 1])
    slab.set_constraint(c)

    for ads in adsorbates:
        for st in sites:

            ads_slab = slab.copy()
            par, chi = ads.split('-')
            print('Now working on system with the following spec:')
            print('metal: {} ads: {}, site: {}\n'.format(i, ads, st))
            h_z = get_height(i, par, st)
            #add_adsorbate(ads_slab, par, h_z, st)
            #if chi == 'H' and not par == 'S':
            #    add_adsorbate(ads_slab, chi, h_z+1.05, st)
            #elif chi == 'H' and par == 'S':
            #    add_adsorbate(ads_slab, chi, h_z+1.30, st)

            if ads == 'C-H':
                molecule = Atoms('CH',
                                 positions=[(0., 0., 0.), (0., 0., 1.1)])
                add_adsorbate(ads_slab, molecule, h_z, st)

            elif ads == 'O-H':
                molecule = Atoms('OH',
                                 positions=[(0., 0., 0.), (0., 0., 0.95)])
                add_adsorbate(ads_slab, molecule, h_z, st)

            elif ads == 'N-H':
                molecule = Atoms('NH',
                                 positions=[(0., 0., 0.), (0., 0., 1.05)])
                add_adsorbate(ads_slab, molecule, h_z, st)

            elif ads == 'S-H':
                molecule = Atoms('SH',
                                 positions=[(0., 0., 0.), (0., 0., 1.25)])
                add_adsorbate(ads_slab, molecule, h_z, st)

            elif ads == 'C-H2':
                molecule = Atoms('CH2',
                                  positions=[(0., 0., 0.),
                                             (0.63, 0.63, 0.63),
                                             (-0.63, -0.63, 0.63)])
                add_adsorbate(ads_slab, molecule, h_z, st)

            elif ads == 'C-H3':
                molecule = Atoms('CH3',
                                  positions=[(0., 0., 0.),
                                             (1.09, 0., 0.63),
                                             (-0.45, 0.8895, 0.63),
                                             (-0.45, -0.8895, 0.63)])
                add_adsorbate(ads_slab, molecule, h_z, st)

            elif ads == 'C-H4':
                molecule = Atoms('CH4',
                                  positions=[(0., 0., 0.),
                                             (0.63, 0.63, 0.63),
                                             (-0.63, -0.63, 0.63),
                                             (0.63, -0.63, -0.63),
                                             (-0.63, 0.63, -0.63)])
                add_adsorbate(ads_slab, molecule, h_z+0.5, st)

            else:
                molecule = Atoms('OH2',
                                  positions=[(0., 0., 0.),
                                             (0., 0.763, 0.59),
                                             (0., -0.763, 0.59)])
                add_adsorbate(ads_slab, molecule, h_z, st)


            ads_slab.wrap()
            parameters['spinpol'] = False
            if i in pol_metal:
                parameters['spinpol'] = True
                ads_slab.set_initial_magnetic_moments(
                        [2.0 if atom.symbol in pol_metal else 0
                            for atom in ads_slab])

            keys['metal'] = i
            keys['site'] = st
            keys['adsorbate'] = ads

            db_ab_A.write(ads_slab, key_value_pairs=keys, data=parameters)








