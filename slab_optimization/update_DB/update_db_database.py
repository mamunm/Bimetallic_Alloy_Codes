#!/usr/bin/env python
# coding: utf-8

from qescripts.fwio import encode_to_atoms
import ase
import os
from ase.db import connect
from ase.io import Trajectory, read
from ase.visualize import view
from ase import Atoms
from collections import deque
from subprocess import call

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

db = connect('Surface.db')

DB_dir = ('/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Surface/'
          'A/Calculations')
list_A = os.listdir(DB_dir)
list_A = sorted(list_A)
bad_lat = ['Al', 'Cu', 'Co', 'Cr', 'Sc', 'Tc','V', 'Zr', 'Zn', 'Y', 
           'W', 'Ta', 'Pd', 'Nb', 'La', 'Hf', 'Ga', 'Fe', 'Au', 'Ag']
for i in bad_lat:
    if i in list_A:
        list_A.remove(i)

for a in list_A:
    
    print('Processing file:{0}'.format(a))
    in_traj_file = (DB_dir + '/' + a + '/Initial_Slab.traj')
    traj_file = (DB_dir + '/' + a + '/qn.traj')
    atoms0 = read(in_traj_file)
    atoms = read(traj_file)
    keys['metal'] = a
    keys['unique_slab_formula'] = keys['metal']
    keys['SBS_symbol'] = 'A1'
    keys_ini =  keys.copy()
    with open(DB_dir + '/' + a + '/out_scf.log') as f:
        for line in f:
            if 'BFGSLineSearch:' in line:
                keys_ini['PE'] = float(line.split()[-2].strip('*'))
                break
    f = deque(open(DB_dir + '/' + a + '/out_scf.log'), 1)
    for line in f:
        keys['PE'] = float(line.split()[-2].strip('*'))
    print(keys_ini['PE'], keys['PE'])
    db.write(atoms0, key_value_pairs=keys_ini, data=parameters) 
    db.write(atoms, key_value_pairs=keys, data=parameters) 

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Surface/A3B'
list_A3B = os.listdir(DB_dir)
list_A3B = sorted(list_A3B)

for a in list_A3B:
        
    print('Processing file:{0}'.format(a))
    in_traj_file = (DB_dir + '/' + a + '/Initial_Slab.traj')
    traj_file = (DB_dir + '/' + a + '/qn.traj')
    atoms0 = read(in_traj_file)
    atoms = read(traj_file)
    keys['metal'] = a
    keys['unique_slab_formula'] = keys['metal']
    keys['SBS_symbol'] = 'L12'
    keys_ini =  keys.copy()
    with open(DB_dir + '/' + a + '/out_scf.log') as f:
        for line in f:
            if 'BFGSLineSearch:' in line:
                keys_ini['PE'] = float(line.split()[-2].strip('*'))
                break
    f = deque(open(DB_dir + '/' + a + '/out_scf.log'), 1)
    for line in f:
        keys['PE'] = float(line.split()[-2].strip('*')) 
    print(keys_ini['PE'], keys['PE'])
    db.write(atoms0, key_value_pairs=keys_ini, data=parameters) 
    db.write(atoms, key_value_pairs=keys, data=parameters) 

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Surface/AB'
list_AB = os.listdir(DB_dir)
list_AB = sorted(list_AB)

for a in list_AB:
        
    print('Processing file:{0}'.format(a))
    in_traj_file = (DB_dir + '/' + a + '/Initial_Slab.traj')
    traj_file = (DB_dir + '/' + a + '/qn.traj')
    try:
        atoms0 = read(in_traj_file)
    except DeprecationWarning:
        call(["python", "-m", "ase.io.trajectory", in_traj_file])
        atoms0 = read(in_traj_file)
    except IOError:
        try:
            atoms0 = read(DB_dir + '/' + a + 
                    '/Initial_Slab_with_constraint.traj')
        except DeprecationWarning:
            call(["python", "-m", "ase.io.trajectory", DB_dir + '/' + a + 
                    '/Initial_Slab_with_constraint.traj'])
            atoms0 = read(DB_dir + '/' + a +                                    
                                        '/Initial_Slab_with_constraint.traj')

    try:
        atoms = read(traj_file)
    except DeprecationWarning:
        call(["python", "-m", "ase.io.trajectory", traj_file])
        atoms0 = read(traj_file)
    keys['metal'] = a
    m1, m2, _ = a.split('2')
    if m1 < m2:
        AB_name = m1 + '2' + m2 + '2'
    else:
        AB_name = m2 + '2' + m1 + '2'
        
    keys['unique_slab_formula'] = AB_name
    keys['SBS_symbol'] = 'L10'
    keys_ini =  keys.copy()
    with open(DB_dir + '/' + a + '/out_scf.log') as f:
        for line in f:
            if 'BFGSLineSearch:' in line:
                keys_ini['PE'] = float(line.split()[-2].strip('*'))
                break
    f = deque(open(DB_dir + '/' + a + '/out_scf.log'), 1)
    for line in f:
        keys['PE'] = float(line.split()[-2].strip('*')) 
    print(keys_ini['PE'], keys['PE'])
    db.write(atoms0, key_value_pairs=keys_ini, data=parameters) 
    db.write(atoms, key_value_pairs=keys, data=parameters) 



