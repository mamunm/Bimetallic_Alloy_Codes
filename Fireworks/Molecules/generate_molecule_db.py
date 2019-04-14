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
from ase.db import connect
from six.moves import urllib
from six.moves.urllib.error import HTTPError
from ase.io import read, write

# Change any of the desired keywords here.
parameters = {
    'mode': 'relax',
    'opt_algorithm': 'bfgs',
    'xc': 'BEEF-vdW',
    'kpts': (1, 1, 1),
    'pw': 500,
    'dw': 5000,
    'sigma': 0.1,
    'nbands': -20,
    'fmax': 0.05,
    'beefensemble': True,
    'spinpol': True,
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
            'energy': 1e-6 * 13.6,
            'mixing': 0.2,
            'nmix': 10,
            'mix': 4,
            'maxsteps': 2000,
            'diag': 'david'
    },
}


molecules = ['acetylene', 'ethane', 'methane', 'CO2', 'H2', 'H2O2', 
             'formic acid', 'ethanoic acid', 'N2O', 'NO', 'O2', 'ethene', 
             'acetone', 'ethanol', 'ethanal', 'CO', 'CS2', 'H2O', 
             'hydrogen cyanide', 'N2', 'NH3', 'NO2', 'SO2', 'methanol', 
             'H2S', 'SO3', 'SO4', 'benzene', 'cyclohexane', 'cyclohexene', 
             'methanethiol', 'ethanethiol', '1-propanethiol', '2-propanethiol', 
             'butanethiol', 'pentanethiol', '2-propenethiol', 'thiophenol', 
             'propane', 'propanol', 'propene', 'propyne', 'butane', '1-butene', 
             'cis-2-butene', 'trans-2-butene', 'isobutene', '1,3-butadiene', 
             '1-butyne', '2-butyne', 
             'cyclobutadiene', 'cyclobutene', 'cyclobutane', 'cyclobutyne', 
             'squaric acid', '1-butanol', '2-butanol', 'pentane', '1-pentene', 
             'cis-2-pentene', 'trans-2-pentene', '1-pentyne', '2-pentyne', 
             'benzoic acid', 'toluene']
molecules = sorted(molecules)

db = connect('molecules.db')

def fetch_molecule(mol_name):

    print('Processing: {0}'.format(mol_name))

    if ' ' in mol_name:
        mol_name = mol_name.split()[0] + '%20' + mol_name.split()[1]

    link = 'https://cactus.nci.nih.gov/chemical/structure/' + mol_name + '/sdf'
    
    try:
        f = urllib.request.urlopen(link)
        data = f.read()
        data = data.decode('utf-8')
        
        if '%20' in mol_name:
            mol_name = mol_name.replace('%20', '_')
        
        '''
        with open(mol_name + '.sdf', 'w') as ff:
            ff.write(data)
        '''

        data = data.split('\n')
    
    except HTTPError:
        print('The molecule you are looking for does not' 
              ' exist in the database') 
        return None
    
    else:
        try:
            int(data[3].split()[0])
        except ValueError:
            xyz_data = data
        else:
            n_atom = data[3].split()[0]
            xyz_data = [n_atom]
            xyz_data.append(data[0])
        
            for i, dt in enumerate(data):
                if 3 < i < 4 + int(n_atom):
                
                    words = dt.split()
                    temp = [words[3], float(words[0]), 
                            float(words[1]), float(words[2])]
                    xyz_data.append(temp)
        
        if not os.path.exists('mol_structure'):
            os.makedirs('mol_structure')
        
        with open('mol_structure/' + mol_name + '.xyz', 'w') as f:
            for dt in xyz_data:
                if isinstance(dt, list):
                    
                    f.write('%2s %12.4f %12.4f %12.4f\n'%(dt[0], 
                        dt[1], dt[2], dt[3]))
                
                else:
                    
                    f.write('%s\n'%(dt))
        
        print('Processing of {0} done successfully. Congratualltions!'.format(
               mol_name))

        return read('mol_structure/' + mol_name + '.xyz')

for mol in molecules:

    atoms = fetch_molecule(mol)

    atoms.set_cell([10, 10, 10])

    atoms.set_pbc(True)

    atoms.center()
    
    atoms.set_initial_magnetic_moments([2.0 for atom in atoms])
            
    keys = {'molecule': mol}

    db.write(atoms, key_value_pairs=keys, data=parameters)

 
   





