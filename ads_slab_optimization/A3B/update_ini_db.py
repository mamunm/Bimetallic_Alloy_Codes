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
#IDS = [i for i in IDS if i < 1000]
IDS = sorted(IDS)

if os.path.exists('A3B_ini_ads_slab.db'):
    os.remove('A3B_ini_ads_slab.db')

db = connect('A3B_ini_ads_slab.db')

#Processing A3B calculations
for ID in IDS:
    print('Processing id: {0}'.format(ID))
    launch = launchpad.get_fw_dict_by_id(ID)
    
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
    atoms = encode_to_atoms(encoding)
    
    try: 
        keys['symbol'] = launch['name']['calc']['symbol']
        keys['adsorbate'] = launch['name']['calc']['adsorbate']
        keys['site'] = 'hollow'
        keys['SBS_symbol'] = 'L12'
    
    except KeyError:
        continue
        
    db.write(atoms[0], key_value_pairs=keys, data=parameters)

if os.path.exists('input.traj'):
    os.remove('input.traj')

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/A3B'
list_A3B_hol = os.listdir(DB_dir + '/Hollow')
list_A3B_hol = sorted(list_A3B_hol)

for a in list_A3B_hol:
    s_files = os.listdir(DB_dir + '/Hollow/' + a + '/Calculations')
    for b in s_files:
        
        print('Processing file:{0}'.format(a + '_' + b))
        traj_file = (DB_dir + '/Hollow/' + a + '/Calculations/' + b + 
                '/Initial_slab_with_adsorbate.traj')
        atoms = read(traj_file)
        #atoms = atoms[-1]
        keys['symbol'] = a
        ads, st, _, _ = b.split('_')
        keys['adsorbate'] = ads
        keys['site'] = st
        keys['SBS_Symbol'] = 'L12'

        db.write(atoms, key_value_pairs=keys, data=parameters) 


DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/A3B'
list_A3B_top = os.listdir(DB_dir + '/Top')
list_A3B_top = sorted(list_A3B_top)

for a in list_A3B_top:
    s_files = os.listdir(DB_dir + '/Top/' + a + '/Calculations')
    for b in s_files:
        
        print('Processing file:{0}'.format(a + '_' + b))
        traj_file = (DB_dir + '/Top/' + a + '/Calculations/' + b + 
                '/Initial_slab_with_adsorbate.traj')
        atoms = read(traj_file)
        #atoms = atoms[-1]
        keys['symbol'] = a
        ads, st, _, _ = b.split('_')
        keys['adsorbate'] = ads
        keys['site'] = st
        keys['SBS_Symbol'] = 'L12'

        db.write(atoms, key_value_pairs=keys, data=parameters) 



#Bi calculations                                                 
for ID in IDS:                                                                  
    print('Processing id: {0}'.format(ID))                                      
    launch = launchpad.get_fw_dict_by_id(ID)                                    
                                                                                
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']    
    atoms = encode_to_atoms(encoding)                                           
                                                                                
    try:                                                                        
        keys['symbol'] = launch['name']['calc']['metal']                       
        keys['adsorbate'] = launch['name']['calc']['adsorbate']                 
        keys['SBS_symbol'] = launch['name']['calc']['SBS_symbol']              
        keys['site'] = launch['name']['calc']['site']                 
                                                                                
    except KeyError:                                                            
        continue                                                                
                                                                                
                                                                                
    keys['unique_name'] = keys['symbol'] + keys['adsorbate'] + keys['site']     
    if db.count('unique_name={}'.format(ID)) > 0:                               
        print('ALLREADY IN DB???')                                              
        continue                                                                
    if keys['SBS_symbol'] == 'L12' and len(atoms) == 13:                    
        db.write(atoms[0], key_value_pairs=keys, data=parameters)   



DB_dir = '/scratch/users/mamunm/SHERLOCK/Binary_Alloy/Adsorbate-Surface/A3B/'
DB_dir += 'Calculations'
list_A3B_scratch = os.listdir(DB_dir)
list_A3B_scrtch = sorted(list_A3B_scratch)

for a in list_A3B_scratch:
    s_files = os.listdir(DB_dir + '/' + a + '/Calculations')
    for b in s_files:
        c, d = a.split('3')
        fin_traj_path = (DB_dir + '/' + a + '/Calculations/' + b + '/final_'
                + b + '_' + c + '9' + d + '3.traj')
        if os.path.exists(fin_traj_path):
            print('Processing file:{0}'.format(a + '_' + b))
            traj_file = (DB_dir + '/' + a + '/Calculations/' + b + 
                    '/Initial_slab_with_adsorbate.traj')
            atoms = read(traj_file)
            keys['symbol'] = a
            ads, st, _, _ = b.split('_')
            keys['adsorbate'] = ads
            keys['site'] = st
            keys['SBS_Symbol'] = 'L12'

            db.write(atoms, key_value_pairs=keys, data=parameters) 



#Bi calculations                                                 
for ID in IDS:                                                                  
    print('Processing id: {0}'.format(ID))                                      
    launch = launchpad.get_fw_dict_by_id(ID)                                    
                                                                                
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']    
    atoms = encode_to_atoms(encoding)                                           
                                                                                
