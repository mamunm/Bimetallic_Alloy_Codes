#!/usr/bin/env python
# -*-coding: utf-8 -*-

#script.py
#Osman Mamun
#LAST UPDATED: 07-27-2018

from fireworks import LaunchPad
from catkit.flow.fwio import encode_to_atoms
import os
import numpy as np
import re
import subprocess

host = 'suncatls2.slac.stanford.edu'
username, name = 'mamun', 'mamunm_db'
password = os.environ['FW_PASS']

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)

IDS = launchpad.get_fw_ids(query={'state': 'COMPLETED'})


pol_metal = ['Fe', 'Ni', 'Co', 'Mn']



def check_initial_magnet(filename, atom, cluster):
    if cluster == 'SLAC':
        ssh = subprocess.Popen(['ssh', 'mamunm@suncatls1.slac.stanford.edu', 
            'cat', filename], stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE)
        input = ssh.stdout
    else:
        input = file(filename)

    found = False
    for il, line in enumerate(input):
        if found:
            if re.search(atom + '1', line):
                magline = line
                break
        if re.search('atomic species   magnetization', line):
            found = True

    mag = magline.split('   ')[-1].split('(')[0].strip(' ')
    return mag



for ID in IDS[::-1]:
    launch = launchpad.get_fw_dict_by_id(ID)

    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
    atoms = encode_to_atoms(encoding)

    formula = atoms[-1].get_chemical_formula()

    if 'Fe' in formula or 'Mn' in formula or 'Co' in formula or 'Ni' in formula:
        if 'Fe' in formula:
            atom = 'Fe'
        elif 'Mn' in formula:
            atom = 'Mn'
        elif 'Co' in formula:
            atom = 'Co'
        elif 'Ni' in formula:
            atom = 'Ni'
        imag = atoms[-1].get_initial_magnetic_moments()
        if imag is not None:
            if np.all(np.array(imag) == 0.0):
                print(ID, 'Initial magmoms are zero!')
            else:
                try:
                    mag = check_initial_magnet(
                            launch['launches'][-1]['launch_dir'] + 
                            '/log', atom, 
                            launch['launches'][-1]['fworker']['name'])
                    try:                        
                        if float(mag) == 0:
                            print(ID, 
                                  'Initial magmoms are zero in logfile!', 
                                  launch['launches'][-1]['launch_dir'] + 
                                  '/log' )
                        else:
                            print(ID, 'OK')
                    except:
                        print(ID, 'weird logfile!', 
                                launch['launches'][-1]['launch_dir'] + '/log')
                except:
                    print(ID, 'No logfile', 
                            launch['launches'][-1]['launch_dir'] + '/log')
        else:
            print(ID, 'No initial magmoms set!')

