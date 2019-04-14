#!/usr/bin/env python
# -*-coding: utf-8 -*-

#Resubmit_pol.py
#Osman Mamun
#LAST UPDATED: 05-22-2018

from fireworks import LaunchPad, Firework, Workflow, PyTask
from qescripts.fwio import atoms_to_encode, encode_to_atoms
import numpy as np
import re

host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)


ids = sorted(launchpad.get_fw_ids(
    query={'state': 'COMPLETED'}))

ids = [603]
culprit = []
halar_po = []
slac_id = []

for i in ids:
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_dict_by_id(i)
    launchdir = launch['launches'][-1]['launch_dir']
    cluster = launch['launches'][-1]['fworker']['name']
    print(cluster)
    if cluster == 'SLAC':
        print('Dis calculation was done in SLAC cluster!')
        slac_id += [i]
        continue
    
    encoding = launch['launches'][-1]['action']['stored_data']['trajectory']
    atoms = encode_to_atoms(encoding)

    if len(atoms[-1]) != 13:
        continue

    if len(set(atoms[-1].get_chemical_symbols()[-5:-1])) == 1:
        m1, m2 = atoms[-1].get_chemical_symbols()[-1], None
    
    if len(set(atoms[-1].get_chemical_symbols()[-5:-1])) == 2:                  
        m1, m2 = set(atoms[-1].get_chemical_symbols()[-5:-1])  
    
    print(m1, m2)
    pol_metal = ['Fe', 'Co', 'Ni', 'Mn']

    if m1 in pol_metal or m2 in pol_metal:
        if atoms[-1].info['spinpol'] == False:
            print('Dis iz a culprit. Kill it before it lay eggs.') 
            print(m1, m2)
            culprit += [i]
  
        if atoms[-1].info['spinpol'] == True:
            logdir = launchdir + '/log'
            with open(logdir, 'r') as f:
                for line in f:
                    if re.search('total magnet', line):
                        magline = line
            mag = magline.split('=')[-1].split('Bohr')[0].strip(' ')
            print 'spinpol: {}'.format(mag)
            if np.isclose(float(mag), 0):
                print 'mag is zero!'
                print launchdir + '/log'
                halar_po += [i]

np.save('culprit.npy', np.asarray(culprit))
np.save('slac_id.npy', np.asarray(slac_id))
np.save('halar_po.npy', np.asarray(halar_po))
