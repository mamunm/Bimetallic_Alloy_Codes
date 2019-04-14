#!/usr/bin/env python
# -*-coding: utf-8 -*-

#Resubmit_pol.py
#Osman Mamun
#LAST UPDATED: 05-22-2018

from fireworks import LaunchPad, Firework, Workflow, PyTask
from catkit.flow.fwio import atoms_to_encode, encode_to_atoms
import numpy as np
import re
import os
from ase.constraints import FixAtoms

host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)

#IDS = sorted(launchpad.get_fw_ids(query={"state": 'COMPLETED'}))
IDS = [i for i in range(29700, 30700)]
rs = 0

for ID in IDS:
    print('Processing: {0}'.format(ID))
    launch = launchpad.get_fw_by_id(ID)
    atoms = encode_to_atoms(launch.spec['_tasks'][0]['args'][0])[0]
    if len(atoms) > 10 and len(atoms.constraints[0].get_indices()) < 8:
        print(len(atoms))
        print(atoms.constraints[0].get_indices())
        constraint = FixAtoms(indices=[i for i in range(8)])
        atoms.set_constraint(constraint)
        rs += 1
        print(atoms.constraints[0].get_indices())
        atoms._calc = None
        encoding = atoms_to_encode(atoms)
        print(encoding)
        launchpad.update_spec([ID], spec_document={'_tasks.0.args.0': encoding})  
        launchpad.rerun_fw(ID)

if os.path.exists('input.traj'):
    os.remove('input.traj')

print('Total number of resubmitted jobs: {}'.format(rs))
