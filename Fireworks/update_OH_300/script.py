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

IDS = np.arange(51775, 52111)

for i in IDS:
   print('Processing: {0}'.format(i))
   '''
   launch = launchpad.get_fw_by_id(i)
   atoms = encode_to_atoms(launch.spec['_tasks'][0]['args'][0])[0]
   param = atoms.info['calculator_parameters']
   atoms.info = param

   atoms._calc = None
   encoding = atoms_to_encode(atoms)

   launchpad.update_spec([i], spec_document={'_tasks.0.args.0': encoding})
   '''
   launchpad.delete_wf(i)
if os.path.exists('input.traj'):
    os.remove('input.traj')


