#!/usr/bin/env python
# -*-coding: utf-8 -*-

#Resubmit_pol.py
#Osman Mamun
#LAST UPDATED: 05-22-2018

from fireworks import LaunchPad, Firework, Workflow, PyTask
from qescripts.fwio import atoms_to_encode, encode_to_atoms
import numpy as np
import re
import os

host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)

ids = np.load('../identify_problem/problems.npy')[()]
ids = [int(i) for i in ids]

for ID in ids:
    print('Processing: {0}'.format(ID))
    launch = launchpad.get_fw_dict_by_id(ID)
    atoms = encode_to_atoms(launch['spec']['_tasks'][0]['args'][0])[0]
    launchpad.rerun_fw(ID)


if os.path.exists('input.traj'):
    os.remove('input.traj')
