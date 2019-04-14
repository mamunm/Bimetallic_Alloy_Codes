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
from ase.db import connect

host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)

IDS = [i for i in range(28670, 29706)]

for i in IDS:
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_by_id(i)
    info = launch.to_dict()['name']['calc']
    db = connect('new_ab_A.db')
    rerun = False
    for a in db.select():
        kvp = a.key_value_pairs
        cond1 = info['metal'] == kvp['metal']
        cond2 = info['adsorbate'] == kvp['adsorbate']
        cond3 = info['site'] == kvp['site']
        if cond1 and cond2 and cond3:
            atoms = a.toatoms()
            rerun = True
    if rerun:
        atoms._calc = None
        encoding = atoms_to_encode(atoms)
        print('Reruning FW job: {} \n with spec: {}'.format(i, info))
        launchpad.rerun_fw(i)
        launchpad.update_spec([i], spec_document={'_tasks.0.args.0': encoding})
if os.path.exists('input.traj'):
    os.remove('input.traj')

