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

IDS = sorted(launchpad.get_fw_ids(query={"state": 'READY'}))

for i in IDS:
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_by_id(i)
    atoms = encode_to_atoms(launch.spec['_tasks'][0]['args'][0])[0]

    atoms.info['convergence']['energy'] = 1e-4
    atoms.info['convergence']['mixing'] = 0.2
    atoms.info['fmax'] = 0.07
    atoms.info['sigma'] = 0.25

    atoms._calc = None
    encoding = atoms_to_encode(atoms)

    launchpad.update_spec([i], spec_document={'_tasks.0.args.0': encoding})

if os.path.exists('input.traj'):
    os.remove('input.traj')

