#!/usr/bin/env python
# -*-coding: utf-8 -*-

#script.py
#Osman Mamun
#DATE CREATED: 12-03-2018

from fireworks import LaunchPad, Firework, Workflow, PyTask
from ase.io import read
from catkit.flow.fwio import atoms_to_encode, encode_to_atoms
from ase.db import connect
from ase.visualize import view
import numpy as np

host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
        host=host,
        name=name,
        username=username,
        password=password)

#IDS = launchpad.get_fw_ids(query={'state': 'READY'})
#IDS = [i for i in IDS if i < 52000]
IDS = np.load('run_later.npy')[()]
for i in IDS:
    print('Processing: {0}'.format(i))
    launchpad.reignite_fw(i)
