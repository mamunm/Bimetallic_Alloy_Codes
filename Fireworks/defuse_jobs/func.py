#!/usr/bin/env python
from fireworks import LaunchPad, Firework, Workflow, PyTask
from ase.io import read
from catkit.flow.fwio import atoms_to_encode, encode_to_atoms
from ase.db import connect
from ase.visualize import view

host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)




IDS = launchpad.get_fw_ids(query={'state': 'READY'})
IDS = [i for i in IDS if i < 30390]
for i in IDS:
    print('Processing: {0}'.format(i))
    launchpad.defuse_fw(i)

