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




IDS = launchpad.get_fw_ids(query={'state': 'COMPLETED'})
IDS = [i for i in IDS if i > 51775]
#IDS = [31032]
l_dir = []
for i in IDS:
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_dict_by_id(i)
    l_dir += [launch['launches'][-1]['launch_dir']]

print(sum([1 for _ in l_dir if 'nfs' in _]))
print(len(l_dir))

