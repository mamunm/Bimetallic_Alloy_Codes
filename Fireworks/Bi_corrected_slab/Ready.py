#!/usr/bin/env python
from fireworks import LaunchPad, Firework, Workflow, PyTask
from ase.io import read
from qescripts.fwio import atoms_to_encode, encode_to_atoms
from ase.db import connect
from ase.visualize import view

host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)


ids = launchpad.get_fw_ids(
    query={'state': 'READY'}
)

iddd = []
for i in ids:
    if i > 24200:
        print('Processing: {0}'.format(i))
        launch = launchpad.get_fw_by_id(i)
        atoms = encode_to_atoms(launch.spec['_tasks'][0]['args'][0])[0]
        number_atoms = len(atoms.get_atomic_numbers()) 

        if number_atoms == 12:
            iddd.append(i)

print(iddd)
print(len(iddd))
import numpy as np
np.save('iidd.npy', iddd)
