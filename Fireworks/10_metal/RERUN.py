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
    query={'state': 'DEFUSED'}
)

adsorbate = ['C', 'O', 'H', 'S', 'N']
ten_metal = ['Pt', 'Pd', 'Rh', 'Ru', 'Re', 'Au', 'Bi', 'Co', 'Ag', 'Cu']

for i in ids:
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_by_id(i)
    atoms = encode_to_atoms(launch.spec['_tasks'][0]['args'][0])[0]
    
    try:
        metal = launch.name['calc']['symbol']
        ads = launch.name['calc']['adsorbate']
    except KeyError:
        continue



    if ads in adsorbate:
        m1, m2 = metal.split('3')
        if all([m1 in ten_metal, m2 in ten_metal]):
            launchpad.reignite_fw(i)

