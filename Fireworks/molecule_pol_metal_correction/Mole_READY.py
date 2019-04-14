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

ids = [i for i in ids if i > 26000]

for i in ids:
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_by_id(i)
    atoms = encode_to_atoms(launch.spec['_tasks'][0]['args'][0])[0]
    
    try:
        mol = launch.name['calc']['molecule']
    except KeyError:
        continue


    launchpad.defuse_fw(i)
    atoms.set_initial_magnetic_moments([2.0 for atom in atoms])

    atoms._calc = None
    encoding = atoms_to_encode(atoms)

    launchpad.reignite_fw(i)
    launchpad.update_spec([i], spec_document={'_tasks.0.args.0': encoding})
