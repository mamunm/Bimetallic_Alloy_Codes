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

for i in ids:
    if i > 26470:
        print('Processing: {0}'.format(i))
        launch = launchpad.get_fw_by_id(i)
        atoms = encode_to_atoms(launch.spec['_tasks'][0]['args'][0])[0]
    
        try:
            metal = launch.name['metal']
            print(metal)
            if '2' in metal:
                a, b, _ = metal.split('2')
            elif '3' in metal:
                a, b = metal.split('3')
            else:
                continue
        except KeyError:
            continue

        pol_metal = ['Fe', 'Co', 'Ni', 'Mn']

        if a in pol_metal or b in pol_metal:
            if atoms.info['spinpol'] == True:
                atoms.set_initial_magnetic_moments([2.0 if atom.symbol in pol_metal
                                              else 0 for atom in atoms])

                atoms._calc = None
                encoding = atoms_to_encode(atoms)

                launchpad.update_spec([i], spec_document={'_tasks.0.args.0': encoding})
