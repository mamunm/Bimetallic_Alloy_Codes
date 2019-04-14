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




ids = [26486, 26487]

for i in ids:
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_by_id(i)
    atoms = encode_to_atoms(launch.spec['_tasks'][0]['args'][0])[0]

    launchpad.defuse_fw(i)
    
    atoms.info['convergence']['energy'] = 1e-4
    atoms.info['convergence']['mixing'] = 0.2
    atoms.info['fmax'] = 0.07
    atoms.info['sigma'] = 0.25
    atoms.set_initial_magnetic_moments([2.0 if i.symbol == 'Mn' else 0 
                                        for i in atoms])

    atoms._calc = None
    encoding = atoms_to_encode(atoms)

    launchpad.reignite_fw(i)
    launchpad.update_spec([i], spec_document={'_tasks.0.args.0': encoding})
