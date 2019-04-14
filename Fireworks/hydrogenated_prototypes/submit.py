#!/usr/bin/env python
# -*-coding: utf-8 -*-

#submit.py
#Osman Mamun
#LAST UPDATED: 06-11-2018
from fireworks import LaunchPad, Firework, Workflow, PyTask
from ase.db import connect
from catkit.flow.fwio import atoms_to_encode

host = 'suncatls2.slac.stanford.edu'  
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)

db = connect('hydrogenated-prototypes-mamunm.db')

for d in db.select():
    atoms = d.toatoms()
    atoms.info = d.data.parameters
    search_keys = d.key_value_pairs

    encoding = atoms_to_encode(atoms)

    t0 = PyTask(
        func='catkit.flow.fwio.encode_to_atoms',
        args=[encoding])

    t1 = PyTask(
        func='catkit.flow.fwespresso.get_potential_energy',
        stored_data_varname='trajectory')

    firework = Firework(
        [t0, t1],
        spec={'tags': search_keys,
              '_priority': 3,
              'connectivity': d.data.connectivity,
              'surface_atoms': d.data.surface_atoms},
        name='hydrogenated-species-1')

    workflow = Workflow([firework])
    launchpad.add_wf(workflow)
