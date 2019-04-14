#!/usr/bin/env python
from fireworks import LaunchPad, Firework, Workflow, PyTask
from ase.io import read
from qescripts.fwio import atoms_to_encode, encode_to_atoms
from ase.db import connect
from ase.visualize import view
import numpy as np

def get_info(A):
    db = connect('slab.db')
    for i, d in enumerate(db.select(['site=hollow'])):
        slab_an  = d.toatoms().get_atomic_numbers()
        slab_pos  = d.toatoms().positions
        kvp = d.key_value_pairs
        an = A.get_atomic_numbers()
        pos = A.positions
        if all(an == slab_an):
            if all(pos[-1] == slab_pos[-1]):
                ss = kvp['site_id']
                aa = kvp['a']
                cc = kvp['c']
                b_id = kvp['bulk_id']

    print(ss, aa, cc, b_id)    
    return (ss, aa, cc, b_id)

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
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_by_id(i)
    
    
    try:
        symbol = launch.name['calc']['symbol']
        ads = launch.name['calc']['adsorbate']
    except KeyError:
        continue

    else:
        print('Changing spec for calc: {0}'.format(symbol))
        launch.name['calc']['facet'] = '(1, 1, 1)'
        launch.name['calc']['layer'] = 3
        launch.name['calc']['crystal'] = 'fcc'
        launch.name['calc']['vacuum'] = 10
        launch.name['calc']['SBS_symbol'] = 'L12'
        launch.name['calc']['PS_symbol'] = 'cP4'
        site_id, a, c, BID = get_info(atoms)
        launch.name['calc']['site_id'] = site_id
        launch.name['calc']['a'] = a
        launch.name['calc']['c'] = c


        launchpad.update_spec([i], spec_document={'calc': {'symbol': symbol,
            'adsorbate': ads,
            'facet': '(1, 1, 1)',
            'layer': 3,
            'crystal': 'fcc',
            'vacuum': 10,
            'SBS_symbol': 'L12',
            'PS_symbol': 'cP4',
            'site_id': site_id,
            'a': a,
            'c': c,
            'bulk_id': BID}})
        #launchpad.rerun_fw(i)
