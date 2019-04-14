#!/usr/bin/env python
from fireworks import LaunchPad, Firework, Workflow, PyTask
from ase.io import read
from qescripts.fwio import atoms_to_encode, encode_to_atoms
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



IDS = launchpad.get_fw_ids(query={'state': 'COMPLETED'})
IDS = [i for i in IDS if i > 4610]
IDS = IDS[::-1]
pol_metal = ['Fe', 'Co', 'Ni', 'Mn']
weak_pol = ['Co']
strong_pol = ['Fe', 'Ni', 'Mn']
count = 0
j = 1
for ID in IDS:    #print ID
#for db in db_error.select():
    #ID = db.fw_id

    launch = launchpad.get_fw_dict_by_id(ID)

    #atoms = encode_to_atoms(launch['spec']['_tasks'][0]['args'][0])[-1]
    atoms = encode_to_atoms(launch['launches'][-1]['action']['stored_data']['trajectory'])[-1]
    formula = atoms.get_chemical_formula()
    
    if not np.any([p in formula for p in pol_metal]):
        continue
    
    data = atoms.info
    if not data['spinpol']:
        data['spinpol'] = True
    
    if not np.all(atoms.get_initial_magnetic_moments()==0):
        continue
    

    atoms.info = data
    atoms.set_initial_magnetic_moments([2.2 if atom.symbol == 'Fe' else 1.7 if atom.symbol == 'Co' else 0.6 if atom.symbol == 'Ni' else 1.7 if atom.symbol == 'Mn' else 0 for atom in atoms])
                   
    atoms._calc = None

    encoding = atoms_to_encode(atoms)


    print(ID)
    #continue

    launchpad.rerun_fw(ID)
    
    launchpad.update_spec([ID], spec_document={'_tasks.0.args.0': encoding})
