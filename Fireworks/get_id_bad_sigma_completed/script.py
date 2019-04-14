from fireworks import LaunchPad
import ase
from ase.io import Trajectory, read
from ase.visualize import view
from ase import Atoms
import json
import os
import numpy as np
import datetime
import sys 
from pprint import pprint
import re

# Read credentials from a secure location
host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)

# Select the ID of the first completed calcualtion
IDS = launchpad.get_fw_ids(query={'state': 'COMPLETED'})
IDS = sorted(IDS)

a = 0
b = 0
FWID = []

for ID in IDS:
    print(str(ID) + '\n')
    launch = launchpad.get_fw_dict_by_id(ID)
    field = launch['spec']['_tasks'][0]['args'][0]
    m = re.search('"sigma": 0.\d\d', field)
    print(m.group(0)) 
    
    with open('out.txt', 'a') as f:
        f.write('ID: {}   {}\n'.format(ID, m.group(0))) 
