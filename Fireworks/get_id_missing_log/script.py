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

time = datetime.datetime(2018, 1,1)
# Read credentials from a secure location
host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)

# Select the ID of the first completed calcualtion
IDS = launchpad.get_fw_ids(query={'state': 'COMPLETED'})#, 'time_end': {"$gt" : time}})
IDS = sorted(IDS)

a = 0
b = 0
FWID = []

for ID in IDS:
    print(str(ID) + '\n')
    launch = launchpad.get_fw_dict_by_id(ID)
    ldir = launch['launches'][-1]['launch_dir'].encode('utf-8')
    cluster = launch['launches'][-1]['fworker']['name'].encode('utf-8')
    if cluster == "SLAC":
        a += 1
        print(ldir, cluster)
        if os.path.exists(ldir):
            print('Path exists')
        else:
            print("path don't exist.")
            b += 1
            FWID += [ID]

print(a, b)
np.save('fwid.npy', FWID)
