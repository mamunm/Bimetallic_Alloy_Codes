#!/usr/bin/env python
# -*-coding: utf-8 -*-

import numpy as np
from ase.io import read
from ase.visualize import view

coord = []
sym = []
lattice = []
readn = False
read = False
readlat = False
with open('pw.pwo', 'r') as f:
    for line in f:
        if readlat and len(lattice) < 3:
            lattice += [[float(i) for i in line.split()[3:6]]]
        if read:
            if not 'End' in line:
                temp = line.split()[:4]
                if '1' in temp[0]:
                    temp[0] = temp[0][:-1]
                coord += [[float(i) for i in temp[1:]]]
                sym += [temp[0]]
        if 'Begin final coordinates' in line:
            readn = True
        if 'ATOMIC_POSITIONS' in line and readn:
            read = True
        if 'End final coordinates' in line:
            read = False
        if 'crystal axes: (cart. coord. in units of alat)' in line:
            readlat = True

coord = np.array(coord)
lattice = np.array(lattice)
print(lattice)
print(coord)
real_coord = np.matmul(coord, lattice.T)
print(real_coord)
'''
with open('a.xyz', 'w') as f:
    f.write(str(len(sym)))
    f.write('\n \n')
    for i, s in enumerate(sym):
        stringg = s + ' ' + ' '.join([str(j) for j in real_coord[i].tolist()])
        f.write(stringg)
        f.write('\n')

from ase.io import read
ml = read('a.xyz')
view(ml)
'''
