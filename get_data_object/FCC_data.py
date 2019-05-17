#!/usr/bin/env python
# -*-coding: utf-8 -*-

import numpy as np
from collections import namedtuple
from collections import Counter
from ase.db import connect
from collections import defaultdict
from catkit.mydb import FingerprintDB
import pickle

data = {'C': 'Carbon', 
        'O': 'Oxygen', 
        'N': 'Nitrogen', 
        'H': 'Huydrogen', 
        'S': 'Sulfur', 
        'CH': 'CH', 
        'CH2': 'CH2', 
        'CH3': 'CH3', 
        'OH': 'OH', 
        'NH': 'NH', 
        'SH': 'SH'}
pmetals = ['Y', 'La', 'Sc', 'Zr', 'Hf', 'Ti', 'Ta', 'Nb', 'V',
           'Cr', 'Mo', 'W', 'Re', 'Tc', 'Os', 'Ru', 'Ir','Rh',
           'Ni', 'Co', 'Fe', 'Mn', 'Pt', 'Pd', 'Au', 'Ag', 'Cu',
           'Zn',  'Cd', 'Hg', 'Al',  'Ga',  'In', 'Tl', 'Pb', 'Sn', 'Bi']


with open('pure_data.pkl', 'rb') as f:
    Adata = pickle.load(f)

def get_E(A, B, C, D):
    'metal1, metal2, metal3, adsorbate'
    try:
        q1 = 0.33 * (Adata[D][A][3] + 
                     Adata[D][B][3] + 
                     Adata[D][C][3])
    except TypeError:
        q1 = None
    try:
        q2 = 0.33 * (Adata[D[0]][A][3] + 
                     Adata[D[0]][B][3] + 
                     Adata[D[0]][C][3])
    except TypeError:
        q2 = None
    return q1, q2

data = ['CH', 'CH2', 'CH3', 'OH', 'NH', 'SH']
jacobsen_data = {d: None for d in data}
for d in data:
    #f_loc = '/Users/osmanmamun/Delta_Learning_Paper/get_db_fp/{}_fp.db'
    with FingerprintDB('{}_fp.db'.format(d)) as f:
        metadata = f.get_metadata(params=np.arange(1, 6))
        energies = f.get_fingerprints(params=np.array([107]))
    E = []
    EA = []
    EAx = []
    md = []
    for i, m in enumerate(metadata, 1):
        if m[-1] == 'A1':
            continue
        if not 'FCC' in m[3]:
            continue
        A, B, C = m[3].split('|')[0].split('_')
        E += [energies[i]]
        A, B = get_E(A, B, C, d)
        EA += [A]
        EAx += [B]
        md += [m]
        print(m)
        print(energies[i][0], A, B)
    E = np.array(E).reshape(-1)
    EA = np.array(EA)
    EAx = np.array(EAx)
    jacobsen_data[d] = {'E': E, 'EA': EA, 'EAx': EAx, 'md': md}

np.save('jacobsen_data_FCC.npy', jacobsen_data)
        
        




