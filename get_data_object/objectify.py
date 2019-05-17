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

for i in Adata:
    for j, k in Adata[i].items():
        d = [m for m in k[1:] if m is not None]
        if d:
            Adata[i][j] = min(d)
        else:
            Adata[i][j] = None

def get_E(A):
    if A[2] == 'top':
        q1 = Adata[A[0]][A[3]] if Adata[A[0]][A[3]] is not None else 'Nan'
        q2 = Adata[A[0][0]][A[3]] if Adata[A[0][0]][A[3]] is not None else 'Nan'
        return q1, q2
    if A[2] == 'bridge':
        m1, m2 = A[3].split('|')[0].split('_')
        try:
            q1 = 0.5 * (Adata[A[0]][m1] + Adata[A[0]][m2])
        except TypeError:
            q1 = 'Nan'
        try:
            q2 = 0.5 * (Adata[A[0][0]][m1] + Adata[A[0][0]][m2])
        except TypeError:
            q2 = 'Nan'
        return q1, q2
    if A[2] == 'hollow':
        m1, hol = A[3].split('|')
        m1, m2, m3 = m1.split('_')
        if hol == 'FCC':
            try:
                q1 = 0.33 * (Adata[A[0]][m1] + 
                             Adata[A[0]][m2] + 
                             Adata[A[0]][m3])
            except TypeError:
                q1 = 'Nan'
            try:
                q2 = 0.33 * (Adata[A[0][0]][m1] + 
                             Adata[A[0][0]][m2] + 
                             Adata[A[0][0]][m3])
            except TypeError:
                q2 = 'Nan'
            return q1, q2
        if hol == 'HCP':
            try:
                q1 = 0.33 * (Adata[A[0]][m1] + 
                             Adata[A[0]][m2] + 
                             Adata[A[0]][m3])
            except TypeError:
                q1 = 'Nan'
            try:
                q2 = 0.33 * (Adata[A[0][0]][m1] + 
                             Adata[A[0][0]][m2] + 
                             Adata[A[0][0]][m3])
            except TypeError:
                q2 = 'Nan'
            return q1, q2

data = ['CH', 'CH2', 'CH3', 'OH', 'NH', 'SH']
for d in data:
    #f_loc = '/Users/osmanmamun/Delta_Learning_Paper/get_db_fp/{}_fp.db'
    with FingerprintDB('{}_fp.db'.format(d)) as f:
        metadata = f.get_metadata(params=np.arange(1, 6))
        f.parameter_entry('EAP', 'EA based on pure metal')
        f.parameter_entry('EAxP', 'EAx based on pure metal')
        for i, m in enumerate(metadata, 1):
            EA, EAx = get_E(m)
            print(m, EA, EAx)
            f.fingerprint_entry(i, 108, EA)
            f.fingerprint_entry(i, 109, EAx)




