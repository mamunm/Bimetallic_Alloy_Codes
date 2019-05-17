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
#pdata = namedtuple('PD', 'name top bridge FCC HCP')
Adata = {k: {pm: defaultdict(tuple) for pm in pmetals} for k in data}

for d in data:
    f_loc = '/Users/osmanmamun/Delta_Learning_Paper/get_db_fp'
    db = connect(f_loc + '/{}_atoms.db'.format(d))
    for pm in pmetals:
        name = pm
        top, bridge, FCC, HCP = None, None, None, None
        for i in  db.select(sb_symbol="A1", adsorbate=d):
            if not Counter(i.formula) == Counter(pm + '12' + d):
                continue
            if i.site == 'top':
                top = i.reaction_energy
            if i.site == 'bridge':
                bridge = i.reaction_energy
            if i.site == 'hollow':
                if 'FCC' in i.site_type:
                    FCC = i.reaction_energy
                if 'HCP' in i.site_type:
                    HCP = i.reaction_energy
        Adata[d][pm] = (name, top, bridge, FCC, HCP)

with open('pure_data.pkl', 'wb') as f:
    pickle.dump(Adata, f, pickle.HIGHEST_PROTOCOL)
