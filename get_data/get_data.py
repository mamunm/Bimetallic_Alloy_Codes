#!/usr/bin/env python
# -*-coding: utf-8 -*-

#script.py
#Osman Mamun
#DATE CREATED: 10-19-2018
import sys
import json
import re
import numpy as np
import matplotlib
import pylab as p
from ase.db import connect
from ase.visualize import view

from catkit.hub.query import *

def update_data(OLD, NEW):
    R_D = {'reactions': {'edges': [],
                         'totalCount': NEW['reactions']['totalCount']}}

    if OLD['reactions']['totalCount'] == NEW['reactions']['totalCount']:
        return NEW
    else:
        temp = OLD['reactions']['edges']
        for i in NEW['reactions']['edges']:
            if i not in OLD['reactions']['edges']:
                temp += [i]

        R_D['reactions']['edges'] = temp
        return R_D


ordered_metals = ['Y', 'La', 'Sc', 'Zr', 'Hf', 'Ti', 'Ta', 'Nb', 'V',  'Cr',
                  'Mo', 'W', 'Re', 'Tc', 'Os', 'Ru', 'Ir','Rh', 'Ni', 'Co',
                  'Fe', 'Mn', 'Pt', 'Pd', 'Au', 'Ag', 'Cu', 'Zn',  'Cd', 'Hg',
                  'Al',  'Ga',  'In', 'Tl', 'Pb', 'Sn', 'Bi']

references = {}
references['C'] = 'CH4+H2'
references['H'] = 'H2'
references['N'] = 'N2'
references['O'] = 'H2O+H2'
references['S'] = 'H2S'
references['CH'] = 'CH4+H2'
references['CH2'] = 'CH4+H2'
references['CH3'] = 'CH4+H2'
references['NH'] = 'N2+H2'
references['SH'] = 'H2S+H2'
references['OH'] = 'H2O+H2'
references['OH2'] = 'H2O'

adsorbates = ['C', 'H', 'N', 'O', 'S',
              'CH', 'CH2', 'CH3', 'NH', 'SH', 'OH']

for adsorbate in adsorbates:
   data = get_reactions(n_results='all',
                        pubId='MamunHighT2019',
                        #sites=site_choice,
                        reactants=references[adsorbate],
                        products=adsorbate,
                        columns=['surfaceComposition, reactionEnergy',
                                 'sites', 'reactants', 'products'])


   if not os.path.exists('{}_data.npy'.format(adsorbate)):
       print('Number of count for {}: {}'.format(adsorbate, 
             data['reactions']['totalCount']))
       np.save('{}_data.npy'.format(adsorbate), data)
   else:
       old_data = np.load('{}_data.npy'.format(adsorbate))[()]
       data = update_data(old_data, data)
       print('Number of count for {}: {}'.format(adsorbate, 
             data['reactions']['totalCount']))
       np.save('{}_data.npy'.format(adsorbate), data)


##data structue and how to retrieve it:
#import numpy as np
#data = np.load('C_data.npy')[()]
##check the length of the data
#print(data['reactions']['totalCount'])
##To loop over structures
#for d in data['reactions']['edges']:
#    print(d)


