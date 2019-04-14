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

from cathub.query import *

references = {}
references['O'] = 'H2O+H2'
references['OH'] = 'H2O+H2'

adsorbates = ['O', 'OH']

for adsorbate in adsorbates:
   data = get_reactions(n_results='all',
                        pubId='MamunHighT2019',
                        #sites=site_choice,
                        reactants=references[adsorbate],
                        products=adsorbate,
                        columns=['surfaceComposition, reactionEnergy',
                                 'sites', 'reactants', 'products'],
                        subtables=['reactionSystems'])


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


