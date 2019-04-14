#!/usr/bin/env python
# -*-coding: utf-8 -*-

#script.py
#Osman Mamun
#DATE CREATED: 10-22-2018

import numpy as np
from catkit.mydb import FingerprintDB
from catdelta.catscaling import LinearScaling 
import os
from collections import Counter

def is_pol(A):
    pol_metal = ['Fe', 'Co', 'Ni', 'Mn']
    if '3' in A:
        a, b = A.split('3')
        if a in pol_metal or b in pol_metal:
            return True
    elif '2' in A:
        a, b, _ = A.split('2')
        if a in pol_metal or b in pol_metal:
            return True
    else:
        if A in pol_metal:
            return True
    return False

def get_refined_data(X, Y, M, S, I, MAE):
    x = []
    y = []
    m = []
    diagnostic = []
    for xx, yy, mm in zip(X, Y, M):
        if abs(yy - S * xx - I) > 3 * MAE:
            if not is_pol(mm[1]):
                x += [xx]
                y += [yy]
                m += [mm]
                diagnostic += [mm]
        else:
            x += [xx]
            y += [yy]
            m += [mm]
    return x, y, m, diagnostic

def get_data(A, l_A, C, l_C):

    def get_cond(A, B):
        if not set(A) == set(B):
            return False
        if '|' in A:
            if not A.split('|')[1] == B.split('|')[1]:
                return False
            if not Counter(A.split('|')[0].split('_')) == \
                    Counter(B.split('|')[0].split('_')):
                return False
        else:
            if not Counter(A.split('_')) == Counter(B.split('_')):
                return False
        return True


    with FingerprintDB(A) as fp:
        data_A = fp.get_fingerprints(params=np.arange(1, l_A))
        data_A = data_A[:, -1]
        metadata_A = fp.get_metadata(params=np.arange(1, 6))
    with FingerprintDB(C) as fp:
        data_C = fp.get_fingerprints(params=np.arange(1, l_C))
        data_C = data_C[:, -1]
        metadata_C = fp.get_metadata(params=np.arange(1, 6))

    x = []
    y = []
    metadata = []

    for a, b in zip(data_A, metadata_A):
        for c, d in zip(data_C, metadata_C):
            if all([b[1]==d[1],
                    b[2]==d[2],
                    get_cond(b[3], d[3]),
                    b[4]==d[4]]):
                x += [a]
                y += [c]
                metadata += [[b[0], b[1], b[2], b[3], b[4],
                              d[0], d[1], d[2], d[3], d[4]]]

    return x, y, metadata

adsorbates_pair = [('C', 'CH'),
                   ('C', 'CH2'),
                   ('C', 'CH3'),
                   ('O', 'OH'),
                   ('S', 'SH'),
                   ('N', 'NH')]

for ads1, ads2 in adsorbates_pair:
    print('Working on adsorbates: {} and {}'.format(ads1, ads2))
    l_dict = {a : 97 if a in ['C', 'H', 'N', 'S', 'O'] else 108
              for a in ['C', 'H', 'N', 'S', 'O', 'CH', 'CH2', 'CH3',
                        'OH', 'SH', 'NH']}
    l_ads1 = l_dict[ads1]
    l_ads2 = l_dict[ads2]

    x, y, metadata = get_data('../get_db_fp/{}_fp.db'.format(ads1),
                              l_ads1,
                              '../get_db_fp/{}_fp.db'.format(ads2),
                              l_ads2)
    sc = LinearScaling(np.array(x), np.array(y), ads1, ads2, 'eV')
    sc.get_coeff()
    a = sc.plot_scaling()

    if not os.path.exists('scaling_{}_{}'.format(ads1, ads2)):
        os.mkdir('scaling_{}_{}'.format(ads1, ads2))

    scaling_dict = sc.__dict__
    scaling_dict['metadata'] = metadata
    a.savefig('scaling_{}_{}/scaling_plot_all_data.png'.format(ads1, ads2))
    np.save('scaling_{}_{}/all_data.npy'.format(ads1, ads2), scaling_dict)

    x, y, metadata, diagnostic = get_refined_data(x, y, metadata,
                                      sc.slope, sc.intercept, sc.MAE)
    sc = LinearScaling(np.array(x), np.array(y), ads1, ads2, 'eV')
    sc.get_coeff()
    a = sc.plot_scaling()

    scaling_dict = sc.__dict__
    scaling_dict['metadata'] = metadata
    a.savefig('scaling_{}_{}/scaling_plot.png'.format(ads1, ads2))
    np.save('scaling_{}_{}/data.npy'.format(ads1, ads2), scaling_dict)
    np.save('scaling_{}_{}/diagnostic.npy'.format(ads1, ads2), diagnostic)
