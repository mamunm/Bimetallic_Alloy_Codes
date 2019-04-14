#!/usr/bin/env python
# -*-coding: utf-8 -*-

#script.py
#Osman Mamun
#DATE CREATED: 10-22-2018

import numpy as np
from catkit.mydb import FingerprintDB
from scaling import Linear_Scaling
import os
from collections import Counter
from catGP import OMGP, preprocess_data
from sklearn.model_selection import train_test_split
import sys

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
                y += [list(c)]
                metadata += [[b[0], b[1], b[2], b[3], b[4],
                              d[0], d[1], d[2], d[3], d[4]]]
    return np.array(x), np.array(y), metadata


adsorbates_pair = [(sys.argv[1], sys.argv[2])]
tr_ratio = [float(i) for i in sys.argv[3:]]

for ads1, ads2 in adsorbates_pair:
    print('Working on adsorbates: {} and {}'.format(ads1, ads2))
    l_dict = {a : 370 if a in ['C', 'H', 'N', 'S', 'O'] else 381
              for a in ['C', 'H', 'N', 'S', 'O', 'CH', 'CH2', 'CH3',
                        'OH', 'SH', 'NH']}
    l_ads1 = l_dict[ads1]
    l_ads2 = l_dict[ads2]

    s_x, s_y, metadata = get_data('../../../get_db_fp/{}_fp.db'.format(ads1),
                              l_ads1,
                              '../../../get_db_fp/{}_fp.db'.format(ads2),
                              l_ads2)

    if not os.path.exists('gp_scaling_{}_{}'.format(ads1, ads2)):
        os.mkdir('gp_scaling_{}_{}'.format(ads1, ads2))

    min_gp_scaling = {'x' : s_x, 'y' : s_y, 'metadata' : metadata}
    np.save('gp_scaling_{}_{}/min_scaling_data.npy'.format(ads1, ads2),
            min_gp_scaling)
    del min_gp_scaling

    sc = Linear_Scaling(s_x, s_y[:, -1], ads1, ads2, 'eV')
    sc.get_coeff()
    a = sc.plot_scaling()

    scaling_dict = sc.__dict__
    scaling_dict['metadata'] = metadata
    a.savefig('gp_scaling_{}_{}/scaling_plot.png'.format(ads1, ads2))
    np.save('gp_scaling_{}_{}/data.npy'.format(ads1, ads2), scaling_dict)

    print('Working on gp of adsorbate: {}'.format(ads2))

    X = s_y[:, :-1]
    y = s_y[:, -1]
    data = preprocess_data(X, y)
    data.clean_data()
    X, y = data.get_data()

    kernel_recipe = {'ConstantKernel' : [{'RBF' : [1.0,
                                                   {'length_scale' : 1.0}]},
                                      {'constant_value' : 1.0,
                                       'constant_value_bounds' : (3e-7, 3e7)}],
                     'WhiteKernel' : {'noise_level' : 0.1,
                                      'noise_level_bounds' : (1e-5, 1e5)}}

    tr_ratio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    r_state = [10, 20, 42]
    trr_data = {tr: {rs : {} for rs in r_state} for tr in tr_ratio}

    for trr in tr_ratio:
        print('    Working on train ratio: {}'.format(trr))
        for rs in r_state:
            print('    Working on random state: {}'.format(rs))
            X_train, X_test, y_train, y_test, scaling_y_train, \
            scaling_y_test, metadata_train, metadata_test = \
            train_test_split(X, y, s_x, metadata, \
            train_size=trr, random_state=rs)

            MLGP = OMGP(X_train=X_train,
                        X_test=X_test,
                        y_train=y_train,
                        y_test=y_test,
                        kernel_recipe=kernel_recipe,
                        scaling=True,
                        scaling_params={'alpha' : sc.slope,
                                        'gamma' : sc.intercept},
                        scaling_y_train=scaling_y_train,
                        scaling_y_test=scaling_y_test)

            MLGP.run_GP()
            trr_data[trr][rs] = MLGP.__dict__
            trr_data[trr][rs]['metadata_train'] = metadata_train
            trr_data[trr][rs]['metadata_test'] = metadata_test

ads1, ads2 = adsorbates_pair[0]
prev_data = np.load('gp_scaling_{0}_{1}/{1}_rsdata.npy'.format(ads1, ads2))[()]
for trr in tr_ratio:
    prev_data[trr] = trr_data[trr]
np.save('gp_scaling_{0}_{1}/{1}_rsdata.npy'.format(ads1, ads2), prev_data)
