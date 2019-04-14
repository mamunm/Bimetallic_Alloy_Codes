#!/usr/bin/env python
# -*-coding: utf-8 -*-

#script.py
#Osman Mamun
#DATE CREATED: 10-22-2018

import numpy as np
from catkit.mydb import FingerprintDB
from catGP import OMGP, preprocess_data
from scaling import Linear_Scaling
from sklearn.model_selection import train_test_split
import os

def get_data(A):
    def sanity_check(M, N):
        for m, n in zip(M, N):
            if any([mm not in n for mm in m]):
                return False
            else:
                return True

    dpath = '../../../get_db_fp/{}_fp.db'.format(A)
    with FingerprintDB(dpath) as fp:
        data = fp.get_fingerprints(params=np.arange(1, 108))
        metadata = fp.get_metadata(params=np.arange(1, 6))
    dpath = '../../../get_db_fp_catlearn/{}_fp_catlearn.db'.format(A)
    with FingerprintDB(dpath) as fp:
        data_cl = fp.get_fingerprints(params=np.arange(1, 1062))
        metadata2 = fp.get_metadata(params=np.arange(1, 6))
    target = data[:, -1]
    data = data[:, :-1]
    print(data.shape)
    print(data_cl.shape)
    features = np.concatenate((data, data_cl), axis=1)
    if sanity_check(metadata2, metadata):
        return features, target, metadata
    else:
        return None, None, None

adsorbates = ['CH', 'CH2', 'CH3', 'NH', 'SH', 'OH']

for ads in adsorbates:
    print('Working on adsorbate: {}'.format(ads))
    X, y, metadata = get_data(ads)
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
    #r_state = [42]
    trr_data = {tr: {rs : {} for rs in r_state} for tr in tr_ratio}

    for trr in tr_ratio:
        print('    Working on train ratio: {}'.format(trr))
        for rs in r_state:
            print('        Working on random state: {}'.format(rs))
            X_train, X_test, y_train, y_test, metadata_train, metadata_test = \
            train_test_split(X, y, metadata, train_size=trr, random_state=rs)

            MLGP = OMGP(X_train=X_train,
                        X_test=X_test,
                        y_train=y_train,
                        y_test=y_test,
                        kernel_recipe=kernel_recipe)

            MLGP.run_GP()
            trr_data[trr][rs] = MLGP.__dict__
            trr_data[trr][rs]['metadata_train'] = metadata_train
            trr_data[trr][rs]['metadata_test'] = metadata_test

    if not os.path.exists('run_{}'.format(ads)):
        os.mkdir('run_{}'.format(ads))

    np.save('run_{0}/{0}_rsdata.npy'.format(ads), trr_data)

