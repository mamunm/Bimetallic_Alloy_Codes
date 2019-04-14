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
import sys

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

adsorbate = sys.argv[1]
tr_ratio = [float(i) for i in sys.argv[2:]]

X, y, metadata = get_data(adsorbate)
data = preprocess_data(X, y)
data.clean_data()
X, y = data.get_data()

kernel_recipe = {'ConstantKernel' : [{'RBF' : [1.0,
                                              {'length_scale' : 1.0}]},
                                      {'constant_value' : 1.0,
                                       'constant_value_bounds' : (3e-7, 3e7)}],
                 'WhiteKernel' : {'noise_level' : 0.1,
                                  'noise_level_bounds' : (1e-5, 1e5)}}

r_state = [15, 25, 45]
trr_data = {tr: {rs : {} for rs in r_state} for tr in tr_ratio}

for trr in tr_ratio:
    for rs in r_state:
        print('    Working on random state: {}'.format(rs))
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


prev_data = np.load('run_{0}/{0}_rsdata.npy'.format(adsorbate))[()]
for trr in tr_ratio:
    prev_data[trr] = trr_data[trr]
np.save('run_{0}/{0}_rsdata.npy'.format(adsorbate), prev_data)

