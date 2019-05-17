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

FCC_data = np.load('../Jacobsen_data_FCC.npy')[()]
Bulk_data = np.load('/Users/osmanmamun/FireWorks/'
                    'JOB_PREP/enumerated_OH/bulk_data.npy')[()]

def get_data(A):
    y = FCC_data[A]['E']
    md = FCC_data[A]['md']
    lat_data = []
    for m in md:
        if m[-1] == 'L12':
            a, b = m[1].split('3')
            temp = 0.75 * Bulk_data[a]['a'] + 0.25 * Bulk_data[b]['a']
        if m[-1] == 'L10':
            a, b, _ = m[1].split('2')
            temp = 0.5 * Bulk_data[a]['a'] + 0.5 * Bulk_data[b]['a']
        lat_data += [temp]
    lat_data = np.array(lat_data)
    X = np.column_stack([FCC_data[A]['EA'], FCC_data[A]['EAx'], lat_data])
    return X, y

#run OH later
adsorbates = ['CH', 'CH2', 'CH3', 'SH', 'NH', 'OH']

for ads in adsorbates:
    print('Working on adsorbate: {}'.format(ads))
    X, y = get_data(ads)
    data = preprocess_data(X, y)
    data.clean_data()
    X, y = data.get_data()

    kernel_recipe = {'ConstantKernel' : [{'RBF' : [1.0,
                                                   {'length_scale' : 1.0}]},
                                      {'constant_value' : 1.0,
                                       'constant_value_bounds' : (3e-7, 3e7)}],
                     'WhiteKernel' : {'noise_level' : 0.1,
                                      'noise_level_bounds' : (1e-5, 1e5)}}

    r_state = [10, 20, 42, 80, 150, 200, 300, 400]
    #r_state = [42]
    r_data = {rs : {} for rs in r_state}
    for rs in r_state:
        print('    Working on random state: {}'.format(rs))
        X_train, X_test, y_train, y_test = \
            train_test_split(X, y, train_size=0.80, random_state=rs)

        MLGP = OMGP(X_train=X_train,
                    X_test=X_test,
                    y_train=y_train,
                    y_test=y_test,
                    kernel_recipe=kernel_recipe)

        MLGP.run_GP()
        r_data[rs] = MLGP.__dict__
        r_data[rs]['metadata_train'] = metadata_train
        r_data[rs]['metadata_test'] = metadata_test

    if not os.path.exists('run_{}'.format(ads)):
        os.mkdir('run_{}'.format(ads))

    np.save('run_{0}/{0}_rsdata.npy'.format(ads), r_data)
    #LC = MLGP.plot_learning_curve()
    #LC.savefig('run_{0}/{0}_learning_curve.png'.format(ads))
    PP = MLGP.parity_plot(data='train')
    PP.savefig('run_{0}/{0}_parity_plot_train.png'.format(ads))
    PP = MLGP.parity_plot(data='test')
    PP.savefig('run_{0}/{0}_parity_plot_test.png'.format(ads))
    PP = MLGP.parity_plot(data='train', err_bar=True)
    PP.savefig('run_{0}/{0}_parity_plot_train_err_bar.png'.format(ads))
    PP = MLGP.parity_plot(data='test', err_bar=True)
    PP.savefig('run_{0}/{0}_parity_plot_test_err_bar.png'.format(ads))
    #PP = MLGP.parity_plot_fancy(data='train')
    #PP.savefig('run_{0}/{0}_parity_plot_train.png'.format(ads))
    #PP = MLGP.parity_plot_fancy(data='test')
    #PP.savefig('run_{0}/{0}_parity_plot_test.png'.format(ads))

