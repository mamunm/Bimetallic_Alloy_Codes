#!/usr/bin/env python
# -*-coding: utf-8 -*-

#script.py
#Osman Mamun
#DATE CREATED: 10-22-2018

import numpy as np
import re
from catkit.mydb import FingerprintDB
from scaling import Linear_Scaling
import os
from collections import Counter
from catGP import OMGP, preprocess_data
from sklearn.model_selection import train_test_split
import sys 

def get_data(A, B):

    def get_sym(A):
        ordered_metals = ['Y', 'La', 'Sc', 'Zr', 'Hf', 'Ti', 'Ta', 'Nb', 'V',
                          'Cr', 'Mo', 'W', 'Re', 'Tc', 'Os', 'Ru', 'Ir','Rh',
                          'Ni', 'Co', 'Fe', 'Mn', 'Pt', 'Pd', 'Au', 'Ag', 'Cu',
                          'Zn',  'Cd', 'Hg', 'Al',  'Ga',  'In', 'Tl', 'Pb',
                          'Sn', 'Bi']
        if '3' in A:
            return A, 'L12'
        elif A in ordered_metals:
            return A, 'A1'
        else:
            sym = []
            for i in ordered_metals:
                if i in A:
                    sym += [i]

            if not len(sym) == 2:
                A = 'ERROR!!'
            else:
                return '2'.join(sym) + '2', 'L10'

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

    def sanity_check(M, N):
        for m, n in zip(M, N):
            if any([mm not in n for mm in m]):
                return False
            else:
                return True

    dpath = '../../get_db_fp/{}_fp.db'.format(B)
    with FingerprintDB(dpath) as fp:
        data = fp.get_fingerprints(params=np.arange(1, 108))
        md1 = fp.get_metadata(params=np.arange(1, 6))
    dpath = '../../get_db_fp_catlearn/{}_fp_catlearn.db'.format(B)
    with FingerprintDB(dpath) as fp:
        data_cl = fp.get_fingerprints(params=np.arange(1, 1062))
        md2 = fp.get_metadata(params=np.arange(1, 6))
    temp_target = data[:, -1]
    data = data[:, :-1]
    temp_features = np.concatenate((data, data_cl), axis=1)
    del data, data_cl 
    if not sanity_check(md2, md1):
        sys.exit('Something went wrong in sanity check.')


    scaling_data_processed = [] 
    scaling_data = np.load('../../get_data_cathub/{}_data.npy'.format(A))[()]
    scaling_data = scaling_data['reactions']['edges']
    for d in scaling_data:
        if A in ['N', 'H']:
            if not '0.5' in d['node']['reactants']:
                continue
        if A == 'S':
            if 'Hstar' in d['node']['products']:
                continue

        adsorbate = re.findall(r'"([A-Z0-4]+)":', d['node']['sites'])[0]
        reactionenergy = d['node']['reactionEnergy']
        symbol, sb_symbol = get_sym(d['node']['surfaceComposition'])
        site, site_type = re.findall(r'(?<=: )"(.+)"', 
                d['node']['sites'])[0].split('|', 1)
        if site not in ['top', 'bridge', 'hollow']:
            print('Site {} not a valid site_type.'.format(site))
            continue
        if sb_symbol == 'A1':
            site_type = site_type.replace('A', symbol)
        elif sb_symbol == 'L12':
            a, b = symbol.split('3')
            site_type = site_type.replace('A', 'aaa').replace('B', 'bbb')
            site_type = site_type.replace('aaa', a).replace('bbb', b)
        else:
            a, b, _ = symbol.split('2')
            if a > b:
                a, b = b, a
            site_type = site_type.replace('A', 'aaa').replace('B', 'bbb')
            site_type = site_type.replace('aaa', a).replace('bbb', b)
        if ':' in site_type:
            print('This is a problematic site type.\n\n\n')
            continue 
        scaling_data_processed += [[adsorbate, symbol, site, 
            site_type, sb_symbol, reactionenergy]]

    scaling_target = []
    metadata = []
    features= []
    target = []
    
    for a, b, c in zip(md1, temp_features, temp_target):
        for d in scaling_data_processed:
            if all([a[1]==d[1],
                   a[2]==d[2],
                   get_cond(a[3], d[3]),
                   a[4]==d[4]]):
                features += [b]
                target += [c]
                scaling_target += [d[5]]
                metadata += [[a[0], a[1], a[2], a[3], a[4],
                              d[0], d[1], d[2], d[3], d[4]]]
    
    features = np.array(features)
    target = np.array(target)
    scaling_target = np.array(scaling_target)
    
    return features, target, scaling_target, metadata


adsorbates_pair = [('C', 'CH'),
                   ('C', 'CH2'),
                   ('C', 'CH3'),
                   ('O', 'OH'),
                   ('S', 'SH'),
                   ('N', 'NH')]


for ads1, ads2 in adsorbates_pair:
    print('Working on adsorbates: {} and {}'.format(ads1, ads2))

    X, y, sc_y, metadata = get_data(ads1, ads2)
    
    if not os.path.exists('gp_scaling_{}_{}'.format(ads1, ads2)):
        os.mkdir('gp_scaling_{}_{}'.format(ads1, ads2))

    min_gp_scaling = {'x' : sc_y, 'y' : y, 'metadata' : metadata}
    np.save('gp_scaling_{}_{}/min_scaling_data.npy'.format(ads1, ads2),
            min_gp_scaling)
    del min_gp_scaling

    sc = Linear_Scaling(sc_y, y, ads1, ads2, 'eV')
    sc.get_coeff()
    a = sc.plot_scaling()

    scaling_dict = sc.__dict__
    scaling_dict['metadata'] = metadata
    a.savefig('gp_scaling_{}_{}/scaling_plot.png'.format(ads1, ads2))
    np.save('gp_scaling_{}_{}/data.npy'.format(ads1, ads2), scaling_dict)

    print('Working on gp of adsorbate: {}'.format(ads2))

    data = preprocess_data(X, y)
    data.clean_data()
    X, y = data.get_data()

    kernel_recipe = {'ConstantKernel' : [{'RBF' : [1.0,
                                        {'length_scale' : 1.0}]},
                                        {'constant_value' : 1.0,
                                     'constant_value_bounds' : (3e-7, 3e7)}],
                     'WhiteKernel' : {'noise_level' : 0.1,
                                     'noise_level_bounds' : (1e-5, 1e5)}, 
                     'ConstantKernel' : [{'DotProduct' : [1.0,
                                        {'sigma_0' : 1.0}]},
                                     {'constant_value' : 1.0,
                                     'constant_value_bounds' : (3e-7, 3e7)}]}

    r_state = [10, 20, 42, 80, 150, 200, 300, 400]
    r_data = {rs : {} for rs in r_state}
    print('Shpae of y: {}'.format(y.shape))
    print('Shpae of y_sc: {}'.format(sc_y.shape))
    for rs in r_state:
        print('    Working on random state: {}'.format(rs))
        X_train, X_test, y_train, y_test, scaling_y_train, scaling_y_test, \
            metadata_train, metadata_test = \
            train_test_split(X, y, sc_y, metadata, \
            train_size=0.80, random_state=rs)

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
        r_data[rs] = MLGP.__dict__
        r_data[rs]['metadata_train'] = metadata_train
        r_data[rs]['metadata_test'] = metadata_test

    np.save('gp_scaling_{0}_{1}/{1}_rsdata.npy'.format(ads1, ads2), r_data)
    #LC = MLGP.plot_learning_curve()
    #LC.savefig('gp_scaling_{0}_{1}/{1}_learning_curve.png'.format(ads1, ads2))
    PP = MLGP.parity_plot(data='train')
    PP.savefig('gp_scaling_{0}_{1}/{1}_parity_plot_train.png'.format(ads1,
                                                                     ads2))
    PP = MLGP.parity_plot(data='test')
    PP.savefig('gp_scaling_{0}_{1}/{1}_parity_plot_test.png'.format(ads1,
                                                                    ads2))
    PP = MLGP.parity_plot(data='train', err_bar=True)
    string = 'gp_scaling_{0}_{1}/{1}_parity_plot_train_err_bar.png'
    PP.savefig(string.format(ads1, ads2))
    PP = MLGP.parity_plot(data='test', err_bar=True)
    string = 'gp_scaling_{0}_{1}/{1}_parity_plot_test_err_bar.png'
    PP.savefig(string.format(ads1, ads2))
