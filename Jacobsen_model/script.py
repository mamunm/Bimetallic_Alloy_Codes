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

def get_refined_data(A, B):
    M = A.copy()
    N = B.copy()
    mask = np.invert(np.isnan(M) | np.isnan(N))
    return M[mask], N[mask]

adsorbates = ['CH', 'CH2', 'CH3', 'OH', 'SH', 'NH']

for ads in adsorbates:
    print('Working on adsorbates: {}'.format(ads))
    with FingerprintDB('../get_data_object/{}_fp.db'.format(ads)) as f:
        data = f.get_fingerprints(params=np.arange(107, 110))
    EA = data[:, 1]
    EAx = data[:, 2]
    E = data[:, 0]
    del data
    x, y = get_refined_data(EA, E)
    jac1 = LinearScaling(np.array(x), np.array(y), 'EA', 'E', 'eV')
    x, y = get_refined_data(EAx, E)
    jac2 = LinearScaling(np.array(x), np.array(y), 'EAx', 'E', 'eV')
    jac1.get_coeff()
    jac2.get_coeff()
    a1 = jac1.plot_scaling()
    a2 = jac2.plot_scaling()

    if not os.path.exists('Jacobsen_{}'.format(ads)):
        os.mkdir('Jacobsen_{}'.format(ads))

    jac1_dict = jac1.__dict__
    jac2_dict = jac2.__dict__
    a1.savefig('Jacobsen_{}/jacobsen_1_plot.png'.format(ads))
    a2.savefig('Jacobsen_{}/jacobsen_2_plot.png'.format(ads))
    np.save('Jacobsen_{}/jac1_data.npy'.format(ads), jac1_dict)
    np.save('Jacobsen_{}/jac2_data.npy'.format(ads), jac2_dict)

