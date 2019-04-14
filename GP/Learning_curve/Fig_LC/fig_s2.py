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
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl

adsorbates = ['CH', 'CH2', 'CH3', 'OH', 'SH', 'NH']

f = plt.figure(facecolor='white', figsize=(10, 12))
plt.style.use('classic')
#plt.style.use('seaborn-whitegrid')
for i, ads in enumerate(adsorbates):
    print('Working on adsorbates: {}'.format(ads))
    with open(f"../LC_GP/run_{ads}/out.txt", "r") as f:
        MAE_train = []
        MAE_train_std = []
        MAE_test = []
        MAE_test_std = []
        for line in f:
            if 'MAE_train mean:' in line:
                MAE_train += [float(line.split()[-1])]
            if 'MAE_train std:' in line:
                MAE_train_std += [float(line.split()[-1])]
            if 'MAE_test mean:' in line:
                MAE_test += [float(line.split()[-1])]
            if 'MAE_test std:' in line:
                MAE_test_std += [float(line.split()[-1])]

    t_d = np.load(f'../LC_GP/run_{ads}/{ads}_rsdata.npy')[()]
    num_p = [len(t_d[i][list(t_d[i].keys())[0]]['y_train']) for i in t_d]
    del t_d

    ax = plt.subplot(3, 2, i+1)
    #plt.grid(color='black', linestyle='-', linewidth=0.5)
    #ax.grid(color='grey', linewidth=0.5, alpha=0.01, linestyle='--')
    ax.grid(False)
    plt.xlabel("# of training points")
    plt.ylabel("MAE [eV]")
    plt.plot(num_p[:-1], MAE_train[:-1], color='darkgreen', alpha=0.6,
             label='train data', marker='o')
    plt.errorbar(num_p[:-1], MAE_train[:-1], yerr=MAE_train_std[:-1],
                 fmt='.k', color='g', alpha=0.2)
    plt.plot(num_p[:-1], MAE_test[:-1], color='darkblue', alpha=0.6,
             label='test data', marker='o')
    plt.errorbar(num_p[:-1], MAE_test[:-1], yerr=MAE_test_std[:-1],
                 fmt='.k', color='b', alpha=0.2)
    plt.xticks(rotation=45)
    if i not in [1, 2]:
        t_t = r"${}$".format(ads)
    elif i == 1:
        t_t = r"$CH_2$"
    else:
        t_t = r"$CH_3$"

    if i == 0:
        plt.legend(loc='best')

    if i == 0:
        plt.text(0.8, 0.6, t_t, size=10,
                 ha="center", va="center",
                 bbox=dict(facecolor='c',
                           edgecolor='black',
                           boxstyle='circle',
                           pad=1),
                           transform=ax.transAxes)
    else:
        plt.text(0.8, 0.8, t_t, size=10,
                 ha="center", va="center",
                 bbox=dict(facecolor='c',
                           edgecolor='black',
                           boxstyle='circle',
                           pad=1),
                           transform=ax.transAxes)

plt.subplots_adjust(left=0.1,
                bottom=0.1,
                right=0.9,
                top=0.9,
                wspace=0.3,
                hspace=0.3)
#plt.tight_layout()
#plt.savefig('fig_1_v2.png', dpi=1500)
#plt.savefig('fig_1_v2_low_res.png')
plt.savefig('GP_LC.eps', transparent=True, dpi=1500)
plt.show()
