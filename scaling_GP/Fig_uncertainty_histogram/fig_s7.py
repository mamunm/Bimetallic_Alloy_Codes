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
import matplotlib.mlab as mlab

adsorbates_pair = [('C', 'CH'),
                   ('C', 'CH2'),
                   ('C', 'CH3'),
                   ('O', 'OH'),
                   ('S', 'SH'),
                   ('N', 'NH')]

plt.figure(facecolor='white', figsize=(10, 12))
plt.style.use('classic')
#plt.style.use('seaborn-whitegrid')
for i, (ads1, ads2) in enumerate(adsorbates_pair):
    print('Working on adsorbates: {}'.format(ads2))
    data = np.load(f"../main_MS_sc_GP/gp_scaling_{ads1}_{ads2}/data.npy")[()]
    data = data['y'] - data['y_pred']
    ax = plt.subplot(3, 2, i+1)
    #plt.grid(color='black', linestyle='-', linewidth=0.5)
    #ax.grid(color='grey', linewidth=0.5, alpha=0.01, linestyle='--')
    ax.grid(False)
    n, bins, patches = plt.hist(data, 30, range=(-4, 4), normed=1,
                                facecolor='lightgreen', alpha=0.5)

    '''
    # add a 'best fit' line
    y = mlab.normpdf(bins, mu, sigma)
    plt.plot(bins, y, 'r--')
    '''
    plt.xlabel('Error [eV]')

    if i not in [1, 2]:
        t_t = r"${}$".format(ads2)
    elif i == 1:
        t_t = r"$CH_2$"
    else:
        t_t = r"$CH_3$"

    plt.text(0.8, 0.8, t_t, size=14, fontsize=18,
             ha="center", va="center",
             bbox=dict(facecolor='c',
                       edgecolor='black',
                       boxstyle='circle',
                       pad=1),
                       transform=ax.transAxes)
    plt.xlim(-4, 4)
plt.subplots_adjust(left=0.1,
                bottom=0.1,
                right=0.9,
                top=0.9,
                wspace=0.3,
                hspace=0.3)
#plt.tight_layout()
#plt.savefig('fig_1_v2.png', dpi=1500)
#plt.savefig('fig_1_v2_low_res.png')
plt.savefig('scaling_hist.eps', transparent=True, dpi=1500)
plt.show()

