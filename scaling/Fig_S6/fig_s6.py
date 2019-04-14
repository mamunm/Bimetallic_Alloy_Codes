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
from scipy.stats import linregress

def get_coeff(x, y):
    lr = linregress(x, y)
    return lr[0], lr[1]

adsorbates = ['CH', 'CH2', 'CH3', 'OH', 'SH', 'NH']

f = plt.figure(facecolor='white', figsize=(10, 12))
plt.style.use('seaborn-colorblind')
#plt.style.use('seaborn-whitegrid')
ca = {'CH': 'C', 'CH2': 'C', 'CH3': 'C',
      'OH': 'O', 'SH': 'S', 'NH': 'N'}

for i, ads in enumerate(adsorbates):
    print('Working on adsorbates: {}'.format(ads))
    data = np.load('../scaling_{}_{}/data.npy'.format(ca[ads], ads))[()]
    X = data['x']
    Y = data['y']

    m_s = []
    m_i = []
    s_s = []
    s_i = []
    n_data = []

    for j in range(2, len(X), 2):
        n_data += [j]
        temp = []
        for k in range(100):
            index = np.random.choice(np.arange(len(X)), j)
            x = X[index]
            y = Y[index]
            temp += [get_coeff(x, y)]
        a = np.array([k[0] for k in temp])
        b = np.array([k[1] for k in temp])
        m_s += [a.mean()]
        m_i += [b.mean()]
        s_i += [a.std()]
        s_s += [b.std()]

    ax = plt.subplot(3, 2, i+1)
    #plt.grid(color='black', linestyle='-', linewidth=0.5)
    #ax.grid(color='grey', linewidth=0.5, alpha=0.01, linestyle='--')
    nn = 50
    ax.grid(color='grey', alpha=0.3)
    plt.plot(n_data[:nn], m_s[:nn], '-', linewidth=4,
             c='darkblue', label="slope")
    plt.errorbar(n_data[:nn], m_s[:nn],
                 yerr=s_s[:nn], fmt='.k', color='blue')
    plt.plot(n_data[:nn], m_i[:nn], '-', linewidth=4,
             c='darkgreen', label="intercept")
    plt.errorbar(n_data[:nn], m_i[:nn],
                 yerr=s_i[:nn], fmt='.k', color='green')
    plt.xlabel('# data points')
    plt.ylabel('Coefficient value')
    if i == 0:
        plt.ylim([-2, 2])
    elif i in [1, 4]:
        plt.ylim([-0.5, 1.5])
    else:
        plt.ylim([-1.5, 1.5])
    if i == 0:
        plt.legend(loc='best')
    if i not in [1, 2]:
        t_t = r"${}$".format(ads)
    elif i == 1:
        t_t = r"$CH_2$"
    else:
        t_t = r"$CH_3$"
    
    plt.xlim(-4, 4)
    if i in [1, 4]:
        plt.text(0.8, 0.8, t_t, size=10,
                 ha="center", va="center",
                 bbox=dict(facecolor='c',
                           edgecolor='black',
                           boxstyle='circle',
                           pad=1),
                           transform=ax.transAxes)
    elif i in [2, 3, 5]:
        plt.text(0.8, 0.2, t_t, size=10,
                 ha="center", va="center",
                 bbox=dict(facecolor='c',
                           edgecolor='black',
                           boxstyle='circle',
                           pad=1),
                           transform=ax.transAxes)
    else:
        plt.text(0.8, 0.5, t_t, size=10,
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
plt.savefig('scaling_LC.eps', dpi=1500)
plt.show()
