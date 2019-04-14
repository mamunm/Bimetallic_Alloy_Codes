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
    with open(f"../GP_Dot_Product/gp_scaling_{ads1}_{ads2}/out.txt", "r") as f:
        for line in f:
            if 'MAE_train mean:' in line:
                MAE_train = float(line.split()[-1])
            if 'MAE_train std:' in line:
                MAE_train_std = float(line.split()[-1])
            if 'MAE_test mean:' in line:
                MAE_test = float(line.split()[-1])
            if 'MAE_test std:' in line:
                MAE_test_std = float(line.split()[-1])
    data = np.load(f"../GP_Dot_Product/gp_scaling_{ads1}_{ads2}/{ads2}_rsdata.npy")[()][42]
    ax = plt.subplot(3, 2, i+1)
    #plt.grid(color='black', linestyle='-', linewidth=0.5)
    #ax.grid(color='grey', linewidth=0.5, alpha=0.01, linestyle='--')
    ax.grid(False)
    plt.xlabel(r"$\Delta E_{DFT} [eV]$")
    plt.ylabel(r"$\Delta E_{GP} [eV]$")
    xx = data['y_test']
    yy = data['y_pred_test']
    yy_un = data['y_pred_test_uncertainty']
    if ads2 == 'CH3':
        xx = [a for a in xx if not a > 4]
        yy = [b for a, b in zip(xx, yy) if not a > 4]
        yy_un = [b for a, b in zip(xx, yy_un) if not a > 4]
    plt.scatter(xx, yy,
                c=yy_un,
                alpha=0.7, cmap='viridis')
    plt.clim(0.2, 0.9)
    plt.plot(xx, xx, '-', color='black')
    if i not in [1, 2]:
        t_t = r"${}$".format(ads2)
    elif i == 1:
        t_t = r"$CH_2$"
    else:
        t_t = r"$CH_3$"

    plt.text(0.15, 0.8, t_t, size=10,
             ha="center", va="center",
             bbox=dict(facecolor='c',
                       edgecolor='black',
                       boxstyle='circle',
                       pad=1),
                       transform=ax.transAxes)

    text = r"$MAE_{0}$ = {1:0.2f} $\pm$ {2:0.3f}{3}"
    text += "$MAE_{4}$ = {5:0.2f} $\pm$ {6:0.3f}"
    text = text.format('{train}',
                       MAE_train,
                       MAE_train_std,
                       '\n',
                       '{test}',
                       MAE_test,
                       MAE_test_std)
    plt.text(0.62, 0.15, text,
             ha='center', va='center', fontsize=10,
             bbox=dict(facecolor='w',
             edgecolor='black',
             boxstyle='round',
             pad=1),
             transform=ax.transAxes)
    if i == 4:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.xaxis.set_ticks(np.arange(-2.5, 2, 1))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_ticks(np.arange(-2.5, 2, 1))
    if i == 2:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.xaxis.set_ticks(np.arange(-2, 5, 1))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_ticks(np.arange(-2, 5, 1))
        plt.plot([-2, 5], [-2, 5], '-', color='black')
    if i == 3:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.xaxis.set_ticks(np.arange(-3, 1.5, 1))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_ticks(np.arange(-3, 1.5, 1))

    plt.colorbar(label='uncertainty of prediction [eV]',
                 format='%0.2f')
plt.subplots_adjust(left=0.1,
                bottom=0.1,
                right=0.9,
                top=0.9,
                wspace=0.3,
                hspace=0.3)
#plt.tight_layout()
#plt.savefig('fig_1_v2.png', dpi=1500)
#plt.savefig('fig_1_v2_low_res.png')
plt.savefig('fig_3_test.eps', transparent=True, dpi=1500)
plt.show()

