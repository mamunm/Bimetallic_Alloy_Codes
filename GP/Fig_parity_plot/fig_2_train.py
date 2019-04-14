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

adsorbates = ['CH', 'CH2', 'CH3', 'OH', 'SH', 'NH']

plt.figure(facecolor='white', figsize=(10, 12))
plt.style.use('classic')
#plt.style.use('seaborn-whitegrid')
for i, ads in enumerate(adsorbates):
    print('Working on adsorbates: {}'.format(ads))
    data = np.load(f"../main_MS_GP/run_{ads}/{ads}_rsdata.npy")[()][42]
    ax = plt.subplot(3, 2, i+1)
    #plt.grid(color='black', linestyle='-', linewidth=0.5)
    #ax.grid(color='grey', linewidth=0.5, alpha=0.01, linestyle='--')
    ax.grid(False)
    plt.xlabel(r"$\Delta E_{DFT} [eV]$")
    plt.ylabel(r"$\Delta E_{GP} [eV]$")
    plt.scatter(data['y_train'], data['y_pred_train'],
                c=data['y_pred_train_uncertainty'],
                alpha=0.7, cmap='viridis')
    plt.colorbar(label='unc. of prediction [eV]')
    plt.clim(0.2, 0.9)
    plt.plot(data['y_train'], data['y_train'], '-', color='black')
    if i not in [1, 2]:
        t_t = r"${}$".format(ads)
    elif i == 1:
        t_t = r"$CH_2$"
    else:
        t_t = r"$CH_3$"

    plt.text(0.15, 0.8, t_t, size=10,
             ha="center", va="center",
             bbox=dict(facecolor='orange',
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
plt.savefig('fig_2_train.eps', transparent=True, dpi=1500)
plt.show()

