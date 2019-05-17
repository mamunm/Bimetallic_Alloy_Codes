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
plt.rcParams['axes.axisbelow'] = True

adsorbates = ['CH', 'CH2', 'CH3', 'OH', 'NH', 'SH']
MAE_SC = [0.19, 0.26, 0.36, 0.25, 0.23, 0.26]
MAE_GP = [0.18, 0.27, 0.18, 0.19, 0.24, 0.21]
MAE_GP_std = [0.012, 0.022, 0.007, 0.013, 0.006, 0.006]
MAE_DL = [0.11, 0.19, 0.14, 0.14, 0.18, 0.12]
MAE_DL_std = [0.001, 0.014, 0.007, 0.008, 0.014, 0.009]


plt.figure(facecolor='white', figsize=(10, 6))
#plt.style.use('seaborn-white')

index = np.arange(len(adsorbates))
bar_width = 0.2
opacity = 1

plt.grid(c='grey', linestyle='--', alpha=0.2, axis='y')


rects1 = plt.bar(index, MAE_SC, bar_width,
                 alpha=0.4,
                 color='red',
                 label='Scaling')

rects2 = plt.bar(index+bar_width, MAE_GP, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='Gaussian Process',
                 yerr=MAE_GP_std,
                 ecolor='darkblue',
                 capsize=4,
                 error_kw={'elinewidth': 3})

rects3 = plt.bar(index+2*bar_width, MAE_DL, bar_width,
                 alpha=opacity,
                 color='green',
                 label='Delta Learning',
                 yerr=MAE_DL_std,
                 ecolor='darkgreen',
                 capsize=4,
                 error_kw={'elinewidth': 3})

plt.axhline(y=0.10, linewidth=3, color='grey', linestyle='-')
plt.xlabel('Adsorbates', fontsize=16)
plt.ylabel('Mean Absolute Error [eV]', fontsize=16)
r_ads = []
for ads in adsorbates:
    if ads not in ['CH2', 'CH3']:
        r_ads += [r"${}*$".format(ads)]
    elif ads == 'CH2':
        r_ads += [r"$CH_2*$"]
    else:
        r_ads += [r"$CH_3*$"]
plt.xticks(index + 1 * bar_width, r_ads,
           fontsize=16) # rotation=45)
plt.yticks(fontsize=16)
plt.legend(loc='upper left')
#plt.tight_layout()
plt.savefig('fig_MAE_SC_GP_DL.eps', dpi=1500)
plt.show()

