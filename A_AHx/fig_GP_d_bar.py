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
MAE_d0 = []
MAE_d0_std = []
MAE_d1 = []
MAE_d1_std = []
MAE_d2 = []
MAE_d2_std = []

for ads in adsorbates:
    for i in range(3):
        with open(f"GP_d{i}/run_{ads}/out.txt", "r") as f:
            for line in f:
                if 'MAE_test mean:' in line:
                    globals()[f'MAE_d{i}'] += [float(line.split()[-1])]
                if 'MAE_test std:' in line:
                    globals()[f'MAE_d{i}_std'] += [float(line.split()[-1])]



plt.figure(facecolor='white', figsize=(10, 6))
#plt.style.use('seaborn-white')

index = np.arange(len(adsorbates))
bar_width = 0.2
opacity = 1

plt.grid(c='grey', linestyle='--', alpha=0.2, axis='y')


rects1 = plt.bar(index, MAE_d0, bar_width,
                 alpha=0.4,
                 color='red',
                 label=r'$1^{st}\:NN$',
                 yerr=MAE_d0_std,
                 ecolor='darkred',
                 capsize=4,
                 error_kw={'elinewidth': 3})

rects2 = plt.bar(index+bar_width, MAE_d1, bar_width,
                 alpha=opacity,
                 color='blue',
                 label=r'$2^{nd}\:NN$',
                 yerr=MAE_d1_std,
                 ecolor='darkblue',
                 capsize=4,
                 error_kw={'elinewidth': 3})

rects3 = plt.bar(index+2*bar_width, MAE_d2, bar_width,
                 alpha=opacity,
                 color='green',
                 label=r'$3^{rd}\:NN$',
                 yerr=MAE_d2_std,
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
plt.savefig('fig_MAE_d_bar.eps', dpi=1500)
plt.show()

