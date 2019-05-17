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
    with open(f"run_{ads}/out.txt", "r") as f:
        for line in f:
            if 'MAE_train mean:' in line:
                MAE_train_GP += [float(line.split()[-1])]
            if 'MAE_train std:' in line:
                MAE_train_std_GP += [float(line.split()[-1])]
            if 'MAE_test mean:' in line:
                MAE_test_GP += [float(line.split()[-1])]
            if 'MAE_test std:' in line:
                MAE_test_std_GP += [float(line.split()[-1])]

for ads in adsorbates:
    ads_pref = ads[0]
    with open(f"../scaling_GP/main_MS_sc_GP/gp_scaling_{ads_pref}_{ads}/out.txt",
            "r") as f:
        for line in f:
            if 'MAE_train mean:' in line:
                MAE_train_SCGP += [float(line.split()[-1])]
            if 'MAE_train std:' in line:
                MAE_train_std_SCGP += [float(line.split()[-1])]
            if 'MAE_test mean:' in line:
                MAE_test_SCGP += [float(line.split()[-1])]
            if 'MAE_test std:' in line:
                MAE_test_std_SCGP += [float(line.split()[-1])]

MAE_scaling = [0.19, 0.26, 0.36, 0.25, 0.26, 0.23]

plt.figure(facecolor='white', figsize=(10, 6))
#plt.style.use('seaborn-white')

index = np.arange(len(adsorbates))
bar_width = 0.2
opacity = 1

plt.grid(c='grey', linestyle='--', alpha=0.2, axis='y')

rects1 = plt.bar(index, MAE_scaling, bar_width,
                 alpha=opacity,
                 color='red',
                 label='scaling relation')

rects2 = plt.bar(index+bar_width, MAE_test_GP, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='machine learning',
                 yerr=MAE_test_std_GP,
                 ecolor='darkblue',
                 capsize=4,
                 error_kw={'elinewidth': 3})

rects3 = plt.bar(index+2*bar_width, MAE_test_SCGP, bar_width,
                 alpha=opacity,
                 color='green',
                 label='delta learning',
                 yerr=MAE_test_std_SCGP,
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
plt.savefig('fig_MAE.eps', transparent=True, dpi=1500)
plt.show()

