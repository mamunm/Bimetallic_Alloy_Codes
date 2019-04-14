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

adsorbates = ['CH', 'CH2', 'CH3', 'OH', 'SH', 'NH']
MAE_train = []
MAE_test = []
MAE_train_std = []
MAE_test_std = []

for ads in adsorbates:
    with open(f"../main_MS_GP/run_{ads}/out.txt", "r") as f:
        for line in f:
            if 'MAE_train mean:' in line:
                MAE_train += [float(line.split()[-1])]
            if 'MAE_train std:' in line:
                MAE_train_std += [float(line.split()[-1])]
            if 'MAE_test mean:' in line:
                MAE_test += [float(line.split()[-1])]
            if 'MAE_test std:' in line:
                MAE_test_std += [float(line.split()[-1])]

plt.figure(facecolor='white', figsize=(8, 6))
#plt.style.use('seaborn-white')

index = np.arange(len(adsorbates))
bar_width = 0.3
opacity = 1

plt.grid(c='grey', linestyle='--', alpha=0.2, axis='y')
rects1 = plt.bar(index, MAE_train, bar_width,
                 alpha=opacity,
                 color='darkred',
                 label='training data',
                 yerr=MAE_train_std,
                 ecolor='darkgreen',
                 capsize=4,
                 error_kw={'elinewidth': 3})

rects2 = plt.bar(index+bar_width, MAE_test, bar_width,
                 alpha=opacity,
                 color='darkblue',
                 label='testing data',
                 yerr=MAE_test_std,
                 ecolor='darkgreen',
                 capsize=4,
                 error_kw={'elinewidth': 3})

plt.axhline(y=0.15, linewidth=2, color='grey', linestyle='-')
plt.axhline(y=0.30, linewidth=2, color='grey', linestyle='-')
plt.xlabel('Adsorbates', fontsize=12)
plt.ylabel('Mean Absolute Error [eV]', fontsize=12)
r_ads = []
for ads in adsorbates:
    if ads not in ['CH2', 'CH3']:
        r_ads += [r"${}*$".format(ads)]
    elif ads == 'CH2':
        r_ads += [r"$CH_2*$"]
    else:
        r_ads += [r"$CH_3*$"]
plt.xticks(index + 0.5 * bar_width, r_ads,
           fontsize=12) # rotation=45)
plt.yticks(fontsize=12)
plt.legend(loc='best')
#plt.tight_layout()
plt.savefig('fig_3.eps', transparent=True, dpi=1500)
plt.show()

