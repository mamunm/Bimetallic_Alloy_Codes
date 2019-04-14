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
adsorbates_pair = [('C', 'CH'),
                   ('C', 'CH2'),
                   ('C', 'CH3'),
                   ('O', 'OH'),
                   ('S', 'SH'),
                   ('N', 'NH')]

plt.figure(facecolor='white', figsize=(15, 8))
plt.style.use('classic')
for i, (ads1, ads2) in enumerate(adsorbates_pair):
    print('Working on adsorbates: {} and {}'.format(ads1, ads2))
    data = np.load(f'../scaling_{ads1}_{ads2}/data.npy')[()]
    ax = plt.subplot(2, 3, i+1)
    #plt.grid(color='black', linestyle='-', linewidth=0.5)
    ax.grid(False)
    title = "{0}* vs. {1}*".format(data['y_name'], data['x_name'])
    plt.title(title)
    plt.xlabel(r"$E_{0}\;[{1}]$".format('{' + data['x_name'] + '}',
                                            data['e_unit']))
    plt.ylabel(r"$E_{0}\;[{1}]$".format('{' + data['y_name'] + '}',
                                            data['e_unit']))
    plt.scatter(data['x'],
                data['y'],
                color='darkred',
                alpha=0.6,
                label='data points')
    plt.plot(data['x'],
             data['y_pred'],
             color='black',
             alpha=0.9,
             label="Scaling line")
    text = "slope = {0:0.2f}\nintercept = {1:0.2f}\n"
    text += "r^2 = {2:0.2f}\nMAE={3:0.2f}"
    text = text.format(data['slope'],
                       data['intercept'],
                       data['r_squared'],
                       data['MAE'])
    plt.text(0.75, 0.2, text,
             ha='center', va='center', fontsize=10,
             bbox=dict(facecolor='w',
             edgecolor='black',
             boxstyle='round',
             pad=1),
             transform=ax.transAxes)
    if i == 0:
        ax.legend(loc=2, prop={'size': 10})

plt.subplots_adjust(left=0.1,
                bottom=0.1,
                right=0.9,
                top=0.9,
                wspace=0.3,
                hspace=0.3)
#plt.tight_layout()
#plt.savefig('fig_1.png', dpi=1500)
plt.savefig('fig_1_low_res.png')
plt.show()

