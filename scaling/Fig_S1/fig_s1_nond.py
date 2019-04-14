#/usr/bin/env python
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

def get_min_data(D):
    slab = list(set([i[1] for i in D['metadata']]))
    nond = ['Sc', 'Y', 'Al', 'Ga', 'Sn', 'Bi',
            'Pb', 'In', 'Ti', 'Zn', 'Cd', 'Hg',
            'Cu', 'Ag', 'Au']
    slab = [i for i in slab if not any([j in i for j in nond])]
    print(slab)
    x = []
    y = []
    sl = []
    for s in slab:
        x_temp = [i for i, j in zip(D['x'], D['metadata']) if j[1] == s]
        y_temp = [i for i, j in zip(D['y'], D['metadata']) if j[1] == s]
        if x_temp and y_temp:
            x += [min(x_temp)]
            y += [min(y_temp)]
            sl += [s]
    return x, y, sl

def encircle(x,y, ax=None, **kw):
    if not ax: ax=plt.gca()
    p = np.c_[x,y]
    mean = np.mean(p, axis=0)
    d = p-mean
    r = np.max(np.sqrt(d[:,0]**2+d[:,1]**2 ))
    circ = plt.Circle(mean, radius=1.05*r,**kw)
    ax.add_patch(circ)

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
    print('Working on adsorbates: {} and {}'.format(ads1, ads2))
    data = np.load(f'../scaling_{ads1}_{ads2}/data.npy')[()]
    ax = plt.subplot(3, 2, i+1)
    #plt.grid(color='black', linestyle='-', linewidth=0.5)
    ax.grid(True)
    title = "{0}* vs. {1}*".format(data['y_name'], data['x_name'])
    plt.title(title)
    plt.xlabel(r"$E_{0}\;[{1}]$".format('{' + data['x_name'] + '}',
                                            data['e_unit']))
    plt.ylabel(r"$E_{0}\;[{1}]$".format('{' + data['y_name'] + '}',
                                            data['e_unit']))
    x_min, y_min, metal_type = get_min_data(data)
    sc = Linear_Scaling(np.array(x_min), np.array(y_min), ads1, ads2, 'eV')
    sc.get_coeff()
    if i < 3:
        plt.ylim([-1, 7])
    plt.scatter(sc.x, sc.y,
                facecolors='none', edgecolors='darkred',
                alpha=0.4)
    plt.plot(sc.x, sc.y_pred,
             color='black',
             alpha=0.9)
             #label="Scaling line")
    p_slope = {'CH': 0.75, 'CH2': 0.5, 'CH3': 0.25,
               'OH': 0.5, 'NH': 0.67, 'SH': 0.5}
    text = r"$\gamma$ = {0:0.2f}{1}$\xi$ = {2:0.2f}{3}"
    text += r"$r^2$ = {4:0.2f}{5}MAE={6:0.2f}{7}$\gamma_{8}$={9:0.2f}"
    text = text.format(sc.slope,
                       '\n',
                       sc.intercept,
                       '\n',
                       sc.r_squared,
                       '\n',
                       sc.MAE,
                       '\n',
                       '{pred}',
                       p_slope[ads2])
    if i not in [1, 2]:
        plt.text(0.75, 0.3, text,
                 ha='center', va='center', fontsize=10,
                 bbox=dict(facecolor='w',
                 edgecolor='black',
                 boxstyle='round',
                 pad=1),
                 transform=ax.transAxes)
    else:
        plt.text(0.2, 0.7, text,
                 ha='center', va='center', fontsize=10,
                 bbox=dict(facecolor='w',
                 edgecolor='black',
                 boxstyle='round',
                 pad=1),
                 transform=ax.transAxes)
    if i == 0:
        ax.legend(loc=2,
                  scatterpoints=1,
                  prop={'size': 8})
    '''
    if i in [1, 2]:
        x2 = [y for y, z in zip(xhollow, yhollow) if z < 0.2]
        y2 = [z for y, z in zip(xhollow, yhollow) if z < 0.2]
        encircle(x2, y2, ec="orange", fc="none", linewidth=2)
    if i in [4, 5]:
        x2 = [y for y, z in zip(xhollow, yhollow) if z < -2]
        y2 = [z for y, z in zip(xhollow, yhollow) if z < -2]
        encircle(x2, y2, ec="orange", fc="none", linewidth=2)
    if i in [3]:
        x2 = [y for y, z in zip(xtop, ytop) if all([-2 < y < 0,
                                                          z < -1])]
        y2 = [z for y, z in zip(xtop, ytop) if all([-2 < y < 0,
                                                          z < -1])]
        encircle(x2, y2, ec="orange", fc="none", linewidth=2)
    '''
plt.subplots_adjust(left=0.1,
                bottom=0.1,
                right=0.9,
                top=0.9,
                wspace=0.3,
                hspace=0.3)
#plt.tight_layout()
plt.savefig('fig_s1_nond.eps', transparent=True, dpi=1500)
plt.show()

