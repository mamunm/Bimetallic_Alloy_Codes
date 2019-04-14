#!/usr/bin/env python
# coding: utf-8

import ase
import os
from ase import Atoms

ten_metal = ['Pt', 'Pd', 'Rh', 'Ru', 'Re', 'Au', 'Bi', 'Co', 'Ag', 'Cu'] 
ads = ['C', 'O', 'H', 'N', 'S']
lttm = []

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/A3B'
list_A3B_top = os.listdir(DB_dir + '/Top')
list_A3B_top = sorted(list_A3B_top)

for i in list_A3B_top:
    m1, m2 = i.split('3')
    if all([m1 in ten_metal, m2 in ten_metal]):
        for a in ads:
            lttm += [i + '_' + m1 + '_' + a]
            lttm += [i + '_' + m2 + '_' + a]


DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/AB'
list_AB_top = os.listdir(DB_dir)
list_AB_top = sorted(list_AB_top)

for i in list_AB_top:
    m1, m2, _ = i.split('2')
    if all([m1 in ten_metal, m2 in ten_metal]):
        for a in ads:
            if m1 < m2:
                lttm += [m1 + '2' + m2 + '2' + '_' + m1 + '_' + a]
                lttm += [m1 + '2' + m2 + '2' + '_' + m2 + '_' + a]
            else:
                lttm += [m2 + '2' + m1 + '2' + '_' + m1 + '_' + a]
                lttm += [m2 + '2' + m1 + '2' + '_' + m2 + '_' + a]

DB_dir = '/scratch/users/mamunm/SHERLOCK/Binary_Alloy/Adsorbate-Surface/A3B/'
DB_dir += 'Calculations'
list_A3B_scratch = os.listdir(DB_dir)
list_A3B_scrtch = sorted(list_A3B_scratch)

for a in list_A3B_scratch:
    m1, m2 = a.split('3')
    if not all([m1 in ten_metal, m2 in ten_metal]):
        continue
    try:
        s_files = os.listdir(DB_dir + '/' + a + '/Calculations')
    except OSError:
        continue
    
    for b in s_files:
        fin_traj_path = (DB_dir + '/' + a + '/Calculations/' + b + '/final_'
                + b + '_' + m1 + '9' + m2 + '3.traj')
        if os.path.exists(fin_traj_path):
            k, _, _, l = b.split('_')
            lttm += [a + '_' + l + '_' + k]

DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/'
DB_dir += 'AB_from_SLAC_SUNCAT/Calculations'
list_AB_scratch = os.listdir(DB_dir)
list_AB_scrtch = sorted(list_AB_scratch)

for a in list_AB_scratch:
    m1, m2, _ = a.split('2')
    if not all([m1 in ten_metal, m2 in ten_metal]):
        continue
    try:
        s_files = os.listdir(DB_dir + '/' + a + '/Calculations')
    except OSError:
        continue
    
    for b in s_files:
        fin_traj_path = (DB_dir + '/' + a + '/Calculations/' + b + '/final_'
                + b + '_' + m1 + '6' + m2 + '6.traj')
        if os.path.exists(fin_traj_path):
            k, _, _, l = b.split('_')
            if m1 < m2:
                lttm += [m1 + '2' + m2 + '2' + '_' + l + '_' + k]
            else:
                lttm += [m2 + '2' + m1 + '2' + '_' + l + '_' + k]


with open('top_ten_metal.txt', 'w') as f:
    for i in lttm:
        f.write(i + '\n')
