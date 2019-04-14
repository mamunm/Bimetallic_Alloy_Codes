#!/usr/bin/env python                                                           
# coding: utf-8                                                                 

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
"""
Created on Nov 29, 2017
"""
__author__ = "Osman Mamun"
__version__ = "0.1"
__maintainer__ = "Osman Mamun"
__email__ = "mamunm@stanford.edu"
__date__ = "Nov, 2017"

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
import os
from ase.io import read, write
import glob
import math
import time

def dxy(AAA, BBB):
    return(math.sqrt((AAA[0] - BBB[0])**2 + (AAA[1] - BBB[1])**2))

def dz(AAA, BBB):
    return abs(AAA[2] - BBB[2])

def is_equivalent(AA, BB):
    con_val = []
    dxy_val = []
    dz_val = []
    AA_coord = AA.positions[-4:]
    BB_coord = BB.positions[-4:]
    for i in range(len(AA_coord)):
        #print(AA_coord[i], BB_coord[i])
        dxy_val.append(dxy(AA_coord[i], BB_coord[i]))
        dz_val.append(dz(AA_coord[i], BB_coord[i]))
        cond1 = dxy_val[-1] < 2.5
        cond2 = dz_val[-1] < 2.5
        #print(cond1, cond2)
        if cond1 and cond2:
            con_val.append(1)
        else:
            con_val.append(0)
    return all(aa == 1 for aa in con_val), dxy_val, dz_val

file_names = ['Initial_A3B_trajs.traj', 'Initial_AB_trajs.traj',
              'Final_A3B_trajs.traj', 'Final_AB_trajs.traj',
              'A3B_info.txt', 'AB_info.txt', 'A3B_debug.txt', 'AB_debug.txt',
              'slab_xy-dist.png', 'slab_z-dist.png']
for f in file_names:
    if os.path.exists(f):
        os.remove(f)

s_time = time.time()
DB_Dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Surface/'
list_a3b_files = os.listdir(DB_Dir + 'A3B')
list_ab_files = os.listdir(DB_Dir + 'AB')
bundled_a3b_init_traj = [read(DB_Dir + 'A3B/' + i + '/Initial_Slab.traj') 
                     for i in list_a3b_files]
write('Initial_A3B_trajs.traj', bundled_a3b_init_traj)
bundled_ab_init_traj = [read(DB_Dir + 'AB/' + i + '/Initial_Slab.traj') 
                     for i in list_ab_files]
write('Initial_AB_trajs.traj', bundled_ab_init_traj)
bundled_a3b_fin_traj = [read(glob.glob(DB_Dir + 'A3B/' + i + 
                             '/final_*.traj')[0]) for i in list_a3b_files]
write('Final_A3B_trajs.traj', bundled_a3b_fin_traj)
bundled_ab_fin_traj = [read(glob.glob(DB_Dir + 'AB/' + i + '/final_*.traj')[0]) 
                    for i in list_ab_files]
write('Final_AB_trajs.traj', bundled_ab_fin_traj)

Yes = 0
No = 0
DXY = []
DZ = []
with open('A3B_info.txt', 'w') as f:
    f.write('Information for A3B slab calculations:\n')
    for i, a3b in enumerate(list_a3b_files):
        A = bundled_a3b_init_traj[i]
        B = bundled_a3b_fin_traj[i]
        is_eq, dxyv, dzv = is_equivalent(A, B)
        with open('A3B_debug.txt', 'a') as ff:
            ff.write('{0}. Status of {1}: {2}\n'.format(i+1, a3b, is_eq))
            ff.write('Distance info:\n')
            for m in range(len(dzv)):
                ff.write('{0:2.4f}   {1:2.4f}\n'.format(dxyv[m], dzv[m]))
                DXY += [dxyv[m]]
                DZ += [dzv[m]]
        if is_eq:
            Yes += 1
            f.write('Status of {0:5s}: Successful\n'.format(a3b))
        else:
            No += 1
            f.write('Status of {0:5s}: Unsuccessful\n'.format(a3b))
print('Number of successful A3B jobs: %d' %Yes)
print('Number of unsuccessful A3B jobs: %d' %No)
Yes = 0
No = 0
with open('AB_info.txt', 'w') as f:
    f.write('Information for AB slab calculations:\n')
    for i, ab in enumerate(list_ab_files):
        A = bundled_ab_init_traj[i]
        B = bundled_ab_fin_traj[i]
        is_eq, dxyv, dzv = is_equivalent(A, B)
        with open('AB_debug.txt', 'a') as ff:
            ff.write('{0}. Status of {1}: {2}\n'.format(i+1, ab, is_eq))
            ff.write('Distance info:\n')
            for m in range(len(dzv)):
                ff.write('{0:2.4f}   {1:2.4f}\n'.format(dxyv[m], dzv[m]))
                DXY += [dxyv[m]]
                DZ += [dzv[m]]
        if is_eq:
            Yes += 1
            f.write('{0}. Status of {1}: Successful\n'.format(i, ab))
        else:
            No += 1
            f.write('{0}. Status of {1}: Unsuccessful\n'.format(i, ab))
print('Number of successful AB jobs: %d' %Yes)
print('Number of unsuccessful AB jobs: %d' %No)
print('\nTotal Elapsed Time: %2.4f seconds' %(time.time() - s_time))

import matplotlib.pyplot as plt
plt.figure(figsize=(10, 8))
plt.hist(DXY, bins=100, alpha=0.5, label='top layer xy displacement')
#plt.ylim(0, 350)
#plt.xlim(0, max(DXY))
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('slab_xy-dist.png')

plt.figure(figsize=(10, 8))
plt.hist(DZ, bins=100, alpha=0.5, label='top layer z displacement')
#plt.ylim(0, 150)
#plt.xlim(0, max(DZ))
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('slab_z-dist.png')
