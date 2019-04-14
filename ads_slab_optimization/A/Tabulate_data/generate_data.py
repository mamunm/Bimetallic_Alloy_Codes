#!/usr/bin/env python
# coding: utf-8

import os
from ase.db import connect
import numpy as np
from tabulate import tabulate
import time


def get_slab_energy(M):
    
    for i in db_slab.select():
        if i['formula'] == M + '12':
            
            return i['epot']


def get_ads_nrg_site(METAL, ADSORBATE, SITE):
    
    for i in db_ads_slab.select():
        
        if i['metal'].encode('utf-8') == METAL: 
            if i['adsorbate'].encode('utf-8') == ADSORBATE: 
                if i['site'].encode('utf-8') == SITE:
            
                    return (i['energy'], 
                            get_final_ads_site(METAL, ADSORBATE, SITE, 
                            i.toatoms()))
        
    else: 
            
        return 'N/A', 'N/A'

def compare_slab(A, B):
    
    pos0 = A.positions[-5:-1]
    pos1 = B.positions[-5:-1]
    lat = A.get_cell()[0][0]
    d_xy = np.empty(4)
    d_z = np.empty(4)
    
    for i in range(len(pos0)):
        
        d_xy[i] = np.linalg.norm(pos0[i][:2] - pos1[i][:2]) 
        d_z[i] = np.linalg.norm(pos0[i][2] - pos1[i][2])
    
    cond1 = np.all(np.absolute(d_xy) < 0.05 * lat)
    cond2 = np.all(np.absolute(d_z) < 0.15 * lat)
    
    if cond1 and cond2:
        
        return True
    
    else:
        
        return False

def is_subsurface(A):
    
    pos0 = A.positions[-5:-1][:, 2]
    pos1 = A.positions[-1][2]
    lat = A.get_cell()[0][0]
    
    if np.any((pos1 - pos0) < 0):
        
        return True
    
    else:
        
        return False

def is_ontop(A):

    pos0 = A.positions[-5:-1][:, :2]
    pos1 = A.positions[-1][:2]
    d_xy = np.empty(4)
    
    for i in range(len(pos0)):
        
        d_xy[i] = np.linalg.norm(pos0[i] - pos1) 
    
    if np.any(abs(d_xy) < 1e-2):
        
        return True
    
    else:
        
        return False

def is_bridge(A):

    pos0 = A.positions[-5][:2]
    pos1 = A.positions[-4][:2]
    pos2 = A.positions[-1][:2]
    
    bridge_dist = A.get_cell()[0][0] / 4
    
    d_xy = np.empty(2)
    d_xy[0] = np.linalg.norm(pos0 - pos2) 
    d_xy[1] = np.linalg.norm(pos1 - pos2) 
    
    if np.any(abs(d_xy - bridge_dist) < 8e-2):
        
        return True
    
    else:
        
        return False

def get_final_ads_site(M, A, S, IS):

    in_slab = IS

    for z in db_ini_ads_slab.select():
        
        if z['metal'].encode('utf-8') == M:                         
            if z['adsorbate'].encode('utf-8') == A:                   
                if z['site'].encode('utf-8') == S:
                    
                    fin_slab = z.toatoms()

    if compare_slab(in_slab, fin_slab):
        if not is_subsurface(fin_slab):
            if S == 'fcc' or S == 'hcp':
                
                return S
            
            if is_ontop(fin_slab):
                
                return 'ontop'
            
            else:
                if is_bridge(fin_slab):
                    
                    return 'bridge'
                
                else:
                    
                    return 'hollow_site'
        
        else:
            
            return 'Subsurface'
    
    else:
        
        return 'Reconfigured'



start_time = time.time()

del_file = ['A_ads_data.npy', 'Mol_data.npy']
for f in del_file:
    if os.path.exists(f):
        os.remove(f)


db_slab = connect('DB_files/Surface.db')
db_mol = connect('DB_files/molecules.db')
db_ads_slab = connect('DB_files/A_ads_slab.db')
db_ini_ads_slab = connect('DB_files/initial_A_ads_slab.db')

metals=['Al','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Y', 
        'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','La','Hf', 
        'Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']
metals = sorted(metals)

d = {}
ads_keys = ['*', 'H', 'C', 'O', 'N', 'S']
site_keys = ['ontop', 'bridge', 'fcc', 'hcp']

header = ['Site', 'E* [eV]', 'E_H* [eV]', 'Fin_Ads_Site', 'E_C* [eV]', 
          'Fin_Ads_Site', 'E_O* [eV]', 'Fin_Ads_Site', 'E_N* [eV]', 
          'Fin_Ads_Site', 'E_S* [eV]', 'Fin_Ads_Site']

for m in metals:

    d[m] = {k:{} for k in ads_keys}

    for ak in d[m]:

        if ak == '*':

            d[m][ak]['energy'] = get_slab_energy(m)
        
        else:
            
            for sk in site_keys:
                
                a, b = get_ads_nrg_site(m, ak, sk)
                d[m][ak][sk] = {}
                d[m][ak][sk]['energy'] = a
                d[m][ak][sk]['final_site'] = b


np.save('A_ads_data.npy', d)

Mol_data = {}

for dt in db_mol.select():

    Mol_data[dt['formula']] = dt['epot']  

np.save('Mol_data.npy', Mol_data)
print('\n')
print(tabulate([[a, b] for a, b in Mol_data.items()], 
     ['Molecule', 'Energy [eV]'], tablefmt='grid'))


print('\nTotal Elapsed time: {0} seconds'.format(time.time() - start_time))
