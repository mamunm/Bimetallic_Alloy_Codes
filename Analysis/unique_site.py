#!/usr/bin/env python                                                         
# coding: utf-8                                                               

from __future__ import division, unicode_literals                             
                                                                              
"""                                                                           
Created on Oct 18, 2017 
Identify all the unique sites on a flat surface.
"""                                                                           
__author__ = "Osman Mamun"                                                    
__version__ = "0.1"                                                           
__maintainer__ = "Osman Mamun"                                                
__email__ = "mamunm@stanford.edu"                                             
__date__ = "Oct 18, 2017"                                                     
                                                                              
import os                                                                     
import sys                                                                    
from ase.io import read                                                       
from ase.build import make_supercell                                          
import numpy as np                                                            
from copy import deepcopy                                                     
from spglib import get_symmetry                                              
from ase.build import add_adsorbate
from ase.io import write

def distance(A, B):                                                            
    '''                                                                        
        This function returns the distance between two atoms                   
    '''                                                                        
    return np.sqrt((A[2] - B[2])**2 + (A[3] - B[3])**2 + (A[4] - B[4])**2)     

def top_sites(slab):                                                          
    '''                                                                       
        This function returns all the identified top sites                    
    '''                                                                       
    slab_coord = slab.positions                                               
    ts = np.array(slab_coord[ abs(slab_coord[:,2] -                           
                  max(slab_coord[:,2])) < 1e-1 ],                             
                  dtype = 'double', order = 'C')                              
    symbol = np.array(slab.get_chemical_symbols()[-len(ts):],                 
                                dtype = 'str', order = 'C')                   
    atomic_number = np.array(slab.get_atomic_numbers()[-len(ts):],            
                                dtype = 'str', order = 'C')                   
    x = np.array(ts[:,0], dtype = 'double', order = 'C')                      
    y = np.array(ts[:,1], dtype = 'double', order = 'C')                      
    z = np.array(ts[:,2], dtype = 'double', order = 'C')                      
    return [symbol, atomic_number, x, y, z]                                  

def bridge_sites(ts, supercell):                                           
    '''                                                                    
        This function returns all the identified top sites                 
    '''                                                                    
    symbol = []                                                            
    atomic_number = []                                                     
    x = []                                                                 
    y = []                                                                 
    z = []                                                                 
    bridge_distance = 0.5 * slab.cell[0][0]                                
    for i in range(len(ts[0])):                                            
        tsA = [ts[k][i] for k in range(len(ts))]                           
        for j in range(len(supercell[0])):                                 
            scB = [supercell[l][j] for l in range(len(supercell))]         
            d = distance(tsA, scB)                                         
            c1 = d < 1.1 * bridge_distance                                 
            c2 = abs(ts[4][i] - supercell[4][j]) < 0.1                     
            if c1 and d > 0.1 and c2:                                      
                sym = ts[0][i] + supercell[0][j]                           
                an = str(int(float(ts[1][i]) * 0.5 +                       
                             float(supercell[1][j]) * 0.5))                
                x_t = ts[2][i] * 0.5 + supercell[2][j] * 0.5               
                y_t = ts[3][i] * 0.5 + supercell[3][j] * 0.5               
                z_t = ts[4][i] * 0.5 + supercell[4][j] * 0.5               
                if i == 0:                                                 
                    symbol.append(sym)                                     
                    atomic_number.append(an)                               
                    x.append(x_t)                                          
                    y.append(y_t)                                          
                    z.append(z_t)                                          
                else:                                                      
                    append = []                                            
                    for ll in range(len(x)):                               
                        cond1 = abs(x_t - x[ll]) < 1e-3                    
                        cond2 = abs(y_t - y[ll]) < 1e-3                    
                        if cond1 and cond2:                                
                            append.append(False)                           
                        else:                                              
                            append.append(True)                            
                    if not False in append:                                
                        symbol.append(sym)                                 
                        atomic_number.append(an)                           
                        x.append(x_t)                                      
                        y.append(y_t)                                      
                        z.append(z_t)                                      
    symbol = np.array(symbol, dtype = 'str', order = 'C')                  
    atomic_number = np.array(atomic_number, dtype = 'str', order = 'C')    
    x = np.array(x, dtype = 'double', order = 'C')                         
    y = np.array(y, dtype = 'double', order = 'C')                         
    z = np.array(z, dtype = 'double', order = 'C')                         
    return [symbol, atomic_number, x, y, z]                                

def get_supercell(slab):                                                    
    '''                                                                     
        This function returns a 3 by 3 supercell of the original cell       
    '''                                                                     
    P = [[3, 0, 0], [0, 3, 0], [0, 0, 1]]                                   
    supercell = make_supercell(slab, P)                                     
    slab_cell = slab.cell                                                   
    x = np.array(supercell.positions[:,0] - (slab_cell[0][0] +              
                 slab_cell[1][0]), dtype = 'double', order = 'C')           
    y = np.array(supercell.positions[:,1] - (slab_cell[1][1]),              
                 dtype = 'double', order = 'C')                             
    z = np.array(supercell.positions[:,2],                                  
                 dtype = 'double', order = 'C')                             
    symbol = np.array(supercell.get_chemical_symbols(),                     
                      dtype = 'str', order = 'C')                           
    atomic_number = np.array(supercell.get_atomic_numbers(),                
                             dtype = 'str', order = 'C')                    
    return [symbol, atomic_number, x, y, z]                                 

def hollow_sites(ts, supercell):
    '''
        This function returns all the identified top sites
    '''
    symbol = []
    atomic_number = []
    x = []
    y = []
    z = []
    bridge_distance = 0.5 * slab.cell[0][0]
    for i in range(len(ts[0])):
        tsA = [ts[k][i] for k in range(len(ts))]
        for j, scBB in enumerate(supercell[0]):
            scB = [supercell[l][j] for l in range(len(supercell))]
            d1 = distance(tsA, scB)
            c1 = d1 < 1.1 * bridge_distance
            c2 = d1 > 0.1
            c3 = abs(ts[4][i] - supercell[4][j]) < 0.1
            #c4 = abs(ts[3][i] - supercell[3][j]) < 0.1
            if c1 and c2 and c3:
                for k, scCC in enumerate(supercell[0]):
                    if j != k:
                        scC = [supercell[l][k] for l in range(len(supercell))]
                        d2 = distance(tsA, scC)
                        d3 = distance(scB, scC)
                        c5 = d2 < 1.1 * bridge_distance
                        c6 = d3 < 1.1 * bridge_distance
                        c7 = d2 > 0.1
                        c8 = abs(ts[4][i] -
                                 supercell[4][k]) < 0.1
                        #c9 = abs(ts[3][i] -
                        #         supercell[3][k]) > 0.1
                        if c5 and c6 and c7 and c8:
                            sym = ts[0][i] + supercell[0][j] + supercell[0][k]
                            an = '%0.3f' %((float(ts[1][i]) +
                                        float(supercell[1][j]) +
                                        float(supercell[1][k])) / 3.0)
                            x_t = (ts[2][i] + supercell[2][j]
                                  + supercell[2][k]) / 3.0
                            y_t = (ts[3][i] + supercell[3][j]
                                  + supercell[3][k]) / 3.0
                            z_t = (ts[4][i] + supercell[4][j]
                                  + supercell[4][k]) / 3.0
                            if sym == []:
                                symbol.append(sym)
                                atomic_number.append(an)
                                x.append(x_t)
                                y.append(y_t)
                                z.append(z_t)
                            else:
                                append = []
                                for ll in range(len(x)):
                                    cond1 = abs(x_t - x[ll]) < 1e-3
                                    cond2 = abs(y_t - y[ll]) < 1e-3
                                    if cond1 and cond2:
                                        append.append(False)
                                    else:
                                        append.append(True)
                                if not False in append:
                                    symbol.append(sym)
                                    atomic_number.append(an)
                                    x.append(x_t)
                                    y.append(y_t)
                                    z.append(z_t)
    symbol = np.array(symbol, dtype = 'str', order = 'C')
    atomic_number = np.array(atomic_number, dtype = 'str', order = 'C')
    x = np.array(x, dtype = 'double', order = 'C')
    y = np.array(y, dtype = 'double', order = 'C')
    z = np.array(z, dtype = 'double', order = 'C')
    return [symbol, atomic_number, x, y, z]

def UBHS(A, R, T, CD):
    frac = np.empty([len(A[0]), 3])
    for i in range(len(A[0])):
        for j in range(3):
            frac[i][j] = A[j+2][i]
    fractional = np.linalg.solve(cell_dim.T, frac.T).T
    for j in range(len(fractional)):
        fractional[j, :] %= 1.0
    fractional[fractional > 0.99999] = 0.0
    #print("fractional:\n", fractional)
    ind = []
    temp = np.empty([len(R),3])
    for i in range(len(fractional)):
        for j, rot in enumerate(R):
            temp[j, :] = np.matmul(fractional[i], rot.T) + T[j]
        for j in range(len(temp)):
            temp[j, :] %= 1.0
        temp[temp > 0.99999] = 0.0
        #print("RT for %d:\n" %i, temp)
        ind.append(unique_index(fractional, temp))
    return ind

def unique_index(F, T):
    ret_ind = []
    for i in range(len(F)):
        for j in range(len(T)):
            c1 = abs(F[i][0] - T[j][0]) < 1e-2
            c2 = abs(F[i][1] - T[j][1]) < 1e-2
            c3 = abs(F[i][2] - T[j][2]) < 1e-2
            if c1 and c2 and c3 and not i in ret_ind:
                ret_ind.append(i)
    return ret_ind

#Main
input_traj_file_name = 'Initial_Slab.traj'
ads_slab = read(input_traj_file_name)
slab = ads_slab[:-1]
slab_an = slab.get_chemical_symbols()
slab_coord = slab.positions
slab_frac_coord = slab.get_scaled_positions()
cell_dim = slab.cell
supercell = get_supercell(slab)
topsite = top_sites(slab)
bridgesite = bridge_sites(topsite, supercell)                                  
hollowsite = hollow_sites(topsite, supercell)

print('\nSpecification of the given crystal structure:\n')
print('Lattice parameters:\n')
print(cell_dim, '\n')
print('Absolute coordinates:\n')
for i in range(len(slab_an)):
    print('{0:2s}{1:-8.2f}{2:-8.2f}{3:-8.2f}'.format(slab_an[i], 
          slab_coord[i][0], slab_coord[i][1], slab_coord[i][2]))
print('\nFractional coordinates:\n')
for i in range(len(slab_an)):
    print('{0:2s}{1:-8.2f}{2:-8.2f}{3:-8.2f}'.format(slab_an[i], 
          slab_frac_coord[i][0], slab_frac_coord[i][1], slab_frac_coord[i][2]))

print('\nAll the possible adsorption sites on the top layer::\n')
print('\nTop sites:\n')
print('   %-10s %-20s %-20s %-12s %-12s %-12s' %('Index', 'chemical_symbol',
      'atomic_number', 'x', 'y', 'z'))
for i in range(len(topsite[0])):
    print(' %-15s %-20s %-15s %-12.2f %-12.2f %-12.2f' %('top_site_' + 
          str(i), topsite[0][i], topsite[1][i], topsite[2][i], 
          topsite[3][i], topsite[4][i]))
print('\n')
bridgesite = bridge_sites(topsite, supercell)                                  
print('\nBridge sites:\n')                                                     
print('   %-10s %-20s %-20s %-12s %-12s %-12s' %('Index', 'chemical_symbol',   
      'atomic_number', 'x', 'y', 'z'))                                         
for i in range(len(bridgesite[0])):                                            
    print(' %-15s %-20s %-15s %-12.2f %-12.2f %-12.2f' %('bridge_site_' + 
          str(i), bridgesite[0][i], bridgesite[1][i], bridgesite[2][i], 
          bridgesite[3][i], bridgesite[4][i]))                       
print('\n')
print('\nHollow sites:\n')                                                     
print('   %-10s %-20s %-20s %-12s %-12s %-12s' %('Index', 'chemical_symbol',   
      'atomic_number', 'x', 'y', 'z'))                                         
for i in range(len(hollowsite[0])):                                            
    print(' %-15s %-20s %-15s %-12.2f %-12.2f %-12.2f' %('hollow_site_'        
          + str(i), hollowsite[0][i], hollowsite[1][i],                    
          hollowsite[2][i], hollowsite[3][i], hollowsite[4][i]))
print('\n')
print('\nIdentifying unique adsorption site based on crystal symmetry' 
      ' operation:\n')
#Symmetry related top sites
a = get_symmetry(slab, symprec=1e-5)
equivalent_atoms = a["equivalent_atoms"]
#rotations = a["rotations"]
rotations = [rot for rot in a["rotations"] if rot[2][2] == 1]
#translations = a["translations"]
translations = [trans for i, trans in enumerate(a["translations"]) 
                if a["rotations"][i][2][2] == 1 ]
print('Output from the spglib symmetry operation on this crystal surface:\n')
print("Equivalent atoms: \n", equivalent_atoms, "\n")
print("Rotations: \n")
for i in range(len(rotations)):
    print("Rotation matrix %d :" %i)
    print(rotations[i], "\n")
print("Translations: \n")
for i in range(len(translations)):
    print("Translation vector %d :" %i)
    print(translations[i], "\n")
print('\n')

_ , index, ind_arr = np.unique(a["equivalent_atoms"][-4:], return_index=True, 
                               return_inverse=True)
mk = [[] for _ in range(len(index))]
for i in range(len(index)):
    for j in range(len(ind_arr)):
        if i == ind_arr[j]:
            mk[i].append(j)
top_slab = read(input_traj_file_name)  
print('Unique top sites:\n')
print('   %-10s %-20s %-20s %-12s %-12s %-12s %-15s' %('Index', 
      'Chemical_symbol', 'Atomic_number', 'x', 'y', 'z', 
      'Equivalent site index'))
for i, j in enumerate(index):
    add_adsorbate(top_slab, 'O', 1.0, (topsite[2][j], topsite[3][j])) 
    print(' %-15s %-20s %-15s %-12.2f %-12.2f %-12.2f %-15s' 
          %('unique_top_site_' + str(i), topsite[0][j], 
          topsite[1][j], topsite[2][j], 
          topsite[3][j], topsite[4][j], mk[i]))
print('\n')
write('uts.traj', top_slab)
for i in range(len(topsite[0])):
    add_adsorbate(top_slab, 'O', 1.0, (topsite[2][i], topsite[3][i]))
write('ts.traj', top_slab)

index_bs = UBHS(bridgesite, rotations, translations, cell_dim)
mk_bs = []
for x in index_bs:
    if x not in mk_bs:
        mk_bs.append(x)
mk_bs_sets = [set(l) for l in mk_bs]
mk_bs = [l for l, s in zip(mk_bs, mk_bs_sets) if not any(s < other 
         for other in mk_bs_sets)]
unique_index_bs = []
for i in range(len(mk_bs)):
    unique_index_bs.append(mk_bs[i][0])
bridge_slab = read(input_traj_file_name)
print('Unique bridge sites:\n')
print('   %-20s %-20s %-20s %-12s %-12s %-12s %-15s' %('Index',
      'Chemical_symbol', 'Atomic_number', 'x', 'y', 'z',
      'Equivalent site index'))
for i, j in enumerate(unique_index_bs):
    add_adsorbate(bridge_slab, 'O', 1.0, (bridgesite[2][j], bridgesite[3][j]))
    print(' %-25s %-20s %-15s %-12.2f %-12.2f %-12.2f %-15s' 
          %('unique_bridge_site_' + str(i), bridgesite[0][j], 
          bridgesite[1][j], bridgesite[2][j],
          bridgesite[3][j], bridgesite[4][j], mk_bs[i]))
print('\n')
write('ubs.traj', bridge_slab)
for i in range(len(bridgesite[0])):
    add_adsorbate(bridge_slab, 'O', 1.0, (bridgesite[2][i], bridgesite[3][i]))
write('bs.traj', bridge_slab)


index_hs = UBHS(hollowsite, rotations, translations, cell_dim)
mk_hs = []
for x in index_hs:
    if x not in mk_hs:
        mk_hs.append(x)
mk_hs_sets = [set(l) for l in mk_hs]
mk_hs = [l for l, s in zip(mk_hs, mk_hs_sets) if not any(s < other 
         for other in mk_hs_sets)]
unique_index_hs = []
for i in range(len(mk_hs)):
    unique_index_hs.append(mk_hs[i][0])
hollow_slab = read(input_traj_file_name)
print('Unique hollow sites:\n')
print('   %-20s %-20s %-20s %-12s %-12s %-12s %-15s' %('Index',
      'Chemical_symbol', 'Atomic_number', 'x', 'y', 'z',
      'Equivalent site index'))
for i, j in enumerate(unique_index_hs):
    add_adsorbate(hollow_slab, 'O', 1.0, (hollowsite[2][j], hollowsite[3][j]))
    print(' %-25s %-20s %-15s %-12.2f %-12.2f %-12.2f %-15s' 
          %('unique_hollow_site_' + str(i), hollowsite[0][j], 
          hollowsite[1][j], hollowsite[2][j],
          hollowsite[3][j], hollowsite[4][j], mk_hs[i]))
print('\n')
write('uhs.traj', hollow_slab)
for i in range(len(hollowsite[0])):
    add_adsorbate(hollow_slab, 'O', 1.0, (hollowsite[2][i], hollowsite[3][i]))
write('hs.traj', hollow_slab)


