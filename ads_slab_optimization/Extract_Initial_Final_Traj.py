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
import time
import math

def dxy(AAA, BBB):                                                              
    return(math.sqrt((AAA[0] - BBB[0])**2 + (AAA[1] - BBB[1])**2))              
                                                                                
def dz(AAA, BBB):                                                               
    return abs(AAA[2] - BBB[2])                                                 
                                                                                
def is_equivalent(AA, BB):                                                      
    con_val = []                                                                
    AA_coord = AA.positions[-5:-1]                                            
    BB_coord = BB.positions[-5:-1]                                             
    for i in range(len(AA_coord)):                                              
        #print(AA_coord[i], BB_coord[i])                                        
        cond1 = dxy(AA_coord[i], BB_coord[i]) < 1e-1                            
        cond2 = dz(AA_coord[i], BB_coord[i]) < 9e-1                             
        #print(cond1, cond2)                                                    
        if cond1 and cond2:                                                     
            con_val.append(1)                                                   
        else:                                                                   
            con_val.append(0)                                                   
    return all(aa == 1 for aa in con_val)  

s_time = time.time() 
DB_Dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/'
list_a3b_files = os.listdir(DB_Dir + 'A3B')
bundled_a3b_init_traj = [read('A3B/' + i + '/Calculations/' + j + 
                '/Initial_slab_with_adsorbate.traj') for i in list_a3b_files
                for j in os.listdir(DB_Dir + 'A3B/' + i + '/Calculations/')]
write('Initial_trajs.traj', bundled_a3b_init_traj)

bundled_a3b_fin_traj = [read(glob.glob('A3B/' + i + '/Calculations/' + j + 
                '/final_*.traj')[0]) for i in list_a3b_files
                for j in os.listdir(DB_Dir + 'A3B/' + i + '/Calculations/')]
write('Final_trajs.traj', bundled_a3b_fin_traj)

Yes = 0                                                                         
No = 0                                                                          

with open('A3B_info.txt', 'w') as f:  
    f.write('Information for A3B ads lab calculations:\n')                  
    for i, a3b in enumerate(list_a3b_files): 
        list_temp = os.listdir(DB_Dir + 'A3B/' + a3b + '/Calculations/')
        for j, ads_a3b in enumerate(list_temp):
            A = bundled_a3b_init_traj[i * 10 + j]                         
            B = bundled_a3b_fin_traj[i * 10 + j]                              
            if is_equivalent(A, B):                                         
                Yes += 1                                        
                p_st = 'Status of {0} {1} job: Surface not reconfigured\n'
                f.write(p_st.format(a3b, ads_a3b))
            else:                                                           
                No += 1 
                p_st = 'Status of {0} {1} job: Surface reconfigured\n'
                f.write(p_st.format(a3b, ads_a3b))

print('Number of successful A3B adsorption jobs: %d' %Yes)                      
print('Number of unsuccessful A3B adsorption jobs: %d' %No)                    


print('\nTotal Elapsed Time: %2.4f seconds' %(time.time() - s_time))   


