#!/usr/bin/env python
# -*-coding: utf-8 -*-
#process.py

import sys 
import os 
import numpy as np  
import sqlite3 
from ase.db import connect  
from collections import deque  

slab_db = connect('/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf/Processing/From_Winther/master_classified_new.db') 
con = slab_db._connect()
cur = con.cursor() 
count = 0 
print('Starting data processing...')
for row in slab_db.select('inp_file,relaxed=1'): 
    print('Processing {0}'.format(row.id))
    if row.get('fw_id', None): 
        print('FW!')    
        continue           
    energy = None
    if not row.SB_symbol=='A1':
        dir = os.path.dirname(row.inp_file) 
        dir += '/Calculations/{ads}_{site}_site_{site_id}'.format(ads=row.adsorbate, 
             site=row.initial_site, site_id=row.initial_site_id)              
    else: 
        dir = os.path.dirname(row.inp_file) 
        dir +=  '/Calculations/{ads}_{site}'.format(ads=row.adsorbate, 
                                                    site=row.initial_site)
    try: 
        f = deque(open(dir + '/out_scf.log'), 1)      
    except:
        #slab_db.delete([row.id])       
        print(dir)       
        #sys.exit()     
        continue               
    for line in f:
        energy = float(line.split()[-2].strip('*'))                            
    diff = energy - row.energy 
    #print(diff)  
    if not np.isclose(diff, 0, atol=1e-5): 
        count += 1   
        print(row.formula)      
        print('ENERGY DIFF:  {}'.format(diff))      
        print(count)       
        cur.execute('UPDATE systems set energy={energy} WHERE id={id}'.format(energy=energy, id=row.id))     
con.commit()              
con.close()
