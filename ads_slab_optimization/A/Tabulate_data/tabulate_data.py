#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
from tabulate import tabulate
import time

def get_rxn_list(B):                                                            
    list_B = []                                                                 
    list_B.append('CH4+' + B + '-->' + 'C-' + B + '+2 H2')                      
    list_B.append('CO+' + B + '-->' + 'C-' + B + '+0.5 O2')                     
    list_B.append('CO2+' + B + '-->' + 'C-' + B + '+O2')                        
    list_B.append('CO2+' + B + '-->' + 'O-' + B + '+CO')                        
    list_B.append('0.5 H2+' + B + '-->' + 'H-' + B)                             
    list_B.append('H2O+' + B + '-->' + 'O-' + B + '+H2')                        
    list_B.append('H2O+' + B + '-->' + 'H-' + B + '+0.5 H2+0.5 O2')             
    list_B.append('H2O2+' + B + '-->' + 'O-' + B + '+H2O')                      
    list_B.append('H2O2+' + B + '-->' + 'H-' + B + '+0.5 H2+O2')                
    list_B.append('CH2O2+' + B + '-->' + 'H-' + B + '+CO2+0.5 H2')              
    list_B.append('CH2O2+' + B + '-->' + 'C-' + B + '+O2+H2')                   
    list_B.append('CH2O2+' + B + '-->' + 'O-' + B + '+CO+H2')                   
    list_B.append('0.5 N2+' + B + '-->' + 'N-' + B)                             
    list_B.append('N2O+' + B + '-->' + 'N-' + B  + '+NO')                       
    list_B.append('N2O+' + B + '-->' + 'O-' + B  + '+N2')                       
    list_B.append('H3N+' + B + '-->' + 'N-' + B  + '+1.5 H2')                   
    list_B.append('NO+' + B + '-->' + 'N-' + B  + '+0.5 O2')                    
    list_B.append('NO+' + B + '-->' + 'O-' + B  + '+0.5 N2')                    
    list_B.append('NO2+' + B + '-->' + 'N-' + B  + '+O2')                       
    list_B.append('NO2+' + B + '-->' + 'O-' + B  + '+NO')                       
    list_B.append('NO2+' + B + '-->' + 'O-' + B  + '+0.5 N2+0.5 O2')            
    list_B.append('0.5 O2+' + B + '-->' + 'O-' + B)                             
    list_B.append('O2S+' + B + '-->' + 'S-' + B + '+O2')                        
    return list_B                                            


def get_ads_nrg(A, B):
    
    RC = True

    Ereac = 0                                                                   
    Eprod = 0 

    r1, r2 = A.split('-->')[0].split('+')                                       
    p2, p3 = None, None 

    if A.split('-->')[1].count('+') == 2:                                       
        p1, p2, p3 = A.split('-->')[1].split('+')                               
    
    elif A.split('-->')[1].count('+') == 1:                                     
        p1, p2 = A.split('-->')[1].split('+')                                   
    
    else:                                                                       
        p1 = A.split('-->')[1]                                                  
    
    if '0.5 ' in r1:                                                            
        Ereac += 0.5 * float(Mol_data[r1.split()[1]])        
    
    else:                                                                       
        Ereac += float(Mol_data[r1])                     
    
    Ereac += float(ads_data[r2]['*']['energy'])  

    try:
        Eprod += float(ads_data[r2][p1.split('-')[0]][B]['energy'])    
    except:
        RC = False

    if p2 is not None:                          

        if '0.5 ' in p2:                                                        
            Eprod += 0.5 * float(Mol_data[p2.split()[1]])    
        
        elif '2 ' in p2:                                                        
            Eprod += 2 * float(Mol_data[p2.split()[1]])      
        
        elif '1.5 ' in p2:                                                      
            Eprod += 1.5 * float(Mol_data[p2.split()[1]])    
        
        else:                                                                   
            Eprod += float(Mol_data[p2])


    if p3 is not None:                                                          
        
        if '0.5 ' in p3:                                                        
            Eprod += 0.5 * float(Mol_data[p3.split()[1]])    
        
        elif '2 ' in p3:                                                        
            Eprod += 2 * float(Mol_data[p3.split()[1]])      
        
        elif '1.5 ' in p3:                                                      
            Eprod += 1.5 * float(Mol_data[p3.split()[1]])    
        
        else:                                                                   
            Eprod += float(Mol_data[p3])                     
    

    if RC:

        DE = '%2.4f' %(Eprod - Ereac)                                
        Ereac = '%2.4f' %(Ereac)                                       
        Eprod = '%2.4f' %(Eprod)                                      
    

        with open('debug.txt', 'a') as dd:                              
        
            dd.write('Reaction: {0}\n'.format(A))                           
            dd.write('    Energy of reactants: {0} eV\n'.format(Ereac))  
            dd.write('    Energy of products:  {0} eV\n'.format(Eprod))      
            dd.write('    Adsorption energy:  {0} eV\n'.format(DE))       
    
        return DE

    else:

        return 'N/A'


############## main() ##############
start_time = time.time()

del_file = ['debug.txt', 'Compiled_data.txt', 'reaction_data.npy']
for f in del_file:
    if os.path.exists(f):
        os.remove(f)

ads_data = np.load('A_ads_data.npy')[()]
Mol_data = np.load('Mol_data.npy')[()]


ads_keys = ['H', 'C', 'O', 'N', 'S']
site_keys = ['ontop', 'bridge', 'fcc', 'hcp']


DATA_TOP = []
DATA_BRIDGE = []
DATA_FCC = []
DATA_HCP = []



keys = sorted(ads_data)

for key in keys:

    dt = [key]
    db = [key]
    df = [key]
    dh = [key]

    dt.append(ads_data[key]['*']['energy'])
    db.append(ads_data[key]['*']['energy'])
    df.append(ads_data[key]['*']['energy'])
    dh.append(ads_data[key]['*']['energy'])

    for a in ads_keys:

        try:
            dt.append('{0:0.2f}'.format(float(
                            ads_data[key][a]['ontop']['energy'])))
        except ValueError:
            dt.append('{0:s}'.format(
                            ads_data[key][a]['ontop']['energy']))
            
        try:
            db.append('{0:0.2f}'.format(float(
                            ads_data[key][a]['bridge']['energy'])))
        except ValueError:
            db.append('{0:s}'.format(
                            ads_data[key][a]['bridge']['energy']))
            
        try:
            df.append('{0:0.2f}'.format(float(
                            ads_data[key][a]['fcc']['energy'])))
        except ValueError:
            df.append('{0:s}'.format(
                            ads_data[key][a]['fcc']['energy']))
        
        try:
            dh.append('{0:0.2f}'.format(float(
                            ads_data[key][a]['hcp']['energy'])))
        except ValueError:
            dh.append('{0:s}'.format(
                            ads_data[key][a]['hcp']['energy']))
        
        dt.append('{0:s}'.format(
                            ads_data[key][a]['ontop']['final_site']))
        
        db.append('{0:s}'.format(
                            ads_data[key][a]['bridge']['final_site']))
        
        df.append('{0:s}'.format(
                            ads_data[key][a]['fcc']['final_site']))
        
        dh.append('{0:s}'.format(
                            ads_data[key][a]['hcp']['final_site']))
    

    DATA_TOP.append(dt)
    DATA_BRIDGE.append(db)
    DATA_FCC.append(df)
    DATA_HCP.append(dh)

    dt = []
    db = []
    df = []
    dh = []


rxn_data = []                                                                   
for f in keys:                                                                
    temp = []                                                                   
    rxn_list = get_rxn_list(f)                                                  
    for rl in rxn_list:                                                         
        f_data = []                                                             
        f_data.append(rl)                                                       
        for s in site_keys:                                                     
            f_data.append(get_ads_nrg(rl, s))                                   
        temp.append(f_data)                                                     
    rxn_data.append(temp)                                                       
                                                                                
           
np.save('reaction_data.npy', rxn_data)
                                                                     
header = ['Site', 'E* [eV]', 'E_H* [eV]', 'Fin_Ads_Site', 'E_C* [eV]', 
          'Fin_Ads_Site', 'E_O* [eV]', 'Fin_Ads_Site', 'E_N* [eV]', 
          'Fin_Ads_Site', 'E_S* [eV]', 'Fin_Ads_Site']

pr_st = '\n\nInformation for calculations started on {0} site:\n\n'

pr_st_2 = '\n\nAdsorption energetics on {0} in eV:\n\n'     

header2 = ['Reaction', 'ontop site', 'bridge site', 'fcc site', 'hcp site']     

with open('Compiled_data.txt', 'w') as dd:  

    dd.write(pr_st.format('ontop'))
    dd.write(tabulate(DATA_TOP, header, tablefmt='grid', 
        numalign='center', stralign='center', floatfmt='0.2f'))

    dd.write(pr_st.format('bridge'))
    dd.write(tabulate(DATA_BRIDGE, header, tablefmt='grid',
        numalign='center', stralign='center', floatfmt='0.2f'))

    dd.write(pr_st.format('fcc'))
    dd.write(tabulate(DATA_FCC, header, tablefmt='grid', 
        numalign='center', stralign='center', floatfmt='0.2f'))

    dd.write(pr_st.format('hcp'))
    dd.write(tabulate(DATA_HCP, header, tablefmt='grid',
        numalign='center', stralign='center', floatfmt='0.2f'))

    dd.write('\n\nInformation of gas phase molecules:\n\n')
    dd.write(tabulate([[a, b] for a, b in Mol_data.items()], 
             ['Molecule', 'Energy [eV]'], tablefmt='grid',
             numalign='center', stralign='center', floatfmt='0.2f'))
    
    dd.write('\n')                                                              
    
    for f, data in zip(keys, rxn_data):  

        dd.write(pr_st_2.format(f))                                            
        dd.write(tabulate(data, header2, tablefmt='grid',                       
            numalign='center', stralign='center', floatfmt='0.2f'))  

        dd.write('\n')                             


print('\nTotal Elapsed time: {0} seconds'.format(time.time() - start_time))
