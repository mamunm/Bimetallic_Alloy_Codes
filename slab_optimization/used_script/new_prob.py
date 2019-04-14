#!/usr/bin/env python
# -*-coding: utf-8 -*-

#new_prob.py
#Osman Mamun
#LAST UPDATED: 06-27-2018
from shutil import copyfile
def inplace_change(filename, old_string, new_string):                           
    # Safely read the input filename using 'with'                               
    with open(filename) as f:                                                   
        s = f.read()                                                            
        if old_string not in s:                                                 
            print '"{old_string}" not found in {filename}.'.format(**locals())  
            return                                                              
                                                                                
    # Safely write the changed content, if found in the file                    
    with open(filename, 'w') as f:                                              
        print 'Changing "{old_string}" to "{new_string}" in {filename}'.format(**locals())
        s = s.replace(old_string, new_string)                                   
        f.write(s)    

new_prob = ['Fe2Sc2', 'Mn2Sc2', 'Ni2Sc2', 'Fe2Ti2', 'Mn2Ti2', 'Ni2Ti2',
             'Fe2V2', 'Mn2V2', 'Ni2V2']

for i in new_prob:
    des = i + '/' + i + '.py'
    copyfile('proto.py', des)
    inplace_change(des, '--job-name=Al2Fe2_Slab', 
                   '--job-name={0}_Slab'.format(i)) 
    inplace_change(des, '--job-name=Al2Fe2_surface', 
                   '--job-name={0}_surface'.format(i)) 
    inplace_change(des, 'Fe2Sc2', 
                   '{0}'.format(i)) 


