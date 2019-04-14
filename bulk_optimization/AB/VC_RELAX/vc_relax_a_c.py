#!/bin/env python
# coding: utf-8
from __future__ import division, unicode_literals
"""
Created on Oct 18, 2017
"""
__author__ = "Osman Mamun"
__version__ = "0.1"
__maintainer__ = "Osman Mamun"
__email__ = "mamunm@stanford.edu"
__date__ = "Oct 18, 2017"

import os 
import sys
from shutil import copyfile
import glob
from itertools import islice

a = []
c = []
filename = []
List = os.listdir('vc_relax/.')
List = filter(lambda x: x != 'Co2Mn2', List)
for f in List:
    if (os.path.exists('vc_relax/'+f+'/cell_relax/log')):
        print 'Opening file: vc_relax/'+f+'/cell_relax/log'
        with open('vc_relax/'+f+'/cell_relax/log') as ff:
            for line in ff:
                if 'CELL_PARAMETERS' in line:
                    AA=''.join(islice(ff,3))
        print 'Lattice matrix found in file'+f+'/cell_relax/log is'    
        print AA
        filename.append(f)
        a.append(AA.split('\n')[0].split()[0])
        c.append(AA.split('\n')[2].split()[2])

with open('Output.txt','w') as fff:
    fff.write('Calculated optimized lattice constant using variable cell approach:\n\n')
    fff.write('Structure              a              c                  a/c\n')
    fff.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    for i in range(len(filename)):
        fff.write('%9s %15s %15s %15s\n' %(filename[i],a[i],c[i],str(float(c[i])/float(a[i]))))
 
       
        
     

