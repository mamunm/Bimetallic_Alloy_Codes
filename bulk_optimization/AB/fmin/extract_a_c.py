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
List = os.walk('.').next()[1]
#print 'The file list looks like:'
#print List
#List = os.listdir('.')
#List = filter(lambda x: x != 'Co2Mn2', List)
for f in List:
    if (os.path.exists( f + '/' + f + '_bulk.out')):
        print 'Extracting data from: ' + f + '/' + f + '_bulk.out'
        with open(f + '/' + f + '_bulk.out') as ff:
            for line in ff:
                if 'Lattice constant a:' in line:
                    filename.append(f)
                    a.append(line.split()[3])
                if 'Lattice constant c:' in line:
                    c.append(line.split()[3])
        print 'Optmization data found in folder :' + f 

with open('Output.txt','w') as fff:
    fff.write('Calculated optimized lattice constant using scipy fmin approach:\n\n')
    fff.write('Structure              a              c                  a/c\n')
    fff.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    for i in range(len(filename)):
        fff.write('%9s %15s %15s %15s\n' %(filename[i],a[i],c[i],str(float(c[i])/float(a[i]))))
 
       
        
     

