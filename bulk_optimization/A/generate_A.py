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
from shutil import rmtree
from subprocess import call
from glob import glob

files = glob('*')
files = [i for i in files if i not in ['Cu', 'generate_A.py']]
for f in files:
    rmtree(f)

def inplace_change(filename, old_string, new_string):
    # Safely read the input filename using 'with'
    with open(filename) as f:
        s = f.read()
        if old_string not in s:
            print('"{old_string}" not found in {filename}.'.format(**locals()))
            return

    # Safely write the changed content, if found in the file
    with open(filename, 'w') as f:
        print('Changing "{old_string}" to "{new_string}" in {filename}'.format(**locals()))
        s = s.replace(old_string, new_string)
        f.write(s)

def mkdirp(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory) 

metals=['Al','Sc','Ti','V','Cr','Mn','Fe','Co','Ni', 'Zn','Ga','Y', \
        'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','La','Hf', \
        'Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi']

for i,metal1 in enumerate(metals):
    name=metal1
    mkdirp(name)
    copyfile('Cu/Cu.py',name+'/'+name+'.py')
    inplace_change(name+'/'+name+'.py',"#SBATCH --job-name=Cu_bulk","#SBATCH --job-name=" + metal1 +"_bulk")
    inplace_change(name+'/'+name+'.py',"metal =  'Cu'","metal =  '" + metal1 + "'")
