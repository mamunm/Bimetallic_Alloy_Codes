#!/usr/bin/env python
# -*-coding: utf-8 -*-

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

def mkdirp(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)

def read_lat(A):
    f=open("/home/mamunm/DATABASE/Binary_Alloy_Project/Bulk/A/"+A+"/bulk.out")
    for line in f:
        if "Lattice constant:" in line:
           return line.split()[2]


def lattice_guess_L12(A,B):
    #assert (isinstance(A,str), 'Arguments must be of string type.'
    #assert (isinstance(B,str), 'Arguments must be of string type.'
    lat1=read_lat(A)
    lat2=read_lat(B)
    return float(lat1)*0.75+float(lat2)*0.25


metals=['Al','Sc','Ti','V','Cr','Mn','Fe','Ni','Cu','Zn','Ga', 'La',\
        'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Hf', 'Y', \
        'Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']
pol_metal=['Fe','Ni','Mn']

for i,metal1 in enumerate(metals):
    name=metal1+"3Co"
    mkdirp(name)
    copyfile('../Cu3Re/A3B_bulk.py',name+'/'+name+'.py')
    a=lattice_guess_L12(metal1,'Co')
    inplace_change(name+'/'+name+'.py',"a=3.95","a=" + str(a))
    inplace_change(name+'/'+name+'.py',"metal = 'Cu'","metal =  '"+metal1+"'")
    inplace_change(name+'/'+name+'.py',"metal2 = 'Re'","metal2 =  'Co'")
    inplace_change(name+'/'+name+'.py',"Cu3Re_bulk",name +"_bulk")
    inplace_change(name+'/'+name+'.py',"06:00:00","08:00:00")
    if metal1 in pol_metal:
        inplace_change(name+'/'+name+'.py',"spinpol=False,","spinpol=True,")
        inplace_change(name+'/'+name+'.py',"atoms.set_pbc((1,1,1))","atoms.set_pbc((1,1,1)) \n    for atom in atoms: \n        atom.magmom=2")
