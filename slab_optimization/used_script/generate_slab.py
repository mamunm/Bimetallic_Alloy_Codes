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

prob_metal = ['AlFe', 'AlMn', 'AlNi', 'CrFe', 'CrMn', 'CrNi', 'FeSc', 'MnSc',
              'NiSc', 'FeTi', 'MnTi', 'NiTi', 'FeV', 'MnV', 'NiV', 'Ag3Co',
              'Al3Co', 'Au3Co', 'Cd3Co', 'Cr3Co', 'Cu3Co', 'Ga3Co', 'Hf3Co',
              'Hg3Co', 'In3Co', 'Ir3Co', 'La3Co', 'Mo3Co', 'Nb3Co', 'Os3Co',
              'Pb3Co', 'Pd3Co', 'Pt3Co', 'Re3Co', 'Rh3Co', 'Ru3Co', 'Sc3Co',
              'Sn3Co', 'Ta3Co', 'Tc3Co', 'Ti3Co', 'Tl3Co', 'V3Co', 'W3Co',
              'Y3Co', 'Zn3Co', 'Zr3Co']

for i in prob_metal:
    if '3' in i:
        name = i
        sbs = 'A3B'
        lib_file_loc = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/'
        lib_file_loc += 'Surface/' + sbs + '/' + name + '/' + name  
        lib_file_loc += '_surface.py'
    else:
        name = i[:2] + '2' + i[2:] + '2'
        sbs = 'AB'
        lib_file_loc = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/'
        lib_file_loc += 'Surface/' + sbs + '/' + name + '/' + name + '.py'
        print(name)

    mkdirp(name)
    working_file_loc = name + '/' + name + '.py'
    try:
        copyfile( lib_file_loc, working_file_loc)
        inplace_change(working_file_loc , "#!/home/vossj/suncat/bin/python", 
                   "#!/home/users/vossj/suncat/bin/python_s2.0")
        inplace_change(working_file_loc , "#SBATCH -p iric,owners,normal", 
                   "#SBATCH -p iric,owners,normal, suncat")
        inplace_change(working_file_loc , "/home/vossj/suncat/psp/gbrv1.5pbe", 
                   "/home/users/vossj/suncat/psp/gbrv1.5pbe")
        inplace_change(working_file_loc , "spinpol=False", "spinpol=True")
        inplace_change(working_file_loc , 
                   "#atom.set_initial_magnetic_moments([2.0 for atom in Slab])", 
                   "atom.set_initial_magnetic_moments([2.0 for atom in Slab])") 
    except IOError:
        print(name)
