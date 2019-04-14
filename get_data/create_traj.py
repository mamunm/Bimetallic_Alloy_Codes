#!/usr/bin/env python
# -*-coding: utf-8 -*-

from ase.io import read, write 
from glob import glob
import ase

data = ['O', 'OH']

for datum in data:
    fs = glob('{}/*'.format(datum))
    for f in fs:
        print('Wroking on {}'.format(f))
        try:
            d = read(f + '/pw.pwo')
            write(f + '/opt.traj', d)
        except AssertionError:
            print('Found assertion error.')
            with open('debug_traj.txt', 'a+') as g:
                g.write('Assertion problem found in file: {}\n'.format(f))
        except ase.io.formats.UnknownFileTypeError:
            with open('debug_traj.txt', 'a+') as g:
                g.write('Empty file problem found in file: {}\n'.format(f))
