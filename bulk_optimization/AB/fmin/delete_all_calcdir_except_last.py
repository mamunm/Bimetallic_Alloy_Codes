#!/usr/bin/env python

import os
import glob
import shutil

List = os.walk('.').next()[1]

for files in List:
    with open(files + '/' + files +'_bulk.out') as f:
        for line in f:
            if 'Function evaluations' in line:
                niter = line.split()[2]
    list_calc = glob.glob(files + '/calcdir-*')
    for lc in list_calc:
        if not niter in lc:
            shutil.rmtree(lc)
            
