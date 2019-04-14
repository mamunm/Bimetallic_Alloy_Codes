#!/usr/bin/env python
# -*-coding: utf-8 -*-

#test.py
#Osman Mamun
#DATE CREATED: 10-09-2018

import numpy as np
from scaling import Linear_Scaling

x = np.arange(20)
y = np.arange(20) + np.random.rand(20)

sc = Linear_Scaling(x, y, 'H', 'C', 'eV')
sc.get_coeff()
a = sc.plot_scaling()
a.show()

