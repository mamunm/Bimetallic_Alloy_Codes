#!/usr/bin/env python
# -*-coding: utf-8 -*-

#__init__.py
#Osman Mamun
#LAST UPDATED: 09-10-2018

#from .neuralnet import NN
from .data_preprocess import preprocess_data
from .run_gp import OMGP

__all__ = ['preprocess_data', 'OMGP']
__version__ = '0.0.1'
