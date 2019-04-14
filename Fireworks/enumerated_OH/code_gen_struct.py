#!/usr/bin/env python
# -*-coding: utf-8 -*-

#code.py
#Osman Mamun
#LAST UPDATED: 08-20-2018

from catkit.gen.surface import SlabGenerator
from ase.build import bulk
import numpy as np
from catkit.gen.adsorption import Builder
from ase import Atom
from ase import Atoms
from ase.build import add_adsorbate
from catkit.build import molecule
from collections import Counter
from pprint import pprint
from ase.db import connect

def get_struct(
        symbol='Ag2Au2',
        adsorbate='OH',
        SB_symbol='L10',
        a = 4,
        c = 3):
    """It returns the structure generated with catkit
       using the keyword given."""
    size = (1, 1)
    layers = 3
    miller = (1, 1, 1)
    vacuum = 10
    get_height = dict([('C', 1.0),
                       ('O', 1.0),
                       ('N', 1.0),
                       ('S', 1.0),
                       ('H', 1.0)])
    if adsorbate in ['CH', 'OH', 'NH', 'SH']:
        tags = [-1, -2]
        n_atoms = 14
    elif adsorbate in ['CH2', 'OH2']:
        tags = [-1, -2, -2]
        n_atoms = 15
    elif adsorbate == 'CH3':
        tags = [-1, -2, -2, -2]
        n_atoms = 16
    else:
        tags = [-1]
        n_atoms = 13

    if '2' in symbol:
        m1, m2, _ = symbol.split('2')
        #b = bulk(m1, cubic=True, crystalstructure='fcc', a=a, c=c)
        #b[2].symbol = m2
        #b[3].symbol = m2
        name = m1 + '2' + m2 + '2'
        b=Atoms(name,
                    scaled_positions=[(0, 0, 0),
                                      (0.5, 0.5, 0),
                                      (0.5, 0, 0.5),
                                      (0, 0.5, 0.5)],
                    cell=[a, a, c],
                    pbc=True)
        b.set_pbc((1,1,1))
    elif '3' in symbol:
        m1, m2 = symbol.split('3')
        b = bulk(m1, cubic=True, crystalstructure='fcc', a=a, c=c)
        b[3].symbol = m2
    else:
        b = bulk(symbol, cubic=True, crystalstructure='fcc', a=a, c=c)

    gen = SlabGenerator(b,
            miller_index=miller,
            vacuum=vacuum,
            layers=layers,
            fixed=2,
            layer_type='trim',
            standardize_bulk=False,
            primitive=False)

    terminations = gen.get_unique_terminations()
    if len(terminations) > 1:
        raise IncosistencyError('There are more than one temrination.')

    slab = gen.get_slab(size=size, iterm=-1)
    mm = {'Fe': 3, 'Ni': 1, 'Co': 1, 'Mn': 1}
    slab.set_initial_magnetic_moments([mm.get(i, 0)
                                      for i in slab.get_chemical_symbols()])

    adsorbate = molecule(adsorbate)[0]
    adsorbate.set_tags(tags)
    builder = Builder(slab)
    structures = builder.add_adsorbate(adsorbate, index=-1)
    return structures

if __name__ == '__main__':
    db = connect('OH.db')
    bulks = np.load('bulk_data.npy')[()]

    for b in bulks:
        a = bulks[b]['a']
        c = bulks[b]['c']
        if '3' in b:
            sbs = 'L12'
            m, n = b.split('3')
            mi = (1, 1, 1)
        elif '2' in b:
            sbs = 'L10'
            m, n, _ = b.split('2')
            mi = (0, 1, 1)
        else:
            sbs = 'A1'
            m = b
            n = 'N/A'
            mi = (1, 1, 1)
        pol_metal = ['Fe', 'Ni', 'Mn', 'Co']
        stru = get_struct(symbol=b, a=a, c=c, SB_symbol=sbs)
        for s_ind, s in enumerate(stru):
            if m in pol_metal or n in pol_metal:
                spol = True
            else:
                spol = False
            data = {'connectivity': s.connectivity,
                    'parameters': {'beefensemble': True,
                                   'calcstress': True,
                                   'convergence': {'diag': 'david',
                                                   'energy': 0.000136,
                                                   'maxsteps': 2000,
                                                   'mix': 4,
                                                   'mixing': 0.3,
                                                   'nmix': 10},
                                   'dipole': {'status': True},
                                   'dw': 5000,
                                   'fmax': 0.05,
                                   'kpts': np.array([6, 6, 1]),
                                   'mode': 'relax',
                                   'nbands': -10,
                                   'nosym': True,
                                   'opt_algorithm': 'bfgs',
                                   'outdir': '.',
                                   'output': {'avoidio': False,
                                   'removesave': True,
                                   'removewf': True,
                                   'wf_collect': False},
                                   'pw': 500,
                                   'sigma': 0.15,
                                   'spinpol': spol,
                                   'xc': 'BEEF-vdW'},
                    'surface_atoms': np.array([ 8,  9, 10, 11])}
            kvp = {'SB_symbol': sbs,
                   'ads0': 'OH',
                   'bulk_name': b,
                   'fixed': 2,
                   'layers': 3,
                   'metal0': m,
                   'metal1': n,
                   'miller_index': str(mi),
                   'site_index': s_ind,
                   'vacuum': 10}
            db.write(s, key_value_pairs=kvp, data=data)
