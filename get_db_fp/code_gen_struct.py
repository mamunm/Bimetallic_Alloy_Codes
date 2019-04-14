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
from classification import (get_site_dict, get_under_bridge ,
                            get_under_hollow)
from pprint import pprint

def get_correct_pos(POS, CELL):
    """Returns the normalized cell position of an atom """
    fractional = np.linalg.solve(CELL.T, POS.T).T
    fractional %= 1.0
    fractional %= 1.0
    fractional[fractional > 0.99999] = 0.0
    return np.matmul(CELL.T, fractional.T).T

def get_struct(
        symbol='Ag2Au2',
        adsorbate='C',
        site='hollow',
        site_type='Ag_Ag_Au|FCC',
        SB_symbol='L10'):
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
    data = np.load('bulk_data.npy')[()]

    try:
        a = data[symbol]['a']
        c = data[symbol]['c']
    except KeyError:
        if '2' in symbol:
            m1, m2, _ = symbol.split('2')
            try:
                a = data[m1 + '2' + m2 + '2']['a']
                c = data[m1 + '2' + m2 + '2']['c']
            except KeyError:
                a = data[m2 + '2' + m1 + '2']['a']
                c = data[m2 + '2' + m1 + '2']['c']

    if '2' in symbol:
        m1, m2, _ = symbol.split('2')
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
    ads_slab = slab.copy()
    coords, conn = gen.adsorption_sites(slab)
    if site == 'top':
        si = 1
    if site == 'bridge':
        si = 2
    if site == 'hollow':
        si = 3

    adsorbate = molecule(adsorbate)[0]
    adsorbate.set_tags(tags)
    builder = Builder(slab)
    structures = builder.add_adsorbate(adsorbate, index=-1)
    ads_slab = None
    for struct in structures:
        connectivity = struct.connectivity
        connectivity = connectivity[~(connectivity==0).all(1)]
        if n_atoms == 13:
            ads_conn = np.array(connectivity[:, -1]).reshape(n_atoms).sum()
            ads_site = [i for i, j in zip(struct.get_chemical_symbols(),
                        np.array(connectivity[:, -1]).reshape(n_atoms))
                        if j == 1]
        if n_atoms == 14:
            ads_conn = np.array(connectivity[:, -2]).reshape(n_atoms).sum() - 1
            ads_site = [i for i, j in zip(struct.get_chemical_symbols(),
                        np.array(connectivity[:, -2]).reshape(n_atoms)[:13])
                        if j == 1]
        if n_atoms == 15:
            ads_conn = np.array(connectivity[:, -3]).reshape(n_atoms).sum() - 2
            ads_site = [i for i, j in zip(struct.get_chemical_symbols(),
                        np.array(connectivity[:, -3]).reshape(n_atoms)[:13])
                        if j == 1]
        if n_atoms == 16:
            ads_conn = np.array(connectivity[:, -4]).reshape(n_atoms).sum() - 3
            ads_site = [i for i, j in zip(struct.get_chemical_symbols(),
                        np.array(connectivity[:, -4]).reshape(n_atoms)[:13])
                        if j == 1]

        if si != ads_conn:
            continue
        if SB_symbol == 'A1':
                if si == 1 or si == 2:
                    ads_slab = struct
                else:
                    if get_under_hollow(struct[:13]) == site_type.split('|')[1]:
                        ads_slab = struct
        else:
            if si == 1:
                if ads_site[0] == site_type:
                    ads_slab = struct
            elif si == 3:
                hol = get_under_hollow(struct[:13])
                a, b = site_type.split('|')
                a = a.split('_')
                if Counter(ads_site) == Counter(a) and b == hol:
                    ads_slab = struct
            else:
                if SB_symbol == 'L10':
                    if '|' not in site_type:
                        if len(set(ads_site)) == 1:
                            if ads_site[0] == site_type:
                                ads_slab = struct
                    else:
                        if not len(set(ads_site)) == 2:
                            continue
                        if get_under_bridge(struct[:13]) == site_type.split('|')[1]:
                            ads_slab = struct
                else:
                    if '|' in site_type and len(set(ads_site)) == 1:
                        if get_under_bridge(struct[:13]) == site_type.split('|')[1]:
                            ads_slab = struct
                    else:
                        ads_slab = struct

    return slab, ads_slab

if __name__ == '__main__':
    slab, ads_slab = get_struct()
    from ase.visualize import view
    view(ads_slab)
