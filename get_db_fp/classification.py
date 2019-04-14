#!/usr/bin/env python
# -*-coding: utf-8 -*-

#claccification.py
#Osman Mamun
#LAST UPDATED: 03-13-2018

import os
import numpy as np
from ase.db import connect
from ase.io import write
from ase.build import make_supercell
from ase.data import covalent_radii as cradii
from scipy.spatial import Voronoi

def compare_slab(A, B):

    hollow_dist = np.linalg.norm(B.positions[8, :2] - B.positions[0, :2])
    a = A.positions[-5:-1]
    b = B.positions[-5:-1]
    allowed_z_movement = 0.5 * cradii[A.get_atomic_numbers()[0]]

    d_xy = np.empty(4)
    #d_z = np.empty(4)

    for i in range(len(a)):

        d_xy[i] = np.linalg.norm(a[i][:2] - b[i][:2])
        #d_z[i] = np.linalg.norm(a[i][2] - b[i][2])

    cond1 = np.all(np.absolute(d_xy) < 0.4 * hollow_dist)
    #cond2 = np.all(np.absolute(d_z) < 0.2 * lat)
    cond2 = np.all([abs(a[i][2] - b[i][2]) < allowed_z_movement
                   for i in range(1,4)])

    if cond1 and cond2:

        return True

    else:

        return False

def is_subsurface(A):

    pos0 = A.positions[-5:-1][:, 2]
    pos1 = A.positions[-1][2]
    metal_covalent_radii = cradii[A.get_atomic_numbers()[-5:-1]]
    ads_covalent_radii =cradii[A.get_atomic_numbers()[-1]]

    if np.any([(pos0[i] - pos1) >  0.5 * metal_covalent_radii[i]
               for i in range(4)]):

        return True

    else:

        return False

def is_ontop(A):

    pos0 = A.positions[-5:-1][:, :2]
    pos1 = A.positions[-1][:2]
    d_xy = np.empty(4)

    for i in range(len(pos0)):

        d_xy[i] = np.linalg.norm(pos0[i] - pos1)

    if np.any(abs(d_xy) < 0.15 * cradii[A.get_atomic_numbers()[0]]):

        return True

    else:

        return False

def is_bridge(A):

    def get_angle(A, B, C):
        AB = B - A
        AC = C - A
        return np.arccos(np.dot(AB, AC) / np.linalg.norm(AB) /
                np.linalg.norm(AC)) * 180 / np.pi

    pos0 = A.positions[-1][:2]
    pos_metal = A.positions[-5:-1, :2]


    cond = False
    for i, pm1 in enumerate(pos_metal):
        for j, pm2 in enumerate(pos_metal):
            if i != j and np.linalg.norm(pm1 - pm2) < 1.1 * A.get_cell()[1][0]:
                if 160 < get_angle(pos0, pm1, pm2) < 200:
                    cond = True

    if cond:

        return True

    else:

        return False

def get_hollow(A):

    from ase.build import make_supercell
    P = [[3, 0, 0], [0, 3, 0], [0, 0, 1]]
    SC = make_supercell(A, P)

    second_layer_Z = A.positions[4][2]
    ads_pos = SC.positions[25][:2]

    second_layer = []
    for pos in SC.positions:
        if abs(pos[2] - second_layer_Z) < 0.05:
            second_layer += [pos[:2]]

    if np.any(np.linalg.norm(second_layer - ads_pos, axis=1) <
            0.25 * cradii[A.get_atomic_numbers()[0]]):

        return 'HCP'

    else:

        return 'FCC'


def get_site_dict(B):

    def apos_not_in_D(A, B):
        for k in D:
            if 'bridge' in k:
                tpos = np.array([D[k]['x'], D[k]['y'], D[k]['z']])
                if np.allclose(tpos, A):
                    return False
        return True

    def vnotinD(A, B):
        for k in D:
            if 'hollow' in k:
                tpos = np.array([D[k]['x'], D[k]['y']])
                if np.allclose(tpos, A):

                    return False

        return True

    C = B.copy()
    C = C[-4:]
    D = {'ontop_site-' + str(i): {'x': pos[0], 'y': pos[1], 'z': pos[2],
         'sym': sym} for i, pos, sym in zip(range(1, 5),
         C.positions, C.get_chemical_symbols())}
    SC = C * (3, 3, 1)
    cell = C.get_cell()

    if cell[1][0] < 0:
        SC.translate([cell[1][0], -cell[1][1], 0])
    else:
        SC.translate([-cell[0][0]-cell[1][0], -cell[1][1], 0])

    top_ind = []
    for i, atom in enumerate(SC):
        for key in D:
            pos = np.array([D[key]['x'], D[key]['y'], D[key]['z']])
            if np.allclose(atom.position, pos, atol=1e-2):
                top_ind += [i]

    t_ind = 5
    for i, atom in enumerate(SC):
        if i not in top_ind:
            k = 'ontop_site-' + str(t_ind)
            D[k] = {}
            D[k]['x'], D[k]['y'], D[k]['z'] = atom.position
            D[k]['sym'] = atom.symbol
            t_ind += 1

    b_dist = 0.5 * cell[0][0]
    b_ind = 1

    for i in top_ind:
        pos = SC.positions[i]
        for j, atom in enumerate(SC):
            apos = atom.position
            cond = np.linalg.norm(apos - pos) < 1.5 * b_dist
            temp = 0.5 * (pos + apos)
            cond2 = apos_not_in_D(temp, D)
            if i != j and cond and cond2:
                k = 'bridge_site-' + str(b_ind)
                D[k] = {}
                D[k]['x'] = temp[0]
                D[k]['y'] = temp[1]
                D[k]['z'] = temp[2]
                D[k]['sym'] = SC.get_chemical_symbols()[i] + '_' + atom.symbol
                b_ind += 1

    h_dist = np.linalg.norm(B.positions[11, :2] - B.positions[4, :2])
    h_ind = 1

    vor = Voronoi(SC.positions[:, :2])

    vertices = vor.vertices

    for v in vertices:
        for t in top_ind:
            cond = np.linalg.norm(v - SC.positions[t, :2]) < 1.5 * h_dist
            cond1 = vnotinD(v, D)
            if cond and cond1:
                k = 'hollow_site-' + str(h_ind)
                D[k] = {}
                D[k]['x'] = v[0]
                D[k]['y'] = v[1]
                D[k]['z'] = max(SC.positions[:, 2])
                min_dist = sorted([np.linalg.norm(v - m[:2])
                    for m in SC.positions])[3]
                symb = [atom.symbol for atom in SC if np.linalg.norm(v -
                    atom.position[:2]) < min_dist]
                D[k]['sym'] = '_'.join(symb)
                h_ind += 1

    return D

def get_under_bridge(A):

    B = A.copy()
    apos = B.positions[-1]
    B = B[np.r_[4:8]]
    SC = B * (3, 3, 1)
    cell = B.get_cell()

    if cell[1][0] < 0:
        SC.translate([cell[1][0], -cell[1][1], 0])
    else:
        SC.translate([-cell[0][0]-cell[1][0], -cell[1][1], 0])

    ret = None
    dis = cell[0][0]

    for ele in SC:
        new_dis = np.linalg.norm(apos - ele.position)
        if new_dis < dis:
            dis = new_dis
            ret = ele.symbol

    return ret

def get_under_hollow(A):

    B = A.copy()
    apos = B.positions[-1]
    B = B[np.r_[4:8]]
    SC = B * (3, 3, 1)
    cell = B.get_cell()

    if cell[1][0] < 0:
        SC.translate([cell[1][0], -cell[1][1], 0])
    else:
        SC.translate([-cell[0][0]-cell[1][0], -cell[1][1], 0])

    ret = 'FCC'

    if np.any([np.linalg.norm(apos[:2] - ele.position[:2]) < 0.5 *
              cradii[ele.number] for ele in SC]):
        ret = 'HCP'

    return ret

def get_site_type(B, S, ):

    Dict = get_site_dict(B)
    apos = B.positions[-1, :2]

    f_a_s = None
    dis = B.get_cell()[0][0]
    Kind = None

    for k in Dict:
        new_dis = np.linalg.norm(apos - np.array([Dict[k]['x'], Dict[k]['y']]))
        if new_dis < dis:
            dis = new_dis
            f_a_s = k.split('_')[0]
            Kind = k

    print(f_a_s)

    if f_a_s == 'ontop':
        s_t = Dict[Kind]['sym']

    if f_a_s == 'bridge':

        a, b = Dict[Kind]['sym'].split('_')

        if S == 'L10':

            if a == b:
                s_t = a
            else:
                s_t = a + '_' + b +'|' + get_under_bridge(B)

        else:

            if a == b:
                s_t = a + '|' + get_under_bridge(B)
            else:
                s_t = a + '_' + b

    if f_a_s == 'hollow':
        s_t = Dict[Kind]['sym'] + '|' + get_under_hollow(B)


    return s_t

def get_ads_energy(A, B, C):

    def get_eslab(M):
        db_S = connect('Surface.db')
        for d in db_S.select():
            if d['unique_formula'] == M:
                return d['energy']

    mol_data = {'CH4': -224.912,
                'H2': -32.920,
                'O2': -882.691,
                'N2': -553.612,
                'H2S': -365.396,
                'H2O': -476.622,
                'CO2': -1044.747,
                'NH3': -326.926}

    Eslab = get_eslab(C)

    ads_nrg = []
    ads_eqn = []

    if B == 'H':

        ads_nrg += [A - Eslab - 0.5 * mol_data['H2']]
        ads_eqn += ['0.5H2+*=>H*']

    elif B == 'C':

        ads_nrg += [A + 2 * mol_data['H2'] - Eslab - mol_data['CH4']]
        ads_eqn += ['CH4+*=>C*+2H2']
        ads_nrg += [A + mol_data['O2'] - Eslab - mol_data['CO2']]
        ads_eqn += ['CO2+*=>C*+O2']

    elif B == 'O':

        ads_nrg += [A - Eslab - 0.5 * mol_data['O2']]
        ads_eqn += ['0.5O2+*=>O*']
        ads_nrg += [A + mol_data['H2'] - Eslab - mol_data['H2O']]
        ads_eqn += ['H2O+*=>O*+H2']

    elif B == 'N':

        ads_nrg += [A - Eslab - 0.5 * mol_data['N2']]
        ads_eqn += ['0.5N2+*=>N*']
        ads_nrg += [A + 1.5 * mol_data['H2'] - Eslab - mol_data['NH3']]
        ads_eqn += ['NH3+*=>N*+1.5H2']

    else:

        ads_nrg += [A + mol_data['H2'] - Eslab - mol_data['H2S']]
        ads_eqn += ['0.5H2S+*=>S*+H2']

    return ads_eqn, ads_nrg


