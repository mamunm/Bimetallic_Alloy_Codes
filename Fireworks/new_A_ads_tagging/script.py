#!/usr/bin/env python
from fireworks import LaunchPad, Firework, Workflow, PyTask
from ase.io import read
from qescripts.fwio import atoms_to_encode, encode_to_atoms
from ase.db import connect
from ase.visualize import view
import numpy as np
import os
from ase.build import make_supercell
from scipy.spatial import Voronoi
from ase.data import covalent_radii as cradii

def get_hollow(B):

    C = B.copy()
    apos = B.positions[-1]
    C = C[np.r_[4:8]]
    SC = C * (3, 3, 1)
    cell = B.get_cell()

    if cell[1][0] < 0:
        SC.translate([cell[1][0], -cell[1][1], 0])
    else:
        SC.translate([-cell[0][0]-cell[1][0], -cell[1][1], 0])

    ret = 'FCC'

    if np.any([np.linalg.norm(apos[:2] - ele.position[:2]) < 0.05 *
              cradii[ele.number] for ele in SC]):
        ret = 'HCP'

    return ret

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
    C = C[-5:-1]
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

def get_site_ads(A):

    Dict = get_site_dict(A)
    apos = A.positions[-1, :2]
    ads = A.get_chemical_symbols()[-1]
    dis = 100    

    for k in Dict:
        new_dis = np.linalg.norm(apos - np.array([Dict[k]['x'], Dict[k]['y']]))
        if new_dis < dis:
            dis = new_dis
            f_a_s = k.split('_')[0]

    if f_a_s == 'hollow':
        f_a_s = get_hollow(A)

    return f_a_s, ads

host = 'suncatls2.slac.stanford.edu'
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)


ids = [i for i in range(30660, 31088)]
dic = {}

db = connect('check.db')

for i in ids:
    print('Processing: {0}'.format(i))
    launch = launchpad.get_fw_by_id(i)
    atoms = encode_to_atoms(launch.spec['_tasks'][0]['args'][0])[0]
   
    try:
        des = launch.name['calc']['description']
        metal = launch.name['calc']['metal']

    except KeyError:
        continue

    if not des == 'New lattice pure metal ads stuffs':
        continue
    
    dic[i] = {}
    site, ads = get_site_ads(atoms)
    dic[i]['site'] = site
    dic[i]['adsorbate'] = ads
    if site == 'ontop':
        atoms.set_tags([1 for _ in atoms])
    elif site == 'bridge':
        atoms.set_tags([2 for _ in atoms])
    elif site == 'FCC':
        atoms.set_tags([3 for _ in atoms])
    else:
        atoms.set_tags([4 for _ in atoms])

    db.write(atoms)

np.save('dict.npy', dic)
