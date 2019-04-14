#!/usr/bin/env python
# -*-coding: utf-8 -*-

#submit.py
#Osman Mamun
#LAST UPDATED: 06-11-2018
from fireworks import LaunchPad, Firework, Workflow, PyTask
from ase.db import connect
from catkit.flow.fwio import atoms_to_encode
from glob import glob

def final_traj_exists(M1, M2, A, AM):
    val = False
    mainloc = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf'
    abloc += mainloc + '/AB/'
    if M1 < M2:
        abloc += '{0}2{1}2'.format(M1, M2)
    else:
        abloc += '{1}2{0}2'.format(M1, M2)
    abloc += '/Calculations/{0}_top_site_{1}/'.format(A, AM)
    ff = glob(abloc + 'final*.traj')
    if len(ff) == 1:
        val = True
    
    if not val:
        abloc += mainloc + '/AB_from_SLAC_SUNCAT/'                   
        if M1 < M2:                                                        
            abloc += '{0}2{1}2'.format(M1, M2)
        else:                                           
            abloc += '{1}2{0}2'.format(M1, M2)                                  
            abloc += '/Calculations/{0}_top_site_{1}/'.format(A, AM)            
            ff = glob(abloc + 'final*.traj')                                    
            if len(ff) == 1:                                         
                val = True                   
    if val:
        print('Omitting this caluclation!')
    return val


print('Starting work...')
host = 'suncatls2.slac.stanford.edu'  
username, name, password = netrc().authenticators(host)

launchpad = LaunchPad(
    host=host,
    name=name,
    username=username,
    password=password)

print('Connected to launchpad!')
db = connect('single-atom-prototypes.db')

parameters = {
              'mode': 'relax',
              'opt_algorithm': 'bfgs',
              'xc': 'BEEF-vdW',
              'kpts': (6, 6, 1),
              'pw': 500,
              'dw': 5000,
              'sigma': 0.15,
              'nbands': -10,
              'fmax': 0.05,
              'beefensemble': True,
              'nosym': True,
              'outdir': '.',
              'calcstress': True,
              'dipole': {'status': True},
              'output': {
                         'removesave': True,
                         'avoidio': False,
                         'removewf': True,
                         'wf_collect': False,
                        },
              'convergence': {
                              'energy': 1e-5 * 13.6,
                              'mixing': 0.3,
                              'nmix': 10,
                              'mix': 4,
                              'maxsteps': 2000,
                              'diag': 'david'},
             }

print('Starting the loop!')
tot = 0
missing = 0

for d in len(db)
    at = db.get(d)
    if not at.SB_symbol == L10:
        continue
    if not at.ads0_coord == 1:  
        continue
    tot += 1
    print('working on the {0} system.'.format(d.formula))
    if not at.ads0 in ['H', 'C', 'N', 'O', 'S']: 
        continue
    if not '6' in at.formula:
        continue
    m1, m2, a = at.formula.split('6')
    for i, b in enumerate(at.positions[-5:-1]):
        if np.allclose(b[:2], a.positions[-1][:2]):
            am = at.symbols[i]

    print('Working on system with m1 = {0}, m2 = {1}, a = {2}'.format(m1, m2,
        a))

    if not final_traj_exists(m1, m2, a, am):
        continue

    missing += 1
    '''
    atoms = at.toatoms()
    atoms.info = parameters

    search_keys = at.key_value_pairs   
    encoding = atoms_to_encode(atoms)        
    fire_name = '{}_top_{}'.format(at.formula,
                                   at.site_index)
    t0 = PyTask(
        func='catkit.flow.fwio.encode_to_atoms', 
        args=[encoding])
        
    t1 = PyTask(func='catkit.flow.fwespresso.get_potential_energy',
                stored_data_varname='trajectory')
        
    firework = Firework([t0, t1],
                        spec={'tags': search_keys,
                              '_priority': 0,
                              'connectivity': d.data.connectivity,
                              'surface_atoms': d.data.surface_atoms},
                        name=fire_name)

    workflow = Workflow([firework], name='L10_top')
    launchpad.add_wf(workflow) 
    '''

print('Number of total jobs: {0}'.format(tot))
print('Number of new jobs: {0}'.format(missing))
