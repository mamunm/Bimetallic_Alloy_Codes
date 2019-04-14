#!/usr/bin/env python
# -*-coding: utf-8 -*-

#!/home/vossj/suncat/bin/python
#SBATCH -p iric,owners,normal
#SBATCH --job-name=Ag3Au_ads_slab
#SBATCH --output=In3Cu_ads_slab.out
#SBATCH --error=In3Cu_ads_slab.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER@stanford.edu
#SBATCH --ntasks-per-node=16
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
"""
Created on Nov 29, 2017
"""
__author__ = "Osman Mamun"
__version__ = "0.1"
__maintainer__ = "Osman Mamun"
__email__ = "mamunm@stanford.edu"
__date__ = "Nov, 2017"

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
########Import section########
import numpy as np
import matplotlib
matplotlib.use('Agg')
from ase.utils.eos import *
from ase.units import kJ
from ase.lattice import surface
from ase import *
from espresso import espresso
from ase.db import connect
from ase.atoms import Atoms
from ase.lattice import surface
from ase.build.supercells import make_supercell
from ase.io import read,write
from ase.optimize import QuasiNewton
from ase.io import read
from ase.constraints import FixAtoms
from ase.dft.bee import BEEF_Ensemble
import os,sys
from math import sqrt
import time
from ase.build import add_adsorbate
from ase.constraints import FixedLine
from ase.data import atomic_numbers
from ase.data import covalent_radii
from copy import deepcopy
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
########Functions############
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#

def __init():
    '''
        This function (main) will find all the unique adsorption site on a
        L10 and L12 catalytic surface. This code is a simplified version
        of a generic algorithm based on the following assumptions:
        1) Given slab contains an irreducible representation of the catalysts slab
        2) Local chemical environment is specified in terms of 1st nearest neighbor
        on the top skin and 2nd neares neighbor in the layer underneath.
    '''

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def top_sites_determination(slab):
    '''
        This function returns all the identified top sites
    '''
    slab_coord = slab.positions
    ts = np.array(slab_coord[ abs(slab_coord[:,2] - max(slab_coord[:,2]))
                  < 1e-1 ], dtype = 'double', order = 'C')
    symbol = np.array(slab.get_chemical_symbols()[-len(ts):],
                                dtype = 'str', order = 'C')
    atomic_number = np.array(slab.get_atomic_numbers()[-len(ts):],
                                dtype = 'str', order = 'C')
    x = np.array(ts[:,0], dtype = 'double', order = 'C')
    y = np.array(ts[:,1], dtype = 'double', order = 'C')
    z = np.array(ts[:,2], dtype = 'double', order = 'C')
    return [symbol, atomic_number, x, y, z]

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def get_supercell(slab):
    '''
        This function returns a 3 by 3 supercell of the original cell
    '''
    P = [[3, 0, 0], [0, 3, 0], [0, 0, 1]]
    supercell = make_supercell(slab, P)
    slab_cell = slab.cell
    x = np.array(supercell.positions[:,0]
                 - (slab_cell[0][0] + slab_cell[1][0]),
                 dtype = 'double', order = 'C')
    y = np.array(supercell.positions[:,1]
                 - (slab_cell[1][1]),
                 dtype = 'double', order = 'C')
    z = np.array(supercell.positions[:,2],
                 dtype = 'double', order = 'C')
    symbol = np.array(supercell.get_chemical_symbols(),
                      dtype = 'str', order = 'C')
    atomic_number = np.array(supercell.get_atomic_numbers(),
                             dtype = 'str', order = 'C')
    return [symbol, atomic_number, x, y, z]

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def distance(A, B):
    '''
        This function returns the distance between two atoms
    '''
    return np.sqrt((A[2] - B[2])**2 + (A[3] - B[3])**2 + (A[4] - B[4])**2)

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def get_neighborlist(A, supercell, cell_dim):
    '''
        This function returns the neighborlist around a particular top sites.
        It scans for neighboring metal atoms within one nearest neighbor in the
        lateral direction and one vertical layer in the negative z direction.
    '''
    neighborlist = [[] for _ in range(5)]
    for i in range(len(supercell[0])):
        B = [supercell[j][i] for j in range(len(supercell))]
        d = distance(A, B)
        if d < 0.05 * cell_dim[0][0]:
            for j in range(5):
                neighborlist[j].append(supercell[j][i])
    for i in range(len(supercell[0])):
        B = [supercell[j][i] for j in range(len(supercell))]
        d = distance(A, B)
        if not d < 0.05 * cell_dim[0][0]:
            if abs(A[4] - B[4]) < 0.1 and d < 0.51 * cell_dim[0][0]:
                for j in range(5):
                    neighborlist[j].append(supercell[j][i])
            layer_dist =  (cell_dim[0][0]/ np.sqrt(6.0))
            condition1 = A[4] - B[4] > 0.5 * layer_dist
            condition2 = d < 1.05 * cell_dim[0][0] / np.sqrt(2.0)
            if condition1 and condition2:
                for j in range(5):
                    neighborlist[j].append(supercell[j][i])

    for i in range(5):
        if i == 0 or i == 1:
            neighborlist[i] = np.array(neighborlist[i],
                                       dtype = 'str', order = 'C')
        else:
            neighborlist[i] = np.array(neighborlist[i],
                                       dtype = 'double', order = 'C')
    return neighborlist

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def write_neighborlist( neighborlist, index):
    '''
        This function writes the xyz coordinate file of a particular neighborlist
    '''
    if not os.path.exists('Neighbor_Cluster'):
        os.makedirs('Neighbor_Cluster')
    filename = 'Neighbor_Cluster/neighborlist-' + str(index + 1) + '.xyz'
    with open(filename, 'w') as ff:
        ff.write('%6s \n' %str(len(neighborlist[0])))
        ff.write('system \n')
        for i in range(len(neighborlist[0])):
            symbol = neighborlist[0][i]
            x = neighborlist[2][i]
            y = neighborlist[3][i]
            z = neighborlist[4][i]
            ff.write('%4s  %10.2E  %10.2E   %10.2E \n' %(symbol, x, y, z))

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def centralized(A):
    '''
        This function returns the translated coordinates of a neighborlist with
        respect to the central atom
    '''
    M = deepcopy(A)
    trans = [M[i][0] for i in range(2,5)]
    for i in range(len(M[0])):
        for j in range(2,5):
            M[j][i] = M[j][i] - trans[j-2]
    return M

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def rotate(A, theta):
    '''
        This function rotates the neighborlist at a given angle
    '''
    N = deepcopy(A)
    cos = np.cos(theta*np.pi/180)
    sin = np.sin(theta*np.pi/180)
    for i in range(len(N[0])):
        N[2][i] = A[2][i] * cos - A[3][i] * sin
        N[3][i] = A[2][i] * sin + A[3][i] * cos
    return N

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def compare_struct(A, B):
    '''
        This function compares two neighborlist irrespective of their order.
        It returns 1 if the structures are same and 0 if not same.
    '''
    O = deepcopy(A)
    P = deepcopy(B)
    returnval = 0
    bin_val = np.zeros((len(O[0])))
    for i in range(len(O[0])):
        for j in range(len(P[0])):
            if O[0][i] == P[0][j]:
               cond2 = abs(O[2][i] - P[2][j]) < 1e-4
               cond3 = abs(O[3][i] - P[3][j]) < 1e-4
               cond4 = abs(O[4][i] - P[4][j]) < 1e-4
               if cond2 and cond3 and cond4:
                   bin_val[i] = 1
    if not 0 in bin_val:
        returnval = 1
    return returnval

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def compare_rotate(A, B):
    '''
        This function compares different rotated positions
        of neighborlist B to unrotated A. If any of the rotated
        B is same as A, it will return False.
    '''
    Q = deepcopy(A)
    R = deepcopy(B)
    returnval_2 = True
    theta = np.arange(0, 365, 5)
    bin_val_2 = np.zeros(len(theta))
    for i in range(len(theta)):
        RR = rotate(R, theta[i])
        bin_val_2[i] = compare_struct(Q, RR)
        del RR
        if 1 in bin_val_2:
            returnval_2 = False
    return returnval_2

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def diff_site(A, B):
    '''
        This function compares neighborlist of A and B, and
        determines wheather they are same or different.
        Returns false, if same, and true, if different.
    '''
    S = deepcopy(A)
    T = deepcopy(B)
    returnval_3 = True
    if S[0][0] == T[0][0] and len(S[0]) == len(T[0]):
        SS = centralized(S)
        TT = centralized(T)
        returnval_3 = compare_rotate(SS, TT)
    return returnval_3

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def unique_top_site(topsite, supercell, cell_dim):
    '''
        This function returns the coordinates of the identified unique top sites.
    '''
    neighborlist = [[] for _ in range(len(topsite[0]))]
    uniquetopsite = [[] for _ in range(len(topsite))]
    for i in range(len(topsite[0])):
        atomos = [topsite[j][i] for j in range(len(topsite))]
        neighborlist[i] = get_neighborlist(atomos, supercell, cell_dim)
        write_neighborlist(neighborlist[i], i)
    omit = []
    for i in range(len(topsite[0])):
        if i not in omit:
            for j in range(len(topsite[0])):
                if i < j:
                    res = diff_site(neighborlist[i], neighborlist[j])
                    if not res:
                        omit.append(j)
    for i in range(len(topsite[0])):
        if i not in omit:
            for j in range(5):
                uniquetopsite[j].append(topsite[j][i])
    for i in range(len(topsite[0])):
        write_neighborlist(neighborlist[i], i+4)
    for key in range(5):
        uniquetopsite[key] = np.array(uniquetopsite[key])
    return uniquetopsite
'''
    cos = np.cos(120*np.pi/180)
    sin = np.sin(120*np.pi/180)
    PPPP = centralized(neighborlist[0])
    QQQQ = centralized(neighborlist[1])
    RRRR = rotate(QQQQ, cos, sin)
    write_neighborlist(PPPP, 100)
    write_neighborlist(QQQQ, 101)
    write_neighborlist(RRRR, 102)
    print('neighborlist PPPP:', PPPP)
    print('neighborlist QQQQ:', QQQQ)
    print('neighborlist RRRR:', RRRR)
    compare_struct(PPPP, RRRR)
'''

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def main_site_determination(input_traj_file_name, WO =False):
    #input_traj_file_name = 'Initial_Slab_with_constraint.traj'
    slab = read(input_traj_file_name)
    slab_coord = slab.positions
    cell_dim = slab.cell
    supercell = get_supercell(slab)
    topsite = top_sites_determination(slab)
    uts = unique_top_site(topsite, supercell, cell_dim)
    if WO:
        f.write('Site specification:\n')
        f.write('\nTop sites:\n')
        f.write('   %-10s %-20s %-25s %-25s %-25s %-25s\n'
                %('Index', 'chemical_symbol', 'atomic_number', 'x', 'y', 'z'))
        for i in range(len(topsite[0])):
            f.write(' %-15s %-20s %-15s %-25s %-25s %-25s\n'
                    %('top_site_' + str(i + 1), topsite[0][i], topsite[1][i],
                    str(topsite[2][i]), str(topsite[3][i]), str(topsite[4][i])))
        f.write('\nUnique top sites:\n')
        f.write('   %-20s %-20s %-28s %-25s %-25s %-25s\n'
                %('Index', 'chemical_symbol', 'atomic_number', 'x', 'y', 'z'))
        for i, sym in enumerate(uts[0]):
            f.write(' %-25s %-20s %-15s %-25s %-25s %-25s\n'
                    %('unique_top_site_' + str(i + 1), uts[0][i], uts[1][i],
                    str(uts[2][i]), str(uts[3][i]), str(uts[4][i])))
        f.write('\n')
    return uts

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def get_height(adsorbate, metal):
    '''
        This function returns the arithmetic average of the covalent radii of an adsorbate-metal system
    '''
    ads_ato_num = atomic_numbers[adsorbate]
    ads_cov_rad = covalent_radii[ads_ato_num]
    metal_ato_num = atomic_numbers[metal]
    metal_cov_rad = covalent_radii[metal_ato_num]
    return ads_cov_rad + metal_cov_rad

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
def run_geo_opt(bulk_name, final_traj_name, site_coord):
    '''
        This function runs the DFT calculations for the given trajectory name. It puts the adsorbate on the
        passed site coordinate.
    '''
    time2= time.time()
    f.write('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    f.write('Starting calculations: %30s\n'
            %(final_traj_name.split('.')[0]))
    f.write('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

    db=connect('/home/mamunm/DATABASE/adsorbate-slab.db')
    ######Calculations parameters#######
    pw=500.
    dw=5000.
    xc='BEEF-vdW'
    kpts=(6,6,1)
    code='Qunatum-Espresso'
    facet=(1,1,1)
    layer=3
    vacuum=10
    a = final_traj_name.split('_')
    adsorbate = a[1]
    site = a[2] + '_' + a[3] + '_' + a[4]
    slab_name = a[5].split('.')[0]
    height = get_height(adsorbate, a[4])
    f.write('\nComputed initial adsorption height: %4.2f Angstrom\n' %(height))
    if '6' in slab_name:
        SBS_symbol = 'L10' #Strukterbericht symbol
        PS_symbol = 'tP2'  #Pearson symbol
    else:
        SBS_symbol = 'L12' #Strukterbericht symbol
        PS_symbol = 'cP4'  #Pearson symbol

    supercell='2x2'
    fmaxout=0.05
    crystal='fcc'
    pol_metal = ['Fe', 'Ni', 'Co', 'Mn']
    spinpol = False
    DB_dir = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Surface/'

    if '6' in slab_name:
        metal1 = slab_name.split('6')[0]
        metal2 = slab_name.split('6')[1]
        path_to_slab = DB_dir + '/AB/'
        if os.path.exists(path_to_slab + bulk_name):
            slab = read(path_to_slab + bulk_name +
            '/Initial_Slab_with_constraint.traj')
        if metal1 in pol_metal or metal2 in pol_metal:
            spinpol = True
            slab.set_initial_magnetic_moments([2.0 if atom.symbol in pol_metal else 0 for atom in slab])
    else:
        metal1 = slab_name.split('9')[0]
        metal2 = slab_name.split('9')[1].split('3')[0]
        path_to_slab = DB_dir + '/A3B/'
        if os.path.exists(path_to_slab + bulk_name):
            slab = read(path_to_slab + bulk_name +
                   '/Initial_Slab_with_constraint.traj')
        if metal1 in pol_metal or metal2 in pol_metal:
            spinpol = True
            slab.set_initial_magnetic_moments([2.0 if atom.symbol in pol_metal else 0 for atom in slab])

    calc = espresso(pw = pw,
                    dw = dw,
                    xc = xc,
                    kpts = kpts,
                    spinpol = spinpol,
                    nbands = -10,
                    occupations = 'smearing', #'smeraing', 'fixed', 'tetrahedra'
                    smearing = 'fd', #Fermi-Dirac
                    sigma = 0.15,
                    calcstress = True,
                    mode = 'relax',
                    nosym = True,
                    beefensemble = True,
                    convergence= {'energy' : 1e-5*13.6,
                                  'mixing' : 0.3,
                                  'nmix' : 10,
                                  'mix' : 4,
                                  'maxsteps' : 2000,
                                  'diag' : 'david'
                                 },
                    psppath='/home/vossj/suncat/psp/gbrv1.5pbe',
                    output = {'avoidio':False,'removewf':True,'wf_collect':False},
                    outdir='Calculations/'+ adsorbate + '_' + site + '/calcdir',
                    dipole={'status':True})

    #Add adsorbate on the surface
    add_adsorbate(slab, adsorbate, height, site_coord)
    slab.write('Initial_slab_with_adsorbate.traj')
    #if qn.traj exists, restart from the last steps
    if os.path.exists('qn.traj') and os.path.getsize('qn.traj') != 0:
        slab=read('qn.traj',index=-1)
    slab.set_calculator(calc)
    dyn = QuasiNewton(slab,logfile='out_scf.log',trajectory='qn.traj')
    dyn.run(fmax=fmaxout)
    energy = slab.get_potential_energy()
    ens=BEEF_Ensemble(calc)
    d_energy=ens.get_ensemble_energies()
    calc.stop()

    f.write('Final SCF energy: %18f [eV]\n'%(energy))
    f.write('Deviation in energy: %18f [eV]\n'%(d_energy.std()))

    slab.write(final_traj_name)
    relaxed_atoms = calc.get_final_structure()
    psppath='/home/vossj/suncat/psp/gbrv1.5pbe'
    pseudos = psppath.split('/')[-1]
    f.write('Potential Energy: %2.10E eV\n' %(energy))
    f.write('Maximum force: %2.6f eV/A\n' %(fmaxout))
    f.write('Final Structure: %-300s\n' %(str(relaxed_atoms)))
    db.write(relaxed_atoms,
             crystal = crystal,
             epot = energy,
             supercell = supercell,
             layer = layer,
             SBS_symbol = SBS_symbol,
             vacuum = vacuum,
             facet = str(facet),
             PS_symbol = PS_symbol,
             slab_name = slab_name,
             adsorbate = adsorbate,
             site = site,
             name=str(slab.get_chemical_formula()),
             fmaxout = fmaxout,
             xc = xc,
             pw = pw,
             dw = dw,
             Ensemble_energy = d_energy.std(),
             code = code,
             psp = pseudos)
    f.write('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    f.write("Elapsed time for %20s :: %s seconds\n"
            %(final_traj_name.split('.')[0], time.time()-time2))
    f.write('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
################Main#################
start_from = None
time_so_far = 0
if os.path.exists('Output.txt'):
    with open('Output.txt', 'r') as f:
        for line in f:
            if 'Total elapsed time' in line:
                sys.exit('This calculation was succesfully completed.')
            else:
                if 'Starting calculations:' in line:
                    start_from = line.split()[2]
                if 'Elapsed time for' in line:
                    time_so_far += float(line.split()[5])

#Change this two parameter only
system_class = 'L12'
bulk_name = 'In3Cu'

if '2' in bulk_name:
    metal1 = bulk_name.split('2')[0]
    metal2 = bulk_name.split('2')[1]
    slab_name = metal1 + '6' + metal2 + '6'
else:
    metal1, metal2 = bulk_name.split('3')
    slab_name = metal1 + '9' + metal2 + '3'
crystal = 'fcc'
if system_class == 'L10':
    DB_dir = ('/scratch/users/mamunm/DATABASE/' +
              'Binary_Alloy_Project/Surface/AB/')
else:
    DB_dir = ('/scratch/users/mamunm/DATABASE/' +
              'Binary_Alloy_Project/Surface/A3B/')

input_traj_dir = DB_dir + bulk_name + '/Initial_Slab_with_constraint.traj'

with open('Output.txt', 'a+', 0) as f:
    if not start_from == None:
        start_time=time.time()
        uts = main_site_determination(input_traj_dir, WO = False)
        root_dir = os.getcwd()
        if not os.path.exists('Calculations'):
            os.makedirs('Calculations/')
        sub_dir = root_dir + '/' + 'Calculations/'
        os.chdir(sub_dir)
        adsorbates = ['H', 'O', 'N', 'C', 'S']
        for i, adsorbate in enumerate(adsorbates):
                for j, site in enumerate(uts[0]):
                    calc_name = adsorbate + "_top_site_" + uts[0][j]
                    if not os.path.exists(calc_name):
                        os.makedirs(calc_name)
                    os.chdir(calc_name)
                    final_traj_name = "final_" + calc_name + "_" + slab_name + ".traj"
                    if not os.path.exists(final_traj_name):
                        site = (uts[2][j], uts[3][j])
                        run_geo_opt(bulk_name, final_traj_name, site)
                    os.chdir(sub_dir)

    else:
        start_time=time.time()
        f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                '+++++++++++++++++\n')
        f.write('Specification of the calculation:\n')
        f.write('Name: Adsorbate on slab calculation\n')
        lt = time.localtime()
        f.write('Start time: {0:02d}/{1:02d}/{2:04d} at {3:02d}:{4:02d}:{5:02d}\n'
                .format(lt[1],lt[2],lt[0],lt[3],lt[4],lt[5]))
        f.write('Initializing calculation with the BEEF-vdw functional using Quantum '
                'Espresso Package.\n')
        f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                '+++++++++++++++++\n')
        nearest_neighbor = 1
        vertical_search_layer = 1
        f.write('\nIdentifying unique top sites based'
                ' on the local chemical environment.\n')
        f.write('Local chemical environment specifications:\n')
        f.write('Nearest neighbor considered in the horizontal direction : {0:1d}\n'
                .format(nearest_neighbor))
        f.write('Vertical depth considered: {0:1d}\n'
                .format(vertical_search_layer))
        f.write('\nNote to self:\n    I used spglib symmetry'
                ' operations on both L10 and L12 and found\n'
                'two unique top sites for both of them.\n\n')
        uts = main_site_determination(input_traj_dir, WO = True)
        root_dir = os.getcwd()
        if not os.path.exists('Calculations'):
            os.makedirs('Calculations/')
        sub_dir = root_dir + '/' + 'Calculations/'
        os.chdir(sub_dir)
        adsorbates = ['H', 'O', 'N', 'C', 'S']
        for i, adsorbate in enumerate(adsorbates):
            for j, site in enumerate(uts[0]):
                calc_name = adsorbate + "_top_site_" + uts[0][j]
                if not os.path.exists(calc_name):
                    os.makedirs(calc_name)
                os.chdir(calc_name)
                final_traj_name = "final_" + calc_name + "_" + slab_name + ".traj"
                if not os.path.exists(final_traj_name):
                    site = (uts[2][j], uts[3][j])
                    run_geo_opt(bulk_name, final_traj_name, site)
                os.chdir(sub_dir)

    f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    f.write("Everything done successfully.\nSee .out file for the results.\n")
    f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    f.write("      Total elapsed time :: %s seconds\n" %(time.time()-start_time + time_so_far))
    f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    lt = time.localtime()
    f.write('End time: {0:02d}/{1:02d}/{2:04d} at {3:02d}:{4:02d}:{5:02d}\n'
                .format(lt[1],lt[2],lt[0],lt[3],lt[4],lt[5]))
    f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
