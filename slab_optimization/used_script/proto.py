#!/usr/bin/env python
#SBATCH -p suncat,owners,normal 
#SBATCH --job-name=Al2Fe2_Slab
#SBATCH --output=Al2Fe2_surface.out
#SBATCH --error=Al2Fe2_surface.err
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER@stanford.edu
#SBATCH --ntasks-per-node=16
#SBATCH -C CPU_GEN:HSW
#SBATCH --qos=normal
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
from ase.io import write
from ase.optimize import QuasiNewton
from ase.io import read
from ase.constraints import FixAtoms
from ase.dft.bee import BEEF_Ensemble
import os
import sys

def get_slab(A):
    db = connect('/scratch/users/winther/binary_alloys/data/slab-mamunm.db')
    for i in db.select():
        if i['slab_name'] == ''.join(A.split('2')):
            return i.toatoms()


pw=500.
dw=5000.
xc='BEEF-vdW'
kpts=(6,6,1)
code='Qunatum-Espresso'
facet=(1,1,1)
layer=3
vacuum=10
SBS='L10' #Strukterbericht symbol
supercell='2x2'
fmaxout=0.05
crystal='fcc'

name = 'Fe2Sc2'
Slab = get_slab(name)
Slab.set_initial_magnetic_moments([2.0 for atom in Slab])

calc = espresso(pw=pw, 
                dw=dw,    
                xc=xc,    
                kpts=kpts,
                spinpol=True,
                nbands=-30,
                occupations = 'smearing', #'smeraing', 'fixed', 'tetrahedra'
                smearing='fd', #Fermi-Dirac
                sigma=0.1,
                calcstress=True,
                mode='relax',
                nosym=True,
                beefensemble=True,
                convergence= {'energy':1e-6*13.6,    
                              'mixing':0.1,
                              'nmix':10,
                              'mix':4,
                              'maxsteps':500,
                              'diag':'david'
                              }, 
                psppath='/home/users/vossj/suncat/psp/gbrv1.5pbe',
                output = {'avoidio':False,'removewf':True,'wf_collect':False},
                outdir='calcdir',
                dipole={'status':True}) 

#if qn.traj exists, restart from the last steps
if os.path.exists('qn.traj') and os.path.getsize('qn.traj') != 0:
    Slab=read('qn.traj',index=-1)
from math import sqrt
Slab.set_calculator(calc)          
write('Initial_Slab_with_constraint.traj', Slab)
dyn = QuasiNewton(Slab,logfile='out_scf.log',trajectory='qn.traj')
dyn.run(fmax=fmaxout)
energy = Slab.get_potential_energy()
ens=BEEF_Ensemble(calc)
d_energy=ens.get_ensemble_energies()
calc.stop()

print 'Final SCF energy: %18f [eV]'%(energy)
print 'Deviation in energy: %18f [eV]'%(d_energy.std())

print '\n'

Slab.write('final_' + Slab.get_chemical_formula() + '.traj')
relaxed_atoms = calc.get_final_structure()
psppath='/home/users/vossj/suncat/psp/gbrv1.5pbe'
pseudos = psppath.split('/')[-1]
print 'Potential Energy:', energy, 'eV'
print 'Maximum force:', fmaxout, 'eV/A'
print 'Final Structure:', relaxed_atoms

