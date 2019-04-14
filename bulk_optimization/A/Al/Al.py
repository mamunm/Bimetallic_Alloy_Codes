#!/home/vossj/suncat/bin/python
#SBATCH -p iric,owners,normal 
#SBATCH --job-name=Al_bulk
#SBATCH --output=bulk.out
#SBATCH --error=bulk.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER@stanford.edu
#SBATCH --ntasks-per-node=16

import numpy as np    
import matplotlib
matplotlib.use('Agg') 
from ase.utils.eos import * 
from ase.units import kJ
from ase.lattice import bulk
from ase import *
from espresso import espresso
from ase.db import connect

db = connect('/home/mamunm/DATABASE/bulk.db')
metal =  'Al'
metal2 =  None

a=3.95     
strains = np.linspace(0.6, 1.40, 50) 
#strains = np.linspace(0.9, 1.1, 3) 
crystal = 'fcc'

pw=800.
dw=8000.
xc='BEEF-vdW'
kpts=(12,12,12)
code='Qunatum-Espresso'

volumes = []  
energies = []
print '\n'
print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print '                   Equation of state data                  '
print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print 'Lattice constant        Energy [eV]          Volume [AA^3] '
print '-----------------------------------------------------------'
for i in strains: 
    if metal2:
      atoms = bulk(metal, crystal, a=a*i, cubic=True)
      atoms.set_chemical_symbols(metal+'2'+metal2+'2')
    else:
      atoms = bulk(metal, crystal, a*i, cubic=True)
    atoms.set_pbc((1,1,1))                
    calc = espresso(pw=pw, 
                dw=dw,    
                xc=xc,    
                kpts=kpts,
                nbands=-20,
                sigma=0.1,
                mode='scf',
                convergence= {'energy':1e-6,    
                              'mixing':0.1,
                              'nmix':10,
                              'mix':4,
                              'maxsteps':500,
                              'diag':'david'
                              }, 
                psppath='/home/vossj/suncat/psp/gbrv1.5pbe',
                output = {'avoidio':True,'removewf':True,'wf_collect':False},
                outdir='calcdir') 

    atoms.set_calculator(calc)          
    volumes.append(atoms.get_volume())    
    energy=atoms.get_potential_energy()   
    energies.append(energy)               
    print '%16f %18f %22f'%(a*i,energy,atoms.get_volume())

print '\n'
eos = EquationOfState(volumes, energies) 
v0, e0, B = eos.fit()                    

if metal2:
  best_a = (v0)**(1./3.) 
  atoms = bulk(metal, crystal, a=best_a, cubic=True)
  atoms.set_chemical_symbols(metal+'2'+metal2+'2')
else:
  best_a = (v0)**(1./3.) 
  atoms = bulk(metal, crystal, best_a, cubic=True)

print 'Lattice constant:', best_a, 'Angstrom'
print 'Bulk modulus:', B / kJ * 1e24, 'GPa'
print '(Fitted) total energy at equilibrium latt. const.:', e0, 'eV'
eos.plot(atoms.get_chemical_formula()+'-eos.png')  



atoms.set_pbc((1,1,1))
calc = espresso(pw=pw,
                dw=dw,
                xc=xc,
                mode='scf',
                opt_algorithm='bfgs',
                beefensemble=True,
                kpts=kpts,
                nbands=-20,
                sigma=0.1,
                convergence= {'energy':1e-7,
                              'mixing':0.1,
                              'nmix':10,
                              'mix':4,
                              'maxsteps':500,
                              'diag':'david',
                              'mixing_mode':'local-TF'
                              },
                psppath='/home/vossj/suncat/psp/gbrv1.5pbe',
                output = {'avoidio':True,'removewf':True,'wf_collect':False},
                outdir='calcdir')

atoms.write('final_' + atoms.get_chemical_formula() + '.traj')
atoms.set_calculator(calc) 
forces = atoms.get_forces(atoms)
epot = atoms.get_potential_energy()
ens = calc.get_nonselfconsistent_energies()
relaxed_atoms = calc.get_final_structure()
fmaxout = np.sqrt((forces**2).sum(axis=1).max())
psppath='/home/vossj/suncat/psp/gbrv1.5pbe'
pseudos = psppath.split('/')[-1]
a=atoms.get_cell_lengths_and_angles()[0],
b=atoms.get_cell_lengths_and_angles()[1],
c=atoms.get_cell_lengths_and_angles()[2],
print 'Forces on each atom:'
print forces
print 'Potential Energy:', epot, 'eV'
print 'Optimized lattice constant (a):', a, 'Angstrom'
print 'Optimized lattice constant (b):', b, 'Angstrom'
print 'Optimized lattice constant (c):', c, 'Angstrom'
print 'Maximum force:', fmaxout, 'eV/A'
print 'Final Structure:', relaxed_atoms
print 'Esemble Energy:', ens, 'eV'
db.write(relaxed_atoms,
         epot=epot,
         species=str(atoms.get_chemical_formula()),
         fmaxout=fmaxout,
         xc=xc,
         pw=pw,
         dw=dw,
         code=code,
         lattice_constant_a=a[0],
         lattice_constant_b=b[0],
         lattice_constant_c=c[0],
         psp=pseudos,
         data={'BEEFens': ens})
