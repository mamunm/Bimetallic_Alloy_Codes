#+BEGIN_SRC python :results output org drawer
#!/usr/bin/env python
#SBATCH --job-name=job
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=stdout
#SBATCH --error=stderr
#SBATCH -p dev

import os
from ase.io import read
from qeio import cd
from ase.db import connect

path = '/scratch/users/mamunm/DATABASE/Binary_Alloy_Project/Ads-Surf'

db = connect('mamun-a3b-top-traj.db')

cnt = 0
for r, d, f in os.walk(path):
    if 'qn.traj' in f:

        with cd(r):
            images = read('qn.traj', ':')

        for i, atoms in enumerate(images):
            tags = {
                'step': i,
                'rstep': len(images) - i - 1,
                'traj': cnt,
                'path': r,
            }

            db.write(atoms, key_value_pairs=tags)

        cnt += 1
#+END_SRC


#+BEGIN_SRC python :results output org drawer
from ase.db import connect
from ase.visualize import view
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm

db = connect('mamun-a3b-top-traj.db')

traj, img0, img1 = [], [], []
for d in db.select('step=0'):
    img0 += [d.toatoms()]
    traj += [d.key_value_pairs['traj']]

for d in db.select('rstep=0'):
    img1 += [d.toatoms()]

reconfig, non_top, moved_site,  = [], [], []
d_xy, d_z, offsite = [], [], []
for i, atoms0 in enumerate(img0):
    atoms1 = img1[i]

    atoms0 *= (5, 5, 1)
    atoms1 *= (5, 5, 1)

    atoms0.wrap()
    atoms1.wrap()

    pos = atoms0.positions
    pos1 = pos[168][:2]

    match = False
    for p in [164, 165, 166, 167]:
        if np.allclose(pos[p][:2], pos1, rtol=1e-2):
            match = p

    if not match:
        non_top += [traj[i]]
        continue

    pos0 = atoms1.positions[match][:2]
    pos1 = atoms1.positions[168][:2]

    d_ads = norm(pos1 - pos0)

    pos0 = atoms0.positions[164:168].T[:2].T
    pos1 = atoms1.positions[164:168].T[:2].T
    dxy = norm(pos1 - pos0)
    d_xy += [dxy]

    pos0 = atoms0.positions[164:168].T[2].T
    pos1 = atoms1.positions[164:168].T[2].T
    dz = norm(pos1 - pos0)
    d_z += [dz]

    offsite += [d_ads]

    # Cutoff values chosen from figure
    if d_ads > 0.9:
        moved_site += [traj[i]]

    if dz > 1.0 or dxy > 1.4:
        reconfig += [traj[i]]

print(non_top)

print('{:.1%} Images not initially on top site'.format(len(non_top) / len(traj)))
print('{:.1%} Adsorbates moved from top site'.format(len(moved_site) / len(traj)))
print('{:.1%} Reconfigurations'.format(len(reconfig) / len(traj)))

plt.figure(figsize=(6, 4))
plt.hist(offsite, bins=200, alpha=0.5, label='adsorbate xy displacement')
plt.hist(d_xy, bins=200, alpha=0.5, label='top layer xy displacement')
plt.hist(d_z, bins=200, alpha=0.5, label='top layer z displacement')
plt.ylim(0, 100)
plt.xlim(0, max(offsite))
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('./images/ads-xy-dist.png')
print('[[./images/ads-xy-dist.png]]')
#+END_SRC

'''
#+RESULTS:
:RESULTS:
[8, 15, 66, 70, 73, 74, 97, 198, 277, 329, 361, 370, 445, 469, 499, 501, 594, 608, 660, 698, 736, 823, 1330, 1350, 1825, 1892, 2249, 2289, 2471, 2599, 2616, 2757, 2839, 2892, 2920, 2953, 3063, 3164, 3455, 3506, 3527, 3535, 3698, 3771]
1.2% Images not initially on top site
36.7% Adsorbates moved from top site
5.5% Reconfigurations
[[./images/ads-xy-dist.png]]
:END:

Some configurations of initial configurations which are greater than 0.01 angstrom  displaced from the top site.
'''
#+BEGIN_SRC python :results output org drawer
from ase.db import connect
from ase.visualize import view
'''
non_top = [8, 15, 66, 70, 73, 74, 97, 198, 277, 329, 361, 370, 445, 469, 499, 501, 594, 608, 660, 698, 736, 823, 1330, 1350, 1825, 1892, 2249, 2289, 2471, 2599, 2616, 2757, 2839, 2892, 2920, 2953, 3063, 3164, 3455, 3506, 3527, 3535, 3698, 3771]
'''
db = connect('mamun-a3b-top-traj.db')

images = []
for i in non_top:
    atoms = db.get_atoms(i)
    images += [d.toatoms()]

write('image.traj', images)
#view(images)
#+END_SRC
