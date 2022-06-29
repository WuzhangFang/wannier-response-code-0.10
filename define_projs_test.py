"""
Creates a file which defines on which atoms are the wannier orbitals centered

Takes as an input the wannier90.wout file

Produces some debug output and wann_projs file, which is then read by the fortran code

The format of wann_projs is:
1. number of atoms
2. For each atom, print :
   number of orbitals,
   orbital index
"""

import numpy as np

with open('wannier90.wout') as f:
    wout = f.readlines()

for i in range(len(wout)):
    if 'Cartesian Coordinate (Ang)' in wout[i]:
        struct_start = i

#This read the atomic positions
atoms = []
i = struct_start+2
while True:
    if '---------------------' in wout[i]:
        break
    atoms.append([float(wout[i].split()[j]) for j in range(7,10)])
    i += 1
natoms = i-struct_start-2
atoms = np.array(atoms)

for i in range(len(wout)):
    if 'Final State' in wout[i]:
        final_start = i

#This read the centers of wannier orbitals in the final state
pos = []
i = final_start+1
while True:
    if 'Sum of centres and spreads' in wout[i]:
        break
    pos.append(wout[i].split()[6:9])
    i += 1

nwan = i - final_start-1

for p in pos:
    p[0] = p[0].replace(',','')
    p[1] = p[1].replace(',','')
    for i in range(3):
        p[i] = float(p[i])

pos = np.array(pos)

#for a wannier orbital it finds the nearest atom
#returns(number of the atom, distance)
def find_closest(pos,atoms):
    for i in range(natoms):
        dist = np.linalg.norm(pos-atoms[i])
        if i == 0:
            minr = (0,dist)
        if dist < minr[1]:
            minr = (i,dist)
    return minr

#for a number of atom, it finds the distance to nearest neighbour
def find_nn(iat,atoms):
    for i in range(natoms):
        dist = np.linalg.norm(atoms[iat]-atoms[i])
        if i == iat:
            continue
        if i == 0 or (iat == 0 and i == 1):
            minr = (0,dist)
        if dist < minr[1]:
            minr = (i,dist)
    return minr

#finds how many orbitals are on each atom
counts = np.zeros((natoms))
for p in pos:
    for i in range(natoms):
        if find_closest(p,atoms)[0] == i:
            counts[i] += 1

#finds average and maximum distances of wann orbitals for each atom
projs = [[counts[j]] for j in range(natoms)]
print(projs)
dists = [[0,float(0)] for j in range(natoms)]
for i in range(len(pos)):
    minr = find_closest(pos[i],atoms)
    projs[minr[0]].append(i+1)
    dists[minr[0]][0] += minr[1]/counts[minr[0]]
    if minr[1] > dists[minr[0]][1]:
        dists[minr[0]][1] = minr[1]

print('Atomic positions:')
at = 0
for atom in atoms:
    at += 1
    print(at, " ".join("%9.5f" % f for f in atom))
print('')

print('Nearest neighbour distance for each atom:')
at = 0
for i in range(natoms):
    at += 1
    print(at, "%9.5f" % find_nn(i,atoms)[1])
print('')

print('Average distance of orbitals from the atom center and the maximum distance:')
at = 0
for dist in dists:
    at = at+1
    print(at, " ".join("%10.6f" % f for f in dist))
print('')

#tests if the spin-up and spin-down works properly:
for at in range(natoms):
    for i in range(int(counts[at]/2)):
        if (projs[at][i+1] + nwan/2 not in projs[at]):
            print('!!! Warning !!!')
            print('orb', projs[at][i+1], 'is on atom',at, 'but orb', projs[at][i+1] + nwan/2, 'is not')

print('Number of orbitals for each atom and orbital numbers:')
at = 0
for proj in projs:
    at = at+1
    print(at, proj[0], proj[1:-1])

maxnorb = int(max(counts))
#creates the output file
with open('wann_projs','w') as f:
    f.write(str(natoms)+'\n')
    for proj in projs:
        f.write(str(int(proj[0]))+'\n')
        for p in range(1,len(proj)):
            f.write(str(proj[p])+' ')
        for i in range(len(proj),maxnorb+1):
            f.write('0 ')
        f.write('\n')



