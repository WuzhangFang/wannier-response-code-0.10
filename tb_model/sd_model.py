"""
This is a program for constructing an sd tight-binding models.

The model is stored in class Model, which also provides function
for writing the Hamiltonian to a tb_hr.dat file in a format that
can be read by the fortran linear response code.

The model is defined in input file model.in, which has an .ini
syntax.

This program when it is run reads the model.in file and writes the
tb_hr.dat file, but it can also be used as a library.
"""


from backports import configparser
import itertools
import sympy as sp
from sympy import sympify as spf
from sympy import I
from sympy.physics.matrices import msigma
from sympy.physics.quantum import TensorProduct
from sympy.functions import re,im,exp

from find_nn import find_nn2

class Model:
    """
    This class stores the model.

    The model can be read from an input file using read_input function.
    """

    def __init__(self,B,atoms,t,J,mags,dim):
        self.B = B
        self.atoms = atoms
        self.t = t
        self.n_atoms = atoms.shape[0]
        self.J = J
        self.mags = mags
        self.dim = dim

        self.debug = False

    def find_nn(self):
        """
        Finds nearest neighbours.

        Only finds the nearest neighbourse, but could be easily modified for
        next to nearest etc.

        Returns:
            nn: A list, which for each atom contains a list of all nearest neighbours.
                Nearest neighbours have format:
                [atom_index,(n0,n1),distance]
        """

        dim = self.dim

        ncells = 2

        B = self.B
        atoms = self.atoms

        prec = 1e-5

        def find_position(n,atom_n):
            """
            Finds position of atom with index atom_n in unit cell n0,n1
            """

            #first find position of the atom in the 1. unit cell
            position = sp.zeros(1,dim)
            for i in range(dim):
                for j in range(dim):
                    position[j] += B[i,j]*atoms[atom_n,i]
            #now move to correct unit cell:
            for j in range(dim):
                for i in range(dim):
                    position[j] += B[i,j]*n[i]
            return position

        def find_distn(neighbs,prev_dist):
            """
            Finds the nearest neighbour distance, when prev_dist=0.
            When prev_dist != 0 then finds the nearest neighbour distance
            that is larger than prev_dist. This could be used for next nearest...
            """

            min_dist = sp.oo #at the beginning we set the min_dist to infinity
            for n in neighbs:
                if n[2] < min_dist and n[2] > prev_dist:
                    min_dist = n[2]

            return min_dist

        nn = []
        iter1 = range(-ncells,ncells+1,1)
        iterlist = []
        for i in range(dim):
            iterlist.append(iter1)
        cells_iter = list(itertools.product(*iterlist))
        for a in range(self.n_atoms):
            neighbs = []
            for n in cells_iter:
                for a2 in range(self.n_atoms):
                    dist = find_position(n,a2)-find_position([0]*dim,a)
                    neighbs.append([a2,n,dist.norm()])

            min_dist = find_distn(neighbs,0)
            nn_a = []
            for ne in neighbs:
                if ne[2] == min_dist:
                    nn_a.append(ne)

            nn.append(nn_a)

        if self.debug:
            for nes in nn:
                for n in nes:
                    print n
                print ''

        return nn

    def create_Hr(self):
        """
        This creates a tight-binding Hamiltonian in a similar format to that of wannier90.

        Returns: Hr: format is a dictionary, where each key has a format:
            (n0,n1,i,j), where n0,n1 are the cell coordinates and i,j are orbital numbers
            in this cases, there are only two orbitals per atom, spin-up orbital on atom i
            has orbital number i, spin-down i+n_atoms
            Hr[n0,n1,i,j] is <0,0,i|H|n0,n1,j>
        """

        nn = self.find_nn()
        Hr = {}
        for i in range(self.n_atoms):
            for neg in nn[i]:
                cell = [0]*self.dim
                for j in range(self.dim):
                    cell[j] = neg[1][j]
                cell = tuple(cell)
                key = (cell,i+1,neg[0]+1)
                if key in Hr:
                    Hr[key] += -self.t
                else:
                    Hr[key] = -self.t
                key2 = (cell,i+1+self.n_atoms,neg[0]+1+self.n_atoms)
                if key2 in Hr:
                    Hr[key2] += -self.t
                else:
                    Hr[key2] = -self.t

        #this constructs the spin-operators projected on individual atoms
        ms = sp.zeros(self.n_atoms*2)
        for at in range(self.n_atoms):
            proj = sp.zeros(self.n_atoms)
            proj[at,at] = 1
            for i in range(3):
                ms += self.mags[at,i]*TensorProduct(msigma(i+1),proj)

        for i in range(self.n_atoms*2):
            for j in range(self.n_atoms*2):
                cell = (0,)*self.dim
                key = (cell,i+1,j+1)
                if key in Hr:
                    Hr[key] += self.J*ms[i,j]
                else:
                    Hr[key] = self.J*ms[i,j]

        return Hr

    def write_wann_Hr(self):
        """
        This prints the Hr Hamiltonian in a file wannier90_hr.dat in the same format
        as wannier90 uses.

        !!!This is not the correct way how to do it and is left here only as a reference.!!!
        """

        Hr = self.create_Hr()
        rs = []
        for key in Hr:
            r = (key[0]) 
            if r not in rs: 
                rs.append(r)
        n_rs = len(rs)

        with open('wannier90_hr.dat','w') as f:
            f.write('\n')
            f.write(str(self.n_atoms*2)+'\n')
            f.write(str(n_rs)+'\n')
            nrlines = n_rs/15+1
            for l in range(nrlines):
                for i in range(15):
                    if l*15 + i < n_rs:
                        f.write('1  ')
                f.write('\n')
            for key in sorted(Hr):
                f.write(str(key[0][0])+'  ')
                try:
                    f.write(str(key[0][1])+'  ')
                except:
                    f.write(str(0)+'  ')
                try:
                    f.write(str(key[0][2])+'  ')
                except:
                    f.write(str(0)+'  ')

                f.write(str(key[1])+'  ')
                f.write(str(key[2])+'  ')

                f.write(str(float(re(Hr[key])))+ '  ')
                f.write(str(float(im(Hr[key]))))
                
                f.write('\n')

        with open('POSCAR','w') as f:
            f.write('\n')
            f.write('1\n')
            Bout = sp.eye(3)
            Bout[0:self.dim,0:self.dim] = self.B[:,:]
            for i in range(3):
                for j in range(3):
                    f.write(str(float(Bout[i,j]))+'  ')
                f.write('\n')

    def write_Hr(self):
        """
        This prints the Hamiltonian in the tight-binding basis in a format
        that can be easily used for transforming to a k-dependent Hamiltonian.
        It writes a tb_hr.dat file.
        """

        Hr = self.create_Hr()

        with open('tb_hr.dat','w') as f:
            f.write('\n')
            f.write(str(self.n_atoms*2)+'\n')
            f.write(str(len(Hr))+'\n')
            for key in sorted(Hr):
                R = sp.zeros(1,3)
                at_a = (key[2]-1) % self.n_atoms
                at_b = (key[1]-1) % self.n_atoms
                for i in range(3):
                    R[i] += key[0][i] - self.atoms[at_b,i]  + self.atoms[at_a,i]
                    f.write(str(float(R[i]))+'  ')
                
                f.write(str(key[1])+'  ')
                f.write(str(key[2])+'  ')

                f.write(str(float(re(Hr[key])))+ '  ')
                f.write(str(float(im(Hr[key]))))
                
                f.write('\n')

        with open('POSCAR','w') as f:
            f.write('\n')
            f.write('1\n')
            Bout = sp.eye(3)
            Bout[0:self.dim,0:self.dim] = self.B[:,:]
            for i in range(3):
                for j in range(3):
                    f.write(str(float(Bout[i,j]))+'  ')
                f.write('\n')
    
    def create_Hk(self,k):

        Hr = self.create_Hr()
        Hk = sp.zeros(self.n_atoms*2)

        for key in Hr:

            R = sp.zeros(1,3)
            for i in range(3):
                R += self.B[i,:]*key[0][i] 

            pa = sp.zeros(1,3)
            pb = sp.zeros(1,3)
            for i in range(3):
                at_a = (key[2]-1) % self.n_atoms
                at_b = (key[1]-1) % self.n_atoms
                pa += self.B[i,:]*self.atoms[at_a,i] 
                pb += self.B[i,:]*self.atoms[at_b,i] 

            rho = R + pa - pb
            f = k.dot(R)
            Hk[key[1]-1,key[2]-1] += exp(I*k.dot(rho))*Hr[key]

        return Hk


def read_input(inp_file):
    """
    Reads an input file and creates a model object.
    """

    config = configparser.ConfigParser()
    config.read(inp_file)

    dim = int(config['structure']['dim'])

    latt = config['structure']['lattice'].split('\n')
    B = sp.zeros(dim)
    for i in range(dim):
        for j in range(dim):
            B[i,j] = spf(latt[i+1].split()[j])

    B = B * spf(config['structure']['lattice_scale'])

    n_atoms = int(config['structure']['n_atoms'])
    atoms = sp.zeros(n_atoms,dim)
    at = config['structure']['atoms'].split('\n')
    for i in range(n_atoms):
        for j in range(dim):
            atoms[i,j] = spf(at[i+1].split()[j])

    t = spf(config['tb']['hopping'])
    J = spf(config['tb']['exchange'])
    mags = sp.zeros(n_atoms,3)
    ms = config['tb']['mag_moms'].split('\n')
    for i in range(n_atoms):
        for j in range(3):
            mags[i,j] = spf(ms[i+1].split()[j])

    #we normalize the magnetic moments
    for at in range(n_atoms):
        mags[at,:] = mags[at,:] / mags[at,:].norm()

    model = Model(B,atoms,t,J,mags,dim)

    return model

model = read_input('model.in')

model.write_Hr()
