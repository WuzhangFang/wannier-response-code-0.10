"""
This is used for reading the output from the calculation in a nympy format.

Can be used interactively in ipython, allows for easy analysis and plotting...

The main subroutine is read_result, which reads the results and stores it as 
a class result.

Examples:
    To get a tensor with gamma=0.001, formula number 2, projection 1 and nk=100**3 do:

        res.get_result((0.001,2,1,100**3))

    For a spin-Hall effect there are no projections so do:

        res.get_result((0.001,2,100**3))
"""
import numpy as np
class result:
    def __init__(self,npars):
        self.dat = {}
        self.pars = []
        self.npars = npars

    def add_tensor(self,pars,tensor):
        if  (type(pars) != tuple) or (len(pars) != self.npars):
            raise TypeError
        self.dat[pars] = tensor
        if pars in self.pars:
            print 'Tensor for this pars is already present!'
        else:
            self.pars.append(pars)

    def get_tensor(self,pars):
        return self.dat[pars]

def read_output(dire='.',typ='cisp'):

    if typ == 'cisp':
        with open(dire+'/linres_out') as f:
            lines = f.readlines()
    if typ == 'she':
        with open(dire+'/she_out') as f:
            lines = f.readlines()

    for i in range(len(lines)):
        if 'Unformatted output:' in lines[i]:
            start = i+1

    for i in range(start,len(lines)):
        linc = lines[i].split()
         
        if typ == 'cisp':
            for j in [0,1,3,4]:
                linc[j] = int(linc[j])
            for j in [2,5]:
                linc[j] = float(linc[j])
        if typ == 'she':
            for j in [0,1,2,4]:
                linc[j] = int(linc[j])
            for j in [3,5]:
                linc[j] = float(linc[j])

        lines[i]=linc

    return lines, start

def read_cisp_tensor(gamma,nfor,nproj,dire='.'):
    X = np.zeros((3,3))
    lines, start = read_output(dire,typ='cisp')
    
    for i in range(start,len(lines)):
        if abs(lines[i][2]-gamma) < 1e-10 and lines[i][3] == nfor and lines[i][4] == nproj:
            X[lines[i][0]-1,lines[i][1]-1] = lines[i][5]

    return X

def read_she_tensor(gamma,nfor,dire='.'):
    X = np.zeros((3,3,3))
    lines, start = read_output(dire,typ='she')

    for i in range(start,len(lines)):
        if abs(lines[i][3]-gamma) < 1e-10 and lines[i][4] == nfor:
            X[lines[i][0]-1,lines[i][1]-1,lines[i][2]-1] = lines[i][5]

    return X

def find_pars(dire='.',typ='cisp'):

    lines, start = read_output(dire,typ=typ)

    gams = []
    fors = []
    projs = []

    if typ == 'cisp':
        gam_ind = 2
        fors_ind = 3
    if typ == 'she':
        gam_ind = 3
        fors_ind = 4

    for i in range(start,len(lines)):
        #I believe there should be no problems with rounding errors here
        if lines[i][gam_ind] not in gams:
            gams.append(lines[i][gam_ind])
        if lines[i][fors_ind] not in fors:
            fors.append(lines[i][fors_ind])
        if typ == 'cisp':
            if lines[i][4] not in projs:
                projs.append(lines[i][4])

    return gams,fors,projs

def read_result(dire='.',typ='cisp',res=None):
    """
    Reads the result.

    Args:
        dire(string): the directory where the result is located
        typ(string): type of calculation, 'cisp' for 2 operator 
            linear response, 'she' for spin-hall effect
        res(class(result)): add tensors to a previous result
    """

    if res == None:
        if typ == 'cisp':
            res = result(4)
        if typ == 'she':
            res = result(3)

    gams, fors, projs = find_pars(dire,typ=typ)

    with open(dire+'/input') as f:
        inp = f.readlines()
    nk = inp[1].split()[0:3]
    for i,nki in enumerate(nk):
        nk[i] = int(nki)
    nkt = nk[0]*nk[1]*nk[2]

    for gam in gams:
        for f in fors:
            if typ == 'cisp':
                for proj in projs:
                    X = read_cisp_tensor(gam,f,proj,dire)
                    res.add_tensor((gam,f,proj,nkt),X)
            if typ == 'she':
                X = read_she_tensor(gam,f,dire)
                res.add_tensor((gam,f,nkt),X)

    res.gams = gams
    res.fors = fors
    if typ == 'cisp':
        res.projs = projs
    try:
        res.nks.append(nkt)
    except:
        res.nks = [nkt]

    return res




