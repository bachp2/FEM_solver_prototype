import numpy as np
import math
import enum

class ForceStatus(enum.IntEnum):
    UNKNOWN = 0
    KNOWN   = 1
    SOLVED  = 2 

class Truss:
    def __init__(self,nodxy,elenod, elemat, elefab):
        self.nodxy  = nodxy
        self.elenod = elenod
        self.elemat = elemat
        self.elefab = elefab
        print('create instance')

    def _makeElemStiffness(self,n1xy,n2xy,E,A):
        x21 = n2xy[0] - n1xy[0]
        y21 = n2xy[1] - n1xy[1]
        L  = math.sqrt(x21**2+y21**2)
        c   = x21/L*1.0
        s   = y21/L*1.0
        return (E*A/L)*np.matrix([[c**2, c*s,-c**2,-c*s],[c**2, c*s,-c**2,-c*s],
                              [-c**2,-s*c, c**2, s*c],[-s*c,-s**2, s*c, s**2]
                         ])
    
    def _getEFT(self,ni,nj):
        eft = [2*ni,2*ni+1,2*nj,2*nj+1]
        return eft
    
    def makeForceTag(self, nodxyknown, numnodes):
        a = np.empty(numnodes, dtype=int)
        for x in range(numnodes):
            if x in nodxyknown:
                a[x] = int(ForceStatus.KNOWN)
            else:
                a[x] = int(ForceStatus.UNKNOWN)
        return a

    def modForceTag(self,ftag):
        for tag in ftag:
            if tag == int(ForceStatus.UNKNOWN):
                tag = int(ForceStatus.SOLVED)

    def makeModMasterStiffness(self, ftag, K):
        Kmod = K.copy()
        c = 0
        for i in ftag:
            if i == int(ForceStatus.UNKNOWN):
                Kmod[c]   = 0
                Kmod[:,c] = 0
                Kmod[c,c] = 1
            c += 1
        return Kmod
        
    #the definition is not technically correct since this doesn't apply to zero constraints
    def makeModMasterStiffnessConstraints(self, ftag, f, K, constraints):
        Kmod = K.copy()
        f = np.add(f,constraints)
        c = 0
        for i in ftag:
            if i == int(ForceStatus.UNKNOWN):
                Kmod[c]   = 0
                Kmod[c,c] = 1
            c += 1
        #print(Kmod)
        c = 0
        for i in ftag:
            if i == int(ForceStatus.KNOWN):
                f[c] -= np.dot(Kmod[c],constraints)
            c += 1
        #print(f)
        c = 0
        for i in ftag:
            if i == int(ForceStatus.UNKNOWN):
                Kmod[:,c] = 0
                Kmod[c,c] = 1
            c += 1
        
        return [f,Kmod]
        
    def makeMasterStiffness(self, nodxyz, elenod, elemat, elefab):
        numnod = len(nodxyz)
        numele = len(elenod)
        K = np.zeros((numnod*2,numnod*2))
        for e in range(numele):
            (ni,nj) = elenod[e]
            eft = self._getEFT(ni,nj)
            E   = elemat[e]
            A   = elefab[e]
            Ke  = self._makeElemStiffness(nodxyz[ni],nodxyz[nj],E,A)
            for i in range(4):
                ii = eft[i]
                for j in range(i,4):
                    jj = eft[j]
                    K[ii,jj] += Ke[i,j]
                    K[jj,ii] = K[ii,jj]
        return K

    def _calcEleInternForce(self,n1xy,n2xy,E,A,ue):
        x21 = n2xy[0] - n1xy[0]
        y21 = n2xy[1] - n1xy[1]
        LLe = x21**2+y21**2
        dux = ue[2]-ue[0]
        duy = ue[3]-ue[1]
        return (E*A/LLe)*(x21*dux + y21*duy)

    def retrieveInternForces(self, nodxyz, elenod, elemat, elefab, u):
        numele = len(elenod)
        internforces = []
        for e in range(numele):
            (ni,nj) = elenod[e]
            E   = elemat[e]
            A   = elefab[e]
            ue  = [u[ni*2], u[2*ni+1], u[nj*2], u[nj*2+1]]
            internforces.append(self._calcEleInternForce(nodxyz[ni],nodxyz[nj],E,A,ue))
        return internforces
    
    def retrieveInternStresses(self,internf, elefab):
        return np.divide(internf.copy(),elefab)