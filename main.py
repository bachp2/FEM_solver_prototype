import numpy as np
import math
from Truss import Truss, ForceStatus
if __name__ == '__main__':
    #constriants
    constraints = [0,-0.5,0,0.4,0,0]
    #apply external force on nodes
    f      = [0,0,0,0,2,1]
    #node positions
    nodxyz = [[0,0],[10,0],[10,10]]
    #link with nodes at both ends
    elenod = [[0,1],[2,1],[2,0]]
    #Young's modulus array of elements 
    elemat = [100,100,100]
    #cross-sectional area array of elements  
    elefab = [1,.5,2*math.sqrt(2)]
    utag   = []
    truss = Truss(nodxyz,elenod,elemat,elefab)
    ftag   = truss.makeForceTag([2,4,5],len(f))
    K    = truss.makeMasterStiffness(nodxyz,elenod,elemat,elefab)
    [f,Kmod] = truss.makeModMasterStiffnessConstraints(ftag, f, K, constraints)
    u    = np.linalg.solve(Kmod,f)
    fmod = np.dot(K,u)
    internf  = truss.retrieveInternForces(nodxyz,elenod,elemat,elefab,u)
    internst = truss.retrieveInternStresses(internf,elefab)
    print('displacement vector')
    print(u)
    print('force vector')
    print(f)
    print('Mod Master matrix')
    print(Kmod)

