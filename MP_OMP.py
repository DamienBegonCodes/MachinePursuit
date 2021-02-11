import numpy as np
import math as m

def machinePursuit(eps,D,x):
    D = np.array(D)
    x = np.array(x)
    K = len(D[0])
    N = len(D)
    R = x
    alpha = np.zeros(K)
    xChap = np.zeros(N)
    k = 0

    while (np.linalg.norm(R)>eps and k<K):
        m = np.zeros(K)
        for j in range(K):
            m[j] = abs(np.vdot(D[:,j],R))/np.linalg.norm(D[:,j])
        ind = np.argmax(m)
        zOpt = np.vdot(D[:,ind],R)/np.linalg.norm(D[:,ind])**2
        xChap += zOpt*D[:,ind]
        R -= zOpt*D[:,ind]
        alpha[ind] = zOpt
        print("Indice = ", ind+1)
        k += 1
    print("k = ", k, "\nxChap = ", xChap, "\nR = ", R, "\naplha = ", alpha)
    return alpha

D = [[1/2*m.sqrt(2), 1/3*m.sqrt(3), 1/3*m.sqrt(6), 2/3, -1/3], 
[-1/2*m.sqrt(2), -1/3*m.sqrt(3), -1/6*m.sqrt(6), 2/3, -2/3],
[0, -1/3*m.sqrt(3), 1/6*m.sqrt(6), 1/3, 2/3]]

x = [4/3 - 1/2*m.sqrt(2), 4/3 + 1/2*m.sqrt(2), 2/3]

eps = 1e-4

# machinePursuit(eps,D,x)

"""
L approximation parcimonieuse alpha=
  - 1.
    0.
    0.
    2.
    0.
Nombre d iérations nécessaire=
    2.
Les positions des atomes actifs
    4.    1.

    4.154D-16
"""


D = [[1,1,2,5,0,0,3,-2,1,2,2,2],
[0,-1,-1,1,0,0,5,0,2,2,7,-1],
[1,1,1,5,1,2,2,1,1,1,1,5],
[1,5,2,2,5,0,-4,5,1,5,0,0],
[0,2,2,1,1,0,0,0,0,4,-1,-2],
[-1,2,2,2,-2,-3,-4,1,1,1,1,0]]

x = [-10,-10,1,21,0,9]

eps = 1e-1

# machinePursuit(eps,D,x)

def OMP(eps,D,x):
    D = np.array(D)
    x = np.array(x)
    K = len(D[0])
    N = len(D)
    R = x
    alpha = np.zeros(K)
    
    k = 0
    Inds = [i for i in range(K)]

    while (np.linalg.norm(R)>eps and k<K):
        m = np.zeros(K)
        for j in Inds:
            m[j] = abs(np.vdot(D[:,j],R))/np.linalg.norm(D[:,j])
        ind = np.argmax(m)
        print("Indice = ", ind+1)
        
        if (len(Inds) == K):
            Atomes = np.array([D[:,ind]])
        else:
            Atomes = np.append(Atomes,[D[:,ind]],axis=0)
        Inds.remove(ind)
        

        Z = np.zeros(len(Atomes))
        Phi = np.linalg.pinv(Atomes)
        xChap = np.zeros(N)
        R = x
        
        for i in range(len(Atomes)):
            Z[i] = np.vdot(Phi[:,i],x)
            xChap += Z[i]*Atomes[i]
            R = np.add(R, -Z[i]*Atomes[i])
        
        alpha[ind] = Z[len(Atomes)-1]
        
        k += 1
        
        
    print("k = ", k, "\nxChap = ", xChap, "\nR = ", R, "\naplha = ", alpha)
    return alpha


D = [[1/2*m.sqrt(2), 1/3*m.sqrt(3), 1/3*m.sqrt(6), 2/3, -1/3], 
[-1/2*m.sqrt(2), -1/3*m.sqrt(3), -1/6*m.sqrt(6), 2/3, -2/3],
[0, -1/3*m.sqrt(3), 1/6*m.sqrt(6), 1/3, 2/3]]

x = [4/3 - 1/2*m.sqrt(2), 4/3 + 1/2*m.sqrt(2), 2/3]

eps = 1e-4

# OMP(eps,D,x)

D = [[1,1,2,5,0,0,3,-2,1,2,2,2],
[0,-1,-1,1,0,0,5,0,2,2,7,-1],
[1,1,1,5,1,2,2,1,1,1,1,5],
[1,5,2,2,5,0,-4,5,1,5,0,0],
[0,2,2,1,1,0,0,0,0,4,-1,-2],
[-1,2,2,2,-2,-3,-4,1,1,1,1,0]]

x = [-10,-10,1,21,0,9]

eps = 1e-4

OMP(eps,D,x)
