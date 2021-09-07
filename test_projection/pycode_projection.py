###############################################################################
" pycode_projection "
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------
import pylib_lgn as lg
###############################################################################
"""

This script can be used to test the projection procedure presented in
Section 3.5.1 The numerical flux function

The corresponding results are shown in Fig. 5 of Section 3.5.1.

The notations used are the same as in Section 3.5.1.

"""
###############################################################################
"""

Script Overview:

    n: the number of Gauss nodes

    t: 1D array that defines a uniform mesh on [-1,1]

    f1,f2: the vectors V1D[u1], V1D[u2]

    g1: the vector given by g1=W11*f1+W12*f2
    g2: the vector given by g2=W21*f1+W22*f2

    u1: the polynomial u1 defined on the mesh t1=t-1
    u2: the polynomial u2 defined on the mesh t2=t+1

    v1: the polynomial v1 defined on the mesh t1=t-1
    v2: the polynomial v2 defined on the mesh t2=t+1

"""
###############################################################################

n=4

lgn=lg.lgn1d(n)

x1=lgn.x-1.
x2=lgn.x+1

f1=np.sin(x1)
f2=np.sin(x2)

t=np.linspace(-1.,1.,51)
T=lg.getlgn(t,n).dot(lgn.C)

t1=t-1
t2=t+1

###############################################################################

Pa=lgn.Pa
Pb=lgn.Pb

Ua=np.dot(lgn.Pa,lgn.U)
Ub=np.dot(lgn.Pb,lgn.U)

Q=np.vstack([np.hstack([Pb,-Pa]),
             np.hstack([Ub,-Ua])])

u,s,v=np.linalg.svd(Q)

s=np.hstack([np.diag(s),np.zeros((2,2*n-2))])

v=v.T

Q=(v[:,2::]).dot(v[:,2::].T)

g=Q.dot(np.hstack([f1,f2]))

g=np.split(g,2)

g1=g[0]
g2=g[1]

###############################################################################

u1=T.dot(f1)
u2=T.dot(f2)

v1=T.dot(g1)
v2=T.dot(g2)

###############################################################################
" Here is the draft of Fig. 5 "

plt.plot(t1,v1,'g-')
plt.plot(t2,v2,'g-')

plt.plot(t1[::4],u1[::4],'ob')
plt.plot(t2[::4],u2[::4],'or')

###############################################################################
