###############################################################################
" pycode_hps_a=0 "
###############################################################################
import pylib_hps as hp
from pylib_hps import np
###############################################################################
"""

This script implements the HPS scheme for solving the Poisson equation (C.1)

The problem formulated in Cartesian coordinates (alpha=0) is considered.

"""
###############################################################################
def mask_add(root):
    return True

def getorder(err):

    k=np.sort(err.keys())

    odr=np.zeros(len(k)-1)

    for i in range(len(k)-1):
        odr[i]=np.log2(err[k[i]]/err[k[i+1]])

    return odr

###############################################################################
"""

Input parameters:

    k: the refinement level

    n: the number of Gauss nodes per direction for each element

"""

k=3
n=6

#------------------------------------------------------------------------------

lgn=hp.lg.lgn2d(n,2)

grn=hp.lg.elm2d(lgn)

hps=hp.unit([4.,0.],[4.,4.])

for i in range(k-1):
    hps.add(mask_add)

hps.getelm()

###############################################################################

hps.setnod_lgn(lgn)

hps.getrib()
hps.setedg()

hps.setpsn(grn)
hps.psndtn(grn)

###############################################################################

for elm in hps.elm:

    nod=elm.nod_lgn

    psn=elm.psn

    psn.r[:,:]=-2.*np.sin(nod.y+nod.x)

hps.psnrhs()

###############################################################################

for i in range(4):
    
    (x,y)=hps.edg[i]

    hps.top.rhs.f[i]=np.sin(x+y)

hps.psnsol()

###############################################################################

err=[]

err_x=[]
err_y=[]

for elm in hps.elm:

    x=elm.nod_lgn.x
    y=elm.nod_lgn.y
    
    g=np.sin(x+y)

    gx=np.cos(x+y)
    gy=np.cos(x+y)

    f=elm.psn.f[0]

    fx=elm.psn.f[1]
    fy=elm.psn.f[2]

#------------------------------------------------------------------------------

    c2=f-g

    c2x=fx-gx
    c2y=fy-gy

    c2=np.sum(lgn.I*(c2*c2),axis=0)*(elm.nod_lgn.h[0]*elm.nod_lgn.h[0])

    c2x=np.sum(lgn.I*(c2x*c2x),axis=0)*(elm.nod_lgn.h[0]*elm.nod_lgn.h[0])
    c2y=np.sum(lgn.I*(c2y*c2y),axis=0)*(elm.nod_lgn.h[0]*elm.nod_lgn.h[0])

    err.append(c2)

    err_x.append(c2x)
    err_y.append(c2y)

#------------------------------------------------------------------------------
"""

Output values:
    
    err: the global L2 error in the solution

    err_x: the global L2 error in the derivative with respect to x

    err_y: the global L2 error in the derivative with respect to y

"""

err=np.sqrt(np.sum(np.asarray(err)))

err_x=np.sqrt(np.sum(np.asarray(err_x)))
err_y=np.sqrt(np.sum(np.asarray(err_y)))

###############################################################################
