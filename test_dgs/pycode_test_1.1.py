###############################################################################
" pycode_gauss_rib_0_1 "
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
#------------------------------------------------------------------------------
import pylib_lgn as lg
import pylib_dgs as dg
###############################################################################
"""
Introduction:

    This script is used to test the DGSEM scheme.

    The considered test problem is Test 1.1.

    The computation is performed for one period, i.e., until t=8.

    Use plotelm() to plot the grid nodes.

"""
###############################################################################
def gauss(x,y,ax,ay,b,t):

    s=(0.25/b)+t
    q=0.25/(b*s)

    f=q*np.exp(-q*(x-ax*t)**2)*np.exp(-q*(y-ay*t)**2)

    return f

def plotelm():
    
    plt.plot(nod.x,nod.y,'.b')
    plt.axis("equal")

###############################################################################
"""
Input parameters:

    n: the number of Gauss nodes per direction for each finite element
    m: the number of finite elements per direction for each block element

    ax_: the advection velocity in the x-direction
    ay_: the advection velocity in the y-direction

    bx_: the diffusion coefficient in the x-direction
    by_: the diffusion coefficient in the y-direction

    dt: the time step
    nt: the number if time steps

"""

n=4
m=16

ax_=0.0
ay_=2.0

bx_=0.02
by_=0.02

dt=2e-3
nt=4000

###############################################################################

lgn=lg.lgn2d(n,m)

dgs=lg.dgs2d(lgn)

nod=lg.nod2d([0.,0.],[4.,4.],lgn)

fld=dg.dgs_box2d(nod)

fld.f[0]=np.exp(-nod.x*nod.x-nod.y*nod.y)

fld.ax.fill(ax_)
fld.ay.fill(ay_)

fld.bx.fill(bx_*fld.g[0])
fld.by.fill(by_*fld.g[1])

rib=[dg.dgs_rib([fld,fld],0),dg.dgs_rib([fld,fld],1)]

###############################################################################

for j in range(nt):

    fld.edge(0,dgs)
    rib[0].flux(dgs)
    rib[1].flux(dgs)

    fld.delta(0,dgs)
    fld.f[1]=fld.f[0]+dt*fld.df[0]

    fld.edge(1,dgs)
    rib[0].flux(dgs)
    rib[1].flux(dgs)

    fld.delta(1,dgs)
    fld.f[0]=(fld.f[0]+fld.f[1]+dt*fld.df[0])*0.5

###############################################################################
""""
Output values:

    L1: the L-infinity error in the solution
    L2: the global L2-error in the solution

"""

g=gauss(nod.x,nod.y,0.,0.,0.02,dt*nt)

err=fld.f[0]-g

L1,L2=fld.geterr(err,lgn)

print(L1)
print(L2)

x,y=nod.getmsh()
f=fld.getsol(0)

###############################################################################
" 2D plot of the solution "

plt.contourf(x,y,f,levels=np.linspace(0.,1.,21),cmap="jet")
plt.axis("equal")

###############################################################################

