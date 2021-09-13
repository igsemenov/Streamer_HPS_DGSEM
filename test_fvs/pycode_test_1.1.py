###############################################################################
" pycode_gauss_rib_0_1 "
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
#------------------------------------------------------------------------------
import pylib_lgn as lg
import pylib_fvs as fv
###############################################################################
"""
Introduction:

    This script is used to test the FV scheme.

    The considered test problem is Test 1.1.

    The computation is performed for one period, i.e., run till t=8.

    Use plotelm() to plot the grid nodes

"""
###############################################################################
def gauss(x,y,ax,ay,b,t):

    s=(0.25/b)+t
    q=0.25/(b*s)

    f=q*np.exp(-q*(x-ax*t)**2)*np.exp(-q*(y-ay*t)**2)

    return f

def plotelm():
    
    plt.plot(fld.x,fld.y,'.b')
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
    nt: the number of time steps

"""

n=4
m=8

dt=5e-3
nt=1600

ax_=1.
ay_=1.

bx_=0.05
by_=0.05

###############################################################################

c=[0.,0.]
h=[4.,4.]

lgn=lg.lgn2d(n,m)
fvs=lg.fvs2d(lgn)

fld=fv.fvs_box2d(c,h,fvs)

fld.f[0]=np.exp(-fld.x**2-fld.y**2)

fld.ax.fill(ax_)
fld.ay.fill(ay_)

fld.bx.fill(bx_*fld.gx)
fld.by.fill(by_*fld.gy)

rib=[fv.fvs_rib([fld,fld],0),fv.fvs_rib([fld,fld],1)]

###############################################################################
for j in range(nt):

    fld.edge(0)
    rib[0].flux()
    rib[1].flux()

    fld.delta(0)
    fld.f[1]=fld.f[0]+dt*fld.df[0]

    fld.edge(1)
    rib[0].flux()
    rib[1].flux()

    fld.delta(1)
    fld.f[0]=0.5*(fld.f[0]+fld.f[1]+dt*fld.df[0])

###############################################################################
""""
Output values:

    L1: the L-infinity error in the solution

"""

r=gauss(fld.x,fld.y,0.,0.,bx_,dt*nt)

L1=np.amax(abs(r-fld.f[0]))

print(L1)

###############################################################################
" 2D plot of the solution "

fmin=np.amin(fld.f[0])
fmax=np.amax(fld.f[0])

plt.contourf(fld.x,fld.y,fld.f[0],levels=np.linspace(fmin,fmax),cmap="jet")
plt.axis("equal")

###############################################################################
