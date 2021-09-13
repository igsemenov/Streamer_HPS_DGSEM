###############################################################################
" pycode_gauss_rib_4_5 "
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

    The considered test problem is Test 1.3.

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

    rgb=[".r",".b",".g"]
    for i in [0,1,2]:
        plt.plot(fld[i].x,fld[i].y,rgb[i])
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
nt=3200

ax_=-2.0
ay_=0.5

bx_=0.02
by_=0.02

###############################################################################

lgn=lg.lgn2d(n,m)
fvs=lg.fvs2d(lgn)

fld=[0,0,0]

fld[0]=fv.fvs_box2d([-2.,-2.],[2.,2.],fvs)
fld[1]=fv.fvs_box2d([-2.,+2.],[2.,2.],fvs)
fld[2]=fv.fvs_box2d([+4.,+0.],[4.,4.],fvs)

x0=0.
y0=0.

for fld_ in fld:

    fld_.f[0]=np.exp(-(fld_.x-x0)**2-(fld_.y-y0)**2)

    fld_.ax.fill(ax_)
    fld_.ay.fill(ay_)

    fld_.bx.fill(bx_*fld_.gx)
    fld_.by.fill(by_*fld_.gy)

rib=np.zeros(5).tolist()

rib[0]=fv.fvs_rib([fld[0],fld[1],fld[2]],5)
rib[1]=fv.fvs_rib([fld[0],fld[1],fld[2]],4)

rib[2]=fv.fvs_rib([fld[1],fld[0]],0)
rib[3]=fv.fvs_rib([fld[2],fld[2]],0)

rib[4]=fv.fvs_rib([fld[0],fld[1]],0)

###############################################################################

for j in range(nt):

    for fld_ in fld:
        fld_.edge(0)
    for rib_ in rib:
        rib_.flux()

    for fld_ in fld:
        fld_.delta(0)    
        fld_.f[1]=fld_.f[0]+dt*fld_.df[0]

    for fld_ in fld:
        fld_.edge(1)
    for rib_ in rib:
        rib_.flux()

    for fld_ in fld:
        fld_.delta(1)    
        fld_.f[0]=(fld_.f[0]+fld_.f[1]+dt*fld_.df[0])*0.5

###############################################################################
""""
Output values:

    L1: the L-infinity error in the solution for each block element

"""

L1=list()

for fld_ in fld:

    r=gauss(fld_.x-x0,fld_.y-y0,0.25,0.,bx_,nt*dt)
    L1.append(np.amax(abs(r-fld_.f[0])))

print(L1)

###############################################################################
" 2D plot of the solution "

fmin=min([np.amin(fld_.f[0]) for fld_ in fld])
fmax=max([np.amax(fld_.f[0]) for fld_ in fld])

for fld_ in fld:
    plt.contourf(fld_.x,fld_.y,fld_.f[0],
                 levels=np.linspace(fmin,fmax),cmap="jet")
plt.axis("equal")

###############################################################################
