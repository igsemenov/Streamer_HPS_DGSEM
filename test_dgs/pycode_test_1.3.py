###############################################################################
" pycode_gauss_rib_4_5 "
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

    The considered test problem is Test 1.3.

    The computation is performed until t=1.

    Use plotelm() to plot the grid nodes.

"""
###############################################################################
def gauss(x,y,ax,ay,b,t):

    s=(0.25/b)+t
    q=0.25/(b*s)

    f=q*np.exp(-q*(x-ax*t)**2)*np.exp(-q*(y-ay*t)**2)

    return f

def plotelm():

    "Different colors correspond to different block elements"
    
    rgb=[".r",".b",".g"]
    for i in [0,1,2]:
        x,y=nod[i].getmsh()
        plt.plot(x,y,rgb[i])
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

ax_=4.0
ay_=0.0

bx_=0.02
by_=0.02

dt=1e-3
nt=1000

###############################################################################

lgn=lg.lgn2d(n,m)

dgs=lg.dgs2d(lgn)

nod=[0,0,0]
fld=[0,0,0]

nod[0]=lg.nod2d([-2.,-2.],[2.,2.],lgn)
nod[1]=lg.nod2d([-2.,+2.],[2.,2.],lgn)
nod[2]=lg.nod2d([+4.,+0.],[4.,4.],lgn)

fld[0]=dg.dgs_box2d(nod[0])
fld[1]=dg.dgs_box2d(nod[1])
fld[2]=dg.dgs_box2d(nod[2])

for i in [0,1,2]:

    fld[i].f[0]=np.exp(-nod[i].x*nod[i].x-nod[i].y*nod[i].y)


    fld[i].ax.fill(ax_)
    fld[i].ay.fill(ay_)

    fld[i].bx.fill(bx_*fld[i].g[0])
    fld[i].by.fill(by_*fld[i].g[1])

rib=np.zeros(5).tolist()

rib[0]=dg.dgs_rib([fld[0],fld[1],fld[2]],5)
rib[1]=dg.dgs_rib([fld[0],fld[1],fld[2]],4)

rib[2]=dg.dgs_rib([fld[1],fld[0]],0)
rib[3]=dg.dgs_rib([fld[2],fld[2]],0)

rib[4]=dg.dgs_rib([fld[0],fld[1]],0)

###############################################################################

for j in range(nt):

    for fld_ in fld:
        fld_.edge(0,dgs)
    for rib_ in rib:
        rib_.flux(dgs)

    for fld_ in fld:
        fld_.delta(0,dgs)    
        fld_.f[1]=fld_.f[0]+dt*fld_.df[0]

    for fld_ in fld:
        fld_.edge(1,dgs)
    for rib_ in rib:
        rib_.flux(dgs)

    for fld_ in fld:
        fld_.delta(1,dgs)    
        fld_.f[0]=(fld_.f[0]+fld_.f[1]+dt*fld_.df[0])*0.5

###############################################################################
""""
Output values:

    L1: the L-infinity error in the solution for each block element
    L2: the L2-error in the solution for each block element

"""

t_=nt*dt

L1=list()
L2=list()

for i in [0,1,2]:

    x,y=nod[i].getmsh()
    f=fld[i].getsol(0)
    plt.contourf(x,y,f,levels=np.linspace(-0.001,1.,41),cmap="jet")
    plt.axis("equal");

    r=gauss(nod[i].x,nod[i].y,ax_,ay_,bx_,dt*nt)

    L1_,L2_=fld[i].geterr(r-fld[i].f[0],lgn)

    L1.append(L1_)
    L2.append(L2_)

print(L1)
print(L2)

###############################################################################
"2D plot of the solution"

for i in [0,1,2]:

    x,y=nod[i].getmsh()
    f=fld[i].getsol(0)
    plt.contourf(x,y,f,levels=np.linspace(-0.001,1.,41),cmap="jet")

plt.axis("equal")

###############################################################################
