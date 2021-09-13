###############################################################################
" pycode_sin2d_rib_2_3 "
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

    The considered test problem is Test. 2.2.

    The computation is performed until t=0.5.

    Use plotelm() to plot the grid nodes

"""
###############################################################################
def getrho(x,y,t):

    r=-0.9*np.sin(x+y)*np.sin(-t+x+y)+1.1*np.sin(x+y)*np.cos(-t+x+y)\
    +1.1*np.sin(-t+x+y)*np.cos(x+y)+0.2*np.sin(-t+x+y)\
    +0.9*np.cos(x+y)*np.cos(-t+x+y)+1.0*np.cos(-t+x+y)

    return r

def plotelm():

    "Different colors correspond to different block elements"

    rgb=[".r",".b",".g",".k",".m"]
    for i in [0,1,2,3,4]:
        x,y=nod_[i].getmsh()
        plt.plot(x,y,rgb[i])

    plt.axis("equal")

###############################################################################
"""
Input parameters:

    n: the number of Gauss nodes per direction for each finite element
    m: the number of finite elements per direction for each block element

    dt: the time step
    nt: the number of time steps

"""

n=4
m=16

dt=5e-4
nt=1000

###############################################################################

lgn=lg.lgn2d(n,m)

dgs=lg.dgs2d(lgn)

nod_=[0,0,0,0,0]

nod_[0]=lg.nod2d([-np.pi,-np.pi],[np.pi,np.pi],lgn)
nod_[1]=lg.nod2d([-np.pi,+np.pi],[np.pi,np.pi],lgn)
nod_[2]=lg.nod2d([+np.pi,-np.pi],[np.pi,np.pi],lgn)
nod_[3]=lg.nod2d([+np.pi,+np.pi],[np.pi,np.pi],lgn)

nod_[4]=lg.nod2d([0.,-4.*np.pi],[2.*np.pi,2.*np.pi],lgn)

fld_=[0,0,0,0,0]

for i in [0,1,2,3,4]:
    fld_[i]=dg.dgs_box2d(nod_[i])

for i in [0,1,2,3,4]:

    fld=fld_[i]
    nod=nod_[i]

    fld.f[0]=np.sin(nod.x+nod.y)

    fld.ax[:,:]=1.0+np.sin(nod.x+nod.y)
    fld.ay[:,:]=1.0+np.cos(nod.x+nod.y)

    fld.bx[:,:]=0.1+0.1*np.cos(nod.x+nod.y)
    fld.by[:,:]=0.1+0.1*np.sin(nod.x+nod.y)

    fld.bx=fld.bx*fld.g[0]
    fld.by=fld.by*fld.g[1]

rib_=np.zeros(9).tolist()

rib_[0]=dg.dgs_rib([fld_[0],fld_[1]],0)
rib_[1]=dg.dgs_rib([fld_[2],fld_[3]],0)

rib_[2]=dg.dgs_rib([fld_[1],fld_[3]],1)
rib_[3]=dg.dgs_rib([fld_[0],fld_[2]],1)

rib_[4]=dg.dgs_rib([fld_[3],fld_[1]],1)
rib_[5]=dg.dgs_rib([fld_[2],fld_[0]],1)

rib_[6]=dg.dgs_rib([fld_[4],fld_[4]],1)

rib_[7]=dg.dgs_rib([fld_[0],fld_[2],fld_[4]],2)
rib_[8]=dg.dgs_rib([fld_[1],fld_[3],fld_[4]],3)

###############################################################################

for j in range(nt):

    t_=j*dt
    
    for i in  [0,1,2,3,4]:    
        fld_[i].edge(0,dgs)
        fld_[i].df[1]=getrho(nod_[i].x,nod_[i].y,t_)

    for rib in rib_:
        rib.flux(dgs)

    for fld in fld_:
        fld.delta(0,dgs)    
        fld.f[1]=fld.f[0]+dt*fld.df[0]

    t_=(j+1)*dt

    for i in  [0,1,2,3,4]:    
        fld_[i].edge(1,dgs)
        fld_[i].df[1]=getrho(nod_[i].x,nod_[i].y,t_)

    for rib in rib_:
        rib.flux(dgs)

    for fld in fld_:
        fld.delta(1,dgs)    
        fld.f[0]=(fld.f[0]+fld.f[1]+dt*fld.df[0])*0.5

###############################################################################
"""
Output values:

    L1: the L-infinity error in the solution for each block element
    L2: the global L2-error in the solution for each block element

"""

t_=nt*dt

L1=list()
L2=list()

for i in [0,1,2,3,4]:
    
    fld=fld_[i]
    nod=nod_[i]

    x,y=nod.getmsh()
    f=fld.getsol(0)

    r=np.sin(nod.x+nod.y-t_)

    L1_,L2_=fld.geterr(r-fld.f[0],lgn)

    L1.append(L1_)
    L2.append(L2_)

print(L1)
print(L2)

###############################################################################
" 2D plot of the solution "

for i in [0,1,2,3,4]:

    fld=fld_[i]
    nod=nod_[i]

    x,y=nod.getmsh()
    f=fld.getsol(0)
    plt.contourf(x,y,f,levels=np.linspace(-1.01,1.01,41),cmap="jet")

plt.axis("equal")

###############################################################################
