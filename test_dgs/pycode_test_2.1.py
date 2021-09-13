###############################################################################
" pycode_sin2d_rib_0_1 "
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

    The considered test problem is Test. 2.1.

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
    
    plt.plot(nod.x,nod.y,'.b')
    
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

nod=lg.nod2d([0.,0.],[np.pi,np.pi],lgn)

fld=dg.dgs_box2d(nod)

fld.f[0]=np.sin(nod.x+nod.y)

fld.ax[:,:]=1.0+np.sin(nod.x+nod.y)
fld.ay[:,:]=1.0+np.cos(nod.x+nod.y)

fld.bx[:,:]=0.1+0.1*np.cos(nod.x+nod.y)
fld.by[:,:]=0.1+0.1*np.sin(nod.x+nod.y)

fld.bx=fld.bx*fld.g[0]
fld.by=fld.by*fld.g[1]

rib=[dg.dgs_rib([fld,fld],0),dg.dgs_rib([fld,fld],1)]

###############################################################################

for j in range(nt):

    fld.edge(0,dgs)
    rib[0].flux(dgs)
    rib[1].flux(dgs)

    t_=j*dt
    fld.df[1]=getrho(nod.x,nod.y,t_)

    fld.delta(0,dgs)    
    fld.f[1]=fld.f[0]+dt*fld.df[0]

    fld.edge(1,dgs)
    rib[0].flux(dgs)
    rib[1].flux(dgs)

    t_=(j+1)*dt
    fld.df[1]=getrho(nod.x,nod.y,t_)

    fld.delta(1,dgs)    
    fld.f[0]=(fld.f[0]+fld.f[1]+dt*fld.df[0])*0.5

###############################################################################
""""
Output values:

    L1: the L-infinity error in the solution
    L2: the global L2-error in the solution

"""

t_=nt*dt

x,y=nod.getmsh()
f=fld.getsol(0)

r=np.sin(nod.x+nod.y-t_)

L1,L2=fld.geterr(r-fld.f[0],lgn)

print(L1)
print(L2)

###############################################################################
" 2D plot of the solution"

plt.contourf(x,y,f,levels=np.linspace(-1.0,1.0,41),cmap="jet")
plt.axis("equal")

###############################################################################
