###############################################################################
" pycode_sind2d_rib_0_1 "
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

    The considered test problem is Test 2.1.

    The computation is performed for one period, i.e., until  t=1.

    Use plotelm() to plot the grid nodes

"""
###############################################################################
def getrho(x,y,t):

    r=-0.9*np.sin(x+y)*np.sin(-t+x+y)+1.1*np.sin(x+y)*np.cos(-t+x+y)\
    +1.1*np.sin(-t+x+y)*np.cos(x+y)+0.2*np.sin(-t+x+y)\
    +0.9*np.cos(x+y)*np.cos(-t+x+y)+1.0*np.cos(-t+x+y)

    return r

def plotelm():
    
    plt.plot(fld.x,fld.y,'.b')

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
m=8

dt=1e-3
nt=1000

###############################################################################

lgn=lg.lgn2d(n,m)
fvs=lg.fvs2d(lgn)

fld=fv.fvs_box2d([0.,0.],[np.pi,np.pi],fvs)

fld.f[0]=np.sin(fld.x+fld.y)

fld.ax[:,:]=1.0+np.sin(fld.x+fld.y)
fld.ay[:,:]=1.0+np.cos(fld.x+fld.y)

fld.bx[:,:]=0.1+0.1*np.cos(fld.x+fld.y)
fld.by[:,:]=0.1+0.1*np.sin(fld.x+fld.y)

fld.bx=fld.bx*fld.gx
fld.by=fld.by*fld.gy

rib=[fv.fvs_rib([fld,fld],0),fv.fvs_rib([fld,fld],1)]

###############################################################################

for j in range(nt):

    fld.edge(0)
    rib[0].flux()
    rib[1].flux()

    t_=j*dt
    fld.df[1]=getrho(fld.x,fld.y,t_)

    fld.delta(0)    
    fld.f[1]=fld.f[0]+dt*fld.df[0]

    fld.edge(1)
    rib[0].flux()
    rib[1].flux()

    t_=(j+1)*dt
    fld.df[1]=getrho(fld.x,fld.y,t_)

    fld.delta(1)    
    fld.f[0]=(fld.f[0]+fld.f[1]+dt*fld.df[0])*0.5

###############################################################################
""""
Output values:

    L1: the L-infinity error in the solution

"""

t_=nt*dt

r=np.sin(fld.x+fld.y-t_)

L1=np.amax(abs(r-fld.f[0]))

print(L1)

###############################################################################
" 2D plot of the solution "

plt.contourf(fld.x,fld.y,fld.f[0],levels=np.linspace(-1.,1.,41),cmap="jet")
plt.axis("equal")

###############################################################################
