###############################################################################
" pycode_sin2d_rib_4_5 "
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

    The considered test problem is Test 2.3.

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

    rgb=[".r",".b",".g",".k",".m"]
    for i in [0,1,2,3,4]:
        plt.plot(fld_[i].x,fld_[i].y,rgb[i])
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

fld_=[0,0,0,0,0]

fld_[0]=fv.fvs_box2d([-np.pi,-np.pi],[np.pi,np.pi],fvs)
fld_[1]=fv.fvs_box2d([-np.pi,+np.pi],[np.pi,np.pi],fvs)
fld_[2]=fv.fvs_box2d([+np.pi,-np.pi],[np.pi,np.pi],fvs)
fld_[3]=fv.fvs_box2d([+np.pi,+np.pi],[np.pi,np.pi],fvs)

fld_[4]=fv.fvs_box2d([4.*np.pi,0.],[2.*np.pi,2.*np.pi],fvs)

for fld in fld_:

    fld.f[0]=np.sin(fld.x+fld.y)

    fld.ax[:,:]=1.0+np.sin(fld.x+fld.y)
    fld.ay[:,:]=1.0+np.cos(fld.x+fld.y)

    fld.bx[:,:]=0.1+0.1*np.cos(fld.x+fld.y)
    fld.by[:,:]=0.1+0.1*np.sin(fld.x+fld.y)

    fld.bx=fld.bx*fld.gx
    fld.by=fld.by*fld.gy

rib_=np.zeros(9).tolist()

rib_[0]=fv.fvs_rib([fld_[0],fld_[1]],0)
rib_[1]=fv.fvs_rib([fld_[2],fld_[3]],0)

rib_[2]=fv.fvs_rib([fld_[1],fld_[3]],1)
rib_[3]=fv.fvs_rib([fld_[0],fld_[2]],1)

rib_[4]=fv.fvs_rib([fld_[1],fld_[0]],0)
rib_[5]=fv.fvs_rib([fld_[3],fld_[2]],0)
rib_[6]=fv.fvs_rib([fld_[4],fld_[4]],0)

rib_[7]=fv.fvs_rib([fld_[2],fld_[3],fld_[4]],5)
rib_[8]=fv.fvs_rib([fld_[0],fld_[1],fld_[4]],4)


###############################################################################

for j in range(nt):

    t_=j*dt

    for fld in fld_:    
        fld.edge(0)
        fld.df[1]=getrho(fld.x,fld.y,t_)

    for rib in rib_:
        rib.flux()

    for fld in fld_:
        fld.delta(0)    
        fld.f[1]=fld.f[0]+dt*fld.df[0]

    t_=(j+1)*dt

    for fld in  fld_:    
        fld.edge(1)
        fld.df[1]=getrho(fld.x,fld.y,t_)

    for rib in rib_:
        rib.flux()

    for fld in fld_:
        fld.delta(1)    
        fld.f[0]=(fld.f[0]+fld.f[1]+dt*fld.df[0])*0.5

###############################################################################
""""
Output values:

    L1: the L-infinity error in the solution for each block element

"""

t_=nt*dt

L1=list()

for fld in fld_:

    r=np.sin(fld.x+fld.y-t_)

    L1.append(abs(np.amax(r-fld.f[0])))

print(L1)

###############################################################################
" 2D plot of the solution "

for fld in fld_:
    plt.contourf(fld.x,fld.y,fld.f[0],
                 levels=np.linspace(-1.01,1.01,41),cmap="jet")
plt.axis("equal")


###############################################################################
