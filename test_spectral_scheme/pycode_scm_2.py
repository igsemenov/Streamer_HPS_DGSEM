###############################################################################
" pycode_scm_2 "
###############################################################################
import numpy as np
import pylib_chb as ch
import matplotlib.pyplot as plt
from scipy.special import jv
###############################################################################
def test_solve(n):

    chb=ch.chb2d(n)
    
    x,y=np.meshgrid(chb.x,chb.x)
    
    x=x.flatten()
    y=y.flatten()

###############################################################################

    i=np.arange(n*n).reshape(n,n)
    
    i0=i[0,:]
    i1=i[:,0]
    
    i2=i[-1,:]
    i3=i[:,-1]
    
    ie=i[1:-1,1:-1].flatten()
    
    ib=np.hstack([i0,i1,i2,i3])
    
    xb=x[ib]
    yb=y[ib]

###############################################################################

    ax=np.zeros_like(x)
    
    ax[ie]=1./(x[ie]+1.)
    
    A=chb.Lx+chb.Ly+ax[...,None]*chb.Ux
    
    Re=-4.*2.*jv(0,2.+2.*x)*np.sin(2.*y)
    
    Ae=A.copy()
    Ab=A.copy()
    
    Ae=Ae[:,ie]
    Ae=Ae[ie,:]
    
    Ab=Ab[:,ib]
    Ab=Ab[ie,:]

    Ae=np.linalg.solve(Ae,np.eye(Ae.shape[0]))
    
    Ab=-np.matmul(Ae,Ab)
    
    Rb=np.matmul(Ae,Re[ie])

###############################################################################

    f=x*0.
    
    f[ib]=jv(0,2.*xb+2.)*np.sin(2.*yb)
    
    f[ie]=Ab.dot(f[ib])+Rb
    
    gx=chb.Ux.dot(f)
    gy=chb.Uy.dot(f)

###############################################################################

    x_=2.*(1.+x)
    y_=2.*y
    
    f_ref=jv(0,x_)*np.sin(y_)
    
    gx_ref=2.*0.5*(jv(-1,x_)-jv(1,x_))*np.sin(y_)
    gy_ref=2.*jv(0,x_)*np.cos(y_)
    
    err=np.amax(abs(f-f_ref))
    
    err_x=np.amax(abs(gx-gx_ref))
    err_y=np.amax(abs(gy-gy_ref))

    return np.array([n,err,err_x,err_y])

###############################################################################

err=list()

for n in [4,6,8,10,12,14,16,18]:
    err.append(test_solve(n))

err=np.stack(err)

np.savetxt("err_scm_2.txt",err)

###############################################################################
