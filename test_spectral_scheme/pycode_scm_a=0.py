###############################################################################
" pycode_scm_a=0 "
###############################################################################
import numpy as np
import pylib_chb as ch
###############################################################################
def test_solve(n):

    chb=ch.chb2d(n)
    
    x,y=np.meshgrid(chb.x,chb.x)
    
    x=x.flatten()
    y=y.flatten()

#------------------------------------------------------------------------------

    i=np.arange(n*n).reshape(n,n)
    
    i0=i[0,:]
    i1=i[:,0]
    
    i2=i[-1,:]
    i3=i[:,-1]
    
    ie=i[1:-1,1:-1].flatten()
    
    ib=np.hstack([i0,i1,i2,i3])
    
    xb=x[ib]
    yb=y[ib]

#------------------------------------------------------------------------------

    A=chb.Lx+chb.Ly
    
    Re=-2.*4.*np.sin(2.*x+2.*y)
    
    Ae=A.copy()
    Ab=A.copy()
    
    Ae=Ae[:,ie]
    Ae=Ae[ie,:]
    
    Ab=Ab[:,ib]
    Ab=Ab[ie,:]

    Ae=np.linalg.solve(Ae,np.eye(Ae.shape[0]))

    Ab=-np.matmul(Ae,Ab)
    
    Rb=np.matmul(Ae,Re[ie])

#------------------------------------------------------------------------------

    f=x*0.
    
    f[ib]=np.sin(2.*xb+2.*yb)
    
    f[ie]=Ab.dot(f[ib])+Rb
    
    gx=chb.Ux.dot(f)
    gy=chb.Uy.dot(f)
    
    err=np.amax(abs(f-np.sin(2.*x+2.*y)))
    
    err_x=np.amax(abs(gx-2.*np.cos(2.*x+2.*y)))
    err_y=np.amax(abs(gy-2.*np.cos(2.*x+2.*y)))

    return np.array([n,err,err_x,err_y])

###############################################################################
