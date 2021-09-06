###############################################################################
" pycode_int_2 "
###############################################################################
import numpy as np
import pylib_lgn as lg
from scipy.special import jv
###############################################################################
def test_solve(n):

    lgn=lg.lgn1d(n)
    
    nod=lg.nod2d([0.,0.],[1.,1.],lgn)
    
    grn=lg.elm2d(lgn)

#------------------------------------------------------------------------------

    x=nod.x[:,0]
    y=nod.y[:,0]
    
    Re=-4.*2.*jv(0,2.*(1.+x))*np.sin(2.*y)

#------------------------------------------------------------------------------

    ax=1./(1.+x[...,None])

    Ae=grn.ope.Lx+grn.ope.Ly+ax*grn.ope.Ux

    Ae=np.linalg.solve(Ae,np.eye(n*n))

    Ab=-grn.opb.Lx-grn.opb.Ly-ax*grn.opb.Ux

    G=Ae.dot(Ab)
    g=Ae.dot(Re)

#------------------------------------------------------------------------------

    xe,ye=nod.getedg()
    
    xe=np.hstack(xe)
    ye=np.hstack(ye)
    
    F=jv(0,2.*(xe+1.))*np.sin(2.*ye)
    
    f=grn.ope.Q.dot(g)+(grn.ope.Q.dot(G)+grn.opb.Q).dot(F)
    
    gx=grn.ope.Ux.dot(g)+(grn.ope.Ux.dot(G)+grn.opb.Ux).dot(F)
    gy=grn.ope.Uy.dot(g)+(grn.ope.Uy.dot(G)+grn.opb.Uy).dot(F)

#------------------------------------------------------------------------------

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
