###############################################################################
" pycode_int_a=0 "
###############################################################################
import numpy as np
import pylib_lgn as lg
###############################################################################
def test_solve(n):

    lgn=lg.lgn1d(n)

    nod=lg.nod2d([0.,0.],[1.,1.],lgn)

    grn=lg.elm2d(lgn)

#------------------------------------------------------------------------------

    x=nod.x[:,0]
    y=nod.y[:,0]
    
    Re=-2.*4.*np.sin(2.*x+2.*y)

#------------------------------------------------------------------------------

    Ae=grn.ope.Lx+grn.ope.Ly

    Ae=np.linalg.solve(Ae,np.eye(n*n))

    Ab=-grn.opb.Lx-grn.opb.Ly

    G=Ae.dot(Ab)
    g=Ae.dot(Re)

#------------------------------------------------------------------------------

    xe,ye=nod.getedg()
    
    xe=np.hstack(xe)
    ye=np.hstack(ye)
    
    F=np.sin(2.*xe+2.*ye)
    
    f=grn.ope.Q.dot(g)+(grn.ope.Q.dot(G)+grn.opb.Q).dot(F)
    
    gx=grn.ope.Ux.dot(g)+(grn.ope.Ux.dot(G)+grn.opb.Ux).dot(F)
    gy=grn.ope.Uy.dot(g)+(grn.ope.Uy.dot(G)+grn.opb.Uy).dot(F)

#------------------------------------------------------------------------------

    err=np.amax(abs(f-np.sin(2.*x+2.*y)))

    err_x=np.amax(abs(gx-2.*np.cos(2.*x+2.*y)))
    err_y=np.amax(abs(gy-2.*np.cos(2.*x+2.*y)))

    return n,err,err_x,err_y

###############################################################################
