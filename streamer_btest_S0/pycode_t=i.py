###############################################################################
" pycode: initial, t=0 "
###############################################################################
import pylib_hps as hp
from pylib_hps import np
###############################################################################
"""

This script initializes the mesh and the solution at t=0

The data are saved to pycode/time_0

"""
###############################################################################
def mask_add(root):

    if root.p[0]==0:
        return True
    else:
        return False

###############################################################################

lgn=hp.lg.lgn2d(6,4)
fvs=hp.lg.fvs2d(lgn)
dgs=hp.lg.dgs2d(lgn)
grn=hp.lg.elm2d(lgn)

hps=hp.unit([50.,50.],[50.,50.])

hps.add(mask_add)
hps.add(mask_add)
hps.add(mask_add)
hps.add(mask_add)
hps.add(mask_add)
hps.add(mask_add)

hps.getelm()

###############################################################################

hps.setnod_lgn(lgn)

hps.setfld_fvs(fvs)
hps.setfld_dgs(dgs)

for elm in hps.elm:

    fld=elm.fld_fvs

    fld.f[0][:,:]=1e-6
    fld.f[2][:,:]=1e-6

    fld.f[0]=fld.f[0]+np.exp(-(fld.x/2.1)**2-((fld.y-50.0)/2.7)**2)
    fld.f[2]=fld.f[2]+np.exp(-(fld.x/2.1)**2-((fld.y-50.0)/2.7)**2)

    nod=elm.nod_lgn
    fld=elm.fld_dgs

    fld.f[0][:,:]=1e-6
    fld.f[2][:,:]=1e-6

    fld.f[0]=fld.f[0]+np.exp(-(nod.x/2.1)**2-((nod.y-50.0)/2.7)**2)
    fld.f[2]=fld.f[2]+np.exp(-(nod.x/2.1)**2-((nod.y-50.0)/2.7)**2)

###############################################################################

hps.save("pydata/time_0")

###############################################################################
