###############################################################################
" pycode_t=i: t=0, initial state "
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

lgn=hp.lg.lgn2d(6,8)
fvs=hp.lg.fvs2d(lgn)
dgs=hp.lg.dgs2d(lgn)
grn=hp.lg.elm2d(lgn)

hps=hp.unit([62.5,62.5],[62.5,62.5])

hps.add(mask_add)
hps.add(mask_add)
hps.add(mask_add)

hps.getelm()

hps.elm[11].split()
hps.elm[14].split()
hps.elm[15].split()

hps.getelm()

hps.elm[17].split()
hps.elm[18].split()

hps.getelm()

###############################################################################

hps.setnod_lgn(lgn)

hps.setfld_fvs(fvs)
hps.setfld_dgs(dgs)

for elm in hps.elm:

    fld=elm.fld_fvs

    fld.f[0][:,:]=1e-7
    fld.f[2][:,:]=1e-7

    fld.f[2]=fld.f[2]+0.05*np.exp(-(fld.x/4.)**2-((fld.y-100.0)/4.)**2)

    nod=elm.nod_lgn
    fld=elm.fld_dgs

    fld.f[0][:,:]=1e-7
    fld.f[2][:,:]=1e-7

    fld.f[2]=fld.f[2]+0.05*np.exp(-(nod.x/4.)**2-((nod.y-100.0)/4.)**2)

###############################################################################

hps.save("pydata/time_0")

###############################################################################
