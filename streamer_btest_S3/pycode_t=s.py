###############################################################################
" pycode: step, t>0 "
###############################################################################
import pylib_hps as hp
#------------------------------------------------------------------------------
import time
###############################################################################
"""

This script simulates the streamer propagation on the time interval of 0.5 ns

The starting state is loaded from pydata\time_#, where #=time_moment

The computed results are saved to pydata\time_#, where #=time_moment+1

The time_moment can be 0,1,2,3,4

"""
###############################################################################
time_moment=4
#------------------------------------------------------------------------------
dt=2e-3
mt=5
nt=50
###############################################################################
def mask_sub(root):

    eps=root.eps1

    flag=False

    if root.p[2]>=4:
        if eps<0.2:
            flag=True

    return flag

def mask_add(root):

    eps=root.eps1

    flag=False

    if eps>1.0:
        flag=True

    if root.p[2]==6:
        flag=False

    if root.p[0]>0:
        flag=False

    return flag

###############################################################################

lgn=hp.lg.lgn2d(6,4)
fvs=hp.lg.fvs2d(lgn)
dgs=hp.lg.dgs2d(lgn)
grn=hp.lg.elm2d(lgn)

path="pydata//time_"+str(time_moment)

hps=hp.pickle.load(open(path+"/hps.bin"))

hps.load(path,dgs,fvs,lgn)

hps.getelm()

hps.rib_dgs_fvs()

hps.setpsn(grn)
hps.psndtn(grn)

###############################################################################

tic=time.time()

for s in range(nt):

    print("s="+str(s+1))
    hps.step_dgs_fvs(dt,mt,dgs,fvs,lgn)

    hps.top.geteps(hps.top)

    hps.top.add_amr(hps.top,mask_add,fvs,dgs,grn,lgn)

    hps.top.keyup(hps.top)
    hps.top.sub_amr(hps.top,mask_sub,fvs,dgs,grn,lgn)

    flag_fix=True
    while flag_fix:
        flag_fix=hps.fix_amr(fvs,dgs,grn,lgn)

    hps.rib_dgs_fvs()
    hps.psndtn(grn)

toc=time.time()

print("\n run time: "+str(toc-tic))

###############################################################################

hps.step_psn(lgn)

sol=hps.getsol(4,lgn,fvs)

hps.save("pydata\\time_"+str(time_moment+1))

###############################################################################
