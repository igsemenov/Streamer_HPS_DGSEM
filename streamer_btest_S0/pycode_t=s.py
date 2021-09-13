###############################################################################
" pycode: step, t>0 "
###############################################################################
import pylib_hps as hp
#------------------------------------------------------------------------------
import time
###############################################################################
"""

This script simulates the streamer propagation on the time interval 0.5 ns

The initial state is loaded from pydata\time_#, where #=time_moment

The computed results are saved to pydata\time_#, where #=time_moment+1

"""
###############################################################################
time_moment=4
#------------------------------------------------------------------------------
dt=2e-3
mt=5
nt=50
###############################################################################

lgn=hp.lg.lgn2d(6,4)
fvs=hp.lg.fvs2d(lgn)
dgs=hp.lg.dgs2d(lgn)
grn=hp.lg.elm2d(lgn)

path="pydata/time_"+str(time_moment)

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

toc=time.time()
print("run time: "+str(toc-tic))

###############################################################################

hps.step_psn(lgn)

sol=hps.getsol(4,lgn,fvs)

hps.save("pydata\\time_"+str(time_moment+1))

###############################################################################
