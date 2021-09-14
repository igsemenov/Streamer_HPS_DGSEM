# Introduction
This folder contains the scripts used to perform the simulation **S0** described in **Section 4.1**. 

The scripts must be placed in the same directory. This directory must contain the folder pydata that includes the subfolders time_#, with #=0,1,2,3,4,5.

The simulation steps are as follows:

1. Run pycode_t=i.py. This script initiales the mesh and the solution at t=0.
2. Run pycode_t=s.py sequentially using time_moment=0,1,2,3,4.

At the end of the simulation, the folders pydata\time_#, with #=0,1,2,3,4,5, contain the simulation results.

The results in different folders correspond to different time moments, as indicated below:
pydata\time_0: t=0, initial state
pydata\time_1: t=0.5ns
pydata\time_2: t=1.0ns
pydata\time_3: t=1.5ns
pydata\time_4: t=2.0ns
pydata\time_5: t=2.5ns

