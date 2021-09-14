# Introduction
This folder contains the scripts used to perform the simulation **S0** described in **Section 4.1**. 

The scripts must be placed in the same directory. This directory must contain the folder pydata that includes the subfolders time_#, with #=0,1,2,3,4,5.

The simulation steps are as follows:

1. Run pycode_t=i.py. This script initiales the mesh and the solution at t=0.
2. Run pycode_t=s.py sequentially using time_moment=0,1,2,3,4.

At the end of the simulation, the folders pydata\time_#, with #=0,1,2,3,4,5, contain the simulation results corresponding to the time moments 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 ns, respectively.

The simulation results can be retrieved using the script pycode_sol.py.
