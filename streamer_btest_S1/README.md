# Introduction
This folder contains the scripts used to perform the simulation **S1** described in **Section 4.1**. 

The scripts must be placed in the same working directory.

The working directory must contain the folder pydata that includes the subfolders time_#, where #=0,1,2,3,4,5.

The simulation steps are as follows:

1. Run pycode_t=i.py. This script initializes the mesh and the solution at t=0.
2. Run pycode_t=s.py sequentially, by setting time_moment=0,...,4 in the script file.

After that, the folders pydata\time_#, with #=0,1,2,3,4,5, contain the simulation results corresponding to the time moments 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 ns, respectively.

# Retrieving the simulation results

See README.md in streamer_btest_S0

# Files Overview
## pycode_t=i.py
This script is used to initialize the simulation.
## pycode_t=s.py
This script is used to perform the simulations steps.
## pycode_sol.py
This script is used to retrieve the simulation results.
## pylib_hps.py
This script contains the basic utilities of the HPS scheme and those of the entire simulation procedure.
## pyib_dgs.py
This script contains the basic utilities of the DGSEM scheme.
## pyib_fvs.py
This script contains the basic utilities of the FV scheme.
## pyib_psn.py
This script contains the basic utilities of the HPS scheme within a single block element.
## pylib_lgn.py
This script contains the basic utilities for working with the polynomial basis functions.
