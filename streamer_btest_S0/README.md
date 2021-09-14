# Introduction
This folder contains the scripts used to perform the simulation **S0** described in **Section 4.1**. 

The scripts must be placed in the same directory. This directory must contain the folder pydata that includes the subfolders time_#, with #=0,1,2,3,4,5.

The simulation steps are as follows:

1. Run pycode_t=i.py. This script initializes the mesh and the solution at t=0.
2. Run pycode_t=s.py sequentially with time_moment=0,1,2,3,4.

At the end of the simulation, the folders pydata\time_#, with #=0,1,2,3,4,5, contain the simulation results corresponding to the time moments 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 ns, respectively.

# Retrieving the simulation results

The simulation results can be loaded using the script pycode_sol.py.

The starting point to get the results is the list **hps.elm** that contains the pointers to the block elements.

The computational mesh can be plotted by calling hps.plot().

Suppose elm=hps.elm[#], where #=0,1,2,3,.... Then the structure of the output data is as follows:

 - elm.c - coordinates of the block element center
 - elm.h - half-length of the block element side
 - elm.p[2] - the refinement level of the block element 
 - elm.nod_lgn.x - 2D array that contains the x-coordinate of the Gauss nodes within the block element
 - elm.nod_lgn.y - 2D array that contains the y-coordinate of the Gauss nodes within the block element
 - elm.fld_fvs.x - 2D array that contains the x-coordinate of the FV mesh nodes within the block element
 - elm.fld_fvs.y - 2D array that contains the y-coordinate of the FV mesh nodes within the block element
 - elm.fld_dgs.f[0] - 2D array that contains the electron density at the Gauss nodes
 - elm.fld_dgs.f[2] - 2D array that contains the ion density at the Gauss nodes
 - elm.fld_fvs.f[0] - 2D array that contains the electron density at the FV mesh nodes
 - elm.fld_fvs.f[2] - 2D array that contains the ion density at the FV mesh nodes
 - elm.psn.f[1] - 2D array that contains the derivative of the potential with respect to x at the Gauss nodes
 - elm.psn.f[2] - 2D array that contains the derivative of the potential with respect to y at the Gauss nodes

Using the derivatives of the potential, the corresponding components of the electric field (in kV/cm) are obtained as

Ex=-1809.0\*elm.psn.f[1]

Ey=-1809.0\*elm.psn.f[2]-E0

Here E0=52.0 is the background electric field.

In addition, the following data are provided in the instance **sol**:
- sol.y_fld_axs - 1D array that contains the FV mesh nodes on the axis of symmetry
- sol.y_lgn_axs - 1D array that contains the Gauss nodes on the axis of symmetry
- sol.Ne_axs - 1D array that contains the electron density at the nodes of sol.y_fld_axs
- sol.Ni_axs - 1D array that contains the ion density at the nodes of sol.y_fld_axs
- sol.E_axs - 1D array that contains the axial electric field at the nodes of sol.y_lgn_axs

These results can be quickly visualized using the following functions:
- sol.rho() - plot the distribution of sol.Ne_axs (red) and sol.Ni_axs (blue)
- sol.rho_log() - the same as sol.rho() but using the logarithmic scale for density.
- sol.field() - plot the distribution of sol.E_axs

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
