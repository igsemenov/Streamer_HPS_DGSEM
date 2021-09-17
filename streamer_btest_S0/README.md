# Introduction
This folder contains the scripts used to perform the simulation **S0** described in **Section 4.1**. 

The scripts must be placed in the same working directory.

The working directory must contain the folder pydata that includes the subfolders time_#, where #=0,1,2,3,4,5.

These folders can be created by running the script pycode_dir.py.

The simulation steps are as follows:

1. Run pycode_t=i.py. This script initializes the mesh and the solution at t=0.
2. Run pycode_t=s.py sequentially, by setting time_moment=0,...,4 in the script file.

Folders pydata\time_# contain the simulation results corresponding to the time moments 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 ns, respectively.

**The units for physical quantities are as follows:**
 - time: 1 ns
 - length: 1e-4 m
 - density: 1e20 1/m3
 - electric field: 1 kV/cm
 - electric potential: 1 kV

# Retrieving the simulation results

The simulation results can be loaded using the script pycode_sol.py.

The starting point to get the results is the list **hps.elm** that contains the pointers to the block elements.

The computational mesh can be plotted by calling hps.plot().

Suppose elm=hps.elm[#], where #=0,1,2,3,.... Then the structure of the output data is as follows:

 - elm.c - coordinates of the block element center
 - elm.h - half-length of the block element side
 - elm.p[2] - the refinement level of the block element 
 - elm.nod_lgn.getmsh() - tuple (x,y), where x,y are 2D arrays that contain the coordinates of the Gauss nodes within the block element
 - elm.fld_dgs.getsol(0) - 2D array that contains the electron density at the Gauss nodes
 - elm.fld_dgs.getsol(2) - 2D array that contains the ion density at the Gauss nodes
 - elm.psn.getsol(1) - 2D array that contains the derivative of the potential with respect to x at the Gauss nodes
 - elm.psn.getsol(2) - 2D array that contains the derivative of the potential with respect to y at the Gauss nodes
 - elm.fld_fvs.x - 2D array that contains the x-coordinate of the FV mesh nodes
 - elm.fld_fvs.y - 2D array that contains the y-coordinate of the FV mesh nodes
 - elm.fld_fvs.f[0] - 2D array that contains the electron density at the FV mesh nodes
 - elm.fld_fvs.f[2] - 2D array that contains the ion density at the FV mesh nodes

Using the derivatives of the potential, the corresponding components of the electric field are obtained as

Ex=-1809.\*elm.psn.getsol(1)

Ey=-1809.\*elm.psn.getsol(2)-E0

Here E0=52. is the background electric field.

In addition, the following data are provided in the instance **sol**:
- sol.y_axs - 1D array that contains the FV mesh nodes on the axis of symmetry
- sol.Ne_axs - 1D array that contains the electron density at the FV mesh nodes on the axis of symmetry
- sol.Ni_axs - 1D array that contains the ion density at the FV mesh nodes on the axis of symmetry
- sol.E_axs - 1D array that contains the axial electric field at the FV mesh nodes on the axis of symmetry

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
## pycode_dir.py
This script can be used to create the folder pydata (including required subfolders).
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
