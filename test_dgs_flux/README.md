# Introduction
This folder contains the scripts used to test the numerical flux function for the DGSEM scheme.

The corresponding results are presented in the **Appendix D. Validation of the numerical flux function**. 

In these experiments, a one-dimensional advection-diffusion equation (D.1) was solved.

# Files Overview
## pycode_dgs_flux.py
This script implements the DGSEM scheme for solving equation (D.1). See the script file for further details.
## pylib_lgn.py
This module contains the basic utilities for working with the polynomial basis functions in 1D and 2D.
## pylib_dgs.py
This module contains the basic utilities for implementing the DGSEM scheme in 1D.
