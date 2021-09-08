# Introduction
This folder contains the scripts used to test the HPS solver.

The corresponding results are presented in the **Appendix C. Validation of the HPS solver**.

In these experiments, the Poisson equation (C.1) was solved.

In particular, this folder contains the scripts used in the computations with **alpha=0** (problem formulated in Cartesian coordinates).

# Files Overview
## pycode_hps_a=0.py
This script implements the HPS scheme for solving the Poisson equation (C.1) with alpha=0. See the script file for further details.
## pylib_hps.py
This script contains the basic utilities for the HPS scheme as well as for the DGSEM and FV schemes.
## pylib_psn.py
This script contains the basic utilities for the HPS scheme on the block elements (second level grid).
## pylib_dgs.py
This script contains the basic utilities for the DGSEM scheme.
## pylib_fvs.py
This script contains the basic utilities for the FV scheme.
## pylib_lgn.py
This script contains the basic utilities for working with the polynomial basis functions.
