# Overview

This folder contains the scripts used to perform the test computations presented in the Appendix B. Validation of the spectral method for the Poisson equation.

In these experiments, the model Poisson equation B1 was solved at alpha=0,1.

The computations were performed using two different methods: the spectral collocation method (SCM) and the method based on the integral reformulation of the problem (IRM).

# pycode_int_a=0.py
This script implements the IRM scheme for solving Eq. B1 at alpha=0.
# pycode_int_a=1.py
This script implements the IRM scheme for solving Eq. B1 at alpha=1.
# pycode_scm_a=0.py
This script implements the SCM scheme for solving Eq. B1 at alpha=0.
# pycode_scm_a=1.py
This script implements the SCM scheme for solving Eq. B1 at alpha=1.
# pylib_lgn.py
This module contains basic utilities required for implementing the IRM scheme.
# pylib_chb.py
This module contains basic utilities required for implementing the SCM scheme.
