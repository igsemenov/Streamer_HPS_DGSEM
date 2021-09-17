# Introduction

This folder contains the computer code that was used to perform the numerical experiments presented in the paper

**I. L. Semenov, "A spectral element method for modelling streamer discharges in low-temperature atmospheric-pressure plasmas"**

The code is written in the Python programming language (Python 2.7.16 was used).

In what follows, citations and references correspond to those of the paper.

The frequently used abbreviations are:
 - **DGSEM** - discontinious Galerkin spectral element method
 - **FV** - finite volume
 - **HPS** - the Hierarchical Poincare-Steklov scheme

**Overview of the subfolders:**
 - **streamer_btest_S0** - the code used for the simulation S0 described in Section 4.1.
 - **streamer_btest_S1** - the code used for the simulation S1 described in Section 4.1.
 - **streamer_btest_S2** - the code used for the simulation S2 described in Section 4.1.
 - **streamer_btest_S3** - the code used for the simulation S3 described in Section 4.1.
 - **streamer_btest_S4** - the code used for the simulation S4 described in Section 4.1.
 - **streamer_etest** - the code used for the simulation described in Section 4.2.
 - **test_dgs** - the code used to test the DGSEM scheme (results are not presented in the paper).
 - **test_dgs_flux** - the code used to test the numerical flux function for the DGSEM. The results are presented in the Appendix D.
 - **test_fvs** - the code used to test the FV scheme (results are not presented in the paper).
 - **test_hps_a=0** - the code used to test the HPS scheme (problem in Cartesian coordinates). The results are presented in the Appendix C.
 - **test_hps_a=1** - the code used to test the HPS scheme (axisymmetric problem). The results are presented in the Appendix C.
 - **test_projection** - the code used to test the projection procedure described in Section 3.5.1.
 - **test_spectral_scheme** - the code used to test the spectral scheme for the Poisson equation. The results are presented in the Appendix B.

