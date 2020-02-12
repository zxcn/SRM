# SRM
"SRM.m" is a MATLAB code implementation of the smoothing regularization method in the submission to Optics and Lasers in Engineering, entitled "Reconstructing Stokes parameters from non-uniform division-of-focal-plane modulation".

Input: 
DoFP image (img), modulation parameters (m0, m1, and m2), and regularization parameters (lmd0, lmd1, and lmd2).

Output: 
Stokes parameters (s0, s1, and s2)

Notice that the regularization parameters needs to be per-processed by the scaling and rotation tranformation given in the Appendix A of the submission.
