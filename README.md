# SRM
Smoothing regularization method for reconstructing Stokes parameters from division-of-focal-plane modulation

"SRM.m" is a MATLAB code implementation of the smoothing regularization method in the submission to Optics Express, entitled "Least-squares and smoothing regularization methods for reconstructing Stokes parameters from division-of-focal-plane modulation".

Notice that the regularization parameters needs to be per-processed by the scaling and rotation tranformation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% img: DoFP image.                                                        %
% m0, m1, and m2: Modulation parameters                                   %
% lmd0, lmd1, and lmd2: Regularization parameters.                        %
% s0, s1, and s2: Stokes parameters.                                      %
% Using pcg to solve A*s = b.                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
