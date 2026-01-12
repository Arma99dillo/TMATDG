# TMATDG
MATLAB package for the solution of multiple scattering problems, coupling Trefftz Discontinuos Galerkin methods for Helmholtz scattering with the T-matrix method. This repository contains the scripts used to carry out the numerical tests of the method described in the paper

_TMATDG: applying TDG methods to multiple scattering via T-matrix approximation (2026)_, by Armando Maria Monforte

# Reproducibility Instructions
All the codes have been tested on MATLAB R2023b and R2024b release. The MATLAB Partial Differential Equation Toolbox™ is needed for the correct working of the code. The `TMATROM_OBJECT_ORIENTED_CORE` directory is taken from the TMATROM package from https://github.com/stuart-hawkins/tmatrom with little modifications.

# Contents and Structure
The TMATDG package is based on the TMATROM package and is used to approximate the T-matrix of ploygonal obstacles for multiple scattering problems.

Examples:
-
The following files contain the code to run all the examples in Section 4 of the paper, allowing to reproduce the figures therein:
* In `ExampleSingleScatterer.m` is shown how to approximate the T-matrix of a single polygonal scatterer. It is possible to change the shape of the obstacle just modifying the vertices of the polygon. It is also shown how to compute the solution of the scattering problem for any incident field and how to change the position and orientation of the scatterer; 
* In `ExampleMultipleScattering.m` we report an example of a multiple scattering problem. The parameters, such as the shape and arrangement of the obstacles and the incident field are easily modified;
* The `ExampleParameterDepending.m` generates the plots of section 4.3.3 of the paper.

Auxiliary functions
-
All the auxiliary files are int the `src` directory. We briefly explain what the main functions of the package do. For a detailed explanation see Section 4 of the paper.

Single Obstacle Scattering
-
Scatterers are defined by polygon vertices, type (penetrable or impenetrable) and refractive index. Here are the functions to approximate the T-matrix of a single scatterer:
* `TDGsolver`: Builds and solves the DtN–TDG linear system for a given incident regular wavefunction and computes the far-field; 
* `ComputeTMatrix.m`: Computes the T-matrix and returns it together with a TDGsolver instance;
* `PlotSolution.m`: Solves a scattering problem using the T-matrix approximation and plots the scattered or total field;
* `RotateTmat.m`: Generates the T-matrix and solver for a rotated scatterer by analytically rotating the T-matrix and the mesh. 

Multiple Scattering
-
Multiple scattering configurations consist of N obstacles with at most N_S distinct shapes. Each obstacle is defined by a shape index, position and rotation angle. T-matrices are computed once per shape and then translated and rotated; the functions are as follows:
* `MultiScatt.m`: Computes T-matrices for all distinct shapes and solves the multiple scattering problem using GMRES through TMATROM. Optionally plots the field near scatterers using DtN–TDG when enclosing circles do not intersect, and can save computed matrices and solvers;
* `MultiTmatSolve.m`: Solves a multiple scattering problem using precomputed T-matrices only, avoiding recomputation and enabling fast simulations with new arrangements or incident fields. 
  
