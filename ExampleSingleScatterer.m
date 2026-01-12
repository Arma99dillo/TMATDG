%% Set main parameters
clear; close all; addpath src; addpath TMATROM_OBJECT_ORIENTED_CORE

k = 5; % wavenumber
h = 0.5; % mesh size
p = 20; % number of plane wave directions

% define the scatterer vertices as a nvx2 matrix
scatt.vertices = [-1, -1; -1, 1; 1, 1; 1, -1]; 
scatt.type = 'trans'; 
scatt.n_in = 3+1i;

%% Compute T-matrix
[tmat,solver]=ComputeTMatrix(k,h,p,scatt);

%% Solve and plot solution
theta = -pi/3; % incident angle
uinc = plane_wave(theta,k); % incident plane wave function
PlotPar.type_plot = 'tot'; PlotPar.inside = true; 
PlotPar.limX = [-5,5]; PlotPar.limY = [-5,5];
PlotSolution(tmat,solver,uinc,PlotPar);

%% Scatterer rotation
rotation_angle = pi/6; 
[rotTmat,rotSolver] = RotateTmat(tmat,solver,rotation_angle); % rotate
rotTmat.setOrigin(1-1i); % translate T-matrix
PlotPar.type_plot = 'tot'; PlotPar.inside = true; 
PlotPar.limX = [-5,5]; PlotPar.limY = [-5,5];
PlotSolution(rotTmat,rotSolver,uinc,PlotPar);

%% Change incident field
uinc_new = point_source(3+2i,k); %circular wave
PlotPar.type_plot = 'tot'; PlotPar.inside = true; 
PlotPar.limX = [-5,5]; PlotPar.limY = [-5,5];
PlotSolution(tmat,solver,uinc_new,PlotPar);

