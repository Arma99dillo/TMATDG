%% Set main parameters
clear; close all; addpath src; addpath TMATROM_OBJECT_ORIENTED_CORE

k = 10; % wavenumber
h = 1; % mesh size
p = 15; % number of plane wave directions

% define the shape vertices as a Nx2 matrix
NShape = 3; ScatShape=cell(NShape,1);
ScatShape{1}.vertices = [1/3, 1/3; 1/3, 1; -1/3, 1; -1/3, 1/3; -1, 1/3;
    -1, -1/3; -1/3, -1/3; -1/3, -1; 1/3, -1; 1/3, -1/3; 1, -1/3; 1, 1/3]; % cross
ScatShape{1}.type = 'dir';
ScatShape{2}.vertices = [0, 1; -sqrt(3)/2, -1/2; sqrt(3)/2, -1/2,]; % triangle
ScatShape{2}.type = 'trans'; ScatShape{2}.n_in = 2.5;
ScatShape{3}.vertices = [1, 0; 0, 1; -1, 0; 0,-1]; % square
ScatShape{3}.type = 'dir';

%define the scatterer type, position and rotation angle
ScatArr.shape = [1; 1; 2; 3; 2]; % scatterer shape
ScatArr.pos = [-4-4i; 4-3.5i; 0; -3+4i; 3.5+3i]; % scatterer center position
ScatArr.rot = [-pi/4; 0; 0; 0; pi]; % rotation angle

%% Compute multiple scattering problem
theta = 3*pi/4; % incident angle
uinc = plane_wave(theta,k); % incident plane wave function
PlotPar.inside = true; PlotPar.limX=[-7,7]; PlotPar.limY=[-7,7];
SavePath.choice=true; SavePath.file = 'MultiTest.mat';
MultiScatt(k,h,p,uinc,ScatShape,ScatArr,PlotPar,SavePath)


%% Solve with a different configuration but same shapes
ScatArrNew.shape = [1; 2; 3; 3; 2; 1; 2; 2; 1; 3];
ScatArrNew.pos = [5+5.5i; 0+5i; -4+4.5i; 3-1.5i; -1+0.5i; -5-5i; -0.5-4.5i; 4-5i; 4+2i; -4]; 
ScatArrNew.rot = [0; pi/4; pi/2; 0; 0; 0; pi; -pi/2; pi/4; pi/4];
tmat=load(SavePath.file).tmat; solver=load(SavePath.file).solver;
PlotPar.inside = false; PlotPar.limX=[-7,7]; PlotPar.limY=[-7,7];
MultiTmatSolve(k,uinc,tmat,solver,ScatArrNew,PlotPar);