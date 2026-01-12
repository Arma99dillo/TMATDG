function [tmat,solver] = ComputeTMatrix(kwave,h,nd,scatt,varargin)

points=scatt.vertices;

% check if we have enough points and compute center
if size(points,1) <=2
    error('Not enough points to define a polygon!')
end

% check if they are anticlockwise and flip if necessary
s=0;
for j=1:(size(points,1)-1)
    s=s+(points(j+1,1)-points(j,1))*(points(j+1,2)+points(j,2));
end
s=s+(points(1,1)-points(end,1))*(points(1,2)+points(end,2));

if s>0
    points(:,1) = flip(points(:,1));
    points(:,2) = flip(points(:,2));
end

% compute polygon center
c = mean(points, 1); center = c(1)+1i*c(2);

% solver parameters
param.h=h; param.alpha=1/2; param.beta=1/2; param.delta=1/2;
param.nd=nd; param.d=zeros(2,param.nd); 
% define PW directions
for l=1:param.nd
    param.d(:,l)=[cos((2*pi*l)/param.nd); sin((2*pi*l)/param.nd)];
end
param.q=8; % gaussian quadrature
[param.nodes,param.w]=gaussquad(param.q);

% compute radius to set T-matrix dimension
points = points - mean(points, 1); % center in zero to easy compute distance
distances = sqrt(sum(points.^2, 2)); % distance from the origin
radius = max(distances);
param.vertices = points; % vertices

% set enclosing circle radius and DtN order of truncation
param.R=2*h+radius;
param.M=ceil(kwave*param.R)+5;

% compute T-matrix order of truncation
if nargin == 4 % impose truncation order
    nmax = suggestedorder(kwave,radius);
else 
    nmax = varargin{1};
end
% setup the solver
if strcmp(scatt.type, 'trans')
    param.epsilon=[1, scatt.n_in];
    scatterer = 'trans';
    solver = TDGsolver(kwave,[],scatterer,param);
    disp('Setup TDG solver')
    tic()
    solver.setup();
    toc()
elseif strcmp(scatt.type, 'dir')
    scatterer = 'dir';
    solver = TDGsolver(kwave,[],scatterer,param);
    disp('Setup TDG solver')
    tic()
    solver.setup();
    toc()
else
    error('Non-existent scatterer type!')
end

% setup the T-matrix
disp('Computing T-matrix')
tic()
tmat = ghtmatrix(nmax,kwave,solver,0); %compute always in zero
toc()

if abs(center) > 1e-10
    %set the origin to the scatterer center if necessary
    tmat.setOrigin(center);
end

% save radius of the scatterer
tmat.saveRadius(radius);

return