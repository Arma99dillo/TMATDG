function [tmat,solver] = MultiTmat(kwave,h,nd,ScatShape,varargin)

% Computes T-matrices for every shape in ScatShape

% get parameters
Ntype=size(ScatShape,1);
n = nargin;
if n == 5
    M = varargin{1};
else
    M = 0;
end

% compute expansion orders for every shape and save the maximum
orders = zeros(Ntype,1);
for j=1:Ntype
    points = ScatShape{j}.vertices;
    points = points - mean(points, 1); % center in zero to easy compute distance
    distances = sqrt(sum(points.^2, 2)); % distance from the origin
    radius = max(distances);
    orders(j) = suggestedorder(kwave,radius);
end
nmax = max(orders);

% for each type of shape setup solver and compute T-matrix
parfor j=1:Ntype
    if n == 4
        [tmat{j},solver{j}]=ComputeTMatrix(kwave,h,nd,ScatShape{j},nmax);
    else
        [tmat{j},solver{j}]=ComputeTMatrix(kwave,h,nd,ScatShape{j},nmax,M);
    end
end

end