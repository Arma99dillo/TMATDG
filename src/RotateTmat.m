function [rotTmat,rotSolver] = RotateTmat(tmat,solver,rotation_angle)

% Rotates T-matrix and solver, giving back the rotated ones

% integrers of the expansion
indices = (-tmat.order:tmat.order)';
D = diag(exp(1i * rotation_angle * indices)); %diagonal matrix

if abs(tmat.origin) > 1e-10 % check if the center is the origin
    center = tmat.origin; % save original center
    tmat.setOrigin(0); % translate Tmatrix to the origin
    rotTmat = copy(tmat);
    rotTmat.matrix = D' * tmat.matrix * D; % rotate T-matrix
    rotTmat.setOrigin(center); tmat.setOrigin(center); % translate back
else % centered in the origin
    rotTmat = copy(tmat);
    rotTmat.matrix = D' * tmat.matrix * D; % rotate T-matrix
end

% rotate mesh and vertices using rotation matrix
Rot = [cos(rotation_angle), -sin(rotation_angle); 
    sin(rotation_angle),  cos(rotation_angle)];
rotSolver = copy(solver);
rotSolver.mesh.p = (Rot*rotSolver.mesh.p(:,:)')';

rotSolver.param.vertices = (Rot*rotSolver.param.vertices(:,:)')';

end
