function mesh = GenerateMesh(self)

% Generate mesh for the DtN-TDG method
% Distinguish between DIrichlet and transition problem
if strcmp(self.scatterer,'dir')

    % build mesh
    param=self.param;
    [mesh.p,mesh.t,mesh.I,mesh.B] = GenerateMeshPolygon(param);
    
    % get boundary and internal edges
    [mesh.LI,mesh.LDtN,mesh.LDir] = MeshBoundDir(mesh,param);

    self.type='dir'; %set scatterer type
else

    % build mesh
    param=self.param;    
    [mesh.p,mesh.t,mesh.I,mesh.B] = GenerateMeshPolygonTrans(param);
    
    % get boundary and internal edges and identify regions
    [mesh.E,mesh.LI,mesh.LDtN,mesh.LTrans] = MeshBoundTrans(mesh,param);

    self.type='trans'; %set scatterer type
end

return