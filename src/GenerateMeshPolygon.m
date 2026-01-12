function [p_ord,t_ord,I,B] = GenerateMeshPolygon(param)

% Generate mesh with a Dirichlet polygonal obstacle

% get parameters
h=param.h; R=param.R; vertices=param.vertices;

% Create a PDE model
model = createpde();

% Define the geometry 
n=size(vertices,1);
g2 = [2;n;vertices(:,1);vertices(:,2)]; % polygon
g1 = [1,0,0,R,zeros(1,size(g2,1)-4)]; % circle

geom_matrix = [g1' g2]; set_formula = 'g1-g2';    

% Add the geometry to the PDE model
geometryFromEdges(model, decsg(geom_matrix, set_formula, char('g1', 'g2')'));

% Generate the mesh with Hmax set to h
mesh = generateMesh(model, 'Hmax', h,'GeometricOrder','linear');
[p1,e,t1] = meshToPet(mesh); p=p1'; t=t1(1:3,:)';

% Plot the mesh
figure()
pdeplot(model);


n_nodes = size(p,1);
n_triangles = size (t,1);

% reorder triangles vertices
for i_t=1:n_triangles
    v=p(t(i_t,:),:);
    if det([v(:,1) v(:,2) ones(3,1)]) < 0
        t(i_t,[1 2 3])= t(i_t,[2 1 3]);
    end
end

% reorder triangle nodes such as the internal nodes are the first
old2new = zeros(n_nodes,1); new2old = zeros(n_nodes,1);
i_int = 1; i_bound = n_nodes;
for j = 1:n_nodes
    % for every node, I check if it is a boundary or internal node
    if ismember(j,e(1,:)) || ismember(j,e(2,:))
        old2new(j) = i_bound;
        new2old(i_bound) = j;
        i_bound = i_bound - 1;
    else
        old2new(j) = i_int;
        new2old(i_int) = j;
        i_int = i_int + 1;
    end
end

% reorder p and t
p_ord = p(new2old,:);
t_ord = old2new(t);
B=(i_bound+1):n_nodes; % boundary nodes
I=1:(i_int-1); % internal nodes


end