function [E,LI,LDtN,LTrans] = MeshBoundTrans(mesh,param)

% Divide edges in internal and boundary edges, saving the triangles
% adjiacent to each edge
% Identify the region where each trinagle is located

% get parameters
p=mesh.p; t=mesh.t; B=mesh.B; vertices=param.vertices;

E=zeros(size(t,1),1); %vector identifying the region
for j=1:size(t,1)
    % compute baricenters
    bar = (p(t(j,1),:) + p(t(j,2),:) + p(t(j,3),:))/3;
    if inpolygon(bar(1),bar(2),vertices(:,1),vertices(:,2)) 
        % check if is inside the polygonal scatterer
        E(j)=2;
    else
        E(j)=1;
    end
end

% Divide in internal and boundaryy edges
control = 0; ni=0; nDtN=0; c=0;
for j = 1:size(t,1)
    
    % triangles edges
    Edges = [ t(j,2) t(j,3); t(j,3) t(j,1); t(j,1) t(j,2)]; %lati triangolo
    
    for l = 1:3   
        if ni+nDtN>=3
            li = [Edges(l,2) Edges(l,1)]; % flip edge
            [~, index]=ismember(LI(:,1:2),li,'rows');
            if sum(index)== 1 % check if I already included this edge
                % in this case, it must be an internal edge
                control = 1;
                % save it in the LI matrix
                v=find(index);
                LI(v,4) = j;
            end
        end

        if control == 0 % if it is not already included
            % check if it is a boundary edge
            if(ismember(Edges(l,1),B) && ismember(Edges(l,2),B))
                c=c+1;
            end
            if c==1
                % in this case all boundary edges are on the circular boundary
                nDtN=nDtN+1;
                LDtN(nDtN,1:2)= Edges(l,:);
                LDtN(nDtN,3)=j;
            else 
                % otherwise, add internal edge to LI
                ni = ni+1;
                LI(ni,1:2) = Edges(l,:);
                LI(ni,3) = j;
            end
        end
        control = 0;
        c=0;
    end
    control=0;
    c=0;
end

% Between internal edges, identify which ones are on the transmission boundary
nTrans=0;
for j=1:size(LI,1)
    t1=LI(j,3); t2=LI(j,4); % adjiacent trinagles
    if E(t1) ~= E(t2)
        % if they triangles belong to different regions the edge is on the
        % transmission boundary
        nTrans = nTrans+1;
        LTrans(nTrans,:)=LI(j,:);
    end
end

return