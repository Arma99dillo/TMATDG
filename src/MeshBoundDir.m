function [LI,LDtN,LDir] = MeshBoundDir(mesh,param)

% Divide edges in internal and boundary edges, saving the triangles
% adjiacent to each edge

% get parameters
p=mesh.p; t=mesh.t; B=mesh.B; R=param.R;

% First divide in internal and boundaryy edges
control = 0; ni=0; nb=0; c=0; toll=1e-3;
for j = 1:size(t,1)

    % triangles edges
    Edges = [ t(j,2) t(j,3); t(j,3) t(j,1); t(j,1) t(j,2)]; 
    
    for l = 1:3   
        if ni+nb>=3
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
                % add boundary edge to LB
                nb=nb+1;
                LB(nb,1:2)= Edges(l,:);
                LB(nb,3)=j;
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


% Between boundary edges, I have to distinguish between the ones that are
% on the Dirichlet boundary and the one that are on the circular boundary
nDtN=0; ndir=0;
for j=1:size(LB,1)
    if abs(p(LB(j,1),1)^2+p(LB(j,1),2)^2-R^2)<=toll 
        % if it is on the circular boundary, add it to LDtN
        nDtN=nDtN+1;
        LDtN(nDtN,:)=LB(j,:);
    else 
        % otherwise, add it to LDir
        ndir=ndir+1;
        LDir(ndir,:)=LB(j,:);
    end
end

return