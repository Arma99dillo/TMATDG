function A = MatrixTDGDirTrans(self)

% Build the TDG matrix with Dirichlet condition on the circular boundary

% get parameters
mesh=self.mesh; param=self.param;
t=mesh.t; p=mesh.p; B=mesh.B; K=param.K; LI=mesh.LI; LDtN=mesh.LDtN; LTrans=mesh.LTrans;
alpha=param.alpha; beta=param.beta; epsilon=param.epsilon;
E=mesh.E; nd=param.nd; d=param.d;

% Define the matrix
m = size(t,1);
A = spalloc(m*nd,m*nd,10*m*nd);

Rot =[0,1;-1,0]; % Rotation matrix
t = [t t(:,1)];

% First cycle on triangles to compute internal edges terms 
for T = 1:m 

    k=K*sqrt(epsilon(E(T))); % wavenumber in the triangle

    % find internal edges
    nB=0; nI=0; L_I=zeros(3,2); L_B=zeros(3,2);
    for u=1:3
        if(ismember(t(T,u),B) && ismember(t(T,u+1),B))
            nB=nB+1;
            L_B(nB,:)=[t(T,u), t(T,u+1)];
        else
            nI=nI+1;
            L_I(nI,:)=[t(T,u), t(T,u+1)];
        end
    end
    
    % cycle on PW directions
    A_aux = zeros(nd,nd); % auxiliary matrix
    for s=1:nI
        li = [L_I(s,1) L_I(s,2)]; lig= [L_I(s,2) L_I(s,1)];
        p1=p(L_I(s,1),:)'; p2=p(L_I(s,2),:)'; % endpoints
        n = (Rot*(p2-p1))/norm(p2-p1); % normal
        % check if the edge lies on the transition boundary
        [~, index]=ismember(LTrans(:,1:2),li,'rows'); [~, indexg]=ismember(LTrans(:,1:2),lig,'rows');
        if  sum(index) == 1 || sum(indexg) == 1 
            % if it is on the internal edge define xi as the average
            xi=0.5*(real(K*sqrt(epsilon(1))) + real(K*sqrt(epsilon(2))));
        else
            % otherwise xi is just k
            xi=real(k);
        end
        for l = 1:nd
            for j= 1:nd
                dl = d(:,l); dj = d(:,j); diff=dl-dj;
                A_aux(j,l) = A_aux(j,l) + (1i*k.*(-(1/2).*(dot(dj,n)+dot(dl,n)))-1i.*(1/xi).*k.^2.*beta.*dot(dj,n).*dot(dl,n)-1i.*xi.*alpha).*phi_int(diff,p1,p2,k);
            end
        end
    end
    A((T-1)*nd+1:T*nd,(T-1)*nd+1:T*nd)=A((T-1)*nd+1:T*nd,(T-1)*nd+1:T*nd)+A_aux;
end

% Cycle on internal edges
for L=1:size(LI,1)

    % adjiacent triangles
    t1=LI(L,3); t2=LI(L,4); 
    p1=p(LI(L,1),:)'; p2=p(LI(L,2),:)'; % endpoints
    n = (Rot*(p2-p1))/norm(p2-p1); % normal

    % wavenumber in the triangles
    e1=epsilon(E(t1)); e2=epsilon(E(t2)); k1=K*sqrt(e1); k2=K*sqrt(e2); 
    xi=0.5*(real(k1)+real(k2)); % value of xi
    
    % cycle on PW directions for the first triangle  
    A_aux =zeros(nd,nd); % auxiliary matrix
    for l=1:nd 
        for j=1:nd
        dl = d(:,l); dj = d(:,j); diff=k1.*dl-k2.*dj;
        A_aux(j,l) = A_aux(j,l) + (0.5*(1i*k2*dot(dj,n)+1i*k1*dot(dl,n))+alpha*1i*xi+beta*1i*(1/xi)*k1*k2*dot(dl,n)*dot(dj,n)).*phi_int(diff,p1,p2,1);           
        end
    end
    A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)=A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)+A_aux;  
    
    % cycle on PW directions for the second triangle
    n=-n; k2=K*sqrt(e1); k1=K*sqrt(e2); % change normal and wavenumbers  
    A_aux =zeros(nd,nd); % auxiliary matrix
    for l=1:nd 
        for j=1:nd
        dl = d(:,l); dj = d(:,j); diff=k1.*dl-k2.*dj;
        A_aux(j,l) = A_aux(j,l) + (0.5*(1i*k2*dot(dj,n)+1i*k1*dot(dl,n))+alpha*1i*xi+beta*1i*(1/xi)*k1*k2*dot(dl,n)*dot(dj,n)).*phi_int(diff,p1,p2,1);           
        end
    end
    A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)=A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)+A_aux;   

end

% Cycle on DtN edges (circular boundary) for Dirichlet BC
k = K*sqrt(epsilon(1)); %wavenumber outside the obstacle
for L=1:size(LDtN,1)
    t1=LDtN(L,3); % adjiacent triangle
    p1=p(LDtN(L,1),:)'; p2=p(LDtN(L,2),:)'; % endpoints
    n = (Rot*(p2-p1))/norm(p2-p1); %normal

    % cycle on PW directions
    A_aux =zeros(nd,nd); % auxiliary matrix
    for l=1:nd 
        for j=1:nd
        dl = d(:,l); dj = d(:,j); diff=dl-dj;
        A_aux(j,l) = A_aux(j,l) + 1i*k.*(-alpha-dot(dl,n)).*phi_int(diff,p1,p2,k);           
        end
    end
    A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)=A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)+A_aux; 
end

return