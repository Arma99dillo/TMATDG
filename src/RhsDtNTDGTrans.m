function b = RhsDtNTDGTrans(self,gd,gn)

% Build rhs for the DtN-TDG method, transmission case

% get parameters
mesh=self.mesh; param=self.param;
t=mesh.t; p=mesh.p; K=param.K; alpha=param.alpha; beta=param.beta;
LTrans=mesh.LTrans; nd=param.nd; d=param.d;  E=mesh.E;
q=param.q; nodes=param.nodes; w=param.w; epsilon=param.epsilon;

% Define rhs
m = size(t,1);
b = zeros(m*nd,1);
Rot =[0,1;-1,0]; % rotation matrix

% Cycle on transmission edges
for L=1:size(LTrans,1)
    % adjiacent triangles and wavenumbers
    t1=LTrans(L,4); t2=LTrans(L,3);
    k1=K*sqrt(epsilon(E(t1))); k2=K*sqrt(epsilon(E(t2)));
    xi=0.5*(real(k1)+real(k2));   
    p2=p(LTrans(L,1),:)'; p1=p(LTrans(L,2),:)'; % endpoints
    n_gamma = (Rot*(p2-p1))/norm(p2-p1); % normal pointing outward of the polygon  
    n_K = n_gamma; % normal pointing outward of the mesh element
    % cycle on PW directions for the first triangle
    for j=1:nd 
        dj = d(:,j);
        f = @(t) exp(-1i*k1*(p1(1).*dj(1)+t.*(p2(1)-p1(1)).*dj(1)+p1(2)*dj(2)+t.*(p2(2)-p1(2)).*dj(2))).*gd(p1(1)+t.*(p2(1)-p1(1)) + 1i*(p1(2)+t.*(p2(2)-p1(2))));
        g = @(t) exp(-1i*k1*(p1(1).*dj(1)+t.*(p2(1)-p1(1)).*dj(1)+p1(2)*dj(2)+t.*(p2(2)-p1(2)).*dj(2))).*sum((gn(p1(1)+t.*(p2(1)-p1(1)) + 1i*(p1(2)+t.*(p2(2)-p1(2)))).*n_gamma));
        int1 = norm(p2-p1).*gaussint(f,0,1,q,nodes,w); 
        int2 = norm(p2-p1).*gaussint(g,0,1,q,nodes,w);
        b((t1-1)*nd+j)=b((t1-1)*nd+j) - dot(n_K,n_gamma).*(-0.5.*1i.*k1.*dot(dj,n_K)+1i.*alpha.*xi).*int1 + ...
            -(beta.*k1.*(1/xi).*dot(dj,n_K)-0.5)*int2;
    end 
    
    p1=p(LTrans(L,1),:)'; p2=p(LTrans(L,2),:)'; % endpoints
    n_K = -n_K; % flip normal to the element
    % cycle on PW directions for the second triangle
    for j=1:nd 
        dj = d(:,j);
        f = @(t) exp(-1i*k2*(p1(1).*dj(1)+t.*(p2(1)-p1(1)).*dj(1)+p1(2)*dj(2)+t.*(p2(2)-p1(2)).*dj(2))).*gd(p1(1)+t.*(p2(1)-p1(1)) + 1i*(p1(2)+t.*(p2(2)-p1(2))));
        g = @(t) exp(-1i*k2*(p1(1).*dj(1)+t.*(p2(1)-p1(1)).*dj(1)+p1(2)*dj(2)+t.*(p2(2)-p1(2)).*dj(2))).*sum((gn(p1(1)+t.*(p2(1)-p1(1)) + 1i*(p1(2)+t.*(p2(2)-p1(2)))).*n_gamma));
        int1 = norm(p2-p1).*gaussint(f,0,1,q,nodes,w); 
        int2 = norm(p2-p1).*gaussint(g,0,1,q,nodes,w);
        b((t2-1)*nd+j)=b((t2-1)*nd+j) - dot(n_K,n_gamma).*(-0.5.*1i.*k2.*dot(dj,n_K)+1i.*alpha.*xi).*int1 + ...
            -(beta.*k2.*(1/xi).*dot(dj,n_K)-0.5)*int2;
    end 

end

return