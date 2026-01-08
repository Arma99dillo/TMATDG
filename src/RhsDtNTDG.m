function b = RhsDtNTDG(self,gd)

% Build rhs for the DtN-TDG method, Dirichlet case

% get parameters
mesh=self.mesh; param=self.param;
t=mesh.t; p=mesh.p; K=param.K; alpha=param.alpha; 
LDir=mesh.LDir; nd=param.nd; d=param.d;
q=param.q; nodes=param.nodes; w=param.w;

% Define rhs
m = size(t,1);
b = zeros(m*nd,1);
Rot =[0,1;-1,0]; % rotation matrix

% Cycle on Dirichlet edges
for L=1:size(LDir,1)
    t1=LDir(L,3); % adjiacent triangle
    p1=p(LDir(L,1),:)'; p2=p(LDir(L,2),:)'; % endpoints
    n = (Rot*(p2-p1))/norm(p2-p1); % normal
    % cycle on PW directions
    for j=1:nd 
        dj = d(:,j);
        f = @(t) exp(-1i*K*(p1(1).*dj(1)+t.*(p2(1)-p1(1)).*dj(1)+p1(2)*dj(2)+t.*(p2(2)-p1(2)).*dj(2))).*gd(p1(1)+t.*(p2(1)-p1(1)) + 1i*(p1(2)+t.*(p2(2)-p1(2))));
        int = norm(p2-p1).*gaussint(f,0,1,q,nodes,w);
        b((t1-1)*nd+j)=b((t1-1)*nd+j) + 1i*K.*(dot(dj,n)-alpha).*int;
    end 

end


return