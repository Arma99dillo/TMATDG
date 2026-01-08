function b = RhsTDGDirMulti(self,c,numscat)

% Build rhs for the TDG method with Dirichlet BCs - multiple scattering

% get parameters
mesh=self.mesh; param=self.param;
t=mesh.t; p=mesh.p; K=param.K; alpha=param.alpha; 
LDir=mesh.LDtN; nd=param.nd; d=param.d;
q=param.q; nodes=param.nodes; w=param.w;

% Define rhs
m = size(t,1);
b = zeros(m*nd,1);
Rot =[0,1;-1,0]; % rotation matrix

% Cycle on DtN edges to apply Dirichlet BC
for L=1:size(LDir,1)
    t1=LDir(L,3); % adjiacent triangle
    p1=p(LDir(L,1),:)'; p2=p(LDir(L,2),:)'; % endpoints
    n = (Rot*(p2-p1))/norm(p2-p1); % normal

    % cycle on PW directions
    for j=1:nd 
        dj = d(:,j);

        % get the incident field as sum of all the fields scattered by
        % other obstacles
        for l=1:numscat
            gd = @(x) c{l}.evaluate(x);
            f = @(t) exp(-1i*K*(p1(1).*dj(1)+t.*(p2(1)-p1(1)).*dj(1)+p1(2)*dj(2)+t.*(p2(2)-p1(2)).*dj(2))).*gd(p1(1)+t.*(p2(1)-p1(1)) + 1i*(p1(2)+t.*(p2(2)-p1(2))));
            int = norm(p2-p1).*gaussint(f,0,1,q,nodes,w);
            b((t1-1)*nd+j)=b((t1-1)*nd+j) + 1i*K.*(dot(dj,n)-alpha).*int;
        end
        % sum global incident field
        gd = @(x) self.uinc.evaluate(x);
        f = @(t) exp(-1i*K*(p1(1).*dj(1)+t.*(p2(1)-p1(1)).*dj(1)+p1(2)*dj(2)+t.*(p2(2)-p1(2)).*dj(2))).*gd(p1(1)+t.*(p2(1)-p1(1)) + 1i*(p1(2)+t.*(p2(2)-p1(2))));
        int = norm(p2-p1).*gaussint(f,0,1,q,nodes,w);
        b((t1-1)*nd+j)=b((t1-1)*nd+j) + 1i*K.*(dot(dj,n)-alpha).*int;
    end 

end


return