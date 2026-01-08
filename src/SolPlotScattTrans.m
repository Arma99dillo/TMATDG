function SolPlotScattTrans(self,u)

% Plots the scattered field inside the circle

% get parameters
ui = @(x) self.uinc.evaluate(x); mesh = self.mesh; param = self.param;
t1=mesh.t; p1=mesh.p; nd=param.nd; d=param.d; K=param.K; 
E=mesh.E; epsilon=param.epsilon;
phi = @(x1,x2,d,k) exp(1i*k.*(x1.*d(1)+x2.*d(2))); %basis functions

% Cycle on triangles
for T=1:size(t1,1)

    k=K*sqrt(epsilon(E(T))); % wavenumber inside the triangle
    
    % get solution coefficients on the triangle
    coeff=u((T-1)*nd+1:T*nd);
    
    % mesh on the triangle
    pv = [p1(t1(T,1),:); p1(t1(T,2),:); p1(t1(T,3),:)];
    g = [2;3;pv(:,1);pv(:,2)]; 
    model = createpde();
    geometryFromEdges(model, decsg(g));
    mesh = generateMesh(model, 'Hmax', 0.01,'GeometricOrder','linear');
    p = mesh.Nodes'; t = mesh.Elements'; 

    % compute numerical solution as a linear combination of plane waves
    u_app=zeros(size(p,1),1);
    for v=1:size(p,1)
        for j=1:nd
            dj=d(:,j);
            u_app(v)=u_app(v)+coeff(j).*phi(p(v,1),p(v,2),dj,k);
        end
        if E(T) == 1
            %subtract incident wave only outside the obstacle
            u_app(v) = u_app(v)+ui(p(v,1)+1i*p(v,2));
        end
    end

    % plot solution
    trisurf(t,p(:,1),p(:,2),real(u_app)); hold on
    
end

shading flat; colorbar; view(2)

end