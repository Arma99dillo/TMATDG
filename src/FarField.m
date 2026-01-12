function val = FarField(self,vec,index)

% Evaluate far-field in the directions specified in vec for the incident
% wave specified by index

% get parameters
mesh=self.mesh; param=self.param; u=self.coeffs(:,index);
LDtN=mesh.LDtN; nd=param.nd; d=param.d; p=mesh.p; K=param.K; R=param.R;
nodes=param.nodes; w=param.w; q=param.q;

% Cycle on the evaluation points
val=zeros(length(vec),1);
for j=1:length(vec)
    
    ang=vec(j); d1=cos(ang); d2=sin(ang); %direction angle
    % Cycle on cirular boundary to compute integral
    for L=1:size(LDtN,1)
        t1=LDtN(L,3); % adjiacent triangle
        coeff=u((t1-1)*nd+1:t1*nd); % PW expansion coefficients

        p1=p(LDtN(L,1),:)'; p2=p(LDtN(L,2),:)'; % endpoints
        theta1=mod(atan2(p1(2),p1(1)),2*pi); theta2=mod(atan2(p2(2),p2(1)),2*pi);
        if theta2 < theta1 
            theta2=2*pi+theta2;
        end
        % Cycle on PW directions
        for l = 1:nd
            dl = d(:,l);
            f = @(theta) (cos(theta)*(d1+dl(1))+sin(theta)*(d2+dl(2)))*exp(1i*K*(R*cos(theta)*(dl(1)-d1)+R*sin(theta)*(dl(2)-d2)));
            val(j) = val(j) - (exp(1i*pi/4)/sqrt(8*pi*K))*1i*K*R*coeff(l)*gaussint(f,theta1,theta2,q,nodes,w);
        end
    end

end

return