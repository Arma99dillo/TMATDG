function coeff=FourierCoeff(p,LDtN,k,d,nd,R,M,q,nodes,w)

% Compute Fourier Coefficients of basis functions on circular boundary

toll=1e-10;

inters=(-M:1:M)';
coeff=zeros((2*M+1)*nd*size(LDtN,1),1); v=1;
for L=1:size(LDtN,1)
    p1=p(LDtN(L,1),:)'; p2=p(LDtN(L,2),:)';
    theta1=mod(atan2(p1(2),p1(1)),2*pi); theta2=mod(atan2(p2(2),p2(1)),2*pi);
    if(abs(theta2)<=toll) 
        theta2=2*pi; 
    end
    for l=1:nd
        dl = d(:,l);
        f = @(theta) exp(1i.*k.*(R.*cos(theta).*dl(1)+R.*sin(theta).*dl(2))-1i.*inters.*theta);
        coeff(v:v+2*M)=(1/(2*pi))*gaussint(f,theta1,theta2,q,nodes,w);
        v=v+2*M+1;
    end
end

return