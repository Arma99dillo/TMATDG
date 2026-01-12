function [int1,int2,int3]=Dir_to_Neum(p,L1,L2,d,nd,l,j,k,l1,l2,M,beta_n,coeff,delta,R,q,nodes,w)

% Compute integral of basis functions on the circular bounday with DtN

dl = d(:,l); dj = d(:,j); % PW directions
integers=(-M:1:M)'; % Fourier modes

% get coefficients of basis functions phi_l e phi_j
theta1=mod(atan2(p(l1(1),2),p(l1(1),1)),2*pi); theta2=mod(atan2(p(l1(2),2),p(l1(2),1)),2*pi); % integration extremes
if theta2 < theta1  
    theta2=2*pi+theta2; 
end
theta1_t=mod(atan2(p(l2(1),2),p(l2(1),1)),2*pi); theta2_t=mod(atan2(p(l2(2),2),p(l2(2),1)),2*pi); % integration extremes
if theta2_t < theta1_t  
    theta2_t=2*pi+theta2_t; 
end
coeffl=coeff((2*M+1)*nd*(L1-1)+(l-1)*(2*M+1)+1 : (2*M+1)*nd*(L1-1)+l*(2*M+1));
coeffj=coeff((2*M+1)*nd*(L2-1)+(j-1)*(2*M+1)+1 : (2*M+1)*nd*(L2-1)+j*(2*M+1));

% compute integral of T(phi_l)*phi_j'
f = @(theta) (1- delta.*dj(1).*cos(theta)-delta.*dj(2).*sin(theta)).*exp(1i.*integers.*theta-1i.*k.*(R.*cos(theta).*dj(1)+R.*sin(theta).*dj(2)));
int1=sum(k.*beta_n.*R.*coeffl.*gaussint(f,theta1_t,theta2_t,q,nodes,w));

% compute integral of phi_l*T(phi_j)
g = @(theta) (dl(1).*cos(theta)+dl(2).*sin(theta)).*exp(1i.*k.*(R.*cos(theta).*dl(1)+R.*sin(theta).*dl(2))-1i.*integers.*theta);
int2=sum(k.*beta_n.*R.*(reshape(coeffj',[2*M+1,1])).*gaussint(g,theta1,theta2,q,nodes,w));

% compute integral of T(phi_l)*T(phi_j)'
int3 = sum(k.^2.*abs(beta_n).^2.*coeffl.*(reshape(coeffj',[2*M+1,1])).*2.*pi.*R);


return