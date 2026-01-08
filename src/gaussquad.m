function [x,w]=gaussquad(p)
% Computes weights and nodes for Gaussian quadrature on [0,1],
% exact for polynomials up to degree 2p-1
%
% Input:  p : number of quadrature nodes
% Output: x : nodes in [0,1] (column vector)
%         w : weights        (column vector)

b = (1:(p-1)) ./ sqrt(4*(1:(p-1)).^2-1) ;
J = diag(b,-1) + diag(b,1);
[ev,ew] = eig(J);
x = (diag(ew)+1)/2;
w = ((ev(1,:).*ev(1,:)))';
end