function val = gaussint(f,a,b,q,nodes,w)

% Computes Gaussian quadrature
val = 0;
for l=1:q
    val = val + w(l)*f(a+(b-a).*nodes(l));
end

val = (b-a).*val;

return