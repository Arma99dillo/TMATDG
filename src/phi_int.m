function val = phi_int(d,p1,p2,k)

% Compute the integral of the function exp(i*k*dot(x,d)) on the segment [p1,p2]

if norm(d)<=1e-10
    val = norm(p2-p1);
elseif abs(dot(d,p2-p1)) <= 1e-10
    val= norm(p2-p1).*exp(1i*k*dot(p1,d));
else
    val = (norm(p2-p1)/(1i*k*dot((p2-p1),d))).*exp(1i*k*dot(p1,d)).*expm1(1i*k*dot(p2-p1,d));
end

return