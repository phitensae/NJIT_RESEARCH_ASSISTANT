function [s, fhat] = hankel_trapezoidal(r,f,betamax)
%assumes r is evenly spaced and f = f(r) 

N = length(r); 
h = r(2) - r(1); 
s = linspace(0,betamax,N); 
fhat = zeros(1,N);

for ii = 1:N
    beta = s(ii); 
    fint = f.*besselj(0,beta*r); %this is the kernel of the hankel transform, we just need to integrate it  
    fa = fint(1);
    fb = fint(end); 
    fsum = fint(2:end-1); 
    hankeltransform = h/2 * (fa + 2*sum(fsum) + fb); 
    fhat(ii) = hankeltransform; 
end