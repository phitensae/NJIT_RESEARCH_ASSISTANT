function [s,fhat] = hankel_riemann(r,f,betamax)
%we assume that r is an equally spaced vector and that f = f(r). betamax is
%the maximum value of s that matters 

if length(r) == 1
    s = NaN;
    fhat = NaN;
    return 
end

N = length(r);
s = linspace(0,betamax,N);
fhat = zeros(1,N); 
dr = r(2) - r(1); 

%s, r, fhat and f are all things that have length N. 

for ii = 1:N 
    beta = s(ii);
    g = f.*besselj(0,beta*r).*r.*dr; 
    fhat(ii) = sum(g);
end