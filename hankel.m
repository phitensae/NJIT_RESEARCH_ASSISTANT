function [s,fhat] = hankel(r,f) 
%%suppose that you have a vector r, which contains points at which a
%%function f is sampled radially. i.e. f = f(r) and we sampled it at the
%%points contained in the vector r and its values are at f. then in theory,
%%hankel(r,f) will return a vector s and fhat, corresponding to the hankel
%%transform of f = f(r); fhat = fhat(s). it is assumed that r is in order,
%%and that r has uniform spacing 

%obtain the length of r 
N = length(r); 
s = linspace(0,90,N); 
fhat = zeros(1,N); 
dr = r(2) - r(1); 

%s, r, fhat and f are all things that have length N. 

for ii = 1:N 
    beta = s(ii);
    g = f.*besselj(0,beta*r).*r.*dr; 
    fhat(ii) = sum(g);
end