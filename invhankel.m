function [r,f] = invhankel(s,f0)
%% this function seeks to take in f0 = f0(s) and return the inverse hankel transform of the pair
%% assume equal spacing in s, 

N = length(s); 
r = linspace(0,100,N); 
ds = s(2) - s(1); 
f = zeros(1,N); 

for ii = 1:N
    beta = r(ii); 
    g = f0.*s.*besselj(0,beta*s).*ds;
    f(ii) = sum(g); 
end