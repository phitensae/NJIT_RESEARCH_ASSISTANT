% it is important to remember in all of this that the hankel transform is
% like the fourier transform and so when you do hankel.m it will only
% capture the first 100 frequencies (morally speaking). invhankel.m seeks
% to capture the radial values so, again, it should go out to something
% like 200 m (or something like that). 
%% this section is to test the function hankel.m

clc 
clear all
close all

h = @(x) 1./x; 
r = linspace(.25,500,1000);
f = h(r); 

[beta, fbeta] = hankel(r,f); 
subplot(2,2,1)
error1 = abs(fbeta - h(beta)); 
plot(beta,error1,'b')
title('error of hankel transform using 1/x')
xlabel('\beta')
ylabel('error')

g = @(r) exp(-r.^2); 
gsample = g(r); 
[s, fs] = hankel(r,gsample); 
error2 = abs(fs - 1/2*g(s/2));
subplot(2,2,2)
plot(s,error2,'k')
title('error of hankel transform using gaussian')
xlabel('\beta')
ylabel('error')


e_r = exp(-r); 
[omega, eomega] = hankel(r,e_r); 
t = @(s) 1./(s.^2 + 1^2).^(3/2); %the hankel transform of e^(-ar) is a/(s^2 + a^2)^(3/2)
error3 = abs(eomega - t(omega)); 
subplot(2,2,3)
plot(omega,error3,'r')
title('error of hankel transform using exponential')
xlabel('\beta')
ylabel('error')


%% this section is to test the function invhankel.m
clc
close all
clear all

h = @(x) 1./x; 
s = linspace(.25,500,3000); 
f0 = h(s); 

[r, f] = invhankel(s, f0); 
error = abs(f - h(r)); 
subplot(2,2,1)
plot(r,error)
title('error of hankel transform using 1/x')

g = @(s) 1/2*exp(-s.^2/4);
g0 = g(s); 
[rho, frho] = invhankel(s,g0);
q = @(r) exp(-r.^2);
error1 = abs(frho - q(r)); 
subplot(2,2,2)
plot(r,error1)
title('error of hankel transform using gaussian')

lol = @(s) 1./(s.^2 + 1^2).^(3/2); 
lol0 = lol(s); 
[rawr, frawr] = invhankel(s,lol0); 
error3 = abs(frawr - exp(-rawr)); 
subplot(2,2,3)
plot(rawr,error3)


%% testing the error of the code (relative)
clear all
clc
close all

g = @(r) exp(-r.^2); 
h = @(r) 1/2*g(r/2); %h is the hankel transform of g 
Nmax = 1000;
SAVED = zeros(1,Nmax-1); %create a matrix to hold the maximum error 
relative = zeros(1,Nmax-1); %create a matrix to hold the relative error

for N = 2:Nmax
    r = linspace(0,100,N); 
    f = g(r); 
    [s,fs] = hankel(r,f);
    true = h(s); 
    error = true - fs; %note that error, true and fs are all vectors 
    SAVED(N-1) = norm(error,Inf); %at n-1 because we went from 2 to Nmax, note that saved also has the maximum error
    relative(N-1) = norm(error,Inf)/norm(true,Inf); 
end

n = 2:Nmax; 
plot(n,SAVED) 
xlabel('number of gridpoints, N')
ylabel('max error')
figure
h = 100./(n-1); 
plot(h,relative)
ylabel('relative error = max(error)/max(true)')
xlabel('step size, h = 1/N')
figure 
plot(n,relative)
xlabel('N')
ylabel('relative error')

%% testing the function hankel_trapezoidal.m 
close all
clear all
clc

h = @(x) 1./x; 
r = linspace(.25,500,3000);
f = h(r); 

[beta, fbeta] = hankel_trapezoidal(r,f,10); 
error = abs(fbeta - h(beta)); 
plot(beta,error,'b')
title('error in hankel transform of 1/x')
xlabel('\beta')
ylabel('error')

g = @(r) exp(-r.^2); 
gsample = g(r); 
[beta, gbeta] = hankel_trapezoidal(r,gsample,10); 
error = abs(gbeta - 1/2*g(beta/2)); 
figure 
plot(beta,error,'b')
title('error in the hankel transform of e^{-r^2}')
xlabel('\beta')
ylabel('error')

%% find out how bad things get on hankel_trapezoidal.m
close all
clear all
clc

N = 100:10:1000; 
l = length(N); 

for jj = 1:l
    r = linspace(.25,500,N(jj)); 
    g = @(r) exp(-r.^2); 
    gsample = g(r); 
    [beta, gbeta] = hankel_trapezoidal(r,gsample,10); 
    error = abs(gbeta - 1/2*g(beta/2)); 
    plot(beta,error,'b')
    title('error in the hankel transform of e^{-r^2}')
    xlabel('\beta')
    ylabel('error')
    hold on
end

maxerror = zeros(1,l);
for jj = 1:l 
    r = linspace(.25,500,N(jj));
    g = @(r) exp(-r.^2); 
    gsample = g(r); 
    [beta, gbeta] = hankel_trapezoidal(r,gsample,10); 
    error = abs(gbeta - 1/2*g(beta/2)); 
    maxerror(jj) = max(error);
end

figure
loglog(N,maxerror)
xlabel('N')
ylabel('maximum error at each iteration of N')