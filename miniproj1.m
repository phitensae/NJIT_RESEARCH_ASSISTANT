% this is an investigation into how does noise affect the hankel transform
% and it is done in two phases. the first phase is to see how the hankel
% transform is affected by noise of 0 mean and different variances. the
% second is with regards to the problem outlined in A Direct Method for the
% Estimation of Sediment Sound Speed with a Horizontal Array in Shallow
% Water

%% hankel transform on 0 and then we use some additive white gaussian noise

clc
close all
clear all

r = 0:.25:200; 
c = zeros(size(r)); 
xcorrupt1 = awgn(c,1/100);
xcorrupt2 = awgn(c,1/10); 
xcorrupt3 = awgn(c,1); 
xcorrupt4 = awgn(c,10); 
xcorrupt5 = awgn(c,100); 

[beta1,f1]  = hankel(r,xcorrupt1); 
[beta2, f2] = hankel(r,xcorrupt2); 
[beta3, f3] = hankel(r,xcorrupt3); 
[beta4, f4] = hankel(r,xcorrupt4); 
[beta5,f5] = hankel(r,xcorrupt5); 

plot(beta1,f1,'r-',beta2,f2,'k:',beta3,f3,'b--',beta4,f4,'g-.',beta5,f5,'m-')
legend('SNR = 1/100','SNR = 1/10','SNR = 1','SNR = 10','SNR = 100')

figure 

subplot(2,3,1)
plot(beta1,f1)
title('SNR = 1/100')

subplot(2,3,2)
plot(beta2,f2)
title('SNR = 1/10')

subplot(2,3,3)
plot(beta3,f3)
title('SNR = 1')

subplot(2,3,4)
plot(beta4, f4)
title('SNR = 10')

subplot(2,3,5)
plot(beta5,f5)
title('SNR = 100')

%% hankel transform on a gaussian and effect of noise 
h = @(r) exp(-r.^2); 
r = linspace(0,500,5000);
f = h(r); %hence f is a signal of known gaussians

[s0, fs] = hankel(r,f); 

fcorrupt1 = awgn(f,100);
fcorrupt2 = awgn(f,10); 
fcorrupt3 = awgn(f,1); 
fcorrupt4 = awgn(f,1/10); 
fcorrupt5 = awgn(f,1/100); 

[s1, fscrewed1] = hankel(r,fcorrupt1); 
[s2, fscrewed2] = hankel(r,fcorrupt2); 
[s3, fscrewed3] = hankel(r,fcorrupt3); 
[s4, fscrewed4] = hankel(r,fcorrupt4);
[s5, fscrewed5] = hankel(r,fcorrupt5);

figure

plot(s0,fs,'k-',s1,fscrewed1,'b-',s2,fscrewed2,'g-',s3,fscrewed3,'r-',s4,fscrewed4,'m-',s5,fscrewed5,'y-')
legend('true hankel transform','AWGN = 100','AWGN = 10','AWGN = 1','AWGN = 0.1','AWGN = .01')

figure
subplot(2,3,1) 
plot(s0,fs,'k-',s1,fscrewed1,'b+')
legend('true hankel transform','AWGN = 100')

subplot(2,3,2)
plot(s0,fs,'k-',s2,fscrewed2,'g+')
legend('true hankel transform','AWGN = 10')

subplot(2,3,3)
plot(s0,fs,'k-',s3,fscrewed3,'r+')
legend('true hankel transform','AWGN = 1')

subplot(2,3,4)
plot(s0,fs,'k-',s4,fscrewed4,'m+')
legend('true hankel transform','AWGN = 0.1')

subplot(2,3,5)
plot(s0,fs,'k-',s5,fscrewed5,'y+')
legend('true hankel transform','AWGN = .01')


e1 = abs(fs - fscrewed1); 
e2 = abs(fs - fscrewed2); 
e3 = abs(fs - fscrewed3); 
e4 = abs(fs - fscrewed4); 
e5 = abs(fs - fscrewed5); 

figure 
subplot(2,3,1)
plot(s0,e1)
title('error when SNR = 100')

subplot(2,3,2)
plot(s0,e2)
title('error when SNR = 10')

subplot(2,3,3)
plot(s0,e3)
title('error when SNR = 1')

subplot(2,3,4)
plot(s0,e4)
title('error when SNR = 0.1')

subplot(2,3,5)
plot(s0,e5)
title('error when SNR = 0.01')

subplot(2,3,6)
plot(s0,e1,s0,e2,s0,e3,s0,e4,s0,e5)
legend('SNR = 100','SNR = 10','SNR = 1','SNR = 0.1','SNR = 0.01')

%% solving the pressure field when sound speed profile given is constant
clear all
clc
close all

c = @(z) 150*z./z; %speed profile is a constant 
%the reason why it's written like this is to make sure that the code runs 
%and is dependent on the size of the vector z
SSP = c; %ssp stands for sound speed profile 
f = 10; %using a frequency f = 10 Hz 
w = 2*pi*f; 

k = @(z) w./c(z);

% beta_max = 90; 
z_max = 30; 
% M = 20; %for the number of beta/r entries
% N = 1000; %for the number of z entries 
% 
% b = linspace(0,beta_max,M);
% z = linspace(0,z_max,N); 
% dz = z(2) - z(1); %step size in z
% db = b(2) - b(1); %step size in beta 
% b = linspace(0,beta_max,M);
% z = linspace(0,z_max,N); 
% as it turns out, by the CFL condition, this code will never produce the
% correct solution. we need to edit this such that it does. 
dz = 0.1; %this number is arbitrary for now

v = zeros(N,M);

%solving for the v part
for jj = 1:M
    beta = b(jj); 
    q = @(z) k(z).^2 - beta^2; 
    zsetup = z(2:end-1); 
    qsetup = q(zsetup); 
    d1 = [1,dz^2*qsetup-2,1];
    d2 = [ones(1,N-2),0]; 
    d3 = [0, ones(1,N-2)];
    A = diag(d1) + diag(d2,-1) + diag(d3,1); 
    delta = zeros(N,1);
    delta(1) = 1; 
    v(:,jj) = A\delta; 
end

%getting the u part back 
u = zeros(N,M); 


for ii = 1:N
    [r,u(ii,:)] = invhankel(b,v(ii,:)); 
end
plot(r,u(1,:))
figure 
 
%adding noise to the u 
ucorrupt1 = awgn(u,100); 
ucorrupt2 = awgn(u,10); 
ucorrupt3 = awgn(u,1); 
ucorrupt4 = awgn(u,.1); 

vcorrupt1 = zeros(size(ucorrupt1)); 
vcorrupt2 = zeros(size(ucorrupt2));
vcorrupt3 = zeros(size(ucorrupt3)); 
vcorrupt4 = zeros(size(ucorrupt4)); 

%taking the hankel transform of ucorrupt to get vcorrupt 
for ii = 1:N
    [s1,vcorrupt1(ii,:)] = hankel(r,ucorrupt1(ii,:)); 
    [s2,vcorrupt2(ii,:)] = hankel(r,ucorrupt2(ii,:));
    [s3,vcorrupt3(ii,:)] = hankel(r,ucorrupt3(ii,:));
    [s4,vcorrupt4(ii,:)] = hankel(r,ucorrupt4(ii,:)); 
end
%comparing v to vcorrupt

error1 = abs(vcorrupt1 - v);
error2 = abs(vcorrupt2 - v); 
error3 = abs(vcorrupt3 - v);
error4 = abs(vcorrupt4 - v); 

subplot(2,2,1)
plot(b,error1(1,:),'b',b,error1(end,:),'k')
legend('error at z = 0','error at z = 50')
title('error plot for SNR = 100')
xlabel('\beta or s - transformed variable')

subplot(2,2,2)
plot(b,error2(1,:),'b',b,error2(end,:),'k')
legend('error at z = 0','error at z = 50')
title('error plot for SNR = 10')
xlabel('\beta or s - transformed variable')

subplot(2,2,3)
plot(b,error3(1,:),'b',b,error3(end,:),'k')
legend('error at z = 0','error at z = 50')
title('error plot for SNR = 1')
xlabel('\beta or s - transformed variable')

subplot(2,2,4)
plot(b,error4(1,:),'b',b,error4(end,:),'k')
legend('error at z = 0','error at z = 50')
title('error plot for SNR = 0.1')
xlabel('\beta or s - transformed variable')
