%% solve the forward problem - direct path 
clc
clear all
close all

r = 100; 
z_s = 20; 
z_r = 28; 

 c = @(z) soundspeedprofile(z); %this creates the sound speed profile 
% c = @(z) 1000 + 0.*z; 
p = rayparameter(r,z_s,z_r,c); %finds the ray parameter p 

tau = @(z) 1./(c(z).* (1 - p^2*c(z).^2).^(1/2)); 

arrivaltime = quadgk(tau,z_s,z_r);
fprintf('The arrival time for the direct path is %2.5f \n',arrivaltime)

[t_direct, t_surface] = calculation_arrivaltime(r,z_s,z_r,c);
fprintf('The arrival time (as calculate with by using simple geometry) for the direct path is %2.6f \n',t_direct) 