function [t_direct, t_surface,t_bottom] = calculation_arrivaltime(r,z_s,z_r,WD,c)
%this tries to create an isovelocity profile and uses that to calculate the
%arrival times

z_M = max(max(z_r,z_s),WD); 

if z_M ~= WD
    fprintf('Note: the water column depth is not the maximum. This is a problem')
    t_direct = NaN; 
    t_surface = NaN; 
    t_bottom = NaN;
    return 
end
c_mean = quadgk(c,0,z_M)/z_M; %approximates sound speed profile with isovelocity 

d_direct = sqrt(r^2 + (z_s - z_r)^2); 
t_direct = d_direct/c_mean; 

d_surface = sqrt(r^2 + (z_s + z_r)^2); 
t_surface = d_surface/c_mean; 

d_bottom = sqrt(r^2 + (2*WD - z_r - z_s)^2); 
t_bottom = d_bottom/c_mean;
end