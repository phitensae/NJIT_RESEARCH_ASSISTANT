clear all
close all
clc

c = @(z) soundspeedprofile(z); 

zz = 0:.01:100; 
ff = c(zz); 

plot(zz,ff,'b')

load littledorrit.mat

hold on 

z = ssp(:,1); 
ssprofile = ssp(:,2); 
mean_profile = mean(ssprofile); 

plot(z,ssprofile,'kx',z,mean_profile*ones(size(z)),'r-')

legend('interpolated data','actual data','mean profile')
xlabel('z')
ylabel('c(z)')

reciever_depth = [20 23 28 31 35 38 42 46 50 54 58 61 65 68];
R = 223; 
source_depth = 26; 
watercolumn = 100; 

direct = zeros(1,length(reciever_depth)); 
surface = zeros(1,length(reciever_depth)); 
bottombounce = zeros(1,length(reciever_depth)); 

figure 
ii = 1; 
for z_r = reciever_depth
    [t1 t2 t3] = calculation_arrivaltime(R,source_depth,z_r,watercolumn,c); 
    hydrophone = -z_r * ones(1,3); 
    arrivaltimes = [t1 t2 t3]; 
    plot(arrivaltimes,hydrophone,'x')
    hold on
    direct(ii) = t1; 
    surface(ii) = t2;
    bottombounce(ii) = t3; 
    ii = ii + 1; 
end

figure 

plot(direct,-reciever_depth,'kx',surface,-reciever_depth,'bx',bottombounce,-reciever_depth,'rx')
legend('direct','surface','bottom bounce')
xlabel('time')
ylabel('depth')