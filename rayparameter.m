function p = rayparameter(r,z_s,z_r,c)
%calculates the ray parameter from the depth of the source and the
%receiever, the range and the sound speed profile

c_H = (z_s - z_r)*quadgk(c,z_s,z_r);
p_0 = (c_H *1/(1+(z_r - z_s)/r)^2)^(-1); 

p = p_0; 
temp = p_0;

for ii = 1:100
    f = @(z) (c(z)).^2./(1-p^2.*(c(z)).^2).^(3/2); 
    range_fn = @(z) p.*c(z)./(1-p.^2.*(c(z)).^2).^(1/2);
    r_p = quadgk(range_fn,z_s,z_r); %integration performed in formula 
    drbydp = quadgk(f,z_s,z_r); %integration in formula
    p = p_0 + 1/drbydp * (r - r_p); %p has been updated
    temp = abs(quadgk(range_fn,z_s,z_r));
    if abs(temp - r) < 1e-3
        disp('Succesfully found ray parameter p')
        break
    end
end