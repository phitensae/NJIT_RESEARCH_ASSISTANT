function cc = soundspeedprofile(zz) 
%interpolates the sound speed profile data and evaluates it at zz. need to
%check that zz is in the range of sound speed profile 
load littledorrit.mat; 
cc = interp1(z,ssp(:,2),zz); 
c = ssp(:,2); %the sound speed profile data 


for ii = 1:length(cc)
    if isnan(cc(ii)) %check to see if the sound speed profile interpolant yields nan
        if zz(ii) < min(z) %if where we are interpolating is less than data
            cc(ii) = min(ssp(:,2)); 
        else 
            cc(ii) = max(c); %if where we are interpolating is greater than data
    end
end

end