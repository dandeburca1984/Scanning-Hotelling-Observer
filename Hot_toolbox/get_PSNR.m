function P_SNR = get_PSNR(image,planet_loc,image_centre)

%num_theta_values = 300;

index=1;

radius = norm(planet_loc - image_centre);

rad_conversion = pi/180;

for jj=1:1:361
    
    %theta = (jj-1 ) * (2*pi / 360);
    
    theta = (jj-1) * rad_conversion;
    
    noise_annulus(index,1) = planet_loc(2)+  round( radius * sin(theta));
    
    noise_annulus(index,2) = planet_loc(1) + round(radius * cos(theta));
    
    index = index + 1;
    
end

index = 1;

tmp = image;

for kk=1:1:360
    
    noise_array_full(kk) = image(noise_annulus(kk,1),noise_annulus(kk,2));
    
end

for ii=1:1:360
     
    noise_array(index) = image(noise_annulus(ii,1),noise_annulus(ii,2));
    
    tmp(noise_annulus(ii,1),noise_annulus(ii,2)) = max(max(image));
    
    index = index + 1;
    
end


 I_peak = max(noise_array_full);
 I_mean = mean(noise_array);
 I_std  =  std(noise_array);
 
 P_SNR = (I_peak - I_mean) / I_std;
end