function [x_hat,y_hat,int,t] = ML_Hotelling_estimator(data,psf,options,start_point,r_covariance)

%noise_ring = round(length(data) /2)-5;
%[noise sigma] = get_noise(data,noise_ring);  
%covariance = sigma^2 + noise;
%r_covariance = get_raster_image(covariance);
residual = get_positive_values_only(data);
r_hotelling_image = get_raster_image(residual);

f_prime = @(r_pl)hot_function_position_and_int(r_pl,r_covariance,r_hotelling_image,psf);

[spot_position f_at_position] = fminunc(f_prime,start_point,options);

x_hat =  spot_position(1,1);

y_hat = spot_position(1,2);

int = spot_position(1,3);

t = -f_at_position;

end