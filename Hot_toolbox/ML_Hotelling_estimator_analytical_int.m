function [x_hat,y_hat,int,t] = ML_Hotelling_estimator_analytical_int(data,psf,options,start_point,r_covariance)


residual = get_positive_values_only(data);
r_hotelling_image = get_raster_image(residual);

f_prime = @(r_pl)hot_function_position_analytic_int(r_pl,r_covariance,r_hotelling_image,psf);

[spot_position f_at_position] = fminunc(f_prime,start_point,options);

x_hat =  spot_position(1,1);

y_hat = spot_position(1,2);

psf_tmp = shift_planet(psf,x_hat,y_hat);
r_psf = get_raster_image(psf_tmp);

int = ((r_psf ./ r_covariance) .*  r_hotelling_image) / ...
    ((r_psf ./ r_covariance) .* r_psf);

t = -f_at_position;

end