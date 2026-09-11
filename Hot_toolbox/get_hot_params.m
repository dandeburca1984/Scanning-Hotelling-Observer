function f = get_hot_params(par,data,psf)

A_star = par(1,1);
a_pl = par(1,2);
b_m = par(1,3);
sigma = par(1,4);
r_pl(1,1) = par(1,5);
r_pl(1,2) = par(1,6);

covariance = A_star*psf + b_m + sigma^2;
r_cov = get_raster_image(covariance);

residual = get_positive_values_only( data - psf*A_star);
r_residual = get_raster_image(residual);

psf_planet = shift_planet(psf,r_pl(1,1),r_pl(1,2),1);
r_psf_planet = get_raster_image(psf_planet);

 f = -( sum( ( (a_pl*r_psf_planet) ./ r_cov) ...
        .*(r_residual - 0.5*a_pl*r_psf_planet)));
    
end
