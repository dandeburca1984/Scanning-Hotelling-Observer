function f  = hot_function_position_analytic_int(r_pl,r_covariance,residual_H1,psf)

    xshift = r_pl(1,1);
   
    yshift = r_pl(1,2);
 
    temp = shift_planet(psf,xshift,yshift,1);
    
    r_temp = get_raster_image(temp);
    
    a_pl  =( (r_temp ./ r_covariance .*  residual_H1) / ...
                (r_temp ./ r_covariance .* r_temp) );
                
    f = -( sum( ( (a_pl*r_temp) ./ r_covariance) ...
        .*(residual_H1 - 0.5*a_pl*r_temp)));
    
end
    
    