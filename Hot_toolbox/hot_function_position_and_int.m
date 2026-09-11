function f  = hot_function_position_and_int(r_pl,r_covariance,residual_H1,psf)

    xshift = r_pl(1,1);
   
    yshift = r_pl(1,2);
 
    temp = shift_planet(psf,xshift,yshift,1);
    
    r_temp = get_raster_image(temp);
                
    f = -( sum( ( (r_pl(1,3)*r_temp) ./ r_covariance) ...
        .*(residual_H1 - 0.5*r_pl(1,3)*r_temp)));
    
end
    
    