function f  = hot_function_known_a_pl(r_pl,r_covariance,residual_H1,psf,star_position,a_pl)

    %xshift = round(possible_planet_locations(jj,1));
    xshift = (star_position(1) - r_pl(1));
    %xshift = r_pl(1);
    
    %yshift = round(possible_planet_locations(jj,2));
    yshift = (star_position(2) - r_pl(2));
    %yshift = r_pl(2);
    
    %Shift image to planet location
    temp = shift_planet(psf,xshift,yshift,1);
    
    r_temp = get_raster_image(temp);
    
    %a_pl =( (r_temp ./ r_covariance .*  residual_H1) / ...
       %         (r_temp ./ r_covariance .* r_temp) );
                
    %f = ( sum( ( (a_pl*r_temp) ./ r_covariance) ...
       % .*(residual_H1 - 0.5*a_pl*r_temp)));
    
    f = ( sum( ( (a_pl*r_temp) ./ r_covariance) ...
        .*(residual_H1)));% - 0.5*a_pl*r_temp)));
    
end
    
    