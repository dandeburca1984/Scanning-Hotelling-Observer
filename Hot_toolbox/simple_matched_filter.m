function f = simple_matched_filter(r_pl,psf,residual,a_pl)

 xshift = r_pl(1);
 yshift = r_pl(2);
 
 temp = shift_planet(psf,xshift,yshift,1);
 r_temp = get_raster_image(temp);
 
  
    f = ( sum( ( (a_pl*r_temp)) ...
        .*(residual)));
    
end