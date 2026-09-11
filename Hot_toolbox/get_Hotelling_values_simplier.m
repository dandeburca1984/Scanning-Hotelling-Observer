%new Hot function

function data_map = get_Hotelling_values_simplier(Hot_img,psf,K,star_position)

image_size = size(Hot_img);
    
data_map = zeros(image_size);
    
intensity_map = zeros(image_size);
    
vector_length = length(data_map);
    
residual = get_raster_image(Hot_img);

r_covariance = get_raster_image(K);
 
 %Define R the set of test locations
  index=1;
    
    for xx=90:1:170
        
            for yy=90:1:170
                
                %possible_planet_locations(index,1) =  xx;
                
                %possible_planet_locations(index,2) = yy;
                xshift = (star_position(1) - xx);
                yshift = (star_position(2) - yy);
                
                %Shift image to planet location
                temp = shift_planet(psf,xshift,yshift,1);
     
                %**************PSF at Loc of Companion*********************************
                r_temp = get_raster_image(temp); 
                
                
                 %***********Opt intensity at Loc of Companion**************************
                a_pl(index) =( (r_temp ./ r_covariance .*  residual) / ...
                        (r_temp ./ r_covariance .* r_temp) );
                
                
                %***********Hotelling observer at Loc of Companion*********************
                Hot(index) = ( sum( ( (a_pl(index)*r_temp) ./ r_covariance) ...
                    .*(residual - 0.5*a_pl(index)*r_temp)));
                
                data_map(yy,xx) = Hot(index);
                
                intensity_map(yy,xx) = a_pl(index);
                
                index = index + 1;
            end
    end
    
    save('Hot_maps.mat','data_map','intensity_map');