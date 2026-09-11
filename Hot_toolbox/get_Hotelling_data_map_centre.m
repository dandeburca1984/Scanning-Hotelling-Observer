function data_map = get_Hotelling_data_map_centre(hotelling_image,psf,r_covariance,...
    star_position)

    %**********************************************************************
    %***Function to calculate the Hotelling observer at every pixel in the
    %image frame and then compute the mean Hotelling noise and standard
    %deviation
   
    image_size = size(hotelling_image);
    
    data_map = zeros(126,126);
    %data_map_beta=data_map;
    
    intensity_map = zeros(image_size);
    
    vector_length = length(data_map);
    
    r_hotelling_image = get_raster_image(hotelling_image);
    
    ii=1;
    jj=1;
    
    %sample_pos = [503 513];
    
    test_shifts_map = zeros(vector_length,2);
    
    for aa=(star_position-15):0.1:(star_position+15)
        
        X_position=aa
        
        ii=1;
        
        for bb=(star_position-15):0.1:(star_position+15)
            
            xshift = (aa - star_position(1));
            
            yshift = (bb - star_position(2));
            
            test_shifts_map(ii,1) = xshift;
            test_shifts_map(ii,2) = yshift;
            
          
            
            %Shift image to planet location
            temp = shift_planet(psf,-xshift,-yshift,1);
    
            r_temp = get_raster_image(temp);
    
            a_pl =( (r_temp ./ r_covariance .*  r_hotelling_image) / ...
                (r_temp ./ r_covariance .* r_temp) );
                
            f = ( sum( ( (a_pl*r_temp) ./ r_covariance) ...
                .*(r_hotelling_image - 0.5*a_pl*r_temp)));
            
            data_map(ii,jj) = f;
            
            %data_map_beta(aa,bb) = f;
            
            intensity_map(ii,jj) = a_pl;
            
            %intensity_map_beta(aa,bb) = a_pl;
            
            ii=ii+1;
            
        end
        
        jj=jj+1;
    end
    
    noise_radius = round(vector_length/2) - round(vector_length/8);
    
    [hot_noise hot_sigma] = get_noise(data_map,noise_radius);
    
    save('test_map.mat','data_map','intensity_map','hot_noise','hot_sigma','test_shifts_map');
    
    save('test_shifts_map.mat','test_shifts_map');