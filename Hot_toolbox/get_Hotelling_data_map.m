function data_map = get_Hotelling_data_map(hotelling_image,psf,r_covariance,...
    star_position,sampling,a_pl)

    %**********************************************************************
    %***Function to calculate the Hotelling observer at every pixel in the
    %image frame and then compute the mean Hotelling noise and standard
    %deviation
   
    image_size = size(hotelling_image);
    
    data_map = zeros(image_size/sampling);
   % data_map_beta=data_map;
    
    %intensity_map = zeros(image_size);
    
    vector_length = length(hotelling_image);
    
    r_hotelling_image = get_raster_image(hotelling_image);
    
    %ii=1;
    jj=1;
    
    test_shifts_map = zeros(vector_length,2);
    
    %r_PSF = get_raster_image(psf);
    
    %xsize = length(psf);
    
    for aa=1:sampling:vector_length
        
        XX = aa
        
        ii=1;
        
        xshift = (aa - star_position(1));
        
        for bb=1:sampling:vector_length
            
            yshift = (bb - star_position(2));
            
            test_shifts_map(ii,1) = xshift;
            test_shifts_map(ii,2) = yshift;
       
            %Shift image to planet location
            temp = shift_planet(psf,-xshift,-yshift,1);
    
            r_temp = get_raster_image(temp);
            %r_temp = shift_planet_raster(r_PSF,xsize,-xshift,-yshift);
            %temp2 = get_un_raster_image_two(r_temp);
    
            a_pl_test =( (r_temp ./ r_covariance .*  r_hotelling_image) / ...
                (r_temp ./ r_covariance .* r_temp) );
                
            f = ( sum( ( (a_pl_test*r_temp) ./ r_covariance) .* r_hotelling_image));
                %.*(r_hotelling_image - 0.5*a_pl*r_temp)));
            
            data_map(ii,jj) = f;
            
            %ttt(ii) = f;
            
            %data_map_beta(aa,bb) = f;
            
            %intensity_map(ii,jj) = a_pl;
            
            %intensity_map_beta(aa,bb) = a_pl;
            
            ii=ii+1;
            
        end
        
        jj=jj+1;
    end
    
    %noise_radius = round(vector_length/2) - round(vector_length/8);
    
    %[hot_noise hot_sigma] = get_noise(data_map,noise_radius);
    
    %save('test_map.mat','data_map','intensity_map','hot_noise','hot_sigma','test_shifts_map');
    
    %save('test_shifts_map.mat','test_shifts_map');
    
end