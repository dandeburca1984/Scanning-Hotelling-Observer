function map = get_Hotelling_values(Hot_img,psf,K,star_position)

image_size = size(Hot_img);
    
data_map = zeros(image_size);
    
intensity_map = zeros(image_size);
    
vector_length = length(data_map);
    
residual = get_raster_image(Hot_img);

r_covariance = get_raster_image(K);
 
 %Define R the set of test locations
  index=1;
    
    for xx=1:1:vector_length
        
            for yy=1:1:vector_length
                
                possible_planet_locations(index,1) =  xx;
                
                possible_planet_locations(index,2) = yy;
                
                index = index + 1;
            end
    end
    
    
    %Number of possible planet locations
    num_planet_locations = length(possible_planet_locations);
    
    %*********************************************************
    
    %Define Hotelling variables
    a_pl = zeros(1,num_planet_locations);
    Hot = zeros(1,num_planet_locations);
    
    for jj=1:num_planet_locations
    
        xshift = (star_position(1) - possible_planet_locations(jj,1));
    
        yshift = (star_position(2) - possible_planet_locations(jj,2));
    
        %Shift image to planet location
        temp = shift_planet(psf,xshift,yshift,1);
     
        %**************PSF at Loc of Companion*********************************
        r_temp = get_raster_image(temp); 
    
        %***********Opt intensity at Loc of Companion**************************
        a_pl(jj) =( (r_temp ./ r_covariance .*  residual) / ...
                (r_temp ./ r_covariance .* r_temp) );
    
        %***********Hotelling observer at Loc of Companion*********************
        Hot(jj) = ( sum( ( (a_pl(jj)*r_temp) ./ r_covariance) ...
            .*(residual - 0.5*a_pl(jj)*r_temp)));
    
    end
    
    x = possible_planet_locations(:,1);
    y = possible_planet_locations(:,2);
    xlin = linspace(min(x),max(x),image_size(1));
    ylin = linspace(min(y),max(y),image_size(2));
    [X,Y] = meshgrid(xlin,ylin);
    Z_Hot_test_stats = griddata(x,y,Hot,X,Y,'cubic');
    %Z_Hot_test_stats(125:135,125:135)=0;
    Z_intensity_est = griddata(x,y,a_pl,X,Y,'cubic');
    %Z_intensity_est(125:135,125:135)=0;
    save('Hot_maps.mat','Z_Hot_test_stats','Z_intensity_est');
    
    map = Z_Hot_test_stats;
    
 