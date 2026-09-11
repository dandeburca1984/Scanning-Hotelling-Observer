function final_image = get_image_stack_using_hotelling(base_image, image_sequence,...
    r_covariance,psf,star_position,A_hat,base_point)


running_some = zeros(size(image_sequence{1,1}));
number_of_images = length(image_sequence);


%Mask to supress the central maxima due to the primary star

mask_image = ones(size(base_image));

for ii=star_position(1)-3:1:star_position(1)+3
    
      for  kk = star_position(1)-3:1:star_position(1)+3
    
        mask_image(ii,kk) = 0;
    
      end
    
end

options=optimset('LargeScale','off');
options = optimset(options,'TolX',1e-4);
options = optimset(options,'Display','off');

%r_base_image = get_raster_image(base_image);

%**************FIND STAR POSITION******************************************
%Function handle for Hotelling observer
%f_prime = @(r_pl)hot_function_find_star(r_pl,r_covariance,r_base_image,psf,star_position,A_hat);

%Minimization function
%[star_position] = fminunc(f_prime,star_position,options);

%**************************************************************************
%**************************************************************************


%***************FIND COMPANION POSITON*************************************
%Initial Position 
base_image = get_positive_values_only( base_image.* mask_image) ;

%[max_of_image index]=(max( base_image )); 

%[max_of_image location] = max(max_of_image);

%base_point = [location index(location)];

residual_H1 = get_raster_image ( get_positive_values_only (base_image));

%Function handle for Hotelling observer
f = @(r_pl)hot_function(r_pl,r_covariance,residual_H1,psf,base_point);

%Minimization function
[base_point] = fminunc(f,base_point-1,options);

%**************************************************************************

for kk=1:1:number_of_images
    
    tmp_image_sequence = get_positive_values_only( image_sequence{1,kk}.* mask_image );
    
    [tmp_max_of_image tmp_index]=(max( get_positive_values_only( tmp_image_sequence) )); 

    [tmp_max_of_image tmp_location] = max(tmp_max_of_image);

    input_point = [tmp_location tmp_index(tmp_location)]; 
    
    %**********************************************************************
    %************FIND COMPANION POSITION USING HOTELLING
    
    a_hat = get_estimate_of_A(tmp_image_sequence,psf);
    
    residual_H1 = get_raster_image ( get_positive_values_only ( tmp_image_sequence));
    
    %Function handle for Hotelling observer
    f = @(r_pl)hot_function_find_star(r_pl,r_covariance,residual_H1,psf,input_point,a_hat);

    %Minimization function
    [input_point] = fminunc(f,input_point-1,options);
    %**********************************************************************
   
    xshift = (base_point(1) - input_point(1));
    
    yshift = (base_point(2) - input_point(2));
    
    hot_shifts(kk,1) = xshift;
    hot_shifts(kk,2) = yshift;
    
    shifted_img = shift_planet(image_sequence{1,kk},-xshift,-yshift,1);
   
    output_point(kk,1) = input_point(1) + xshift;
    output_point(kk,2) = input_point(2) + yshift;
    
    
    running_some = running_some + shifted_img;
    
    %shifted_image_seq{1,kk} = shifted_img;
    
    
end

save('hot_shifts.mat','hot_shifts','output_point');%,'shifted_image_seq');

final_image = running_some / number_of_images;


