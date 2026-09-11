function final_image = get_image_stack_knowing_rotation_beta(image_sequence,...
    rotation_angle,star_position,radius)

    number_of_images = length(image_sequence);
    running_some = zeros(size(image_sequence{1,1}));
    base_image = get_positive_values_only(image_sequence{1,1});
    
    %**********************************************************************
    mask_image = ones(size(base_image));
    
    num_theta_values = 100;
    num_rho_values = 5;
    
    for ii=1:num_rho_values
    
            num_theta_values = num_theta_values + 1;
    
            for jj=1:num_theta_values
        
                %Radii at which planet can be located 
                radius = (ii-1); %+ diffraction_radius;
        
                %Angles at which planets can be located
                theta = (jj - 1) * (2 * pi / num_theta_values);
                
                x = round(star_position(1) + radius * sin(theta));
        
                y = round(star_position(2) + radius * cos(theta));
        
                mask_image(x,y) = 0;
        
                %index = index + 1;
        
            end
    end

    %for ii=star_position(1)-9:1:star_position(1)+9
        %for  kk = star_position(1)-9:1:star_position(1)+9
        %mask_image(ii,kk) = 0;
        %end
    %end
    
    
    %**********************************************************************
    
   
    for ii=1:1:number_of_images
        
    angle  = -rotation_angle(ii);

    xshift =  (radius * cos( (angle*pi) / 180)) -radius;
    
    yshift =  (radius*sin( (angle*pi) / 180));
    
    test_shifts(ii,1) = xshift;
    test_shifts(ii,2) = yshift;
    
    shift_image = image_sequence{1,ii} .* mask_image;
    
    shifted_img = shift_planet(shift_image,xshift,yshift,1);
    
    running_some = running_some + shifted_img;
      
    end
    
    final_image = running_some / number_of_images;
    
    save('test_shifts.mat','test_shifts');
    