function final_image = get_image_stack_knowing_rotation(image_sequence,...
    rotation_angle,star_position,psf,base_point,radius,noise,sigma)

    number_of_images = length(image_sequence);
    running_some = zeros(size(image_sequence{1,1}));
    base_image = get_positive_values_only(image_sequence{1,1});
    
    %**********************************************************************
    mask_image = ones(size(base_image));
    
    num_theta_values = 100;
    num_rho_values = 1;
    
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

    for ii=star_position(1)-3:1:star_position(1)+3
        for  kk = star_position(1)-3:1:star_position(1)+3
        mask_image(ii,kk) = 0;
        end
    end
    
    
    %**********************************************************************
    
    
    options=optimset('LargeScale','off');
    options = optimset(options,'TolX',1e-4);
    options = optimset(options,'Display','off');
    
    
    %**********************************************************************

    a_hat_test = get_estimate_of_A(base_image,psf);
    
    covariance = psf * a_hat_test  + sigma^2 + noise;
    
    r_covariance = get_raster_image(covariance);

    residual_H1 = get_raster_image ((base_image));

    %Function handle for Hotelling observer
    f = @(r_pl)hot_function(r_pl,r_covariance,residual_H1,psf,star_position);

    %Minimization function
    [base_point_test] = fminunc(f,base_point,options);

    radius = (norm(star_position - base_point_test));
    %**********************************************************************
    
    for ii=1:1:number_of_images

    xshift =  (radius * cos( (rotation_angle(ii)*pi) / 180)) -radius;
    
    yshift =  (radius*sin( (rotation_angle(ii)*pi) / 180));
    
    test_shifts(ii,1) = xshift;
    test_shifts(ii,2) = yshift;
    
    shift_image = image_sequence{1,ii} .* mask_image;
    
    %shifted_img = shift_planet(shift_image,-xshift,-yshift,1);
    
    shifted_img = shift_planet(shift_image,-xshift,-yshift,1);
    %shifted_img = get_positive_values_only( imrotate(shift_image,rotation_angle(ii),'bicubic','crop'));
    
    running_some = running_some + shifted_img;
    
    
    end
    
    final_image = running_some / number_of_images;
    
    save('test_shifts.mat','test_shifts');
    