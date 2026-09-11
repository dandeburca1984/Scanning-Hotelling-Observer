function [x_hat,y_hat,t] = ML_Hotelling(data,psf,a_pl,options,start_point)

h_size = round(length(data)/2);
star_position = [h_size h_size];

%[max_of_image index ]=max(data); 
%[max_of_image location] = max(max_of_image);
%star_position_ = [location index(location)];
%start_point = star_position_-0.1;

%A_hat = sum(sum(data));
noise_ring = round(length(data) /2)-5;
[noise sigma] = get_noise(data,noise_ring);  
%covariance = psf * A_hat  + sigma^2 + noise;
covariance=sigma + noise;
r_covariance = get_raster_image(covariance);

%Raster scanned version of the PSF subtracted data
residual = get_positive_values_only(data);
r_hotelling_image = get_raster_image(residual);
    
%Function handle for Hotelling observer
%f_prime =
%@(r_pl)hot_function_beta(r_pl,r_covariance,r_hotelling_image,psf,star_position,a_pl);
%mask = ones(size(psf));
%mask(32:34,32:34) = 0;
%r_mask  = get_raster_image(mask);
r_residual = r_hotelling_image;%.*r_mask;

%f_prime = @(r_pl)hot_function_known_a_pl(r_pl,r_covariance,r_residual,psf,star_position,a_pl);
f_prime = @(r_pl)hot_function_beta(r_pl,r_covariance,r_residual,psf,star_position,a_pl);

%[max_of_image index ]=max(residual_masked); 
%[unused, location] = max(max_of_image);
%start_point = [location index(location)];
%start_point  = start_point -0.9;

%scale_length = round(length(data)/2);
%x1 = [scale_length,scale_length];
x1 = [15,15];
x2=[0 0];

%[spot_position f_at_position] = fmincon(f_prime,start_point,[],[],[],[],x2,x1,[],options);

[spot_position f_at_position] = fminunc(f_prime,start_point,options);

x_hat =  spot_position(1,1);

y_hat = spot_position(1,2);
%t =  hot_noise_function(spot_position,r_covariance,r_hotelling_image,psf,star_position);

%radius_of_noise = round((length(data) /2)-1);
%circumference_of_circle = round(2*pi*radius_of_noise);
%num_theta_values = circumference_of_circle;
%noise_annulus = zeros(num_theta_values,2);
%hot_noise_array = zeros(num_theta_values,1);

 %index=1;

    %    for jj=1:num_theta_values
    
       %             theta = (jj - 1) * (2 * pi / num_theta_values);
    
          %          noise_annulus(index,1) = round(star_position(1) + radius_of_noise * sin(theta));
    
             %       noise_annulus(index,2) = round(star_position(2) + radius_of_noise * cos(theta));
    
                %    r_pl = noise_annulus(index,:);
    
                  %  hot_noise_array(jj) = hot_noise_function(r_pl,r_covariance,r_hotelling_image,psf,star_position,a_pl);
    
                    %index = index + 1;
        %end
        
%hot_noise = mean(hot_noise_array);
%hot_sigma = var(hot_noise_array);
%t1 = -f_at_position;    
%PSNR = (t1 - hot_noise) / sqrt( 0.5*( var(t1) + hot_sigma) );
%t=t1;
t = -f_at_position;

end