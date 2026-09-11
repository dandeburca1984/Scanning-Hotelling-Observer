function [noise, sigma] = get_noise(image,radius)


image_size = size(image);

half_image_size = image_size(1)/2;

image_centre = [image_size(1)/2 image_size(1)/2];

num_theta_values = 300;

noise_annulus = zeros(num_theta_values,2);

noise_array = zeros(num_theta_values,1);

index=1;

for jj=1:num_theta_values
    
    theta = (jj - 1) * (2 * pi / num_theta_values);
    
    noise_annulus(index,1) = round(image_centre(1) + radius * sin(theta));
    
    noise_annulus(index,2) = round(image_centre(2) + radius * cos(theta));
    
    index = index + 1;
    
end

index = 1;

tmp = image;

for ii = 1:num_theta_values
    
    noise_array(ii) = image(noise_annulus(ii,1),noise_annulus(ii,2));
    
    tmp(noise_annulus(ii,1),noise_annulus(ii,2)) = max(max(image));
    
    index = index + 1;
    
    
end

%figure;imagesc(tmp);colormap gray; axis square
%h = gcf;
%saveas(h,'noise_estimation_locations','jpg') ;
%movefile('noise_estimation_locations.jpg','test_results/');
%delete(h);

noise_distribution = fit_ML_normal(noise_array);

noise = noise_distribution.u;

sigma = noise_distribution.sig2;

end
