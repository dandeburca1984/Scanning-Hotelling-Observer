clear;
clc;

load('bright_planet_image_seq.mat');
load('parameters_bright_planet.mat');

bright_planet_image_seq = data_cube;

[max_of_image index]=(max(bright_planet_image_seq{1,1})); 
[max_of_image location] = max(max_of_image);
image_centre = [location location];

rmin=0;
rmax=30;

image = bright_planet_image_seq{1,1};
masked_image = zeros(size(image));


for xx=1:1:length(image)
    
    for yy=1:1:length(image)
        
        x1 = image_centre(1,1) - xx;
        y1 = image_centre(1,2) - yy;
        
        rho = sqrt( x1.^2 + y1.^2);
        
        if (rho >= rmin) && (rho<=rmax)
            
            masked_image(xx,yy) = image(xx,yy);
            
        end
        
    end
    
end

load('optimized_psf_test.mat');

test = get_positive_values_only( masked_image - Sr);

figure; imagesc(Sr); axis square;
title(' \bf Reconstructed PSF','Fontsize',14);
figure;imagesc(masked_image);axis square;
title('\bf Original PSF','Fontsize',14);

r_test = get_raster_image(test);

mean_test = mean(r_test);
variance_test = var(r_test);

figure; imagesc(test); axes square

%surf(test,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');