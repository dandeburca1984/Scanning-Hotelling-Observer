%Perform SDI Analysis

clear;
clc;

load('bright_planet_image_seq.mat');
load('parameters_bright_planet.mat');

bright_planet_image_seq = data_cube;
bright_planet_norm_psf = norm_star_img;

%bright_planet_norm_psf=get_positive_values_only (bright_planet_image_seq{1,1} ./ ...
  %  sum(sum(bright_planet_image_seq{1,1})));

rotation_angles = rotation_angle;
bright_planet_primary_brightness = A_hat;

clear data_cube; clear min_angle; clear rotation_angle;
clear shifted_images; clear shifts; clear A_hat; clear a;
clear angle_of_rotation; clear noise; clear sigma; clear star_position;

load('dim_planet_image_seq.mat');
load('parameters_dim_planet.mat');

dim_planet_image_seq = data_cube;
%dim_planet_norm_psf = norm_star_img;

dim_planet_norm_psf = get_positive_values_only( dim_planet_image_seq{1,1} ./ ...
    sum(sum(dim_planet_image_seq{1,1})));

clear data_cube; clear min_angle; clear rotation_angle;
clear shifted_images; clear shifts; clear A_hat; clear a;
clear angle_of_rotation; clear noise; clear sigma; clear star_position;

derotated_img=zeros(size(bright_planet_image_seq{1,1}));

%*************Subtract the images********************

%***************************SDI***********************
for ii=1:1:number_of_images
    
    sdi_seq{1,ii} = get_positive_values_only ( bright_planet_image_seq{1,ii} - dim_planet_image_seq{1,ii});
    
    %Create covariance matrix for each difference image*********
    A_1_tmp = get_estimate_of_A(bright_planet_image_seq{1,ii},bright_planet_norm_psf);
    A_2_tmp = get_estimate_of_A(dim_planet_image_seq{1,ii},dim_planet_norm_psf);
    
    [noise_tmp_1 sigma_tmp_1] = get_noise(bright_planet_image_seq{1,ii},125);  
    [noise_tmp_2 sigma_tmp_2] = get_noise(dim_planet_image_seq{1,ii},125);
    
    covariaces_from_SDI{1,ii} = (sigma_tmp_1+sigma_tmp_2) + (A_1_tmp*bright_planet_norm_psf) + (A_2_tmp*dim_planet_norm_psf);
    
end

%***********************ADI**************************

for kk=1:1:number_of_images
    
    if (kk<=(number_of_images-1))
    
        mean_subtracted_images{1,kk} = get_positive_values_only( sdi_seq{1,kk} - sdi_seq{1,kk+1});
        
        %Create covariance matrix for each difference image*********
        A_1_tmp = get_estimate_of_A(sdi_seq{1,kk},bright_planet_norm_psf);
        
        A_2_tmp = get_estimate_of_A(sdi_seq{1,kk+1},bright_planet_norm_psf);
        
        [noise_tmp_1 sigma_tmp_1] = get_noise(sdi_seq{1,kk},125);
        
        [noise_tmp_2 sigma_tmp_2] = get_noise(sdi_seq{1,kk+1},125);
        
        covariaces_from_ADI{1,kk} = (sigma_tmp_1+sigma_tmp_2) + (A_1_tmp*bright_planet_norm_psf) + (A_2_tmp*bright_planet_norm_psf);
        
    elseif (kk==number_of_images)
            
        mean_subtracted_images{1,kk} = get_positive_values_only( sdi_seq{1,kk} - sdi_seq{1,kk-1});
        
         %Create covariance matrix for each difference image*********
        A_1_tmp = get_estimate_of_A(sdi_seq{1,kk},bright_planet_norm_psf);
        
        A_2_tmp = get_estimate_of_A(sdi_seq{1,kk-1},bright_planet_norm_psf);
        
        [noise_tmp_1 sigma_tmp_1] = get_noise(sdi_seq{1,kk},125);
        
        [noise_tmp_2 sigma_tmp_2] = get_noise(sdi_seq{1,kk-1},125);
        
        covariaces_from_ADI{1,kk} = (sigma_tmp_1+sigma_tmp_2) + (A_1_tmp*bright_planet_norm_psf) + (A_2_tmp*bright_planet_norm_psf);
        
    end
    
end

for hh=1:1:number_of_images
    
    tmp_img = get_positive_values_only( imrotate(mean_subtracted_images{1,hh},rotation_angles(hh),'bicubic','crop'));
    
    derotated_img = derotated_img + tmp_img;
    
end

derotated_img = derotated_img / number_of_images;
r_derotated_img = get_raster_image(derotated_img);

[noise_in_residual sigma_residual] = get_noise(derotated_img,125); 

save('SDI_ADI_final_image.mat','derotated_img');

%**************Hotelling******************************

starting_guess = [105 129];
avg_cov = zeros(size(sdi_seq{1,1}));

for kk=1:1:number_of_images
    
    avg_cov = avg_cov + covariaces_from_ADI{1,kk};
    
end

avg_cov = avg_cov / number_of_images;
r_covariance = get_raster_image(avg_cov);

[max_of_image index]=(max(bright_planet_image_seq{1,1})); 
[max_of_image location] = max(max_of_image);
star_position = [location location];


%***********Options for Minimization****************
options=optimset('LargeScale','off');
options = optimset(options,'TolX',1e-4);
options = optimset(options,'Display','off');

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
searchdir = [];
step_size=[];

%Function handle for Hotelling observer
f = @(r_pl)hot_function(r_pl,r_covariance,r_derotated_img,bright_planet_norm_psf,star_position);

%Minimization function
[hotelling_location] = fminunc(f,starting_guess,options);
save('optimisation_history.mat','history');

xshift = (star_position(1) - hotelling_location(1));
yshift = (star_position(2) - hotelling_location(2));
temp = shift_planet(bright_planet_norm_psf,xshift,yshift,1);
r_temp = get_raster_image(temp);

planet_brightness = ( (r_temp ./ r_covariance .*  r_derotated_img) / ...
                (r_temp ./ r_covariance .* r_temp) );
 

planet_magnitude = 2.5*log10(bright_planet_primary_brightness/planet_brightness);

tmp_image = derotated_img;
%tmp_image(125:133,125:133)=0;
figure; imagesc(tmp_image);
%annotation('ellipse',[0.4036 0.4905 0.0375 0.04762],'Color',[1 0 0]);


%hotelling_data_cube = get_Hotelling_data_map(derotated_img,bright_planet_norm_psf,r_covariance,...
  %  star_position,bright_planet_primary_brightness);

%save('Hotelling_data_map.mat','hotelling_data_cube');

small_psf = bright_planet_norm_psf(90:170,90:170);
small_cov = avg_cov(90:170,90:170);
small_r_cov = get_raster_image(small_cov);
small_star_position = [39 39];

for aa=1:1:number_of_images
    
tmp_image = sdi_seq{1,aa};
%tmp_image(125:133,125:133)=0;

small_derotated_img = derotated_img(90:170,90:170);

small_hot_cube = get_Hotelling_data_map(small_derotated_img,small_psf,small_r_cov,...
    small_star_position,bright_planet_primary_brightness);

load('test_map.mat');

small_hotelling_img_seq{1,aa} = small_hot_cube;
small_hotelling_intensity_img_seq{1,aa} = intensity_map;

end

save('map_sequence.mat','small_hotelling_img_seq','small_hotelling_intensity_img_seq');

%for aa=1:1:number_of_images
    
  %  imagesc(small_hotelling_img_seq{1,aa});
    
%end

