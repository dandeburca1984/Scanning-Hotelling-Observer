function results = Scanning_Hotelling_Observer


clear;
clc;

addpath('Hot_toolbox/');
addpath('New_toolbox/');
warning off;

options=optimset('LargeScale','on');
options = optimset(options,'TolX',2*(eps));
options = optimset(options,'TolF',eps);
options = optimset(options,'Display','off');

cd('data/HD_235089/');
data = fitsread('best_image_PSF38.fits');
data = get_positive_values_only(data);
data = data ./ sum(sum(data));

psf = fitsread('best_image_PSF39.fits');
psf = get_positive_values_only(psf);
psf = psf ./ sum(sum(psf));

cd('../../');

%% Locate the approxiamte location of the parent star
% Peak of Cross-Corellation between data and PSF, approx location of Parent
%Star
[col_hat, row_hat] = m_filter_2(data,psf);   

%% Find exact location of Parent Star using conjugate gradient minimisation 
start_point = [col_hat, row_hat];                         % Start point of minimisation
covariance_matrix =ones(64,64);                        % Assumed covariance of ones
r_cov = get_raster_image(covariance_matrix);    % raster scan of covariance matrix    
f_prime = @(r_pl)hot_function_position_analytic_int(r_pl,r_cov,get_raster_image(data),psf); 
[spot_position f_at_position] = fminunc(f_prime,start_point,options); %conjugate gradient minimisation  
shifted_psf = shift_planet(psf,spot_position(1,1),spot_position(1,2));  %shift the PSF to the location of the parent star
shifted_psf = shifted_psf / sum(sum(shifted_psf)); %normalised shifted PSF
parent_star_intensity_hat = get_raster_image(data) / get_raster_image(shifted_psf); %estimate the intensity of the Parent star, it should be very close to 1 
data_minus_parent = get_positive_values_only( data - shifted_psf*parent_star_intensity_hat);    %subtract estimated Parent star signal from the data

%% Estimate the `flat' background level  
n_ring = round(length(data_minus_parent)/2)-2;
[background_hat, variance_hat] = get_noise(data_minus_parent,n_ring);

%% Use Gaussian Covariance model
covariance_matrix =shifted_psf*parent_star_intensity_hat+ background_hat+  variance_hat^2;
r_cov = get_raster_image(covariance_matrix);

%% Locate the approximate location of the faint companion
[col_hat, row_hat] = m_filter_2(data_minus_parent,psf); 

%% Find the exact location of the companion star
start_point =round(abs( [col_hat, row_hat])) ;                         % Start point of minimisation
[x_hat,y_hat,int_hat,test_stat] = ML_Hotelling_estimator_analytical_int(data_minus_parent,psf,options,start_point,r_cov);

%% update the background model
companion_psf = shift_planet(psf,x_hat,y_hat);
companion_psf = companion_psf / sum(sum(companion_psf));
updated_background = get_positive_values_only( data_minus_parent - companion_psf*int_hat);
[b2, v2] = get_noise(updated_background,n_ring);
new_cov  = shifted_psf*parent_star_intensity_hat + b2 + v2;
r_new_cov = get_raster_image(new_cov);
[x_hat2,y_hat2,int_hat2,test_stat2] = ML_Hotelling_estimator_analytical_int(data_minus_parent,psf,options,start_point,r_new_cov);
delta_m2 =  2.5*log10(parent_star_intensity_hat / int_hat2);

%% Construct Output variable
results.parent_location = spot_position;
results.covariance_matrix = new_cov;
results.background_level = b2;
results.background_variance = v2;
results.companion_location = [x_hat2,y_hat2];
results.companion_intensity = int_hat2;
results.companion_differential_magnitude = delta_m2;

%% Define subfunctions
function positive_image = get_positive_values_only_beta(image)
image(image<0) = 0;
positive_image = image;
end

function [estimated_col,estimated_row]=m_filter_2(data,psf)
A_=(ifftshift(fft2(fftshift((data))))).*(ifftshift(fft2(fftshift((psf)))));
n= length(data);
pad = 5*n;
A_=padarray(A_,[pad pad],0);
A=(ifftshift(ifft2(fftshift(A_))));
A=sqrt(A.*conj(A));%This is the correlation.
[a,b]=find(A==max(max(A)));
a=min(a);a = round(a);b=min(b);b = round(b);%just in case.
%Parabolic interpolation of the position of the maxima, as described by
%Poyneer (Applied Optics 2003).
a_=a-(n+2*pad+3)/2+0.5.*(A(a-1,b)-A(a+1,b))./(A(a-1,b)+A(a+1,b)-2.*A(a,b));% this formula in particular is valid for uneven sizes of images.
b_=b-(n+2*pad+3)/2+0.5.*(A(a,b-1)-A(a,b+1))./(A(a,b-1)+A(a,b+1)-2.*A(a,b));
estimated_row=a_.*n./(n+2*pad);
estimated_col=b_.*n./(n+2*pad);
end


tx = linspace (-1, 1, length(data));
ty = tx;
[xx, yy] = meshgrid (tx, ty);
r = sqrt (xx .^ 2 + yy .^ 2);
r2 = r(r==1,r==1);
% noise_pixel_indices = (r=radius); 


end