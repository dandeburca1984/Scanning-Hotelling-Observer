function J = phase_diversity_joint_est_n_wavelengths(zernike_coeff,data,wavelengths,S,pupil)

%Written by: Daniel Burke, 12/5/10
%Function to perfrom a least squares fit to the pupil phase
%using n wavelength (n must at least equal 2)

%zernike_coeff = vector of zernike coefficients
%data = data cube containing n images
%wavelengths = vector of wavelength values
%S = structure containing telescope parameters

%Error metric for each wavelength
Image_plane_error = zeros(1,length(wavelengths));

n = length(zernike_coeff);

%For loop to run though each wavelength
for ii=1:1:length(wavelengths);
    
    zernike_coeff_tmp = zernike_coeff(1:(n-1)) * (wavelengths(1) / wavelengths(ii));
    
    %H_PSF = Gen1Image(S.lam,S.D,S.scale,256,real(zernike_coeff_tmp));
    H_PSF = Generate_Image(S.lam(ii),S.D,S.scale,257,zernike_coeff_tmp,pupil);
    H_PSF = H_PSF ./ sum(sum(H_PSF));
    
    Image_plane_error(ii) = sum(sum( ( data(:,:,ii) - zernike_coeff(n)*H_PSF ).^2));
    
end

J = sum(Image_plane_error);