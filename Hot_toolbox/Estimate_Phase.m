%Get Phase of PSF

clear;
clc;

load('ADI_complicated_target_and_psf.mat');
load('parameters_bright_planet.mat');
clear A_hat; clear a, clear a_max; clear angle_of_rotation;
clear delta_m; clear noise; clear planet_locations; clear radius;
clear ref_psf; clear sigma; clear target_image;
clear time;

lick_pupil = fitsread('lick_pupil.fits');

%Fourier Modulus |F(u,v)|^2
%Fourier_Moduls = sqrt(norm_star_img);

Fourier_Moduls = get_positive_values_only (fitsread('short_exp.fits'));
Fourier_Moduls = sqrt(Fourier_Moduls);

%Support Contstraint
%autocorrelation_of_object = ifftshift( ifft2(norm_star_img));
%abs_auto_corr = abs(autocorrelation_of_object);
%support = abs_auto_corr ./ max(max(abs_auto_corr));
support = lick_pupil;

%g_initial_tmp = zeros(size(Fourier_Moduls));

%for r=0:1:60
    
        %for angle=0:0.01:360
            
            %x = star_position(1,1)+ round( r .* cos(angle));
            %y = star_position(1,1)+  round(r.* sin(angle));
            %support(x,y) = 1;%image(x,y);
            %g_initial_tmp(x,y) = 1;
        %end
     
%end


%Phase Retieval Loop

%Generate Random Noise Image
randn('state', sum(100*clock));

g_inital = get_positive_values_only( randn(length(Fourier_Moduls)) );
g_inital = g_inital .* support;

G_inital = ifftshift( fft2(g_inital) );

%G_inital = get_positive_values_only( randn(length(Fourier_Moduls)) );

G_prime = G_inital .* ( Fourier_Moduls ./ abs(G_inital) ) ;

g_new = ifft2(G_prime);

%g_new = get_positive_values_only (real(g_new)) .* support;
g_new = (abs(g_new)) .* support;

Fourier_Error(1) = ( sum(sum( abs(G_prime - Fourier_Moduls)))) / (sum(sum( Fourier_Moduls)));
display('Begining Phase Loop');

for ii=2:1:100
    
    g_old = g_new;
    G_old = ( fft2(g_old) );
    
    G_prime = G_old .* ( Fourier_Moduls ./ abs(G_old) ) ;
    
%    G_prime = get_positive_values_only(G_prime);
    
    g_new = ifft2(G_prime);
    %g_new = get_positive_values_only (real(g_new)) .* support;
    g_new = ((g_new)) .* support;
    
    Fourier_Error(ii) = sum(sum(  ( (abs(G_prime) - Fourier_Moduls).^2) )) ...
        ./  sum(sum((Fourier_Moduls.^2)));
    
    second_fourier_error(ii) = sum(sum(  (abs( g_old - g_new)).^2));
    
    phase_estimates{ii} = g_new;
    
    Image_Error(ii) =  sum(sum( ( (Fourier_Moduls.^2) - ( abs(G_prime) ).^2 ) ));
    
    %if (Fourier_Error(ii) < 1e-30), break, end
      
    
end
    
recontstructed_psf = (abs( (G_prime))).^2;
recontstructed_psf = recontstructed_psf ./ max(max(recontstructed_psf));
Fourier_Moduls = Fourier_Moduls.^2;
Fourier_Moduls = Fourier_Moduls ./ max(max(Fourier_Moduls));

[i a] = min(Image_Error);

best_phase = phase_estimates{a};
%best_phase = best_phase ./ max(max(best_phase));

phase_error = get_positive_values_only ((  (Fourier_Moduls) - (recontstructed_psf)));

best_psf = ( abs( ifftshift( fft2( best_phase)))).^2;
best_psf = best_psf ./ max(max(best_psf));

best_error = get_positive_values_only ((  (Fourier_Moduls) - (best_psf)));

best_error_metric = sum(sum(best_error));

est_image_error = sum(sum(phase_error));

%figure; subplot(2,2,1); imagesc(log10(Fourier_Moduls));
%subplot(2,2,2); imagesc( abs(g_new));
%subplot(2,2,3); imagesc(log10(recontstructed_psf));
%subplot(2,2,4); imagesc(g_inital);

%figure;subplot(2,2,1);  imagesc( abs(G_prime));
%subplot(2,2,2); imagesc( abs(g_new));
%subplot(2,2,3); imagesc(abs(best_phase));

%figure; plot(Fourier_Error(10:length(Fourier_Error)));
%figure; plot(Image_Error(10:length(Image_Error)));

%figure; imagesc(best_phase(90:170,90:170));
figure; mesh(abs( best_phase(90:170,90:170)));

    
        