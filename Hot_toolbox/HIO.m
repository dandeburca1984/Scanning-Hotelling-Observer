%Input Output algorithm


clear;
clc;


lick_pupil = fitsread('lick_pupil.fits');

Fourier_Moduls = get_positive_values_only (fitsread('short_exp.fits'));
Fourier_Moduls = sqrt(Fourier_Moduls);

SR_PSF =  abs( fftshift( fft2(lick_pupil))) .^2;
SR_PSF = SR_PSF ./ sum(sum(SR_PSF));
SR_PSF_Max = max(max(SR_PSF));

support = lick_pupil;
inverse_support = 1 - support;

%support region
[row,col] = find(support);
pixel_locations = [row,col];

beta = 1;

%First Loop


g{1} = get_positive_values_only( randn(length(Fourier_Moduls)) );
g{1} = support;

G_inital =  fftshift( fft2(g{1}));

G_prime = G_inital .* ( Fourier_Moduls ./ abs(G_inital) ) ;

g{2} = (ifft2(  ifftshift(G_prime)));

g_tmp_1 = get_positive_values_only ( g{2} .* support);

%g_tmp_2 = ( get_positive_values_only( g{1} - beta*g{2})) .*  inverse_support;

g{2} = g_tmp_1;% + g_tmp_2;

Fourier_Error(1) = ( sum(sum( abs(G_prime - Fourier_Moduls)))) / (sum(sum( Fourier_Moduls)));

Image_Error(1) =  sum(sum( ( (Fourier_Moduls.^2) - ( abs(G_prime) ).^2 ) ));

for ii=2:1:50
    
    input = g{ii};
    
    G_old = fftshift (fft2(input));
    
    G_prime = G_old .* ( Fourier_Moduls ./ abs(G_old) ) ;
    
    g_new = ifft2( ifftshift( G_prime));
    
    g_tmp_1 = get_positive_values_only ( g_new .* support);
    
    %g_tmp_2 = ( get_positive_values_only( input - beta*g_new)) .*  inverse_support;
    
    index=ii+1;
    
    g{index} = g_tmp_1;% + g_tmp_2;
    
    Fourier_Error(ii) = sum(sum(  ( (abs(G_prime) - Fourier_Moduls).^2) )) ...
        ./  sum(sum((Fourier_Moduls.^2)));
    
    Image_Error(ii) =  sum(sum( ( (Fourier_Moduls.^2) - ( abs(G_prime) ).^2 ) ));
    
    %if (Image_Error(ii) < eps), break, end
    
end

%figure;plot(Fourier_Error(5:1:length(Fourier_Error)));
%figure; plot(Image_Error(5:1:length(Image_Error)));

%figure; imagesc(abs(g{length(g)}));


for k = 2:length(g)
        surf(abs(g{k}(90:170,90:170)),'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
daspect([5 5 1])
axis tight
view(90,45)
camlight left
 M(k) = getframe;
end

%movie(M);


[i a] = min(abs(Image_Error));

best_phase = g{a};

figure;
surf(abs(best_phase(90:170,90:170)),'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
daspect([5 5 1])
axis tight
view(90,45)
camlight left

psf = Fourier_Moduls.^2;
psf = psf ./ sum(sum(psf));

SR_of_real_image = max(max(psf)) / SR_PSF_Max;

G_best = fftshift( fft2(best_phase));

best_psf = ( abs( G_best)) .^2;
best_psf = best_psf ./ sum(sum(best_psf));


%diff = sum(sum( get_positive_values_only ( psf - best_psf)));

