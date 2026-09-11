function [x_hat, y_hat] = Hotelling_Filter(data,psf,n,pad)

%Compute the prewitening filter
%pad=0;
%A_hat = get_estimate_of_A(data,psf); %estimate of the intensity 
A_hat = 1;
noise_ring = (length(data) /2)-2;
[noise sigma] = get_noise(data,noise_ring);  
covariance = psf * A_hat  + sigma^2 + noise;

%data_tmp = get_positive_values_only(data-0.5*psf*A_hat)./covariance;% - 0.5*A_hat*psf);
%data_tmp = get_positive_values_only((data-noise)./covariance);
data_tmp = get_positive_values_only(data./covariance);

mask_tmp = (psf);

fourier_hot_mask=conj(ifftshift(fft2(fftshift(((mask_tmp))))));
fourier_hot_residual = (ifftshift(fft2(fftshift((data_tmp)))));

Hot_ =  fourier_hot_mask.*fourier_hot_residual;
Hot_ = padarray(Hot_,[pad pad],0);
Hot = (ifftshift(ifft2(fftshift(Hot_))));
Hot=sqrt(Hot.*conj(Hot));%This is the correlation.

A = Hot;
    
%[a,b]=find(A==max(max(A)));
[a,b]= find(Hot ==max(max(Hot)));
a=min(a);b=min(b);%just in case.

if (a==1)||(b==1)
    a=a+1;
    b=b+1;
end
 
%Parabolic interpolation of the position of the maxima, as described by
%Poyneer (Applied Optics 2003).

a_=a-(n+2*pad+3)/2+0.5.*(A(a-1,b)-A(a+1,b))./(A(a-1,b)+A(a+1,b)-2.*A(a,b));% this formula in particular is valid for uneven sizes of images.

%dbstop if caught error;

b_=b-(n+2*pad+3)/2+0.5.*(A(a,b-1)-A(a,b+1))./(A(a,b-1)+A(a,b+1)-2.*A(a,b));

x_hat=a_.*n./(n+2*pad);

y_hat=b_.*n./(n+2*pad);

[ring_avg ring_var] = get_noise(Hot,noise_ring);

PSNR_hot = (max(max(Hot)) - ring_avg) / (sqrt(ring_var));

%mesh(Hot,'FaceColor','interp','FaceLighting','phong'); camlight right;

end

