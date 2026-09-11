function [x_hat, y_hat,PSNR_hot,t,CC] = Hotelling_Filter2(data,psf,n,pad,noise,sigma,A_hat,star_hat)

%Function to use the Hotelling Observer via interprelating the maximum of
%its cross correlation with the data

%Compute the prewitening filter
noise_ring = (length(data) /2)-5;
[n_ s_] = get_noise(data,noise_ring);  
%c1 = psf * star_hat;%  + s_^2 + n_;
%c2 =  sigma^2 + noise;
c2 = s_^2;% + n_;
%covariance = (c1 + c2)/2;
%covariance=1;

covariance=c2;

hot_residual = data ./ c2;

hot_mask = get_positive_values_only(A_hat*psf);

%hot_mask = A_hat*psf;
%hot_residual = get_positive_values_only (data./c1);
%hot_residual = hot_residual - psf*star_hat;
%pre_whit = hot_residual ./ covariance;
fourier_hot_mask=conj(ifftshift(fft2(fftshift(((hot_mask))))));
fourier_hot_residual = (ifftshift(fft2(fftshift((hot_residual)))));

Hot_ =  fourier_hot_mask.*fourier_hot_residual;
Hot_ = padarray(Hot_,[pad pad],0);
Hot = (ifftshift(ifft2(fftshift(Hot_))));
Hot=sqrt(Hot.*conj(Hot));%This is the correlation.

A = Hot;

search_radius = round(length(A)/5);
star_centre = [ round(length(A)/2) round(length(A)/2)];

K = A( (star_centre(1)-search_radius):(star_centre(1)+search_radius), (star_centre(2)-search_radius):(star_centre(2)+search_radius));

[a,b]=find(A==max(max(K)));

a=min(a);b=min(b);%just in case.

if (a==1); a = a+1;end
if (b==1); b = b+1; end   

if(a==length(A)); a = a-1; end
if(b==length(A)); b = b-1; end


%Parabolic interpolation of the position of the maxima, as described by
%Poyneer (Applied Optics 2003).

a_=a-(n+2*pad+3)/2+0.5.*(A(a-1,b)-A(a+1,b))./(A(a-1,b)+A(a+1,b)-2.*A(a,b));% this formula in particular is valid for uneven sizes of images.

b_=b-(n+2*pad+3)/2+0.5.*(A(a,b-1)-A(a,b+1))./(A(a,b-1)+A(a,b+1)-2.*A(a,b));
 
x_hat=a_.*n./(n+2*pad);
y_hat=b_.*n./(n+2*pad);

%noise_ring_2 = (length(Hot)/2)-1;
%[ring_avg ring_var] = get_noise(Hot,noise_ring_2);
PSNR_hot =1;% (max(max(Hot)) - ring_avg) / (sqrt(ring_var));
t = max(max(K));

CC = K;

end

