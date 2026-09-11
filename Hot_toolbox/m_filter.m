function [estimated_col,estimated_row,PSNR_A,t]=m_filter(N,f_mask,n,pad)

A_=(ifftshift(fft2(fftshift((N))))).*f_mask;

A_=padarray(A_,[pad pad],0);
A=(ifftshift(ifft2(fftshift(A_))));
A=sqrt(A.*conj(A));%This is the correlation.

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
 
estimated_row=a_.*n./(n+2*pad);
estimated_col=b_.*n./(n+2*pad);

noise_ring = (length(A) /2)-1;
[ring_avg ring_var] = get_noise(A,noise_ring);
PSNR_A = (max(max(A)) - ring_avg) / (sqrt(ring_var));
t = max(max(A));

end
 
