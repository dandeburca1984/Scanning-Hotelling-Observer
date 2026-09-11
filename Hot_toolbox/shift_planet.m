%This is a function to shift a planets position in an x and y dir
%It is based upon applying a linear phase shift in the frequency domain of
%the image


function retval = shift_planet(point_spread_function,xshift,yshift,scale)
%****************Read In PSF***********************************************
psf = point_spread_function;
%********Stretch/Oversample Image******************************************
%psf = imresize(psf,scale,'nearest','Antialiasing',1);

%****************FT the PSF ***********************************************
ft_psf = fft2(psf);

%***********Preallocate Mem for shifed array*******************************
ft_psf_shifted = ft_psf;

%************Size of dimensions of array***********************************
xmax = size(psf,1);

%**************************************************************************
%******Shift X direction****************************************************
for xx=1:1:xmax

%One row at a time is shifted
x = ft_psf(xx,:);

s = xshift;

needtr = 0; 

if size(x,1) == 1; x = x(:); needtr = 1; end;

N = size(x,1); 

r = floor(N/2) + 1; 

f = ((1:N) - r) / (N/2); 

%Linear phase shift
p = exp(-j*s*pi*f)'; 

%Ifftshift is used to apply the shift
y = x .* ifftshift(p) ;

if isreal(x); y = real(y); end;

if needtr; y = y.'; end;

ft_psf_shifted(xx,:) = y;

end

%**************************************************************************
%*****Shift Y direction******************************************************
for xx=1:1:xmax

%One row at a time is shifted
x = ft_psf_shifted(:,xx);

s = yshift;

needtr = 0; 

if size(x,1) == 1; x = x(:); needtr = 1; end;

N = size(x,1); 

r = floor(N/2) + 1; 

f = ((1:N) - r) / (N/2); 

%Linear phase shift
p = exp(-j*s*pi*f)';  

%Ifftshift is used to apply the shift
y = x .* ifftshift(p) ;

if isreal(x); y = real(y); end;

if needtr; y = y.'; end;

ft_psf_shifted(:,xx) = y;

end
%**************************************************************************

%Phase shifted spectrum is Inverse FT to return to space domain
ft_psf_shifted = ifft2(ft_psf_shifted);

shifted_psf = abs(ft_psf_shifted);

%Redude image size to origional
%shifted_psf = imresize(shifted_psf,1/scale,'nearest','Antialiasing',1);

%Return Shifted PSF
retval = shifted_psf;
