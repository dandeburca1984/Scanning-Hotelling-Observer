function return_val = shift_planet_raster(r_image,xsize,x_shift,y_shift)

%r_image = get_raster_image(image);
%xsize = length(image);

%Fourier Transform the image
x = (fft(r_image));

%Preallocate memory for FT image
%x = ft_r_image;

%total shift
s = x_shift;

needtr = 0; 

if size(x,1) == 1; x = x(:); needtr = 1; end;

N = length(x); 

r = floor(N/2) + 1; 

f = ((1:N) - r) / (N/2); 

%Linear phase shift
p = exp(-1i*s*pi*f)'; 

%Ifftshift is used to apply the shift

y = x .* ifftshift(p);

if isreal(x); y = real(y); end;

if needtr; y = y.'; end;

ft_psf_r_shifted = y;

yy=ft_psf_r_shifted';
s=y_shift;
p = exp(-1i*s*pi*f)'; 
y2 = yy .* ifftshift(p);
ft_psf_r_shifted = y2';


%Inverse FT

ft_psf_r_shifted = ifft(ft_psf_r_shifted);

return_val =  abs(ft_psf_r_shifted);

end