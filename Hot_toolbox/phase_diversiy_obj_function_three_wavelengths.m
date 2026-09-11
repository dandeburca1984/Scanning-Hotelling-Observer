function J = phase_diversiy_obj_function_three_wavelengths(alpha,D1,D2,D3,pupil,scale_1,scale_2,sigma,zernike_polynomials)%,x,y,t,r)


n = length(alpha);

Z=0*zernike_polynomials(:,:,1);

for k = 4:n
    
    Z=Z+alpha(k)*zernike_polynomials(:,:,k); 
    
end

phase = zeros(length(pupil),length(pupil));

%zgrid = gridfit(x,y,Z,254,254);
%zgrid = padarray(zgrid,[129 129]);


%load('pixel_list.mat');
    
   % for xx=1:1:length(pixel_list)
       
      %      phase(pixel_list(xx,1),pixel_list(xx,2)) = zgrid(pixel_list(xx,1),pixel_list(xx,2));
          
    %end

phase_1 = Z;    
%phase_150=phase;

phase_2 = phase.*scale_1;

phase_3 = phase.*scale_2;

%************************************************************

%three pupil functions
GPF = pupil.* exp(1i*phase_1);

GPF2 = pupil.* exp(1i*phase_2);

GPF3 = pupil .* exp(1i*phase_3);

%PSFS and OTFS
H1_PSF = (abs(ifftshift (ifft2(GPF))) ).^2;
H1_OTF = fftshift( fft2(H1_PSF));


H2_PSF = (abs( ifftshift(ifft2(GPF2 ))) ).^2;
H2_OTF = fftshift(fft2(H2_PSF));

H3_PSF = (abs( ifftshift(ifft2(GPF3 ))) ).^2;
H3_OTF =  fftshift(fft2(H3_PSF));

J =((D1.*H2_OTF - D2.*H1_OTF).^2+(D2.*H3_OTF - D3.*H2_OTF).^2 +(D1.*H3_OTF - D3.*H1_OTF).^2 )...
    ./ (  (  (abs(H1_OTF)).^2  + (sigma(1,1)/sigma(1,2))*((abs(H2_OTF)).^2) + ...
    (sigma(1,2)/sigma(1,3))*((abs(H3_OTF)).^2) ) );


%J2 =   (( abs( D1.*H2_OTF - D2.*H1_OTF)).^2)  ./ ( (abs(H1_OTF)).^2 + (abs(H2_OTF)).^2);

J(isnan(J)) = 0;

J = sum(sum(J));


