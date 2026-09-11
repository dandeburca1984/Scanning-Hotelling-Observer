function J = phase_diversiy_obj_function_elt_full(alpha,D1,D2,pupil,up,zernike_polynomials,x,y,t,r,gamma)


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

phase = Z;    
%phase_150=phase;

phase_164 = phase*up;

%************************************************************

%two pupil functions
GPF = pupil.* exp(i*phase);

GPF2 = pupil.* exp(i*phase_164);

%PSFS and OTFS
H1_PSF = (abs(ifftshift (ifft2(GPF))) ).^2;
H1_OTF = fftshift( fft2(H1_PSF));


H2_PSF = (abs( ifftshift(ifft2(GPF2 ))) ).^2;
H2_OTF = fftshift(fft2(H2_PSF));


J =   (( abs( D1.*H2_OTF - D2.*H1_OTF)).^2)  ./ ( (abs(H1_OTF)).^2 + gamma*((abs(H2_OTF)).^2)  );

J(isnan(J)) = 0;

J = sum(sum(J));


