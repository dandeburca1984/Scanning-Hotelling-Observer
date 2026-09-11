function J = phase_diversiy_obj_function_elt_centre(alpha,D1,D2,pupil,up,zernike_polynomials,x,y,t,r)


n = length(alpha);

alpha_d = alpha*up;

Z_d = 0*zernike_polynomials(:,:,1);

Z=0*zernike_polynomials(:,:,1);

%make the test zernike polynomials

for k = 1:n
    
    Z=Z+alpha(k)*zernike_polynomials(:,:,k);
    
    Z_d=Z_d+alpha_d(k)*zernike_polynomials(:,:,k);
    
end

%two phases
phase = zeros(50,50);
phase_164 = phase;

zgrid = gridfit(x,y,Z,50,50);
zgrid_2 = gridfit(x,y,Z_d,50,50);

image_centre=[25 25];

%radius = 0:0.1:31
%theta = 0:0.1:360

for radius=0:1:24
    
    for theta=0:0.1:360;
        
        xx =  round(image_centre(1) + radius * sin(theta));
        
        yy = round(image_centre(2) + radius * cos(theta));
        
        phase(xx,yy) = zgrid(xx,yy);
        phase_164(xx,yy) = zgrid_2(xx,yy);
        
    end
    
end

phase_150 = padarray(phase,[26 25],'replicate','post');
phase_150 = padarray(phase_150,[25 26],'replicate','pre');


phase_164 = padarray(phase_164,[26 25],'replicate','post');
phase_164 = padarray(phase_164,[25 26],'replicate','pre');
%************************************************************

%two pupil functions
GPF = pupil.* exp(i*phase_150);

GPF2 = pupil.* exp(i*phase_164);

%PSFS and OTFS
H1_PSF = (abs(ifftshift (ifft2(GPF))) ).^2;
H1_OTF = fftshift( fft2(H1_PSF));

H2_PSF = (abs( ifftshift(ifft2(GPF2 ))) ).^2;
H2_OTF = fftshift(fft2(H2_PSF));

J =   (( abs( D1.*H2_OTF - D2.*H1_OTF)).^2)  ./ ( (abs(H1_OTF)).^2 + (abs(H2_OTF)).^2);

J(isnan(J)) = 0;

J = sum(sum(J));


