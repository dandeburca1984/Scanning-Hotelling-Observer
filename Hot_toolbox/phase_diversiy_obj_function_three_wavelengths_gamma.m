function J = phase_diversiy_obj_function_three_wavelengths_gamma(alpha,D1,D2,D3,pupil,scale_1,scale_2,sigma,zernike_polynomials,S)%,x,y,t,r)

%alpha_full = zeros(1,20);

%alpha_full(1,5:20) = alpha(1,16);

H1_PSF = Gen1Image(S.lam,S.D,S.scale,256,real(alpha));
H1_PSF = H1_PSF ./ sum(sum(H1_PSF));
H1_OTF = fftshift( fft2(H1_PSF));

alpha_beta = real(alpha) * scale_1;
H2_PSF = Gen1Image(S.lam2,S.D,S.scale,256,alpha_beta);
H2_PSF = H2_PSF ./ sum(sum(H2_PSF));
H2_OTF = fftshift(fft2(H2_PSF));

alpha_gamma = real(alpha) *scale_2;
H3_PSF = Gen1Image(S.lam3,S.D,S.scale,256,alpha_gamma);
H3_PSF = H3_PSF ./ sum(sum(H3_PSF));
H3_OTF =  fftshift(fft2(H3_PSF));

%J =((D1.*H2_OTF - D2.*H1_OTF)+(D2.*H3_OTF - D3.*H2_OTF) +(D1.*H3_OTF - D3.*H1_OTF) )...
   % ./ (  sqrt(  (abs(H1_OTF)).^2  + (sigma(1,1)/sigma(1,2))*((abs(H2_OTF)).^2) + ...
    %(sigma(1,1)/sigma(1,3))*((abs(H3_OTF)).^2) ) );

gamma = sigma(1,1)/sigma(1,2);
    
%J =  ( D2.*H1_OTF - D1.*H2_OTF)  ./ sqrt( (abs(H1_OTF)).^2 + gamma*((abs(H2_OTF)).^2));

J=((D1 - H1_PSF)).^2 + ((D2 - H2_PSF)).^2 + ((D3 - H3_PSF)).^2;

J(isnan(J)) = 0;

J = sum(sum(J));


