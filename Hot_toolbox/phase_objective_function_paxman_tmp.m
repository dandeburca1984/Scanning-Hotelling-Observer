function J = phase_objective_function_paxman_tmp(Z_coeffs,D1,D2,pupil_1,pupil_2,lambda_scale_factor,Zernikes,gamma,A_hat,S,pf)
%Gonsalves Least Squares function for phase diversity
%This function computes a metric which is used as an indicator as to the
%fit of the estimated phase

%Z_coeffs = vector of zernike coefficents
%D1 = FT of focal plane intensity measuremtns at lambda 1
%D2 = FT of focal plane intensity measuremtns at lambda 2
%pupil = image of the telescope pupil
%lambda_scale_factor = lambda 1 / lambda 2

%npix = find_pupil_size(pupil_1);  %diameter of the pupil
%pad_factor = (length(pupil_1) - npix)/2;  %factor to pad out the pupil with zeros

phase = S.pup_size;

%phase = ZernPhaseScreen(Z_coeffs,npix); %phase in the pupil

%changing ii=1 tp ii=4

        for ii=4:1:length(Z_coeffs)
            
            phase = phase + Z_coeffs(ii)*Zernikes(:,:,ii);
            
        end

%phase = padarray(phase,[pad_factor pad_factor]);

phase_164 = phase*lambda_scale_factor;  %phase in lambda 2

%two pupil functions
%p1=Phase2psf(phase/S.lam,pf,pupil_1);
p1=Phase2psf_beta(phase/S.lam,pupil_1);
%p1 =need2crop(p1,S,S.pf_lam1);
H1_PSF = p1 ./sum(sum(p1));
H1_OTF = fftshift(fft2(H1_PSF));

%GPF2 = pupil.* exp(2*pi*i*phase_164);
%p2=Phase2psf(phase_164/S.lam2,pf,pupil_1);
p2=Phase2psf_beta(phase_164/S.lam2,pupil_1);
%p2 =need2crop(p2,S,S.pf_lam2);
H2_PSF = p2 ./sum(sum(p2));
H2_OTF = fftshift(fft2(H2_PSF));
%****************************************************

%********Gonsalves Metric****************************
%J =   (( abs( D2.*H1_OTF - D1.*H2_OTF)).^2)  ./ sqrt( (abs(H1_OTF)).^2 + gamma*((abs(H2_OTF)).^2));
J =  ( D2.*H1_OTF - D1.*H2_OTF)  ./ sqrt( (abs(H1_OTF)).^2 + gamma*((abs(H2_OTF)).^2));
%****************************************************

%***Remove Nans************************************
J(isnan(J)) = 0;
%****************************************************

J = (abs(J)).^2;

J = sum(sum(J));

end

