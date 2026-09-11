function [estimated_psf,estimated_zern_coeff,starting_guesses] = get_psf_from_phase_diversity_large_pupil(D1,D2,pupil,lambda_1,lambda_2)

clc;

addpath('Hot_toolbox/');

%Pupil Function
pupil = fix_pupil_levels(pupil);

%The two imaged wavelengths
%lambda_1 = 1.64*10^-6;
%lambda_2 = 1.8*10^-6;

%Phase shift equals the ratio of the wavelengths
up=lambda_1/lambda_2;

%changed from 1/2 to 2/1

scale = ((2*pi) / 1.5);
%****************************************************
%*********Define the Zernike Polynomials*************

N = length(pupil);
%N=512;%dimension of theta-R grid
r=linspace(0,1,N);t=linspace(0,2*pi,N);
[t,r]=meshgrid(t,r);
[x,y] = pol2cart(t,r);
p = 4:55;%indexes of Zernike functions to be used in the decomposition
n=length(p);
z = zernfun2(p,r(:),t(:),'norm');
%z = zernfun2(p,r(:),t(:));
zz=reshape(z,N,N,length(p));
%n = length(alpha);

%Star Position
%[max_of_image index]=(max(D1)); 
%[max_of_image location] = max(max_of_image);
%star_position = [location index(location)];

%Strehl Ratio
%sr = max(max(D1)) / max(max(airy));

%***************************************************

D1_FT = fftshift(fft2(D1));
D2_FT = fftshift(fft2(D2));

%***********Options for Minimization****************
options=optimset('LargeScale','off');
options = optimset(options,'TolX',eps/2);
options = optimset(options,'TolF',eps/2);
options = optimset(options,'Display','on');
options = optimset(options,'outputfcn',@outfun,'display','iter');
options = optimset(options,'MaxFunEvals',10^6);
options = optimset(options,'MaxIter',10^6);

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
searchdir = [];
step_size=[];

%Phase Diversity Objective Function
j = @(alpha)phase_diversiy_obj_function_elt_full(alpha,D1_FT,D2_FT,pupil,up,zz,x,y,t,r);

%Inital Guess at Zernike Coeff
alpha_prime = randn(1,52);%zeros(1,52);
alpha_prime = ( alpha_prime ./ max(alpha_prime)) .* (lambda_1*1e-6);

starting_guesses = alpha_prime;

[alpha_hat,fval,exitflag,output] = fminunc(j,alpha_prime,options);  
    
%Plot the estimated phase
Z=0*zz(:,:,1);

for k = 1:n
    
    Z=Z+alpha_hat(k)*zz(:,:,k);
    
end

figure;surf(x,y,Z,'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong');
%*****************************************************
%Estimated PSF
H_hat = Phase2psf_beta(Z,pupil);
%H_hat = get_image_elt(pupil,alpha_hat,zz,x,y);
H_hat = H_hat ./ sum(sum(H_hat));

estimated_psf=H_hat;
estimated_zern_coeff = alpha_hat;

function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
         % Concatenate current search direction with 
         % searchdir.
           searchdir = [searchdir;... 
                        optimValues.searchdirection'];
           step_size = [step_size;   optimValues.stepsize];
           %plot(x(1),x(2),'o');
         % Label points with iteration number and add title.
           text(x(1)+.15,x(2),... 
                num2str(optimValues.iteration));
           title('Sequence of Points Computed by fmincon');
         case 'done'
             hold off
         otherwise
     end
end
end
