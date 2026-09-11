function results = Scanning_Hotelling_Observer(varargin)

% Scanning_Hotelling_Observer impliments the Hotelling observer as
% described in Burke et al 2009 PASP 121, 767-777
%
% Inputs
% data, the pre-processed data to be scanned, assumed to be composed of a
% bright parent signal and a faint companion signal
%
% psf, the psf of the system used to locate and estimate the differential
% intensity of the faint companion source
%
% Outpust
% results, structure function containing the location and intensity of the
% parent and companion signals, the data covariance matrix and the
% background estimates
%

%% Check the inputs
%The switch and case statements check that the required inputs are met
switch nargin
    case 0
        error('No inputs, need at least two! \n');
    case 1 
        error('Only one input, need at least two! \n');
    case 2
        fprintf('Two inputs given\n');
        data = varargin{1,1};
        psf = varargin{1,2};
end

warning off;

%% Setup options file for the conjugate minimisation
options=optimset('LargeScale','on');
options = optimset(options,'TolX',2*(eps));
options = optimset(options,'TolF',eps);
options = optimset(options,'Display','off');

%% Preliminaries on data
data = get_positive_values_only_beta(data); %set all negative values to zero
data = data ./ sum(sum(data));  %normalise the total intensity of the data to 1

psf = get_positive_values_only_beta(psf); %set all negative values to zero
psf = psf ./ sum(sum(psf)); %normalise the total intensity of the psf to 1

%% Noise radial vector definition
% A radial matrix is defined to index the location of the noise floor
half_data_size = round( length(data)/2);
tx = linspace (-half_data_size, half_data_size, length(data));
[xx, yy] = meshgrid (tx, tx);
r = round( sqrt (xx .^ 2 + yy .^ 2));

%% Locate the approxiamte location of the parent star
% Peak of Cross-Corellation between data and PSF, approx location of Parent
%Star
[col_hat, row_hat] = m_filter_2(data,psf);   

%% Find exact location of Parent Star using conjugate gradient minimisation 
start_point = [col_hat, row_hat];                         % Start point of minimisation
covariance_matrix =ones(64,64);                        % Assumed covariance of ones
r_cov = get_raster_image(covariance_matrix);    % raster scan of covariance matrix    
f_prime = @(r_pl)hotelling_test_stat(r_pl,r_cov,get_raster_image(data),psf); 
[spot_position f_at_position] = fminunc(f_prime,start_point,options); %conjugate gradient minimisation  
shifted_psf = shift_planet(psf,spot_position(1,1),spot_position(1,2));  %shift the PSF to the location of the parent star
parent_star_intensity_hat = get_raster_image(data) / get_raster_image(shifted_psf); %estimate the intensity of the Parent star, it should be very close to 1 
data_minus_parent = get_positive_values_only_beta( data - shifted_psf*parent_star_intensity_hat);    %subtract estimated Parent star signal from the data

%% Estimate the `flat' background level  
n_ring = round(length(data_minus_parent)/2)-2; % at the radius n_ring away from the center of the image, the noise is assumed to be flat
[background_hat, variance_hat] = get_noise(data_minus_parent,r,n_ring); %estimate background parameters

%% Construct the Gaussian Covariance model
covariance_matrix =shifted_psf*parent_star_intensity_hat+ background_hat+  variance_hat^2;
r_cov = get_raster_image(covariance_matrix); %raster scan the data

%% Locate the approximate location of the faint companion
[col_hat, row_hat] = m_filter_2(data_minus_parent,psf); 

%% Find the exact location of the companion star
start_point =round(abs( [col_hat, row_hat])) ;                         % Start point of minimisation
%Conjugate gradient minimisation used to find the sub-pixel location of the
%companion star and its intensity
[x_hat,y_hat,int_hat,test_statistic1] = ML_Hotelling_estimator_analytical_int(data_minus_parent,psf,options,start_point,r_cov);


%% update the background model
% The estimated Parent star signal and the companion signal is subtracted
% from the data and the background model is updated
companion_psf = shift_planet(psf,x_hat,y_hat); %psf moved to estimated companion location
updated_background = get_positive_values_only_beta( data_minus_parent - companion_psf*int_hat); % residual background image
[updated_background, updated_background_variance] = get_noise(updated_background,r,n_ring); %updated backgroudn model
new_cov  = shifted_psf*parent_star_intensity_hat + updated_background + updated_background_variance; %updated covariance matrix
r_new_cov = get_raster_image(new_cov);
%Search for the location (and hence intensity) of the companion again using
%the updated background model, test_statistic2 should be greater than test_statistic1 
[x_hat2,y_hat2,int_hat2,test_statistic2] = ML_Hotelling_estimator_analytical_int(data_minus_parent,psf,options,start_point,r_new_cov);
delta_m2 =  2.5*log10(parent_star_intensity_hat / int_hat2); %differential magnitude

%% Construct Output variable
results.parent_location = spot_position;
results.paretn_intensity =parent_star_intensity_hat;
results.covariance_matrix = new_cov;
results.background_level = updated_background;
results.background_variance = updated_background_variance;
results.companion_location = [x_hat2,y_hat2];
results.companion_intensity = int_hat2;
results.companion_differential_magnitude = delta_m2;

%% Define subfunctions
%the subfunctions are located here at the end of file to facilatate ease of
%use of the code in general, with the subfunctions located here there is no need
%for a folder of subfunctions
    function data = get_positive_values_only_beta(data)
    %this function locates all the negative values in a data array and sets them to zero    
        data(data<0) = 0;
    end

    function [estimated_col,estimated_row]=m_filter_2(data,psf)
    %this function computes the cross-corellation (CC) between the data and the psf using Fourier techniques    
        A_=(ifftshift(fft2(fftshift((data))))).*(ifftshift(fft2(fftshift((psf)))));
        n= length(data);
        pad = 0*n; %padding is used in the Fourier approx to negate biasing the estimates towards the center of the CC
%         A_=padarray(A_,[pad pad],0);
        A=(ifftshift(ifft2(fftshift(A_))));
        A=sqrt(A.*conj(A));%This is the correlation.
        [a,b]=find(A==max(max(A))); %locate the peak (integer values)
        a=min(a);a = round(a);b=min(b);b = round(b);%just in case.
        %Parabolic interpolation of the position of the maxima, as described by
        %Poyneer (Applied Optics 2003).
        a_=a-(n+2*pad+3)/2+0.5.*(A(a-1,b)-A(a+1,b))./(A(a-1,b)+A(a+1,b)-2.*A(a,b));% this formula in particular is valid for uneven sizes of images.
        b_=b-(n+2*pad+3)/2+0.5.*(A(a,b-1)-A(a,b+1))./(A(a,b-1)+A(a,b+1)-2.*A(a,b));
        estimated_row=a_.*n./(n+2*pad);
        estimated_col=b_.*n./(n+2*pad);
    end

    function [noise_hat, var_hat] = get_noise(data,r,radius)
    %this function estimates the background noise level and its variance at a given radius    
        r = r == radius;
        noise_array = data(r);
        noise_hat = mean(noise_array);
        var_hat = (std(noise_array))^2;
    end

    function raster_data = get_raster_image(data)
    %this function creates a 1D raster scanned image from 2D data    
        x_size = length(data);
        raster_data = zeros(1,x_size^2);
            for index=1:1:x_size
               raster_data(1, (x_size*(index-1) +1) : x_size*index) = data(index,:);
            end       
    end

    function [x_hat,y_hat,int,t] = ML_Hotelling_estimator_analytical_int(data,psf,options,start_point,r_covariance)
    %this function implements the conjugate gradient minimisation of the hotelling observer    
        residual = get_positive_values_only_beta(data);
        r_hotelling_image = get_raster_image(residual);
        f_prime = @(r_pl)hotelling_test_stat(r_pl,r_covariance,r_hotelling_image,psf); %function handle for the unconstrained minimisation
        [spot_position1 f_at_position] = fminunc(f_prime,start_point,options); %unconstrained minimisation 
        [spot_position_h f_at_position] = fminunc(f_prime,spot_position1,options); %minimisation re-run to confirm location of the minima
        [t,int] = hotelling_test_stat(spot_position_h,r_covariance,r_hotelling_image,psf); %function rerun to output the test statistic and the estimated intensity
        t = abs(t);
        x_hat = spot_position_h(1,1);
        y_hat = spot_position_h(1,2);
    end

    function [f,a_pl]  = hotelling_test_stat(r_pl,r_covariance,residual_H1,psf)
    %this function computes the hotelling test statistic, using the analytic expression for the intensity of the signal    
        temp = shift_planet(psf,r_pl(1,1),r_pl(1,2)); %move psf to test location
        r_temp = get_raster_image(temp); %raster scan the image
        %compute the analytic estiamte of the signal intensity
        a_pl  =( (r_temp ./ r_covariance .*  residual_H1) / ...
                    (r_temp ./ r_covariance .* r_temp) );
        %calcualte the value of the Hotelling test statistic        
        f = -( sum( ( (a_pl*r_temp) ./ r_covariance) ...
            .*(residual_H1 - 0.5*a_pl*r_temp)));
    end


    function shifted_psf = shift_planet(psf,xshift,yshift)
    %this function moves the location of the psf/data, with sub-pixel
    %accuracy, to any location in the image using Fourier techniques
        ft_psf = fft2(psf); %Fourier transform the image
        %Apply the desired X shift to the data
        xmax = size(psf,1);
        s = xshift;
        x = ft_psf(1,:);
        needtr = 0; 
        if size(x,1) == 1; x = x(:); needtr = 1; end;
        N = size(x,1); 
        r_ = floor(N/2) + 1; 
        f = ((1:N) - r_) / (N/2); 
        p = exp(-1i*s*pi*f); %Linear phase shift 
        p_ = ifftshift(p); %ifftshift use to apply the shift to the data
        p_2 = repmat(p_,[xmax 1]); %Replicate and tile array shift values
        ft_psf_shifted = ft_psf .* p_2; %apply to FT data
      
        %Apply the desired Y shift to the data
        s = yshift;
        x = ft_psf_shifted(:,1);
        needtr = 0; 
        if size(x,1) == 1; x = x(:); needtr = 1; end;
        N = size(x,1); 
        r_ = floor(N/2) + 1; 
        f = ((1:N) - r_) / (N/2); 
        p = exp(-1i*s*pi*f)';  %Linear phase shift
        p_ = ifftshift(p);%ifftshift use to apply the shift to the data
        p_2 = repmat(p_,[1 xmax]);%Replicate and tile array shift values
        ft_psf_shifted = ft_psf_shifted .* p_2;%apply to FT data

        ft_psf_shifted = ifft2(ft_psf_shifted); %Phase shifted spectrum is Inverse FT to return to space domain     
        shifted_psf =  abs(ft_psf_shifted);
        shifted_psf = (shifted_psf /sum(sum(shifted_psf)) ) * sum(sum(psf));% the FT can change the values of the data
    end


end