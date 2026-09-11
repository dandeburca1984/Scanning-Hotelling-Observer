function return_val = get_A_least_squares(image,ref_psf)

normalised_ref_psf = ref_psf ./ sum(sum(ref_psf));

%************Use optimal method as a control***
options=optimset('LargeScale','off');
options = optimset(options,'TolX',1e-4);
options = optimset(options,'Display','off');


%Function handle for Hotelling observer
%f =
%@(r_pl)hot_function(r_pl,r_covariance,r_residual,normalised_star_image,star_position);

A_start = get_estimate_of_A(image,ref_psf);

f = @(A)get_image_residual(image,A,normalised_ref_psf);

A_least_squares = fminunc(f,A_start,options);

return_val = A_least_squares;


