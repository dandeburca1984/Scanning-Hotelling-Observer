function return_val = get_estimate_of_A(binary_psf,ref_psf)

binary_max_value = max(max(binary_psf));

max_of_ref_star = max(max(ref_psf));

normalised_ref_psf = ref_psf .* (1/max_of_ref_star);

%**************************************************************************

%***********Scale the ref psf to match the real psf************************

scaled_ref_psf = normalised_ref_psf .* binary_max_value;

return_val = sum(sum(scaled_ref_psf));