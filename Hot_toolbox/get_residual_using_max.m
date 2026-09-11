function residual = get_residual_using_max(data,psf)

psf_max_normalised = psf ./ max(max(psf));

data_max = max(max(data));

residual = get_positive_values_only(data - data_max*psf_max_normalised);

end