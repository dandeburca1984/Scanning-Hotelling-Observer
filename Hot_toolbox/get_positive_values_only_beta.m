function positive_image = get_positive_values_only_beta(image)

image(image<0) = 0;
positive_image = image;
