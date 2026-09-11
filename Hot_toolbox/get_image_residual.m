function residual = get_image_residual(image,A,psf)


residual = sum(sum(     (image - A*psf)^2  ));