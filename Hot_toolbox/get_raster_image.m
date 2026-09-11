function raster_image = get_raster_image(image)

image_size = size(image);

x_size = image_size(1);

raster_image = image(1,:);

for index=2:1:x_size
    
    tmp_image = image(index,:);
    
    raster_image = cat(2,raster_image,tmp_image);
    
end