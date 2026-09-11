function un_raster_image = get_un_raster_image_two(image)

image_size = size(image);

x_size = sqrt(image_size(2));

un_raster_image = image(1,1:x_size);

for index=2:1:x_size
    
    tmp_image = image(1,((index-1)*x_size+1) :index*x_size);
    
    un_raster_image = cat(1,un_raster_image,tmp_image);
    
end