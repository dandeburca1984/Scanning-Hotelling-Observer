function stacked_image = get_image_stack(base_image, image_sequence,star_position)


%base image is the first image in the sequence upon which all the other
%images will be stacked
%image sequence is a cell array filled with the image sequence

%Code applies a mask to the center of the image to supress the maxima due
%to the central star

running_some = zeros(size(image_sequence{1,1}));

mask_image = ones(size(base_image));

for ii=star_position(1)-2:1:star_position(1)+2
    
      for  kk = star_position(1)-2:1:star_position(1)+2
    
        mask_image(ii,kk) = 0;
    
      end
    
end



[max_of_image index]=(max( base_image.*mask_image )); 

[max_of_image location] = max(max_of_image);

base_point = [location index(location)];

number_of_images = length(image_sequence);

for kk=1:1:number_of_images
    
    [tmp_max_of_image tmp_index]=(max( image_sequence{1,kk}.*mask_image )); 

    [tmp_max_of_image tmp_location] = max(tmp_max_of_image);

    input_point = [tmp_location tmp_index(tmp_location)]; 
   
    xshift = (base_point(1) - input_point(1));
    
    yshift = (base_point(2) - input_point(2));
    
    shifted_img = shift_planet(image_sequence{1,kk},-xshift,-yshift,1);
    
    running_some = running_some + shifted_img;
    
end

stacked_image = running_some / number_of_images;