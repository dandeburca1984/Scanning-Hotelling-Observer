function positive_image = get_positive_values_only(image)

image_size = size(image);

%indices

aa = image_size(1);
bb = image_size(2);

for ii=1:1:aa
       for jj=1:1:bb
            
           if image(ii,jj) <= 0
                
               image(ii,jj) = 0;
               
               image(ii,jj) = ceil(image(ii,jj));
                
           end
        end
end

positive_image = image;