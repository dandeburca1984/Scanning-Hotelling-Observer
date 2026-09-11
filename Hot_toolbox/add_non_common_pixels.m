function New_Sr = add_non_common_pixels(Sr,Sr_tmp,old_pixels,new_pixels)

tmp  = Sr_tmp;

[c, ia, ib] = intersect(old_pixels, new_pixels,'rows');

%tmp_new_pix = new_pixels;

%non_common_pixel_locations = new_pixels( (max(ib)+1) :1: length(new_pixels),:);

%Set common pixels equal to zero

if (isempty(c) == 0)

            for kk=1:1:size(c,1)
            
                tmp(c(kk,1),c(kk,2)) = 0;% max(max(O_T));
            
            end
            
end



New_Sr = Sr + tmp;