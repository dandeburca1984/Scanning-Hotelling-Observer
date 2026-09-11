function ref_psf = get_ref_psf(ref_image_seq,target_image,rmin,rmax,dr,g,W,Na,delta_r,star_location)

rad_conversion = pi/180;

Sr = zeros(size(ref_image_seq{1,1}));

count=1;

for r = rmin:dr:rmax
    
        %current_radius = r
    
        %change in angle as a function od radius
        delta_phi = ( (g/2) + ( (2*r)/W)*sqrt( g/(pi*Na)) )^(-1);
        
        %Change in azimuthal width as a function of radius
        azimuthal_width = round( ( r+ (delta_r/2) )*delta_phi) ;
    
    for angle=0:azimuthal_width:360
        
        angle_rad = angle * rad_conversion;
        azimuthal_width_rad = azimuthal_width * rad_conversion;
       
        %area_of_ROI = delta_r*(r + delta_r/2)*delta_phi;
        
        %Define inner radius of Optimization subsection
        %inner_ring = get_circle(r,image_centre);
        
        %Define outer radius of Optimization subsection
        %outer_ring = get_circle(r+delta_r,image_centre);

        %num_theta_values = round (2*pi*r);

        rmin=r; 
        rmax=r+delta_r;
        theta_min=angle_rad;
        theta_max = theta_min + (azimuthal_width_rad);
        
        [O_T pixel_locations]  = get_image_slice_and_pixels(target_image,rmin,rmax,theta_min,theta_max,star_location);
        
        tmp = O_T;
        
            for kk=1:1:length(pixel_locations)
            
                tmp(pixel_locations(kk,1),pixel_locations(kk,2)) = 1;% max(max(O_T));
            
            end
        
        [row,col] = find(tmp);
        
        pixel_locations = [row,col];
        
        ck=get_coefficients(rmin,rmax,theta_min,theta_max,O_T,ref_image_seq,star_location,pixel_locations);
        
        Sr_tmp = get_subtraction_subsection(ck,ref_image_seq,rmin,theta_min,theta_max,star_location,dr);
        
        %The same pixels are adding onto each other in the centre
        if (count==1)
        New_Sr = Sr + Sr_tmp;
        
        elseif(count>1)
            
        old_pixels= get_pixels(New_Sr);
        
        new_pixels = get_pixels(Sr_tmp);
       
        %Add non-common_pixels
        New_Sr = add_non_common_pixels(New_Sr,Sr_tmp,old_pixels,new_pixels);
        
        end
        count=count+1;
        
    end

     
end

ref_psf = get_positive_values_only(New_Sr);