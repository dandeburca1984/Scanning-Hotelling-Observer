%Function to get image slice bounded by radius and angle

function [returned_image slice_indeices] = get_image_slice_and_pixels(image,R_min,R_max,Theta_min,Theta_max,image_centre)

%masked_image = zeros(size(image));
returned_image = zeros(size(image));

%for xx=1:1:length(image)
    
  %   x1 = image_centre(1,1) - xx;
    
    %for yy=1:1:length(image)
        
      %  y1 = image_centre(1,2) - yy;
        
        %rho = sqrt( x1.^2 + y1.^2);
        
        %if (rho >= R_min) && (rho<=R_max)
            
          %  masked_image(xx,yy) = image(xx,yy);
            
        %end
        
    %end
    
%end

count=1;
ii=1;
for r=R_min:0.5:R_max
    
        for angle=Theta_min:0.01:Theta_max
            
            slice_indeices(ii,1) = image_centre(1,1) + round( r .* cos(angle));
            slice_indeices(ii,2) = image_centre(1,1) +  round(r.* sin(angle));
            
            %returned_image(slice_indeices(ii,1),slice_indeices(ii,2)) = masked_image(slice_indeices(ii,1),slice_indeices(ii,2));
            returned_image(slice_indeices(ii,1),slice_indeices(ii,2)) = image(slice_indeices(ii,1),slice_indeices(ii,2));
            
            %slice_indeices(ii,1) = x;
            %slice_indeices(ii,2) = y;
            ii=ii+1;
            
        end
        
        count=count+1;
end