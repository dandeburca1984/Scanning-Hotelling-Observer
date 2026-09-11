%Function to get image slice bounded by radius and angle

function returned_image = get_image_slice(image,R_min,R_max,Theta_min,Theta_max,image_centre)

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
for r=R_min:0.1:R_max
    
        for angle=0:0.01:360
            
            x = image_centre(1,1)+ round( r .* cos(angle));
            y = image_centre(1,1)+  round(r.* sin(angle));
            
            %returned_image(x,y) = masked_image(x,y);
            
            returned_image(x,y) = image(x,y);
   
            %slice_indeices(ii,1) = x;
            %slice_indeices(ii,2) = y;
            ii=ii+1;
            
        end
        
        count=count+1;
end

