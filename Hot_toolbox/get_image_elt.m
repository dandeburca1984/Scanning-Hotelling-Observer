function image = get_image_elt(pupil,alpha,zz,x,y)


Z_hat=0*zz(:,:,1);

for k = 1:length(alpha)
    
    Z_hat=Z_hat+alpha(k)*zz(:,:,k);
    
end

phase = zeros(50,50);
zgrid = gridfit(x,y,Z_hat,50,50);
image_centre=[25 25];

for radius=0:1:24
    
    for theta=0:0.1:360;
        
        xx =  round(image_centre(1) + radius * sin(theta));
        
        yy = round(image_centre(2) + radius * cos(theta));
        
        phase(xx,yy) = zgrid(xx,yy);
        
    end
    
end

big_phase = padarray(phase,[26 25],'replicate','post');
big_phase = padarray(big_phase,[25 26],'replicate','pre');


GPF_hat = pupil.* exp(i*big_phase);

image = (abs( ifftshift(ifft2(GPF_hat))) ).^2;