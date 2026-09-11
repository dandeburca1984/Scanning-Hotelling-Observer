clc;
clear;

load('Hotelling_data_map.mat');
load('test_map.mat');
%load('sim_data.mat');
load('optimisation_history.mat');
load('planet_loc.mat');
set(0,'Units','pixels') 
scrsz = get(0,'ScreenSize');

%hotelling = zeros(length(possible_planet_locations),1);
%int = hotelling;

%for zz=1:1:length(possible_planet_locations)
            
  %          hotelling(zz,1) =  Hot(zz);
            
    %        int(zz,1) = a_pl(zz);
%end
intensity_map(125:133,125:133)=0;

%stat_image = get_un_raster_image(hotelling);
stat_image = intensity_map;
%int_image = get_un_raster_image(int);
int_image = intensity_map;

fig=figure('Position',[500 scrsz(4)/4 700 700]);
set(fig,'DoubleBuffer','on');
set(gca,'xlim',[-500 500],'ylim',[-500 500],...
    'NextPlot','replace','Visible','off')

mov = avifile('example.avi');

for noise_ring=125:-5:50

[noise_annulus,noise_array,mean_of_dist,sigma]=get_noise_distribution(int_image,noise_ring);

index = 1;

tmp = int_image;

for ii = 1:length(noise_array)
    
    %noise_array(ii) = image(noise_annulus(ii,1),noise_annulus(ii,2));
    
    tmp(noise_annulus(ii,1),noise_annulus(ii,2)) = max(max(int_image));
    
    index = index + 1;
    
    
end

%figure('Position',[550 scrsz(4)/3 500 500]); imagesc(tmp);

[n,xout]= hist(noise_array,30);

%figure('Position',[50 scrsz(4)/3 500 500]);bar(xout,n);

mask = ones(size(int_image));

star_location = [129 129];

for kk=star_location(1)-1:1:star_location(1)+1
    
        for jj=star_location(1)-1:1:star_location(1)+1
    
                mask(kk,jj)=0;
                
        end
                
end

masked_image = int_image .* mask;

thresholded_image = (int_image .* ( int_image > max(xout))) .* mask ;

%figure('Position',[1100 scrsz(4)/3 500 500]);  imagesc(thresholded_image);

h=fig;
subplot(2,2,1); imagesc(tmp); title('\bf Intensity map','Fontsize',14);
subplot(2,2,2); bar(xout,n); title('\bf Distribution of H_0 Values','Fontsize',14);
subplot(2,2,3); imagesc(thresholded_image(90:170,90:170)); title('\bf Thresholded Image','Fontsize',14); 
text(80,20,'\bf Theshold equals the maximum','Fontsize',14); 
text(80,30,'\bf of H_0 distribution','Fontsize',14);

F = getframe(gca);
mov = addframe(mov,F);

end


mov = close(mov);

regional_maxima_position = imregionalmax(thresholded_image,4);

image_of_regional_maxima = thresholded_image .* regional_maxima_position;

A = 3.2483e6;

index=1;

for vector_x=1:1:length(image_of_regional_maxima)
    
    for vector_y=1:1:length(image_of_regional_maxima)
        
        if image_of_regional_maxima(vector_x,vector_y)>0
            
            star_list.position{index} = [vector_x vector_y];
            
            star_list.magnitude{index} = 2.5*log10( A/ image_of_regional_maxima(vector_x,vector_y));
            
            index=index+1;
            
        end
    end
end
    